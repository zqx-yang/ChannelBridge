#!/usr/bin/env python3
"""Stage 1a: find the best available structure for each enzyme, in order.

This is the **input preparation** decision, and it comes before everything
else: every tunnel, patch and pose is a statement about the coordinates it
was computed on, so choosing those coordinates is the most consequential
choice in the pipeline and the one most often made by accident.

The order is deliberate, cheapest and most-evidential first:

1. **An experimental structure of this exact protein** (PDB). Best case: real
   coordinates, usually with the cognate ligand bound, and possibly a
   ChannelsDB annotation saying what its tunnels are *for*.
2. **AlphaFill**, if a predicted model exists with the cofactor transplanted.
   AlphaFill takes an AlphaFold DB model and transplants ligands from
   homologous holo structures, reporting the donor, the local identity, the
   RMSD and any clashes. That matters here out of proportion to its size: 21
   of these 24 enzymes define their channel start point from a bound ligand,
   and a bare prediction has none.
3. **AlphaFold DB**, if a model exists for the UniProt entry. Ligand-free, so
   the start point has to come from catalytic residues instead
   (`chbridge caver config --start-cavity`).
4. **Predict it** with AlphaFold3, with an MSA. Only when nothing above
   exists.

**A homolog of a different protein is not on this list**, and that is the
change this script records. The first version of this set used the nearest
deposited homolog for every enzyme without its own structure, which put six
enzymes on distant homologs (down to 0.41 sequence identity) and silently
made four pairs of enzymes share one entry -- so three of the 19 "pairs" were
the same computation twice. A prediction from the enzyme's own sequence has
neither problem.

The homolog does not become useless, it changes job: it stays the
**annotation donor**. A fresh prediction has no ChannelsDB entry, so the
category that says "this tunnel leads to the cognate substrate" has to be
carried over by superposing the annotated homolog and seeding CAVER from the
site it identifies -- which is exactly what `run_examples.py tunnels` already
does for a fetched chain that stops before the surface.

Why prediction is defensible here, when the same project found AlphaFold
geometry unreliable for iMARS fusions: that failure was **inter-domain**. In
`p6z`, 15 AlphaFold2 models of one fusion spanned 0.36 DO because the
predicted placement of the two domains relative to each other had a PAE of
28 A against a ceiling of 31 -- a draw from the prior, not a measurement.
Intra-domain PAE was 3.0-5.7 A. ChannelBridge never takes the relative
placement from a prediction; it searches for it. One enzyme at a time is the
regime where prediction works.

And the obvious objection -- that an apo model's side chains collapse the
channel -- was tested and refuted: holo crystal chains against apo models of
the same enzymes gave cluster-1 bottlenecks of 1.92 vs 1.55 A (4CL), 1.94 vs
2.04 (STS) and 2.26 vs 2.30 (CHS). They agree to about 0.1 A, and the apo
model is wider in two of three.

Usage:

    ./find_structures.py              # report what exists for all 24
    ./find_structures.py --tsv        # structures.tsv, the decision record
    ./find_structures.py --fetch      # download what the decisions chose
    ./find_structures.py --only AS    # one enzyme

Needs network access. Nothing is downloaded without --fetch: this reports what is available
and what must be predicted.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from run_examples import enzymes  # noqa: E402

#: Accession per enzyme, from iMARS Table S2G ("Accession number"), split out
#: of its paired "X and Y" cells. Recorded here rather than re-derived because
#: the workbook is copyrighted and not redistributed; `--verify-against` can
#: check these against your own copy.
ACCESSION = {
    "4'OMT": "Q9LEL5.1", "4CL": "AY376729", "7OMT": "FJ156103.1",
    "ADS": "AY006482", "AS": "Q38710.1", "BAS": "Q94FV7.1",
    "BBS": "AIG92846.1", "CHS": "AY237728.1", "CTS": "A8CDR5.1",
    "CUS": "AKG06601.1", "ERG20": "NP_012368.1", "GES": "Q6USK1.1",
    "GGPPS": "GMC29149.1", "GPPS": "AAN01134.1", "LIC": "L13459",
    "MYS": "CAC41012", "NES": "APB87285.1", "OCS": "A4FVP2.1",
    "RAS": "AJA31666.1", "REPI": "AKO60175.1", "RNMT": "AOR51552.1",
    "SQS": "ADF40697.1", "STS": "ACC76753", "TS": "AAC49310",
}

#: What a UniProt accession looks like, per their own documented format.
_UNIPROT_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$")

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
AFDB_PREDICTION = "https://alphafold.ebi.ac.uk/api/prediction/%s"
ALPHAFILL_JSON = "https://alphafill.eu/v1/aff/%s/json"
RCSB_UNIPROT = "https://data.rcsb.org/rest/v1/core/uniprot/%s/1"
ALPHAFILL_CIF = "https://alphafill.eu/v1/aff/%s"
AFDB_CIF = "https://alphafold.ebi.ac.uk/files/AF-%s-F1-model_v%s.cif"
RCSB_CIF = "https://files.rcsb.org/download/%s.cif"

#: Where fetched coordinates go. Beside the deposited entries, and ignored by
#: git for the same reason they are: this set commits decisions, not data.
ENTRIES = os.path.join(HERE, "_entries")

#: Clashes per transplanted heavy atom above which a transplant is not
#: treated as marking a real cavity.
#:
#: A raw clash count cannot be compared across ligands -- SAM has 27 heavy
#: atoms and GST 18 -- so it is normalised. Above this the ligand overlaps the
#: model's own side chains badly enough that its centroid is probably not in
#: free space, so it cannot give a channel start point and the point has to be
#: derived from the lining residues instead (`caver config --start-cavity`).
#:
#: This also decides whether AlphaFill WINS as a source: it is preferred over
#: a plain model because a ligand-derived start point survives, and a
#: transplant this bad voids that reason.
CLASH_PER_ATOM_MAX = 2.0

#: Ligands AlphaFill transplants that are never a channel's cargo -- ions and
#: cryoprotectants. A transplant of one of these is not a usable start point.
NOT_CARGO = {"HOH", "GOL", "EDO", "PEG", "SO4", "PO4", "CL", "NA", "K", "MG",
             "MN", "ZN", "CA", "FE", "NI", "CO", "CD", "ACT", "DMS", "TRS"}


def _get(url, *, tries=3, timeout=45):
    """GET with retries. Returns bytes, or None on a 404."""
    last = None
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=timeout) as fh:
                return fh.read()
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return None
            last = exc
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            last = exc
        time.sleep(1.5 * (attempt + 1))
    raise RuntimeError("%s: %s" % (url, last))


def to_uniprot(accession):
    """UniProt accession for a GenBank/RefSeq/UniProt accession, or None.

    A versioned UniProt accession (``Q9LEL5.1``) is used as-is; anything else
    is looked up as a cross-reference, which is how a nucleotide accession
    like ``AY376729`` finds its protein.
    """
    bare = accession.split(".")[0]
    # `accession:` only accepts a real UniProt accession and 400s the whole
    # query otherwise, which takes the `xref:` half down with it.
    query = ("accession:%s" % bare) if _UNIPROT_RE.match(bare) else ("xref:%s" % bare)
    url = "%s?%s" % (UNIPROT_SEARCH, urllib.parse.urlencode({
        "query": query,
        "fields": "accession,id,protein_name,organism_name,length,xref_pdb",
        "format": "tsv", "size": "1"}))
    body = _get(url)
    if not body:
        return None
    lines = body.decode("utf-8").strip().split("\n")
    if len(lines) < 2:
        return None
    f = lines[1].split("\t")
    f += [""] * (6 - len(f))
    return {"uniprot": f[0], "entry": f[1], "name": f[2],
            "organism": f[3], "length": int(f[4]) if f[4] else 0,
            # UniProt's own curated mapping of which PDB entries contain THIS
            # protein. This is the authority for "does an experimental
            # structure of this exact enzyme exist" -- see decide().
            "pdb": tuple(x for x in f[5].replace(";", " ").split() if x)}


def entry_contains(pdb_id):
    """What protein a PDB entry actually holds, per RCSB's own mapping.

    The check that matters most in this whole script, because it is the one
    the first version of this set did not do. A sequence search returns a
    ranked list and a score; it does not tell you that the entry you picked
    is a different enzyme. RCSB maps each polymer entity to a UniProt entry,
    so this answers "what did I actually get" directly.
    """
    try:
        body = _get(RCSB_UNIPROT % pdb_id.upper())
    except RuntimeError:
        return None
    if not body:
        return None
    try:
        doc = json.loads(body)
    except ValueError:
        return None
    rows = doc if isinstance(doc, list) else [doc]
    for r in rows:
        acc = (r.get("rcsb_uniprot_container_identifiers") or {}).get("uniprot_id")
        prot = r.get("rcsb_uniprot_protein") or {}
        name = (prot.get("name") or {}).get("value")
        org = (prot.get("source_organism") or {}).get("scientific_name")
        if acc:
            return {"uniprot": acc, "name": name or "", "organism": org or ""}
    return None


def afdb(uniprot):
    """``{version, plddt}`` if AlphaFold DB has a model, else None."""
    body = _get(AFDB_PREDICTION % uniprot)
    if not body:
        return None
    try:
        rows = json.loads(body)
    except ValueError:
        return None
    if not rows:
        return None
    r = rows[0]
    return {"version": r.get("latestVersion"), "plddt": r.get("globalMetricValue"),
            "model": r.get("modelEntityId", "")}


def alphafill(uniprot):
    """Transplanted ligands AlphaFill has for this entry, best identity first."""
    try:
        body = _get(ALPHAFILL_JSON % uniprot)
    except RuntimeError:
        return None
    if not body:
        return None
    try:
        doc = json.loads(body)
    except ValueError:
        return None
    out = []
    for hit in doc.get("hits") or []:
        # every one of these can be present-and-null in a real response
        ident = float((hit.get("alignment") or {}).get("identity") or 0.0)
        for t in hit.get("transplants") or []:
            lig = str(t.get("analogue_id", "")).strip().upper()
            if not lig or lig in NOT_CARGO:
                continue
            clash = t.get("clash") or {}
            n_at = int(clash.get("transplant_atom_count") or 0)
            n_cl = int(clash.get("clash_count") or 0)
            out.append({"ligand": lig, "donor": hit.get("pdb_id", ""),
                        "identity": ident,
                        "rmsd": float(hit.get("global_rmsd") or 0.0),
                        "clashes": n_cl,
                        "atoms": n_at,
                        # A raw clash count cannot be compared across ligands
                        # (SAM has 27 heavy atoms, GST 18), so normalise it.
                        "per_atom": (n_cl / n_at) if n_at else 0.0,
                        # How far the transplant's own neighbourhood is from
                        # the donor's: a large value means the site it was
                        # placed into is not the site it came from.
                        "local_rmsd": float(t.get("local_rmsd") or 0.0)})
    out.sort(key=lambda h: (-h["identity"], h["clashes"]))
    # one row per ligand, keeping its best-identity donor
    seen, keep = set(), []
    for h in out:
        if h["ligand"] in seen:
            continue
        seen.add(h["ligand"])
        keep.append(h)
    return {"version": doc.get("alphafill_version"), "date": doc.get("date"),
            "ligands": keep}


#: Below this fraction of its reference span observed, or with an unobserved
#: stretch longer than MAX_GAP, an experimental entry is not usable for channel
#: work however right the protein in it is.
#:
#: Holding the right protein is necessary and not sufficient: an entry can
#: hold it and not hold all of it. 3TSY is the only experimental structure of
#: this set's 4CL and it is a 4CL::stilbene-synthase FUSION at 3.1 A in which
#: the last 102 residues of 4CL were never modelled -- and 4CL is an
#: ANL-superfamily adenylate-forming enzyme, where the C-terminal domain is
#: the one that closes over the active site. A tunnel computed from that would
#: run into open solvent through a hole where a domain should be, and would
#: look like a wide channel rather than like a missing domain.
#:
#: The coverage is taken over the entry's ALIGNED span, not over the whole
#: UniProt length: 2ONG covers residues 57-599 of limonene synthase because
#: the chloroplast transit peptide is cleaved, which is correct biology and
#: must not read as 10% missing.
MIN_COVERAGE = 0.90
MAX_GAP = 30

#: Decisions that override the generic order, with the reason. An override
#: lives HERE and not as an edit to structures.tsv, because that file is
#: generated: a hand-edit to it is lost the next time anyone runs --tsv, and
#: silently, which is the worst way to lose a decision.
#:
#: Each entry must say what the generic order would have chosen and why that
#: is wrong for this enzyme. An override with a reason that applies to every
#: enzyme is a rule, and belongs in decide() instead.
OVERRIDE = {
    "4CL": ("af3",
            "the generic order would take the AlphaFold DB model (pLDDT 89.8) "
            "because one exists. 4CL is a two-domain ANL-superfamily "
            "adenylate-forming enzyme whose C-terminal domain closes over the "
            "active site, so the channel depends entirely on the INTER-domain "
            "angle -- which is the one thing this project measured AF2 to be "
            "unreliable at (15 AF2 models of one fusion spanned 0.36 DO on an "
            "inter-domain PAE of 28 A against a 31 A ceiling, while "
            "intra-domain PAE was 3.0-5.7 A). The AlphaFold DB model is "
            "AF2-era; AF3 with an MSA is predicted fresh and its PAE can be "
            "read for that domain pair specifically. 3TSY stays in the record "
            "as a partial experimental cross-check: it is a 4CL::stilbene "
            "synthase fusion missing 4CL's last 102 residues, so it cannot be "
            "the source, but it does give observed coordinates for ref 22-459 "
            "to align the prediction against (user's call, 2026-10-27)"),
}

RCSB_ENTITY = "https://data.rcsb.org/rest/v1/core/polymer_entity/%s/%s"
RCSB_INSTANCE = "https://data.rcsb.org/rest/v1/core/polymer_entity_instance/%s/%s"
RCSB_ENTRY = "https://data.rcsb.org/rest/v1/core/entry/%s"


def _json(url):
    try:
        body = _get(url)
    except RuntimeError:
        return None
    if not body:
        return None
    try:
        return json.loads(body)
    except ValueError:
        return None


#: Sequence identity to a candidate entry at or above which that entry holds
#: THIS protein, whatever UniProt does or does not say.
#:
#: This closes a hole that cost two of these inputs. Identity was resolved
#: only through UniProt's curated PDB cross-references, which is the right
#: authority -- but four of these enzymes have accessions UniProt does not
#: index at all, so the lookup returned nothing and the code read "the
#: authority cannot answer" as "no experimental structure exists". Both were
#: then predicted. Comparing the sequences directly answers it: **CUS is
#: 100.0 % identical to 3ALE over 402 aligned residues and GGPPS is 100.0 %
#: identical to 2E8T over 335** (S. cerevisiae GGPPS, holo with FsPP/IPP/Mg
#: at 2.13 A). Their own structures were in the candidate list the whole time.
#:
#: Only CUS ends up using its entry. GGPPS's is then vetoed on COVERAGE --
#: 2E8T has coordinates for 89 % of the aligned span, with 26 of the last 35
#: residues unobserved -- so GGPPS is predicted after all, for a stated
#: reason rather than for want of an answer. Finding the entry and then
#: rejecting it is the point: the first version did neither.
#:
#: 0.95 and not lower: CHS is 91.8 % identical to 1BI5 over 388 residues,
#: which is a close orthologue from another plant and not the same protein,
#: and its prediction stands.
#: AlphaFold DB model version this set is pinned to. One constant,
#: because three independent copies of "6" is three chances to drift
#: when AFDB publishes v7: the downloader, build_inputs.py and
#: derive_start_sites.py all name the same file.
AFDB_VERSION = "6"

SAME_PROTEIN_IDENTITY = 0.95
#: Over at least this many aligned residues, so a short local match cannot
#: carry the claim.
SAME_PROTEIN_PAIRS = 150

RCSB_ENTITY_SEQ = "https://data.rcsb.org/rest/v1/core/polymer_entity/%s/%s"


EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def protein_sequence(accession):
    """The protein sequence for a published accession, protein or nucleotide.

    Tries NCBI's protein database first; a nucleotide accession fails there,
    so the CDS translation is requested instead -- NCBI's own translation,
    not one computed here from a guessed frame.

    Lives here rather than in `make_af3_inputs.py` because `decide()` needs it
    to answer "does this candidate entry hold this protein", and that script
    imports this one.
    """
    for db, rettype in (("protein", "fasta"), ("nuccore", "fasta_cds_aa")):
        url = "%s?%s" % (EFETCH, urllib.parse.urlencode(
            {"db": db, "id": accession, "rettype": rettype, "retmode": "text"}))
        try:
            body = _get(url)
        except RuntimeError:
            continue
        if not body:
            continue
        text = body.decode("utf-8", "replace")
        if not text.lstrip().startswith(">"):
            continue
        blocks = [b for b in text.split(">") if b.strip()]
        if len(blocks) != 1:
            # Which of several is the enzyme is a decision, not a lookup.
            return None, ("%s gave %d sequences" % (accession, len(blocks)))
        lines = blocks[0].split("\n")
        seq = "".join(lines[1:]).strip().replace("*", "")
        bad = sorted(set(seq) - set("ACDEFGHIKLMNPQRSTVWY"))
        if bad:
            return None, ("%s has non-standard residue(s) %s"
                          % (accession, ",".join(bad)))
        return seq, lines[0]
    return None, "no protein sequence for %s" % accession


def entry_sequence(pdb_id, chain=None):
    """The one-letter sequence of an entry's polymer entities, as a list."""
    entry = _json(RCSB_ENTRY % pdb_id.upper())
    if not entry:
        return []
    out = []
    for eid in ((entry.get("rcsb_entry_container_identifiers") or {}).get(
            "polymer_entity_ids") or []):
        ent = _json(RCSB_ENTITY_SEQ % (pdb_id.upper(), eid))
        if not ent:
            continue
        seq = ((ent.get("entity_poly") or {}).get(
            "pdbx_seq_one_letter_code_can") or "").replace("\n", "").strip()
        if seq:
            out.append(seq)
    return out


def holds_same_protein(pdb_id, own_sequence):
    """``(is_same, note)`` by comparing sequences, not by asking UniProt.

    Used only as a FALLBACK: UniProt's curated mapping is the better
    authority when it has an answer. This is for when it has none.
    """
    if not own_sequence:
        return False, ""
    try:
        from seqalign import align
    except ImportError:
        return False, ""
    best = (0.0, 0)
    for seq in entry_sequence(pdb_id):
        pairs, ident = align(seq, own_sequence)
        if len(pairs) >= SAME_PROTEIN_PAIRS and ident > best[0]:
            best = (ident, len(pairs))
    ident, n = best
    if ident >= SAME_PROTEIN_IDENTITY:
        return True, ("%s is %.1f%% identical to this enzyme over %d aligned "
                      "residues, so it holds this protein (UniProt indexes no "
                      "entry for the published accession, so its curated "
                      "cross-references could not answer)" % (pdb_id, 100 * ident, n))
    return False, ("%s is %.1f%% identical over %d residue(s), below the "
                   "%.0f%% needed to call it the same protein"
                   % (pdb_id, 100 * ident, n, 100 * SAME_PROTEIN_IDENTITY))


def entry_coverage(pdb_id, uniprot):
    """How much of *uniprot* the entry *pdb_id* actually has coordinates for.

    Returns ``{chain, span, observed, coverage, longest_gap, gaps}`` or None
    when RCSB does not describe the entry well enough to tell. A residue
    present in the construct but unobserved in the density is absent as far
    as any distance measurement is concerned, so it counts as missing here.
    """
    entry = _json(RCSB_ENTRY % pdb_id.upper())
    if not entry:
        return None
    eids = (entry.get("rcsb_entry_container_identifiers") or {}).get(
        "polymer_entity_ids") or []
    for eid in eids:
        ent = _json(RCSB_ENTITY % (pdb_id.upper(), eid))
        if not ent:
            continue
        regions = []
        for a in ent.get("rcsb_polymer_entity_align") or []:
            # uniprot=None means "whatever this entity maps to": used when
            # UniProt indexes no entry for the enzyme and identity was
            # established by sequence instead, so there is no accession to
            # match against.
            if uniprot is not None and str(
                    a.get("reference_database_accession", "")).upper() != uniprot.upper():
                continue
            regions.extend(a.get("aligned_regions") or [])
        if not regions:
            continue
        chains = (ent.get("rcsb_polymer_entity_container_identifiers") or {}).get(
            "auth_asym_ids") or []
        if not chains:
            continue
        chain = chains[0]
        # entity position -> is it this protein's, and which ref residue
        ref_of = {}
        for r in regions:
            for k in range(r["length"]):
                ref_of[r["entity_beg_seq_id"] + k] = r["ref_beg_seq_id"] + k
        inst = _json(RCSB_INSTANCE % (pdb_id.upper(), chain))
        unobserved = set()
        for f in ((inst or {}).get("rcsb_polymer_instance_feature") or []):
            if f.get("type") != "UNOBSERVED_RESIDUE_XYZ":
                continue
            for q in f.get("feature_positions") or []:
                beg = q["beg_seq_id"]
                for k in range(beg, q.get("end_seq_id", beg) + 1):
                    unobserved.add(k)
        span = sorted(ref_of.values())
        missing = sorted(ref_of[e] for e in ref_of if e in unobserved)
        gaps, run = [], []
        for n in missing:
            if run and n == run[-1] + 1:
                run.append(n)
            else:
                if run:
                    gaps.append((run[0], run[-1]))
                run = [n]
        if run:
            gaps.append((run[0], run[-1]))
        observed = len(span) - len(missing)
        return {"chain": chain, "entity": eid,
                "span": (span[0], span[-1]) if span else (0, 0),
                "observed": observed, "missing": len(missing),
                "coverage": observed / len(span) if span else 0.0,
                "longest_gap": max((b - a + 1 for a, b in gaps), default=0),
                "gaps": gaps}
    return None


#: Buffer and cryoprotectant components seen bound in these entries that are
#: not anybody's cargo. NOT_CARGO covers the ions; these are the organics.
NOT_CARGO_ORGANIC = {"BTB", "EPE", "MES", "TRS", "GOL", "EDO", "PEG", "PG4",
                     "PGE", "1PE", "2PE", "P6G", "MPD", "DMS", "CIT", "FMT",
                     "ACY", "3CX", "3FX", "B3P", "PIN", "QQ4"}


def entry_facts(pdb_id):
    """``{resolution, ligands, mutation, chains}`` for ranking entries."""
    entry = _json(RCSB_ENTRY % pdb_id.upper())
    if not entry:
        return None
    info = entry.get("rcsb_entry_info") or {}
    res = (info.get("resolution_combined") or [None])[0]
    bound = [c.upper() for c in (info.get("nonpolymer_bound_components") or [])]
    cargo = [c for c in bound
             if c not in NOT_CARGO and c not in NOT_CARGO_ORGANIC]
    mutation = ""
    for eid in ((entry.get("rcsb_entry_container_identifiers") or {}).get(
            "polymer_entity_ids") or []):
        ent = _json(RCSB_ENTITY % (pdb_id.upper(), eid))
        m = ((ent or {}).get("rcsb_polymer_entity") or {}).get("pdbx_mutation")
        if m:
            mutation = str(m)
            break
    return {"resolution": res, "ligands": cargo, "bound": bound,
            "mutation": mutation,
            "chains": info.get("deposited_polymer_entity_instance_count")}


def _pick_entry(candidates, uniprot):
    """Which of several entries of the same protein to take, and why.

    Coverage is a veto, not a ranking; among the entries that pass it, three
    things decide, in this order:

    1. **A bound cargo**, because the channel's start point comes from where
       the channelled species sits. An apo entry forces `--start-cavity`.
    2. **Wild type**, because a point mutation in the active site is a
       different channel -- which is the whole reason anyone made the mutant.
    3. **Resolution.**

    Without this the first cross-reference wins, which for BAS would take
    3A5Q (apo, 1.80 A) over 3A5R (4-coumaric acid bound, wild type, 1.60 A),
    and could take 3A5S -- the I207L/L208F mutant, and the entry this set
    had already mistakenly assigned to STS.
    """
    scored = []
    for pid in candidates:
        ok, note = _usable_entry(pid, uniprot)
        if not ok:
            continue
        f = entry_facts(pid) or {}
        scored.append((
            0 if f.get("ligands") else 1,
            0 if not f.get("mutation") else 1,
            f.get("resolution") if f.get("resolution") is not None else 99.0,
            pid, f, note))
    if not scored:
        return None, None, [_usable_entry(p, uniprot)[1] for p in candidates]
    scored.sort()
    _, _, res, pid, f, note = scored[0]
    why = "%s: %s, %s, %.2f A" % (
        note,
        ("bound " + "+".join(f["ligands"][:3])) if f.get("ligands") else "apo",
        f["mutation"] if f.get("mutation") else "wild type",
        res if res != 99.0 else 0.0)
    if len(scored) > 1:
        why += " (preferred over %s)" % ", ".join(
            "%s (%s%s%.2f A)" % (
                x[3], ("bound " + "+".join(x[4]["ligands"][:2]) + ", ")
                if x[4].get("ligands") else "apo, ",
                (x[4]["mutation"] + ", ") if x[4].get("mutation") else "",
                x[2]) for x in scored[1:4])
    return pid, why, []


def _usable_entry(pdb_id, uniprot):
    """``(ok, note)`` -- whether an entry has enough of the protein to measure.

    A failure here does not reject the entry from the record; it stops it
    being the SOURCE of coordinates, and says why in the same sentence.
    """
    cov = entry_coverage(pdb_id, uniprot)
    if cov is None:
        return True, ""
    gap = cov["longest_gap"]
    if cov["coverage"] >= MIN_COVERAGE and gap <= MAX_GAP:
        return True, "%s covers ref %d-%d, %.0f%% observed" % (
            pdb_id, cov["span"][0], cov["span"][1], 100 * cov["coverage"])
    return False, (
        "%s holds this protein but only %.0f%% of its ref %d-%d span has "
        "coordinates (%d residue(s) unobserved, longest run %d: %s)"
        % (pdb_id, 100 * cov["coverage"], cov["span"][0], cov["span"][1],
           cov["missing"], gap,
           ", ".join("%d-%d" % g for g in cov["gaps"][:4])))


def _usable_transplant(fill):
    """The best AlphaFill transplant that could give a start point, or None.

    AlphaFill outranks a plain model only because a ligand-derived start point
    survives. A transplant that clashes into the model's own side chains gives
    no such point, so it cannot carry that preference: 4CL's model has ATP,
    AMP and CoA at 2.6, 4.1 and 4.4 clashes per atom (and CoA's local RMSD is
    12 A), all from the tobacco homolog. Falling through to a prediction there
    is the honest answer, not a worse one.

    Only the BEST-identity tier of donors counts. Accepting any usable
    transplant at any identity would pick by clash count, which is the trap
    `pick_ligands.py` exists to avoid: 4CL's first usable transplant at any
    identity is SLU from a 0.44-identity donor, and taking it would say the
    model has a cofactor when what it has is a compound from a distant
    relative. Over the eleven candidates here the tier rule changes exactly
    one answer -- 4CL -- and every other enzyme keeps a usable transplant in
    its top tier.

    This is NOT the same question as "which ligand is the cargo", which is
    curated in `pick_ligands.py`. GES, MYS and OCS have usable top-tier
    transplants (so the model is worth using) while their chosen cargo, the
    C10 analogue GST, clashes and needs `--start-cavity` anyway.
    """
    ligs = [h for h in (fill or {}).get("ligands", []) if h.get("atoms")]
    if not ligs:
        return None
    best_identity = max(h["identity"] for h in ligs)
    for h in ligs:
        if h["identity"] >= best_identity - 1e-9 and h["per_atom"] <= CLASH_PER_ATOM_MAX:
            return h
    return None


def decide(enz, info, af, fill, own_seq=None):
    """``(source, why, use_entry)`` -- which source wins, and why.

    ``use_entry`` is the PDB entry to actually prepare, which is not always
    the one this set recorded: for a switch it is the entry chosen by
    `_pick_entry`, and it is empty for a model or a prediction. It is returned
    as a value rather than left inside *why* so that a later stage can read it
    without parsing English out of a sentence.
    """
    if enz["name"] in OVERRIDE:
        source, why = OVERRIDE[enz["name"]]
        return source, "OVERRIDE: " + why, ""
    if info is None:
        # UniProt has no entry for this accession, so its curated PDB
        # cross-references cannot answer. That is not the same as there being
        # no structure: compare the sequences directly before predicting.
        same, note = holds_same_protein(enz["entry"], own_seq)
        if same:
            ok, cov = _usable_entry(enz["entry"], None)
            if ok:
                return "pdb", "%s; %s" % (note, cov), str(enz["entry"])
            return "af3", ("%s -- but %s. Predict it" % (note, cov)), ""
        return "af3", ("no UniProt entry found for %s%s"
                       % (enz["accession"], "; " + note if note else "")), ""
    # 1. An experimental structure of THIS protein, decided by UniProt's
    #    curated PDB cross-references.
    #
    #    NOT by the search score in enzymes.tsv, which does not mean what the
    #    first version of this set took it to mean. Six enzymes there score
    #    1.0000 against an entry UniProt says does not contain them at all:
    #    MYS scores 1.0000 on 5UV1 while Q93X23 has no PDB entry, and STS
    #    scores 1.0000 on 3A5S while B2MWQ8 has none. The score ranks how
    #    well a search HIT matched, which is not the same question.
    if info["pdb"]:
        same = [p for p in info["pdb"] if p.upper() == str(enz["entry"]).upper()]
        if same:
            ok, note = _usable_entry(enz["entry"], info["uniprot"])
            if ok:
                return "pdb", ("%s contains this protein (UniProt %s "
                               "cross-references it); %s"
                               % (enz["entry"], info["uniprot"], note)), \
                    str(enz["entry"])
            unusable = [note]
        else:
            # Which of the cross-referenced entries can actually be measured,
            # and which of those is the right one to measure? Holding the
            # right protein is necessary, not sufficient, in both senses.
            pid, why, unusable = _pick_entry(info["pdb"], info["uniprot"])
            if pid:
                return "pdb-other", (
                    "UniProt %s cross-references %s, but this set uses %s, "
                    "which it does not; SWITCH TO %s -- %s"
                    % (info["uniprot"], "+".join(info["pdb"][:4]),
                       enz["entry"], pid, why)), pid
        # Every experimental structure of this protein is too incomplete to
        # measure a channel on. Fall through to a model, and say so: the entry
        # stays in the record as a partial cross-check, not as the source.
        if unusable:
            reason_prefix = "no usable experimental structure -- " + \
                "; ".join(x for x in unusable if x)
        else:
            reason_prefix = "no usable experimental structure"
        usable = _usable_transplant(fill)
        if usable:
            return "alphafill", ("%s. AlphaFill model with %s transplanted from "
                                 "%s (%.0f%% identity, %.1f clashes per atom)"
                                 % (reason_prefix, usable["ligand"],
                                    usable["donor"], 100 * usable["identity"],
                                    usable["per_atom"])), ""
        if af:
            return "afdb", "%s. AlphaFold DB model, pLDDT %.1f" % (
                reason_prefix, af["plddt"]), ""
        return "af3", "%s. Predict it" % reason_prefix, ""
    usable = _usable_transplant(fill)
    if usable:
        return "alphafill", ("model with %s transplanted from %s (%.0f%% identity, "
                             "%.2f A RMSD, %d clash%s = %.1f per atom)"
                             % (usable["ligand"], usable["donor"],
                                100 * usable["identity"], usable["rmsd"],
                                usable["clashes"],
                                "" if usable["clashes"] == 1 else "es",
                                usable["per_atom"])), ""
    if af:
        return "afdb", "model %s v%s, pLDDT %.1f -- ligand-free, so the start " \
                       "point must come from catalytic residues" % (
                           af["model"], af["version"], af["plddt"]), ""
    return "af3", "nothing predicted for %s" % info["uniprot"], ""


def read_decisions(path):
    """The committed structures.tsv, as a list of dicts."""
    rows = []
    with open(path, encoding="utf-8") as fh:
        header = None
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if header is None:
                header = f
                continue
            f += [""] * (len(header) - len(f))
            rows.append(dict(zip(header, f)))
    return rows


def fetch(rows, *, force=False):
    """Download the coordinates each decision calls for, into `_entries/`.

    Nothing is prepared or measured here: this only puts the right file on
    disk under a name that says where it came from, so that a later prep run
    cannot silently use the wrong one.
    """
    os.makedirs(ENTRIES, exist_ok=True)
    todo, done, failed = [], 0, []
    for r in rows:
        name, source, uni = r["enzyme"], r["structure_source"], r["uniprot"]
        if source == "alphafill" and uni:
            url, out = ALPHAFILL_CIF % uni, "alphafill_%s.cif" % uni
        elif source == "afdb" and uni:
            ver = AFDB_VERSION
            url, out = AFDB_CIF % (uni, ver), "afdb_%s_v%s.cif" % (uni, ver)
        elif source in ("pdb", "pdb-other"):
            # `use_entry`, not `current_entry`: for a switch those differ, and
            # fetching the one this set recorded would fetch the wrong protein
            # -- the error the whole file exists to correct.
            pid = (r.get("use_entry") or "").strip()
            if not pid:
                todo.append("%s (no use_entry recorded; re-run --tsv)" % name)
                continue
            url, out = RCSB_CIF % pid, "%s.cif" % pid
        elif source == "af3":
            todo.append(name)
            continue
        else:
            todo.append("%s (unknown source %r)" % (name, source))
            continue
        path = os.path.join(ENTRIES, out)
        if os.path.exists(path) and not force:
            print("%-7s have %s" % (name, out))
            done += 1
            continue
        try:
            body = _get(url, timeout=120)
        except RuntimeError as exc:
            failed.append((name, str(exc)[:90]))
            continue
        if not body:
            failed.append((name, "404 %s" % url))
            continue
        with open(path, "wb") as fh:
            fh.write(body)
        print("%-7s %s  (%.1f MB)" % (name, out, len(body) / 1e6))
        done += 1
    print("\nfetched or already present: %d" % done)
    if todo:
        print("not fetched: %s" % ", ".join(todo))
    for name, why in failed:
        print("FAILED %-7s %s" % (name, why))
    return 1 if failed else 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fetch", action="store_true",
                    help="download the coordinates structures.tsv calls for, "
                         "into _entries/, instead of querying anything")
    ap.add_argument("--force", action="store_true",
                    help="re-download files that are already there")
    ap.add_argument("--tsv", action="store_true",
                    help="emit structures.tsv instead of a report")
    ap.add_argument("--only", action="append", help="one enzyme name; repeatable")
    args = ap.parse_args(argv)

    if args.fetch:
        table = os.path.join(HERE, "structures.tsv")
        if not os.path.exists(table):
            raise SystemExit("no %s yet; run without --fetch first" % table)
        return fetch(read_decisions(table), force=args.force)

    enz = enzymes()
    wanted = args.only or sorted(enz)
    rows = []
    for name in wanted:
        e = enz.get(name)
        if e is None:
            print("no such enzyme: %s" % name, file=sys.stderr)
            continue
        e = dict(e)
        e["name"] = name
        e["accession"] = ACCESSION.get(name, "")
        if not e["accession"]:
            print("%-7s no accession recorded" % name, file=sys.stderr)
            continue
        info = to_uniprot(e["accession"])
        af = afdb(info["uniprot"]) if info else None
        fill = alphafill(info["uniprot"]) if info else None
        actual = entry_contains(e["entry"])
        own_seq = None
        if info is None:
            own_seq, _hdr = protein_sequence(e["accession"])
        source, why, use_entry = decide(e, info, af, fill, own_seq)
        if actual and info and actual["uniprot"] != info["uniprot"]:
            why += "; NOTE %s actually contains %s (%s, %s)" % (
                e["entry"], actual["uniprot"], actual["name"], actual["organism"])
        rows.append((name, e, info, af, fill, source, why, actual, use_entry))
        if not args.tsv:
            print("%-7s %-13s %-11s %s"
                  % (name, e["accession"], info["uniprot"] if info else "-",
                     (info["organism"][:40] + ", %d aa" % info["length"]) if info else ""))
            print("        first version: %s chain %s (recorded search score %s)"
                  % (e["entry"], e["chain"], e["score"]))
            if actual:
                flag = "" if (info and actual["uniprot"] == info["uniprot"]) \
                    else "   <-- DIFFERENT PROTEIN"
                print("          that entry holds: %s  %s [%s]%s"
                      % (actual["uniprot"], actual["name"][:40],
                         actual["organism"], flag))
            print("        UniProt PDB xrefs: %s"
                  % ("+".join(info["pdb"]) if info and info["pdb"] else "none"))
            print("        AFDB:    %s" % (
                "v%s pLDDT %.1f" % (af["version"], af["plddt"]) if af else "none"))
            if fill and fill["ligands"]:
                print("        fill:    %s" % ", ".join(
                    "%s from %s (%.0f%%, %d clashes)"
                    % (h["ligand"], h["donor"], 100 * h["identity"], h["clashes"])
                    for h in fill["ligands"][:4]))
            else:
                print("        fill:    %s" % (
                    "no usable transplant" if fill else "no entry"))
            print("        -> %-10s %s\n" % (source.upper(), why))

    if args.tsv:
        print("\t".join(("enzyme", "accession", "uniprot", "organism", "length",
                         "afdb_plddt", "alphafill_ligands", "structure_source",
                         "current_entry", "current_entry_holds", "use_entry",
                         "reason")))
        for name, e, info, af, fill, source, why, actual, use_entry in rows:
            ligs = ",".join(h["ligand"] for h in (fill or {}).get("ligands", [])[:4])
            print("\t".join((
                name, e["accession"], info["uniprot"] if info else "",
                info["organism"] if info else "", str(info["length"]) if info else "",
                ("%.1f" % af["plddt"]) if af else "", ligs, source,
                str(e["entry"]),
                ("%s %s" % (actual["uniprot"], actual["name"])) if actual else "",
                # The entry to actually prepare, which for a switch is not
                # `current_entry`. Empty for a model or a prediction.
                use_entry,
                # The reason is the point of a decision record. Without it the
                # table says what was chosen and not why, and the next reader
                # has to re-derive it or trust it.
                " ".join(why.split()))))
        return 0

    # Nine fields, because `decide` grew a third return value and `rows`
    # grew with it. These two unpacks did not, so the default mode -- the
    # one the docstring shows -- died with `too many values to unpack`
    # before printing anything, and the counts summary had never run.
    counts = {}
    for _n, _e, _i, _a, _f, source, _w, _ac, _u in rows:
        counts[source] = counts.get(source, 0) + 1
    print("=" * 66)
    for source in ("pdb", "pdb-other", "alphafill", "afdb", "af3"):
        if source in counts:
            print("  %-10s %2d" % (source, counts[source]))
    need = [n for n, _e, _i, _a, _f, s, _w, _ac, _u in rows if s == "af3"]
    if need:
        print("\nneed an AlphaFold3 run: %s" % " ".join(need))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
