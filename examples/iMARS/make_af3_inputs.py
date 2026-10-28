#!/usr/bin/env python3
"""Build AlphaFold3 inputs for the enzymes nothing has predicted yet.

Only the enzymes `find_structures.py` marks `af3` reach here: everything with
an experimental structure of itself, an AlphaFill model or an AlphaFold DB
model is already covered, and predicting it again would be wasted GPU time
on a worse answer.

Sequences come from the accession iMARS published, fetched from NCBI, not
retyped. A nucleotide accession is translated by asking NCBI for the CDS
protein rather than translating it here.

**Oligomeric state is a real decision, not a default.** A monomer is the
cheap answer and is wrong where the active site is shared across the dimer
interface: a type III polyketide synthase builds its initiation cavity partly
from the *partner* subunit (the Met137 of CHS, Ferrer et al. 1999), so a
monomeric model has an incomplete site and therefore an incomplete tunnel.
The *trans*-prenyltransferase fold is the opposite case -- each subunit
carries a complete active site, so a monomer is enough even though the
protein is a dimer in solution. `COPIES` below records which is which and
why; it is the kind of choice that is invisible in the output once made.

MSAs are run (AF3's data pipeline is left on), because these are plant and
bacterial enzymes where a single-sequence prediction is the weakest mode.

    ./make_af3_inputs.py                 # write the JSONs
    ./make_af3_inputs.py --job           # also write a Slurm script

Outputs go under `_entries/af3/`, which is not in version control: this set
commits decisions, not data.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from find_structures import ACCESSION, read_decisions  # noqa: E402

#: Where the finished models land: beside the deposited and downloaded ones,
#: and out of version control with them.
AF3_DIR = os.path.join(HERE, "_entries", "af3")

#: Where the AF3 JOB runs. Not AF3_DIR, and not by choice: the cluster's AF3
#: wrapper refuses an OUT_DIR outside the user's project area, so a run tree
#: under `programs/` fails instantly with "is not inside your project area".
#: The run tree therefore lives here and `--collect` copies the predicted
#: models into AF3_DIR afterwards, which is where everything else looks.
#: No default: where your AF3 runs live is a property of your cluster, not
#: of this example, and a path from ours would fail on yours in a way that
#: looks like a bug in the script. `--job` refuses without it.
AF3_RUN = os.environ.get("CHANNELBRIDGE_AF3_RUN", "")

#: Scheduler settings of the emitted job script. Overridable because they are
#: site-specific: a wall-clock limit of 0 means "no limit" to Slurm and is a
#: refusal on some queues, and a node list is meaningless off one cluster.
AF3_TIME = os.environ.get("CHANNELBRIDGE_AF3_TIME", "0-0:00:00")
AF3_EXTRA = os.environ.get("CHANNELBRIDGE_AF3_SBATCH", "")
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

#: How many copies of each chain to predict, and why. The number is the
#: biological unit the ACTIVE SITE needs, not the one the protein adopts.
COPIES = {
    # Type III PKS: the initiation/elongation cavity is completed by the
    # partner subunit, so a monomer's site -- and its tunnel -- is wrong.
    "CHS": (2, "type III PKS; the active-site cavity is completed by the "
               "partner subunit, so a monomer has an incomplete site"),
    "CUS": (2, "type III PKS, as CHS"),
    # trans-prenyltransferase fold: each subunit has a complete active site,
    # so one chain is enough even though the enzyme dimerises.
    "GGPPS": (1, "trans-prenyltransferase fold; each subunit carries a "
                 "complete active site"),
    "SQS": (1, "squalene/phytoene synthase; the active site lies within one "
               "subunit"),
    # 4CL is an ANL-superfamily adenylate-forming enzyme and acts as a
    # monomer; 3TSY, the only experimental structure of it, has one chain.
    # Its two DOMAINS are what matter here, not two chains: the C-terminal
    # domain rotates over the active site, which is why this one is predicted
    # rather than taken from the AF2-era database model.
    "4CL": (1, "ANL-superfamily adenylate-forming enzyme; monomeric, and the "
               "active site is closed by its own C-terminal domain"),
}

#: Seeds. One is enough to get a structure; this project's own measurements
#: say prediction spread matters for INTER-domain geometry, which is not what
#: these single-enzyme models are used for. Raise it if a model comes back
#: with poor pLDDT at the active site.
SEEDS = [1]


def _get(url, params, *, tries=3, timeout=90):
    full = "%s?%s" % (url, urllib.parse.urlencode(params))
    last = None
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(full, timeout=timeout) as fh:
                return fh.read().decode("utf-8")
        except (urllib.error.URLError, urllib.error.HTTPError, OSError) as exc:
            last = exc
            time.sleep(2.0 * (attempt + 1))
    raise RuntimeError("%s: %s" % (full, last))


def protein_sequence(accession):
    """The protein sequence for an accession, protein or nucleotide.

    Tries the protein database first. A nucleotide accession fails there, so
    the CDS translation is requested instead -- which is NCBI's own
    translation, not one computed here from a guessed frame.
    """
    fasta = None
    try:
        body = _get(EFETCH, {"db": "protein", "id": accession,
                             "rettype": "fasta", "retmode": "text"})
        if body.lstrip().startswith(">"):
            fasta = body
    except RuntimeError:
        pass
    if fasta is None:
        body = _get(EFETCH, {"db": "nuccore", "id": accession,
                             "rettype": "fasta_cds_aa", "retmode": "text"})
        if not body.lstrip().startswith(">"):
            raise RuntimeError("no protein sequence for %s; NCBI said: %s"
                               % (accession, body.strip()[:120]))
        fasta = body
    blocks = [b for b in fasta.split(">") if b.strip()]
    if len(blocks) > 1:
        raise RuntimeError(
            "%s gave %d sequences, so which one is the enzyme is a decision "
            "rather than a lookup: %s" % (accession, len(blocks),
                                          "; ".join(b.split("\n")[0][:50] for b in blocks)))
    lines = blocks[0].split("\n")
    header, seq = lines[0], "".join(lines[1:]).strip().replace("*", "")
    bad = sorted(set(seq) - set("ACDEFGHIKLMNPQRSTVWY"))
    if bad:
        raise RuntimeError("%s has non-standard residue(s) %s, which AF3 will "
                           "reject" % (accession, ",".join(bad)))
    return header, seq


#: Global CA RMSD between AF3 samples above which the samples are not one
#: answer. Below it, taking the top-ranked model loses nothing.
MULTIMODAL_RMSD = 2.0


def sample_spread(src_dir, stem):
    """``(rmsds, paths)`` -- how far AF3's samples are from the top model.

    AF3 writes five samples and a ranking, and `collect` takes the top one.
    That is right only if the five are one answer. Where they are not,
    reducing them to a single file throws away the fact that the prediction
    found two states -- and does it silently, which is the part that matters:
    the one file left behind looks exactly like a confident answer.

    Measured on these five enzymes, three agree to <= 0.9 A and two do not.
    The check cannot tell you WHICH kind of disagreement you have, and the
    two kinds need opposite handling: a disordered terminus should be trimmed
    (GGPPS: an 18-residue N-terminal tail swinging 9 A, with a structured
    core of 326 of 358), while a domain in two conformations must be kept in
    both (4CL: a 102-residue C-terminal domain in two states 13 A apart, both
    closed, which is the ANL adenylation/thioesterification alternation). So
    this reports and saves; it does not decide.
    """
    import glob as _glob
    try:
        import numpy as np
        sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
        from channelbridge.geometry import superpose
        from channelbridge.io import read_mmcif
    except ImportError:
        return [], []

    paths = sorted(_glob.glob(os.path.join(
        src_dir, "seed-*_sample-*", "%s_*_model.cif" % stem)))
    top = os.path.join(src_dir, "%s_model.cif" % stem)
    if len(paths) < 2 or not os.path.exists(top):
        return [], []

    def ca(path):
        a = read_mmcif(path)
        return {(x.chain, x.resseq): np.asarray(x.xyz)
                for x in a.atoms if x.name == "CA" and not x.hetatm}

    ref = ca(top)
    out = []
    for p in paths:
        c = ca(p)
        keys = sorted(set(ref) & set(c))
        if len(keys) < 3:
            out.append(float("nan"))
            continue
        _, _, r = superpose(np.asarray([c[k] for k in keys]),
                            np.asarray([ref[k] for k in keys]))
        out.append(r)
    return out, paths


def collect():
    """Copy finished AF3 models from the run tree into `_entries/`.

    AF3 writes one folder per input, named after the input's `name`, holding
    `<name>_model.cif` plus ranking and confidence files. Only the top model
    and its confidence summary are copied: the rest is regenerable and large,
    and `_entries/` is read by the prep stage, not archived.

    The copied name says where it came from -- `af3_<ENZYME>.cif` -- so a prep
    run cannot confuse it with a deposited entry or an AlphaFold DB model.
    """
    import glob
    import shutil

    os.makedirs(AF3_DIR, exist_ok=True)
    out_root = os.path.join(AF3_RUN, "output")
    if not os.path.isdir(out_root):
        raise SystemExit("no AF3 output yet at %s" % out_root)
    n = 0
    # The enzymes actually predicted, from the decision table -- not the
    # keys of COPIES. CUS is in COPIES and resolved to a deposited entry, so
    # iterating COPIES hunted an output folder that was never made and made
    # `--collect` exit 1 for a run in which everything succeeded.
    wanted = [r["enzyme"] for r
              in read_decisions(os.path.join(HERE, "structures.tsv"))
              if r["structure_source"] == "af3"]
    for name in wanted:
        stem = "imars_%s" % name.replace("'", "p").lower()
        # AF3 may write <stem>/ or <stem>_<timestamp>/ when rerun; take the newest
        dirs = sorted(glob.glob(os.path.join(out_root, stem + "*")),
                      key=os.path.getmtime)
        if not dirs:
            print("%-7s no output folder yet" % name)
            continue
        src_dir = dirs[-1]
        models = sorted(glob.glob(os.path.join(src_dir, "*_model.cif")))
        if not models:
            print("%-7s %s has no model yet" % (name, os.path.basename(src_dir)))
            continue
        dest = os.path.join(AF3_DIR, "af3_%s.cif" % name.replace("'", "p"))
        shutil.copy2(models[0], dest)
        # AF3 prefixes its own files with the input name
        # (`imars_chs_summary_confidences.json`), so a bare
        # `summary_confidences.json` silently finds nothing and the collect
        # reports no confidence at all for a model that has it.
        #
        # Two different files, and the distinction matters: the SUMMARY holds
        # ptm/iptm, while the PAE matrix -- the only evidence available about
        # how confidently two domains are placed relative to each other -- is
        # in the full one. Copying the summary under the name
        # `_confidences.json` would hand a reader a file that looks like the
        # one they want and is missing the matrix.
        # AF3 prefixes its files with the INPUT name, not the folder name,
        # and the folder may be `<stem>_<timestamp>` on a rerun -- which the
        # comment above says and this line did not honour. Getting it wrong
        # copies no confidence file, so `verify_4cl.py` prints "PART 2
        # skipped" and the inter-domain check silently does not happen.
        stem_name = stem
        extra = os.path.join(src_dir, "%s_summary_confidences.json" % stem_name)
        full = os.path.join(src_dir, "%s_confidences.json" % stem_name)
        if os.path.exists(full):
            shutil.copy2(full, os.path.join(
                AF3_DIR, "af3_%s_confidences.json" % name.replace("'", "p")))
        plddt = ""
        if os.path.exists(extra):
            shutil.copy2(extra, os.path.join(
                AF3_DIR, "af3_%s_summary.json" % name.replace("'", "p")))
            try:
                with open(extra, encoding="utf-8") as fh:
                    d = json.load(fh)
                # iptm is an INTER-chain measure, so a monomer has none.
                # Printing "None" there reads as a missing value rather than
                # as a quantity that does not apply.
                iptm = d.get("iptm")
                plddt = "  ptm %.3f  iptm %s" % (
                    d.get("ptm", float("nan")),
                    "%.3f" % iptm if iptm is not None else "n/a (one chain)")
            except (ValueError, OSError):
                pass
        print("%-7s %s <- %s%s" % (name, os.path.basename(dest),
                                   os.path.basename(src_dir), plddt))
        rmsds, sample_paths = sample_spread(src_dir, stem)
        if rmsds and max(rmsds) > MULTIMODAL_RMSD:
            worst = max(range(len(rmsds)), key=lambda i: rmsds[i])
            alt = os.path.join(AF3_DIR, "af3_%s_alt.cif" % name.replace("'", "p"))
            shutil.copy2(sample_paths[worst], alt)
            print("        SAMPLES DISAGREE: %s from the top model (range "
                  "%.2f-%.2f A over %d samples)"
                  % ("%.2f A" % rmsds[worst], min(rmsds), max(rmsds), len(rmsds)))
            print("        saved the most distant one as %s. Decide which it "
                  "is before using either:" % os.path.basename(alt))
            print("        a disordered TERMINUS should be trimmed; a DOMAIN "
                  "in two conformations needs both.")
        n += 1
    print("\ncollected %d of %d" % (n, len(wanted)))
    return 0 if n == len(wanted) else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--job", action="store_true", help="also write a Slurm script")
    ap.add_argument("--collect", action="store_true",
                    help="copy finished models from the run tree into _entries/")
    ap.add_argument("--only", action="append", help="one enzyme; repeatable")
    args = ap.parse_args(argv)

    if args.collect:
        return collect()

    table = os.path.join(HERE, "structures.tsv")
    rows = [r for r in read_decisions(table) if r["structure_source"] == "af3"]
    if args.only:
        rows = [r for r in rows if r["enzyme"] in args.only]
    if not rows:
        raise SystemExit("nothing marked af3 in structures.tsv")

    in_dir = os.path.join(AF3_RUN, "input")
    os.makedirs(in_dir, exist_ok=True)
    written = []
    for r in rows:
        name = r["enzyme"]
        acc = ACCESSION.get(name) or r["accession"]
        if name not in COPIES:
            # The module docstring says oligomeric state is a real decision,
            # not a default. Defaulting to a monomer here would spend an AF3
            # run on a structure whose active site may be completed across
            # an interface that is not in it.
            raise SystemExit(
                "%s: no oligomeric-state decision recorded. Add it to COPIES "
                "with its reason before predicting -- a monomer is a claim "
                "about the active site, not a safe default." % name)
        copies, why = COPIES[name]
        header, seq = protein_sequence(acc)
        ids = [chr(ord("A") + i) for i in range(copies)]
        doc = {
            "name": "imars_%s" % name.replace("'", "p").lower(),
            "modelSeeds": SEEDS,
            "sequences": [{"protein": {"id": ids, "sequence": seq}}],
            "dialect": "alphafold3",
            "version": 1,
        }
        path = os.path.join(in_dir, "%s.json" % doc["name"])
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(doc, fh, indent=2)
            fh.write("\n")
        written.append((name, acc, len(seq), copies, path))
        print("%-7s %-14s %4d aa  x%d chain(s)  %s"
              % (name, acc, len(seq), copies, os.path.basename(path)))
        print("        %s" % header[:96])
        print("        copies: %s" % why)

    if args.job:
        if not AF3_RUN:
            raise SystemExit(
                "--job needs $CHANNELBRIDGE_AF3_RUN: the directory your "
                "AlphaFold3 runs live in. There is no default because that "
                "is a property of your cluster, not of this example.")
        out_dir = os.path.join(AF3_RUN, "output")
        script_dir = os.path.join(AF3_RUN, "script")
        log_dir = os.path.join(AF3_RUN, "log")
        for d in (out_dir, script_dir, log_dir):
            os.makedirs(d, exist_ok=True)
        job = os.path.join(script_dir, "imars_af3.sh")
        with open(job, "w", encoding="utf-8") as fh:
            fh.write("""#!/bin/bash
#SBATCH --job-name=imars_af3
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=28G
#SBATCH --time=%s
#SBATCH --output=%s/%%x_%%j.log
%s
### SITE-SPECIFIC. The three lines above and the `source` line below name
### one cluster's scheduler, queue limits and AlphaFold3 wrapper. Nothing
### here can run anywhere else unedited; the JSONs are the portable part.
###
### ChannelBridge examples/iMARS: the %d enzymes with no structure and no
### predicted model, all in ONE job (the wrapper runs every top-level *.json).
### MSAs are run (the data pipeline is left on), because a single-sequence
### prediction is the weakest mode for these.
###
### The Slurm log goes to NFS, not /scratch_ssd/$USER/: Slurm does not create
### parent directories, that path does not exist on every node, and a job
### whose log cannot be opened FAILS instantly with no log saying why.
###
### 28G and 8 CPUs are the values AF3 jobs here are known to complete on, not
### a guess: Slurm reserves what you ASK for, so an over-request idles a GPU
### that another job could have used. Measured peak on this cluster is ~18G.

OUT_DIR="%s"
SCRIPT_DIR="%s"
INPUT_DIR_NFS="%s"
DELETE_DATA_JSON=false
AF3_ARGS=()

source /home/programs/scripts/AlphaFold3/AF3_wrapper.sh
""" % (AF3_TIME, log_dir, AF3_EXTRA, len(written), out_dir,
       script_dir, in_dir))
        os.chmod(job, 0o755)
        print("\nwrote %s" % job)
        print("submit with:  sbatch %s" % job)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
