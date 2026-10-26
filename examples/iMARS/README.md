# The iMARS example set

Nineteen enzyme pairs over twenty-four enzymes, run through stages 0–4 on
both tunnel branches. This is the set ChannelBridge is tested against when
the question is "does this work on structures nobody prepared for it", as
opposed to [the worked example](../gdh_fnr/pair.yaml), which is curated.

```bash
export CAVER=/path/to/caver            # needed for 26 of the 38 sides
export PYMOL=pymol                     # for stages 3 and 4
./run_examples.py all
./run_examples.py report
./run_examples.py all --only iMARS_01_4CL-STS      # just one
```

CAVER is needed for more than the three unannotated enzymes: most annotated
sides end up on a seeded CAVER run too, for the reason in
[What the run says](#what-the-run-says). Network access is needed for
`fetch` and for the ChannelsDB reports.

Nothing is bundled. `fetch` downloads twenty deposited entries (~15 MB);
everything under `_entries/` and `iMARS_*/{input,output}/` is generated and
not in version control. What *is* committed is `enzymes.tsv`, `pairs.tsv`,
the nineteen `pair.yaml` files and this script — every decision, and none of
the data.

## Where the set comes from

iMARS (Wang et al., *Cell* 188:1349, 2025) reports 36 enzyme fusion
constructs. A construct is a linker variant — fifteen of the 36 are 4CL–STS
with different linkers — and ChannelBridge does not fuse anything, so the
linker is not a variable here. Collapsed on the enzyme pair, the 36 become
the 19 in `pairs.tsv`, numbered by first appearance in the paper's table.

**None of these structures is an iMARS construct.** The paper's enzymes are
plant and yeast proteins and almost none has its own crystal structure, so
each is represented by a homologous deposited entry, found by RCSB sequence
search on the parent sequence. `enzymes.tsv` carries the score, and it
matters: six entries are distant homologs (ADS 0.41, LIC 0.46, RAS 0.52,
REPI 0.58, BAS 0.67, 4CL 0.70), where the tunnel is evidence about the
family rather than about that enzyme.

## The two branches, and how the split fell out

Stage 0 asks whether ChannelsDB already has tunnels for an entry. The
answer is per **entry**, not per protein, and this set shows how sharply:

| | |
|---|---|
| rank-1 search hits that are annotated | **4 of 24** |
| enzymes with an annotated entry somewhere in their top 12 | **22 of 24** |
| candidate entries checked to establish that | **198** |

So searching down the ranks is not an optimisation, it is the difference
between the branch being available and not. Two enzymes (**STS**, **CUS**)
have nothing annotated in twelve candidates each and go to CAVER. A third
(**7OMT**) has exactly one annotated candidate, 5ICE — which 4'OMT already
uses — so it goes to CAVER too rather than put two "different" enzymes on
one set of coordinates. That makes pairs 01, 05 and 18 **mixed**: one side
starts on each branch, inside a single pair file.

That is where each side *starts*, not where its tunnel comes from. As the
next section shows, most annotated sides end up on a CAVER run anyway —
seeded from the site the database identified — because the fetched chain
does not reach the surface. `pairs.tsv` records the starting branch and
`report.tsv` the one that actually produced each tunnel.

## What the run says

**Most ChannelsDB chains never reach the surface, so most of them are not
measured from.** Only 10 of the 31 database-branch enzyme sides have a
candidate chain whose ends fall below ChannelBridge's bulk-solvent threshold
(burial < 15 protein atoms within 8 Å). For 5ICE — a single-chain entry, so
nothing is lost to prepping — every Procognate chain has a *minimum* burial
of 29 to 71 atoms along its whole length. The chains stop where the
database's own surface definition put them, which is not where
ChannelBridge's is.

That matters because D_IC and O_IC are both measured **at the mouth**. A
chain that stops while still buried has no exit, and its last sphere is not
one; taking it anyway produces a DO about the database's truncation point
rather than about the protein.

### So what is done instead

The two sources are good at different things, and the fix is to use each for
what it is good at.

ChannelsDB knows **what a tunnel is for** — its categories say whether a
route leads to the cognate ligand, a cofactor or the catalytic site. That is
the hard part of choosing a tunnel, and CAVER cannot tell you. What it does
not give is a chain that reaches solvent. CAVER is the converse: give it a
large enough shell radius and the chain leaves the protein, but nothing in
its output says which cluster is the substrate route.

So when a fetched chain does not reach the surface, **CAVER is re-run seeded
from where the database put the site**, and the result is checked back
against the database chain. Seed order is the database chain's buried end
first, then the bound ligand if that fails — and the order is not arbitrary:
seeding from the ligand *first* took the number of enzymes whose CAVER route
matched the annotation from 6 down to 2, because CAVER then finds whatever
leaves the ligand site, annotation or not. The ligand is a rescue, not a
preference.

How the 38 sides resolve:

| how the tunnel was obtained | sides |
|---|---|
| ChannelsDB chain used directly (it reached the surface) | 10 |
| ChannelsDB-seeded CAVER | 23 |
| CAVER branch proper (nothing annotated) | 3 |
| failed — site enclosed, no ligand to fall back to | 2 |

`report.tsv` carries the per-side detail: how each tunnel was obtained
(`resolved`), whether a seeded run recovered the annotated route
(`recovery`), and the normalised geometry with its verdict.

**All 36 successful sides now give a `mouth` verdict**, against 15 if the
fetched chains are measured from as they come. The two failures are GES and
NES, both on entry 3N0F: CAVER finds nothing above the priority floor from
the annotated site at probe 1.4, 1.1 *or* 0.9, and neither has a ligand to
retry from. Their sites are enclosed in that protomer, and the honest
outcome is to report that rather than produce a number. Stage 3 still runs
for them — the residue classification is true without a tunnel, it just
cannot say which surface sits near a channel — and stage 4 is skipped.

### Did the seeded runs find the annotated route?

Not always, and this is the part to read before trusting the mechanism. The
23 seeded sides are 13 distinct enzymes, several of which appear in more than
one pair. Of those 13:

* **6 recovered it** — 4CL, 4'OMT, MYS, OCS, REPI, RNMT, with 80–100 % of
  the database chain lying within 2 Å of ours.
* **7 did not** — ADS, AS, BBS, ERG20, RAS, SQS, TS. For SQS and TS that is
  expected, since they fell through to the ligand seed. For the others the
  annotation genuinely did not carry over, and the tunnel used is CAVER's
  own choice with no category attached to it.

Judge this by `Match.recovery()`, not `Match.verdict()`. Our chain is longer
by design, so the exits cannot coincide and `verdict`'s coverage term counts
our extra spheres as misses — read the wrong way round it called six
recovered routes "different tunnel", one of them an exact match over 100 %
of the database chain. `recovery()` asks how much of *their* chain lies on
*ours*, which is the question actually being asked.

Three smaller findings, all now guarded in code:

* **A fetched tunnel may belong to a different protomer.** A report
  annotates the whole deposited assembly and does not say which chain each
  tunnel is in. Prepping chain A of 3KRF and taking its best Procognate
  tunnel gave a chain sitting **32 Å from the protein at both ends** —
  correct coordinates, in a frame whose protein had been deleted. It passed
  every other check, because a detached chain trivially "leaves the
  protein". Four entries here annotate a protomer other than A (ADS and BBS
  chain B, GPPS chain D, LIC chain B); `enzymes.tsv` records the chain and
  the driver re-checks attachment at run time.

* ChannelsDB's `Distance` field is cumulative arc length in MOLE entries and
  something else in CAVER-derived ones, where it is non-monotonic and
  saturates near 1. Reading it as a length reported every CAVER-derived
  tunnel as ≈1 Å — a plausible number nothing contradicted. `DbTunnel.length`
  now measures coordinates; `DbTunnel.distance_is_arc_length` says whether
  the reported series can be used.
* A catalytic triad's centre of gravity lands **inside an atom**, because
  the residues face each other. CAVER's `starting_point_residue` then finds
  zero tunnels and does not say why: STS and CUS gave nothing at probe 0.9
  and 1.4 and shell radius 3, 5 and 10. `caver config --start-cavity` moves
  the point into the pocket (3.7 Å and 4.0 Å here) and both give four
  clusters.

## Choices made here that are not the package's defaults

* **Start points for the CAVER branch** come from catalytic residues
  transferred from a holo homolog by sequence alignment, with the residue
  identity checked at the target: STS Cys157/His296/Asn329 and CUS
  Cys174/His316/Asn349 from CHS 1BI5 (70 % and 52 % identity — CUS's triad
  independently matches the published *Oryza* numbering), 7OMT's
  substrate-lining set from 5ICE (88 %).
* **The ChannelsDB pick rule** is: at least 5 Å long, attached to the chain
  that was prepared, then a chain that leaves the protein, then the widest
  bottleneck. CAVER's own ranking is by *priority* (throughput) and no
  ChannelsDB report carries it, so a mixed pair can have one tunnel chosen
  by throughput and one by bottleneck. That is a difference between the
  branches, not between the enzymes.
* **The downstream enzyme is renamed to chain B** at prep. Most of these
  entries have only chain A, and `chbridge validate` rejects a pair whose
  targets share a chain. Renaming is not moving, so the ChannelsDB branch is
  unaffected — a database sphere chain carries coordinates and no chain id.
* **`site` and `vent`** are the ligand where there is one and the catalytic
  residues where there is not. Neither is narrowed to a reactive atom,
  because for most of these enzymes nobody here knows which atom that is,
  and a guess written as a selection would read as knowledge.
* **The pose grid is left at the package default** for all 19 pairs. A
  per-pair grid would be a claim about these enzymes that this set does not
  support.

## What this set does not do

Stage 5 (the pose search) is not built, so **there are no poses here**. The
stage-4 hotspot maps place the two enzymes at the closest, best-oriented
corner of the default grid; the D_IC, O_IC and DO written into each session
are that grid corner and must not be quoted as a search result. The maps are
for deciding which epitopes to aim at, which does not need the final pose.

Seventeen of the 19 maps are built; pairs 09 and 16 are skipped because GES
and NES have no tunnel. A map missing one side's channel would show a
classification with no CHANNEL class and read as if that enzyme had none, so
it is not drawn.

Nothing in this set is a prediction about whether any of these 19 pairs
would channel, or about the yields iMARS measured. Six enzymes are
represented by distant homologs, seven of the seeded CAVER runs did not
recover the annotated route, and no pose search has been run — the set
demonstrates that the tunnel stages work on uncurated input and shows where
they stop being trustworthy. It is not a result about these enzymes.
