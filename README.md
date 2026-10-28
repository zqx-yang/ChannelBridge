# ChannelBridge

**Channel-aware pose search for bispecific binder design.**

Given two cascade enzymes and the tunnel each uses for a channelled
metabolite, ChannelBridge finds rigid-body poses in which the upstream
product channel opens toward the downstream substrate channel, scores them
with the iMARS **DO** metric, and writes the inputs a de novo binder design
run needs.

> 🚧 **UNDER CONSTRUCTION — alpha** 🚧
>
> Published early so the method can be read and criticised. Stages 0 to 7 run
> end to end -- structure prep, ChannelsDB lookup, CAVER, normalisation,
> tunnel selection, patch selection, the hotspot map, the pose search and
> the generator input -- and are exercised over 19 enzyme pairs in
> [examples/iMARS/](examples/iMARS/README.md). **CaverDock ranking** (stage
> 2c) is not exposed through `chbridge`: it was measured on this set and
> could not return an answer (see
> [docs/caverdock.md](docs/caverdock.md)). The table under
> [Layout](#layout-and-what-is-built) says which is which, per stage.
>
> Interfaces and defaults may change without notice, no release is tagged,
> and results produced with it should be reported with the commit they came
> from.

## What it does

Two enzymes catalyse consecutive steps, and the product of the first is the
substrate of the second. If you hold them so that the channel carrying that
intermediate out of one enzyme points at the channel carrying it into the
other, the intermediate has a shorter way to go. A designed bispecific binder
can hold them that way — if you can say which relative pose to aim for.

**ChannelBridge searches for that pose.** Given two prepared enzymes and a
tunnel in each, it places the second relative to the first over a grid of
separations and orientations, scores every placement with the iMARS DO
metric, and reports the poses a single binder could actually bridge —
together with which surface patch on each enzyme it would have to grip.

Everything else in the package exists to make that search possible or to let
you judge its output:

| | |
|---|---|
| **tunnels** | get them from ChannelsDB, or configure and normalise a CAVER run; optionally rank them by transport energetics with CaverDock |
| **epitopes** | classify every residue into what a binder may and may not touch, and rank candidate patches |
| **figures** | a PyMOL session and an editable SVG per pose, and a map of every candidate epitope to judge before committing to a search |

### The score, and what it cannot see

`DO = -0.0079·D_IC + 0.0012·O_IC + 0.7865` — Wang et al., *Cell* 188:1349
(2025). It is the default because it is the published metric, **not because
it is complete**:

* **O_IC is one-sided.** It is the angle at the *upstream* exit, so it asks
  only whether the donor channel aims at the acceptor's exit. Which way the
  acceptor opens is not in the score at all: a receiving mouth pointing 90°
  away scores exactly like one facing the donor.
* **The intermediate is never modelled.** The gap between the two exits is
  empty; nothing checks that the molecule could cross it or fit either mouth.
* **It is a linear fit**, so a DO outside [0, 1] means the pose is outside
  the range it was fitted on. Nothing is clamped — that is information.

The search *samples* the acceptor's orientation (`tilt`) and records the
resulting angle, so you can score on it. Supply your own coefficients, or add
terms the published fit does not have:

```python
from channelbridge.geometry import DOCoefficients, do_score

mine = DOCoefficients(-0.0079, 0.0012, 0.7865, name="DO+acceptor",
                      extra={"o_ic_rev": 0.0012})
do_score(17.0, 165.0, mine, o_ic_rev=180.0)
```

A set with an extra term reports `is_published == False`, because the result
is then no longer the iMARS DO score and should not be called one.

### Why the tunnel handling is not an afterthought

The search measures three things off each tunnel — where it leaves the
protein, how far that is from the other exit, and which way the mouth points
— and **none of them is a property of the protein alone.** All three depend
on how finely a CAVER build discretised the path and how far past the surface
it kept going.

Run on the same structure, the same config and the same seed, two CAVER
builds gave cluster 1 of the same route as 16 spheres spanning 12.1 Å and as
145 spheres spanning 19.1 Å. It is the same tunnel — the dense chain passes
0.078 Å from the sparse chain's exit — but the exits sit **8.9 Å apart**,
which is 0.07 in DO, more than the spread across the top ten poses.
Normalising the chain first brings that to 1.4 Å. The full measurement,
including what normalisation does *not* fix, is in
[docs/caver-builds.md](docs/caver-builds.md).

## Install

```bash
git clone https://github.com/zqx-yang/ChannelBridge.git
cd ChannelBridge
pip install -e '.[dev]'
pytest
```

Python 3.9+, numpy. PyYAML is optional (pair files may be JSON). PyMOL is
needed for the prep, patch and figure stages, and is run as a **subprocess**
rather than imported; the pose search itself needs only numpy. If `pymol` is
not on `PATH`, point `$CHBRIDGE_PYMOL` at it — every stage that needs PyMOL
also takes `--pymol`. (A stage that cannot find PyMOL exits 127, which is
the failure that once voided 38 runs while printing `ok`, so check the
counts a stage reports, not just its exit code.) The SVG assembler additionally needs scipy, matplotlib
and pillow: `pip install -e '.[figures]'`.

A conda environment that covers PyMOL and a JRE as well:

```bash
conda env create -f environment.yml
conda activate channelbridge
pip install -e '.[dev,figures,yaml]'
```

Or open the repository in a dev container (VS Code, or GitHub Codespaces with
nothing installed locally) — see
[.devcontainer/README.md](.devcontainer/README.md).

**CAVER is not included.** It is GPL-3 software from Masaryk University;
download it from [caver.cz](https://caver.cz/) (free academic registration)
and point `$CAVER` at the folder holding `caver.jar` and `lib/`. ChannelBridge
runs it as a separate process and does not link against it. Java 8, 11 or 21
all work; see [docs/pitfalls.md](docs/pitfalls.md).

## Use

Everything about one enzyme pair lives in a **pair file**, so a run is
reproducible from a single artefact that can sit beside the results in
version control. See [examples/gdh_fnr/pair.yaml](examples/gdh_fnr/pair.yaml)
for a fully commented one.

```bash
# Check a pair file and print what it says.
chbridge validate examples/gdh_fnr/pair.yaml

# Stage 1: keep the atoms the later stages measure. Coordinates never move.
chbridge prep 1RWB.pdb --chains B --keep NAP,BGC --out GDH.pdb

# Stage 2a: a CAVER config, then CAVER. Where no ligand marks the site,
# --start-cavity derives the point from the catalytic residues -- their own
# centre of gravity lands inside an atom, and CAVER then finds nothing and
# does not say why.
chbridge caver config --out GDH.conf --start 44.890,12.027,-6.941 \
    --exclude HOH,WAT,BGC,NAP
chbridge caver config --out STS.conf --structure STS.pdb \
    --start-cavity A:157,A:296,A:329
chbridge caver run out_GDH --structure GDH.pdb --config GDH.conf

# Normalise a CAVER tunnel and print its geometry.
chbridge tunnel caver/output/out_GDH_primary_sr14 \
    --structure caver/in/GDH/gdh.pdb --exclude HOH,WAT,BGC,NAP

# Compare two CAVER runs. A large raw exit separation that shrinks when
# normalised is a sampling artefact, not a different tunnel.
chbridge tunnel buildA/out_X buildB/out_X \
    --structure gdh.pdb --exclude HOH,WAT,BGC,NAP --compare

# The DO score for one geometry, and how much the constant set matters.
chbridge do --dic 16 --oic 150 --all-sets

# Stage 5: the search. The grid comes from the pair file; --max-ok stops it
# early without biasing the answer, because the grid is ordered best-DO-first.
# Writes poses.tsv, patch_pairs.tsv, patch_pair_best.tsv and a complex per pose.
chbridge pose --pair-file examples/gdh_fnr/pair.yaml \
    --up-pdb GDH.pdb --dn-pdb FNR.pdb \
    --tunnel-up GDH_normalised.pdb --tunnel-dn FNR_normalised.pdb \
    --patches-up GDH_hotspot_patches.tsv --patches-dn FNR_hotspot_patches.tsv \
    --out-dir out-complex_GDH_FNR

# Figures for chosen poses: sessions for the top ten, an editable SVG for the
# best one, and one pose's structure copied out. Ray-traces run one at a time.
chbridge figures out-complex_GDH_FNR --pse top10 --svg rank1 --pdb geom007 \
    --prefix run1

# Or in two steps, so the camera is yours: build the session, open it, turn it
# to the frame you want, save it -- then render THAT (or any session of your
# own). Everything written carries the --prefix you choose.
chbridge session out-complex_GDH_FNR --pse rank1 --prefix run1
chbridge render run1_geom001_render.pse --prefix run1 --variant both \
    --dic 15.0 --oic 120.0 --do 0.812

# Classify one enzyme's residues and rank its epitopes. Writes the two
# tables every later stage reads; hand-edit them or supply your own.
chbridge patches GDH.pdb --chain B --name GDH --ligands NAP,BGC \
    --tunnel GDH_cl1_norm.pdb --avoid "resi 145+158"

# One session showing every candidate epitope of both enzymes. Built BEFORE
# the sampling, so it computes its own reference placement rather than taking
# a pose. Writes run1_hotspot-map.pse.
chbridge hotspots --prefix run1 --up-chains B,C,D,E --dn-chains N \
    --up-name GDH --dn-name FNR --ligands NAP,BGC,FDA,A2P \
    --up-pdb GDH.pdb --dn-pdb FNR.pdb \
    --tunnel-up GDH_cl1_norm.pdb --tunnel-dn FNR_cl2_norm.pdb

# Stage 7: turn a finished pose search into binder-design input. One spec per
# pose per tool, carrying the pose's geometry and both enzymes' hotspots.
chbridge generate out-complex_GDH_FNR --list-tools
chbridge generate out-complex_GDH_FNR --tool RFD4-Proteina,RFdiffusion3 \
    --prefix run1 --min-length 100 --max-length 150 --n-designs 1000

# ...or ask for the easier jobs first. `context` crops the partner to what
# sits near the epitope, which is also what keeps the job on a GPU.
chbridge generate out-complex_GDH_FNR --mode monospecific
chbridge generate out-complex_GDH_FNR --mode context --context-cut 30
```

**Three ways to ask, and they are a ladder rather than alternatives:**

| `--mode` | one job per | the partner is | the pose's geometry is |
|---|---|---|---|
| `bispecific` (default) | pose | a second target | **in the job** |
| `monospecific` | enzyme | absent | not in the job at all |
| `context` | enzyme | present, cropped, not targeted | in the job |

A campaign can walk up it: take the monospecific answer as a starting point,
re-run it in context so the binder is not designed into the partner, and keep
the bispecific job for what survives. `monospecific` is the easiest to design
and the only one that says nothing about the bridge — two binders made that
way have no reason to be compatible with one pose, and `chbridge generate`
says so when it writes them.

**Hotspots are chosen, not copied.** A patch is 7–14 residues; the generator
is handed the five (`--n-hotspots`) that a binder can actually lock onto —
the exposed F/Y/W/I/L/M, spread along the sequence so they span the epitope
rather than one turn of a helix. That is `patches.pick_hotspots`, the worked
example's rule, and it is also BindCraft's advice for choosing a site. A
patch with fewer anchors yields fewer hotspots and says so rather than
padding the list with residues the binder cannot grip. RFdiffusion3 then gets
its own `TIP` keyword for every one of them, so this package is not keeping a
second copy of RFd3's per-residue atom table in sync with theirs.

**What stage 7 pins, and why it matters more than it looks.** The target
chains are conditioned per atom (`C_CRD` for RFD4-Proteina, the input
coordinates for the others) and must stay **fixed relative to each other**:
their relative placement *is* the output of stages 2–5. RFD4-Proteina's
`C_CRDGRP` is the obvious-looking alternative and is wrong here — it
preserves a group's internal geometry and leaves its global pose free, so
with two targets it releases the inter-enzyme placement and discards D_IC and
O_IC without failing. Every generated spec carries that sentence. The same
applies to a context chain: it is there to occupy the space it occupies in
the pose, which it cannot do if it is free to move.

Tools that take no spec file get their job-script argument lines written
instead; two are refused with the reason rather than given a spec they cannot
honour:

| tool | spec | stage 7 |
|---|---|---|
| RFdiffusion3 | JSON | written |
| RFD4-Proteina | YAML + mmCIF | written |
| BindCraft | JSON | written |
| flexcraft | none (CLI args) | job-script lines |
| RFdiffusion | none (Hydra args) | job-script lines |
| protpardelle-1c | YAML | **refused** — no two-chain binder mode |
| RFDpoly | none | **refused** — multi-polymer model, not this task |

### Two example sets

[examples/gdh_fnr/](examples/gdh_fnr/pair.yaml) is the worked example: one
curated pair, with its two prepared structures committed so the first command
above runs on a clean clone.

[examples/iMARS/](examples/iMARS/README.md) is the opposite — nineteen enzyme
pairs from the iMARS paper, over twenty-four enzymes, fetched or predicted
rather than bundled. It exists to exercise the pipeline on structures nobody
prepared for it, and it is written as a **guide to input preparation**,
because that is the part that decides what every later number means and the
part this set got wrong first.

**The instructive failure is upstream of any stage.** The first version
represented each enzyme by the top hit of an RCSB sequence search and read
the search score as closeness to the enzyme. It is not — it ranks how well a
*hit* matched the query — and checking every entry against RCSB's own UniProt
mapping showed that **only 2 of the 24 entries held the enzyme they stood
for**. Nothing downstream was wrong: ChannelsDB annotated exactly the entry
it was asked about, CAVER measured exactly the coordinates it was given. That
is what makes an input-preparation error expensive. Every enzyme now has
coordinates of its own protein — 3 PDB entries, 2 switched to their own
experimental entry, 10 AlphaFill, 5 AlphaFold DB and 4 predicted with
AlphaFold3 — each with the reason recorded in `structures.tsv`, and
**38 of 38 enzyme sides produce a normalised tunnel.**

A second consequence is worth reading before relying on the ChannelsDB
branch. Annotation is per **entry**, so the branch is only available where an
enzyme's own coordinates happen to be an annotated PDB entry: here that is
2 of 24, and the other 22 have the branch withdrawn rather than left pointing
at tunnels of a protein the set no longer uses. Where it does apply, most
database chains stop while still buried — where the database's own surface
definition put them — and since D_IC and O_IC are both measured at the mouth,
such a chain has no exit to measure. So each source is used for what it is
good at: ChannelsDB says what a tunnel is *for*, which CAVER cannot, and
CAVER produces a chain that reaches solvent, which the database does not.
When a fetched chain falls short, CAVER is re-run **seeded from where the
database put the site** and the result is checked back against the database
chain.

**Against iMARS's own predictions.** `examples/iMARS/compare_to_imars.py`
puts our best pose beside the paper's own rank-1 linker fusion for each pair,
on the same DO coefficients. Over all 19 pairs our best feasible pose scores
higher for **19 of 19**, median **+0.1518 DO**. That is not a contest and not
a result about enzymes: DO falls 0.0079 per angstrom of D_IC and almost the
whole difference is D_IC — ours run 6–22 Å, theirs 10.1–55.3 Å. A terminal
fusion gets whatever the construct folds into; a pose search is limited by
the two surfaces instead of by the chain termini, and the comparison measures
how much geometric room that buys. The narrowest row is the informative one:
on 4CL–RAS, where their fusion already lands the exits 10.1 Å apart, the
margin is +0.0129 and our D_IC is the worse of the two. No binder has been
designed, and iMARS measured yields.

Their supplement also publishes the scripts that produced those numbers
(Table S7), and [docs/imars-verification.md](docs/imars-verification.md) §2
reads them line by line. Their code confirms the construction we implement —
including that O_IC is one-sided and measured on the upstream channel — and
shows what it leaves open: CAVER's cluster 1 is hard-coded with no filters,
nothing truncates the chain at the surface, and the CAVER parameter file the
script reads is not in the supplement.

The seven stages, which five are optional and what each reads and writes are
in [docs/pipeline.md](docs/pipeline.md). The figure layer — what is in a
session, what is in an SVG and why the channel overlay is traced in 2D rather
than rendered in 3D — is in [docs/figures.md](docs/figures.md).

As a library:

```python
from channelbridge import Tunnel, do_score
from channelbridge.io.pdb import read_pdb

protein = read_pdb("gdh.pdb").exclude_resnames(["HOH", "WAT", "BGC", "NAP"]).xyz
tun = Tunnel.from_caver_cluster("out_GDH_primary_sr14/data/clusters/tun_cl_001.pdb")
tun = tun.normalised(protein)          # orient, truncate at the surface, resample

print(tun.describe(protein))            # exit, span, bottleneck, burial, history
p1, p2 = tun.axis_point(arm=1.0), tun.exit
```

## The metric

For an upstream enzyme A and a downstream enzyme B:

```
P1 = a point on tunnel A, one O_IC arm inward from its exit
P2 = the exit of tunnel A          (upstream)
P3 = the exit of tunnel B          (downstream)

D_IC = |P2 - P3|                   inter-channel distance, Å
O_IC = angle P1-P2-P3              inter-channel orientation, degrees
DO   = -0.0079*D_IC + 0.0012*O_IC + 0.7865
```

O_IC is **one-sided** — measured at the upstream exit only — so which enzyme
is upstream changes the score, and a short downstream tunnel is not a defect.
At these coefficients 1 Å of D_IC is worth 6.6° of O_IC.

Two published constant sets exist and disagree by up to **0.032 DO** over the
grid scanned here, so the pair file records which was used.

Every DO is computed on a rigid docked pose. It says where the geometry
*could* be, not that a binder achieving it exists; that is what a structure
prediction filter tests.

## Layout, and what is built

Seven pipeline stages, five of them optional
([docs/pipeline.md](docs/pipeline.md)). **Built** means implemented and
tested here; **to build** means the stage currently lives only in the
author's project scripts and has a package home waiting for it.

```
  input prep                       <- a structure OF THAT PROTEIN: its own
      |                               entry, else AlphaFill, else AlphaFold DB,
      |                               else AlphaFold3 with MSAs. A homolog is
      v                               a different protein, not a fallback
      0 ChannelsDB search          <- by enzyme or PDB id. No structure
            |                         needed: this DECIDES which one you prepare
     +------+-----------------------------+
     | annotated entry found              | nothing annotated, or you want
     | -> adopt that entry                | your own parameters
     v                                    v
  1 prep that entry                    1 prep your structure
     | (coordinates UNMOVED, or the        |
     |  tunnels no longer apply)           v
     |                               2a CAVER config -> CAVER run
     v                                    |
  fetch its tunnels                       |
     |                                    |
     +------------------+-----------------+
                        v
             2b NORMALISE   (orient, truncate at the surface,
                             resample)    <- either source
                        |
          optional  2c CaverDock: rank by transport energetics
          optional      validate: ours vs ChannelsDB, if you have both
                        |
                        v
              3 patch selection
                        |
                        v
              4 hotspot map (optional: look, decide)
                        |
                        v
              5 DO sampling   (needs the patches)
                        |
             +----------+----------+
             v                     v
     6 session + SVG        7 generator input
```

Two constraints the diagram encodes, because breaking either silently gives
wrong numbers rather than an error:

- **On the ChannelsDB branch, prep must not move coordinates.** Their tunnels
  are in the deposited entry's frame. Extracting chains and deleting waters is
  fine; expanding a biological assembly is not.
- **Nothing reaches the DO sampling without patches.** The scan places a
  binder on a (patch, patch) pair and rejects a pose with none, so there is no
  route from a tunnel to a pose that skips stage 3.

Why each stage sits where it does — including why the ChannelsDB lookup is a
source of tunnels rather than a step in the chain, and why normalisation runs
on both branches — is in [docs/pipeline.md](docs/pipeline.md).

| # | stage | state | where it lives |
|---|---|---|---|
| — | input prep: find or predict a structure of *that* protein | **built** | `io/mmcif.py`, `examples/iMARS/find_structures.py` |
| 0 | ChannelsDB lookup (a source of tunnels, not a step) | **built** | `channelsdb.py` |
| 1 | structure selection and prep (PDB or mmCIF in, PDB out) | **built** | `prep/` |
| 2a | CAVER: configure and run it | **built** | `caver/` |
| 2b | CAVER: normalise its output | **built** | `tunnel.py` |
| 2b | choose which tunnels to keep | **built** | `caver/` |
| 2b | validate ours against ChannelsDB | **built** | `channelsdb.py` |
| 2c | CaverDock: rank tunnels by transport energetics | **not used** | — (measured and rejected: [docs/caverdock.md](docs/caverdock.md)) |
| 3 | hotspot patch selection | **built** | `patches.py` |
| 4 | hotspot map session | **built** | `render/hotspots.py` |
| 5 | DO sampling: the metric | **built** | `geometry.py` |
| 5 | DO sampling: the pose search | **built** | `pose.py` |
| 6 | session, render, SVG | **built** | `render/`, `figures.py` |
| 7 | generator input (RFdiffusion3, RFD4-Proteina, BindCraft, ...) | **built** | `generators/` |

Stage 4 runs **before** stage 5 and therefore cannot borrow a pose; it
computes its own reference placement (`reference_pose.py`). Stage 6 is split
in two — `chbridge session` builds the `.pse`, `chbridge render` draws from
any `.pse`, including one you rotated and saved yourself.

```
src/channelbridge/
  geometry.py        D_IC, O_IC, the DO score and its coefficient sets
  tunnel.py          CAVER sphere chains: orient, truncate, resample, O_IC arm
  reference_pose.py  places the two enzymes at O_IC 180 without a pose search
  pose.py            the pose search: the grid, the four filters, the tables
  patches.py         residue classes and hotspot patches; the scoring is pure
  caver/             CAVER configuration, and which tunnels to carry forward
  channelsdb.py      precomputed tunnels from ChannelsDB 2.0, and the
                     superpose-and-compare check against your own
  config.py          the pair file, and the render: style block inside it
  figures.py         finds poses in an output folder and plans their figures
  cli.py             the chbridge command line
  io/pdb.py          fixed-column PDB reading; CAVER's radius-in-B-factor
  render/
    scenes.py        the pose session: objects, scenes, the angle frame
    hotspots.py      the epitope map: one object per class, patch and envelope
    svg.py           the editable SVG: traced channel, angle marker, labels
    frame.py         cameras, label separation in the camera plane and in depth
    cut.py           cutting geometry (the cross-section scene was dropped)
    style.py         every colour, transparency and size, all overridable
    annotation.py    the D_IC line, O_IC arm, arc and labels, shared by
                     the pose session and the hotspot map
    run_pse.py  run_png.py  run_hotspots.py  run_patches.py
    run_validate.py  entry points run INSIDE PyMOL
  caver/  prep/                               stages 1 and 2a
  generators/                                stage 7: tool input specs
tests/               including a cross-build regression against real CAVER output
docs/                caver-builds.md, figures.md, imars-verification.md,
                     pipeline.md, pitfalls.md
examples/gdh_fnr/    the worked example, curated and committed
examples/iMARS/      19 enzyme pairs from the iMARS paper, both tunnel
                     branches, fetched rather than bundled
.devcontainer/       a dev container; environment.yml for conda
```

PyMOL is confined to one layer. The `render/run_*.py` entry points and the
three modules they drive (`render/scenes.py`, `render/hotspots.py`,
`render/annotation.py`) are the only code that touches `pymol`; everything
else — the geometry, the metric, the patch selection, the frames, the SVG
assembler — is ordinary Python and is covered by pytest. The PyMOL layer
cannot be reached by pytest because it must run inside PyMOL's own
interpreter, so each entry point carries its own selftest:

```bash
pymol -cq src/channelbridge/render/run_pse.py      -- --selftest
pymol -cq src/channelbridge/render/run_hotspots.py -- --selftest
pymol -cq src/channelbridge/render/run_patches.py  -- --selftest
```

Judge those by the `PASSED`/`FAILED` line they print, not by the exit code:
PyMOL is prone to segfaulting during interpreter teardown *after* a script has
finished. Run them with `PYTHONUNBUFFERED=1` or the marker can be lost in the
buffer when it does.

## Citing

ChannelBridge wraps CAVER and uses the iMARS DO metric. Neither endorses it;
both must be cited by anyone publishing results.

**CAVER** — Chovancova E, Pavelka A, Benes P, Strnad O, Brezovsky J,
Kozlikova B, Gora A, Sustr V, Klvana M, Medek P, Biedermannova L, Sochor J,
Damborsky J (2012). CAVER 3.0: A Tool for the Analysis of Transport Pathways
in Dynamic Protein Structures. *PLoS Comput Biol* 8: e1002708.
doi:[10.1371/journal.pcbi.1002708](https://doi.org/10.1371/journal.pcbi.1002708)

**The DO score** — from the iMARS work; cite the iMARS paper for the metric
and its coefficients.

For ChannelBridge itself see [CITATION.cff](CITATION.cff).

## Known limits

- The pose grid is **rigid-body**. A gating side chain that must move to admit
  the metabolite is not modelled.
- CAVER on one crystal structure is a snapshot; radii from a single conformer
  under-report accessibility. The remedy is MD frames in a `-pdb` directory.
- `binder_gap` is a minimum protein–protein distance, **not** "room for a
  binder", and the occupancy capsule is an envelope, not a design.
- Normalisation reconciles D_IC across CAVER builds. It does **not** make
  O_IC build-exact: a curved tunnel mouth retains an uncertainty of order 15°,
  reported as `oic_arm_spread_deg`.

## License

MIT — see [LICENSE](LICENSE).
