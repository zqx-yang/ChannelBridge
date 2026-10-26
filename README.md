# ChannelBridge

**Channel-aware pose search for bispecific binder design.**

Given two cascade enzymes and the tunnel each uses for a channelled
metabolite, ChannelBridge finds rigid-body poses in which the upstream
product channel opens toward the downstream substrate channel, scores them
with the iMARS **DO** metric, and writes the inputs a de novo binder design
run needs.

> 🚧 **UNDER CONSTRUCTION — alpha** 🚧
>
> Published early so the method can be read and criticised. Stages 0 to 4 run
> end to end -- structure prep, ChannelsDB lookup, CAVER, normalisation,
> tunnel selection, patch selection and the hotspot map -- and are exercised
> over 19 enzyme pairs in [examples/iMARS/](examples/iMARS/README.md). The
> **pose search** (stage 5), CaverDock ranking and the generator-input stage
> are not yet exposed through `chbridge` and still live in the prototype
> scripts. The table under [Layout](#layout-and-what-is-built) says which is
> which, per stage.
>
> One consequence for a first-time reader: the commands that consume a pose
> (`chbridge figures`, `session`, `render`) need a pose-search output folder,
> which this package cannot yet produce. Interfaces and defaults may change
> without notice, no release is tagged, and results produced with it should be
> reported with the commit they came from.

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
needed only for the prep, pose and figure stages, and is run as a subprocess
rather than imported. The SVG assembler additionally needs scipy, matplotlib
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

# Figures for chosen poses: sessions for the top ten, an editable SVG for the
# best one, and one pose's structure copied out. Ray-traces run one at a time.
chbridge figures out-complex_GDH_FNR --pse top10 --svg rank1 --pdb geom007 \
    --prefix run1

# Or in two steps, so the camera is yours: build the session, open it, turn it
# to the frame you want, save it -- then render THAT (or any session of your
# own). Everything written carries the --prefix you choose.
chbridge session out-complex_GDH_FNR --pse rank1 --prefix run1
chbridge render run1_geom001_render.pse --prefix run1 --variant both \
    --dic 17.0 --oic 165.0 --do 0.850

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
```

### Two example sets

[examples/gdh_fnr/](examples/gdh_fnr/pair.yaml) is the worked example: one
curated pair, with its two prepared structures committed so the first command
above runs on a clean clone.

[examples/iMARS/](examples/iMARS/README.md) is the opposite — nineteen enzyme
pairs from the iMARS paper, represented by homologous deposited entries that
nobody chose for this purpose, fetched rather than bundled. It exists to
exercise both tunnel branches on structures that were not prepared for them,
and that is what it does: 21 of its 24 enzymes have a ChannelsDB-annotated
entry and 3 do not, so three pairs fetch one tunnel and compute the other.

It is also where the uncomfortable result lives. **Most ChannelsDB chains
stop before the surface** as ChannelBridge defines it — only 10 of 31
database-branch sides have one that reaches it. The chains end where the
database's own surface definition put them, and since D_IC and O_IC are both
measured at the mouth, a chain that stops while buried has no exit to
measure.

So the branch is used for what it is good at. ChannelsDB says what a tunnel
is *for*, which CAVER cannot; CAVER produces a chain that reaches solvent,
which the database does not. When a fetched chain falls short, CAVER is
re-run **seeded from where the database put the site**, and the result is
checked back against the database chain. That takes the set from 15 usable
mouths to 36 of 38 — and on 6 of the 13 seeded enzymes CAVER recovers the
annotated route, which is worth knowing before relying on the rest. The
per-side detail is in that README.

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
| 0 | ChannelsDB lookup (a source of tunnels, not a step) | **built** | `channelsdb.py` |
| 1 | structure selection and prep | **built** | `prep/` |
| 2a | CAVER: configure and run it | **built** | `caver/` |
| 2b | CAVER: normalise its output | **built** | `tunnel.py` |
| 2b | choose which tunnels to keep | **built** | `caver/` |
| 2b | validate ours against ChannelsDB | **built** | `channelsdb.py` |
| 2c | CaverDock: rank tunnels by transport energetics | **to build** | — (run by hand today) |
| 3 | hotspot patch selection | **built** | `patches.py` |
| 4 | hotspot map session | **built** | `render/hotspots.py` |
| 5 | DO sampling: the metric | **built** | `geometry.py` |
| 5 | DO sampling: the pose search | **to build** | `pose/` |
| 6 | session, render, SVG | **built** | `render/`, `figures.py` |
| 7 | generator input (RFdiffusion3 etc.) | **to build** | `rfd3/` |

Stage 4 runs **before** stage 5 and therefore cannot borrow a pose; it
computes its own reference placement (`reference_pose.py`). Stage 6 is split
in two — `chbridge session` builds the `.pse`, `chbridge render` draws from
any `.pse`, including one you rotated and saved yourself.

```
src/channelbridge/
  geometry.py        D_IC, O_IC, the DO score and its coefficient sets
  tunnel.py          CAVER sphere chains: orient, truncate, resample, O_IC arm
  reference_pose.py  places the two enzymes at O_IC 180 without a pose search
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
  pose/  rfd3/                               stage homes, not yet filled
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
