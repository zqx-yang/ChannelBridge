# ChannelBridge

**Channel-aware pose search for bispecific binder design.**

Given two cascade enzymes and the tunnel each uses for a channelled
metabolite, ChannelBridge finds rigid-body poses in which the upstream
product channel opens toward the downstream substrate channel, scores them
with the iMARS **DO** metric, and writes the inputs a de novo binder design
run needs.

> **Status: alpha.** The geometry core, the pair-file loader and the CLI
> commands listed below are implemented and tested. The prep, pose-scan and
> RFdiffusion3 stages are being ported from the prototype scripts and are not
> yet exposed through `chbridge`.

## Why it exists

A CAVER tunnel is a chain of probe spheres, and the pipeline needs three
things from it: where the tunnel leaves the protein, how far that is from the
other enzyme's exit, and which way the mouth points. None of those three is a
property of the protein alone — all depend on how finely the CAVER build
discretised the path and how far past the surface it kept going.

Run on the same structure, the same config and the same seed, two CAVER
builds gave cluster 1 of the GDH primary route as 16 spheres spanning 12.1 Å
and as 145 spheres spanning 19.1 Å. It is the same tunnel: the dense chain
passes 0.078 Å from the sparse chain's exit. But the exits sit **8.9 Å
apart**, which is 0.07 in DO — more than the spread across the top ten poses.

ChannelBridge normalises the sphere chain before measuring anything, which
brings that 8.9 Å down to 1.4 Å. The full measurement, including what the
normalisation does *not* fix, is in [docs/caver-builds.md](docs/caver-builds.md).

## Install

```bash
git clone <this repo>
cd channelbridge
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
    --prefix p6a

# Or in two steps, so the camera is yours: build the session, open it, turn it
# to the frame you want, save it -- then render THAT (or any session of your
# own). Everything written carries the --prefix you choose.
chbridge session out-complex_GDH_FNR --pse rank1 --prefix p6a
chbridge render p6a_geom001_render.pse --prefix p6a --variant both \
    --dic 17.0 --oic 165.0 --do 0.850

# Classify one enzyme's residues and rank its epitopes. Writes the two
# tables every later stage reads; hand-edit them or supply your own.
chbridge patches GDH.pdb --chain B --name GDH --ligands NAP,BGC \
    --tunnel GDH_cl1_norm.pdb --avoid "resi 145+158"

# One session showing every candidate epitope of both enzymes. Built BEFORE
# the sampling, so it computes its own reference placement rather than taking
# a pose. Writes p6a_hotspot-map.pse.
chbridge hotspots --prefix p6a --up-chains B,C,D,E --dn-chains N \
    --up-name GDH --dn-name FNR --ligands NAP,BGC,FDA,A2P \
    --up-pdb GDH.pdb --dn-pdb FNR.pdb \
    --tunnel-up GDH_cl1_norm.pdb --tunnel-dn FNR_cl2_norm.pdb
```

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
 1 prep ──► 2 CAVER ──► 3 patches ──► 4 hotspot map   (look, decide)
                             │
                             ▼
                         5 DO sampling ──► 6 session + SVG
                             │
                             └──────────► 7 generator input
```

| # | stage | state | where it lives |
|---|---|---|---|
| 1 | structure selection and prep | **to build** | `prep/` |
| 2 | CAVER: run it | **to build** | `caver/` |
| 2 | CAVER: normalise its output | **built** | `tunnel.py` |
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
                     entry points run INSIDE PyMOL
  caver/  prep/  pose/  rfd3/                  stage homes, not yet filled
tests/               including a cross-build regression against real CAVER output
docs/                caver-builds.md, figures.md, imars-verification.md,
                     pipeline.md, pitfalls.md
examples/gdh_fnr/    the worked example
.devcontainer/       a dev container; environment.yml for conda
```

The `render/run_*.py` files are the only ones that import PyMOL. Everything
else is ordinary Python and is tested with pytest; those three are checked by
their own selftests, because they must run inside PyMOL's interpreter:

```bash
pymol -cq src/channelbridge/render/run_pse.py -- --selftest
pymol -cq src/channelbridge/render/run_hotspots.py -- --selftest
```

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
