# The pipeline

Seven stages. **Five of the seven are optional**, and each optional stage can
be replaced by the user's own answer rather than skipped — the point is that
nothing downstream has to care which route the input took.

| # | stage | optional | you can instead |
|---|---|---|---|
| 1 | structure selection and preparation | yes | hand in prepared PDBs |
| 2 | CAVER: tunnels, channel residues, shadow | no | — |
| 3 | hotspot patch selection | yes | name your own patches, or edit the suggested ones |
| 4 | hotspot map session | yes | — (a figure, nothing depends on it) |
| 5 | DO sampling | no | — |
| 6 | complex session and SVG | yes | — (figures) |
| 7 | generator input prep (RFdiffusion3 or other) | yes | — |

```
 1 prep ──► 2 CAVER ──► 3 patches ──► 4 hotspot map (look, decide)
                             │
                             ▼
                         5 DO sampling ──► 6 session + SVG
                             │
                             └──────────► 7 generator input
```

## Why stage 4 sits where it does

The hotspot map is shown **before** the DO sampling, so that the epitopes can
be judged — and changed — before committing to a search over them. That one
fact settles the design of the map:

* it cannot take the relative position of the two enzymes from the search's
  output, because at that point there is none;
* so it computes its own reference placement from the two channels alone
  (see below), and
* it must make it easy to say "not that patch, this one", which is why every
  patch is drawn, labelled, and carries the envelope of the shortest binder
  that could sit on it.

Stage 3 is where that answer goes back in. A patch list is just residue
numbers, so accepting the user's own costs nothing downstream.

## The reference placement

O_IC is one-sided — the angle P1-P2-P3 at the **upstream** exit — so
`O_IC = 180°` says only that the downstream exit lies on the ray leaving the
upstream channel. It fixes the downstream exit to `P2 + d·u` and says nothing
about the downstream enzyme's orientation.

So the placement is: put the downstream exit on that ray, scan the downstream
orientation over the pose search's own tilt/azimuth/spin grid, and raise `d`
from the minimum until an orientation exists that does not interpenetrate.
DO falls monotonically with `d`, so the smallest clash-free `d` is the
highest-DO geometry available at O_IC 180.

On the p6a GDH–FNR pair that is **D_IC 9.0 Å, O_IC 180°, DO 0.9314** (tilt 30,
azimuth 240, spin 240; closest approach 3.38 Å). The pose search, run later
over its full grid, independently reaches 8.0 Å at O_IC 180 — one grid step
better, which is the expected size of the difference: the placement's clash
test is a voxel grid, conservative to within one voxel, and its spin grid is
coarser.

Two things this placement deliberately does **not** do:

* it applies no bridgeable-pair, binder-gap or occupancy test. Those are what
  the search is for, and a map that refused to be drawn because no binder
  fits would hide exactly the fact worth seeing. On p6a the ideal geometry is
  one where no patch pair has room for a binder — which is why the search's
  best real pose is at D_IC 17 Å, not 9 Å.
* it is not a pose. It is a frame to look at epitopes in.

## Stage 3: which surface a binder may use

`chbridge patches` classifies every residue of one enzyme into exactly one
of seven classes and then ranks clusters of the usable ones as candidate
epitopes. The classes are a partition, and the ORDER of the tests is the
specification, not an accident — a residue can be in the channel *and* in
the shadow, and which label it gets decides whether the search may put a
binder there:

| | |
|---|---|
| `LIGAND_CORE` | holds the cofactor; touching it breaks the enzyme |
| `CHANNEL` | lines the tunnel, and is exposed |
| `SHADOW` | in the cylinder in front of a mouth, where a binder would cap the channel it is meant to leave open |
| `AVOID` | named by you, or near enough to a named **exposed** residue that a binder centred there would cover it |
| `IFACE` | buried by another chain of the assembly |
| `BURIED` | not exposed enough to bind |
| `ALLOWED` | the usable surface |

Three details that are easy to get wrong and are not guessable from the
output:

* **SHADOW is tested at every opening**, via `any()`, while the table
  reports the primary. An enzyme with two mouths has residues that are
  SHADOW on the second one while sitting far off the primary axis.
* **The classification axis is a fit over the outer spheres**, not the last
  two. Consecutive CAVER spheres are ~1 Å apart, so a pair is dominated by
  placement jitter; the two differ by 9-12°, enough to move residues across
  the SHADOW boundary. O_IC itself still uses the short arm, because it is
  *defined* on one (`Tunnel.axis` vs `Tunnel.fitted_axis`).
* **An AVOID exclusion is grown only from exposed residues.** A catalytic
  residue with 5 Å² of exposed area sits under the surface where no epitope
  can reach it; growing a 12 Å ball from it would delete a third of a usable
  face to protect something already unreachable.

A patch is a cluster of ALLOWED residues around an exposed anchor
(F/Y/W/I/L/M) in a helix or strand, scored on anchors, secondary structure,
hydrophobicity, reach from the mouth and exposure, less a side-chain entropy
penalty. Patches are then taken best-first, skipping any whose centroid is
within `patch_sep` of one already kept, so the result is spread over the
surface instead of six views of one spot.

**Every threshold is a flag** (`--patch-r`, `--shadow-r`, `--exposed`, …)
and every default is in `PatchParams` with its reason. The stage is optional
in the sense that matters: the two tables are plain TSV, so editing them or
writing your own is a supported route, and nothing downstream asks who wrote
them.

## What the score does not see, and what to do about it

`DO = -0.0079 D_IC + 0.0012 O_IC + 0.7865` has two terms, and **O_IC is
one-sided**: it is the angle P1-P2-P3 measured at the *upstream* exit, so it
asks only whether the donor channel is aimed at the acceptor's exit. Which
way the **acceptor** channel opens is not in the score at all. A receiving
mouth pointing 90 degrees sideways, or straight back into its own protein,
scores exactly like one facing the donor squarely.

The pose search samples that orientation anyway -- it is the `tilt`
parameter, gridded `0,30,60,90` by default -- and records the resulting
angle in `poses.tsv` as **`O_IC_rev`**:

| tilt | `O_IC_rev` | |
|---|---|---|
| 0 | ~180 | the two mouths face each other |
| 30 | ~150 | |
| 60 | ~120 | |
| 90 | ~90 | the acceptor opens sideways to the line joining the exits |

So every pose carries the third quantity; nothing scores it. On the worked
GDH-FNR example this matters: of 62,309 poses scanned, 50 are feasible, and
**all 14 best by DO have the acceptor pointing 90-120 degrees away**. The
best feasible pose whose acceptor faces the donor is rank 20, at a cost of
0.0597 DO -- which, at 0.0079 DO per angstrom, is the same as 7.6 A of extra
exit separation. Facing geometries are not unreachable here; they are
unrewarded, and they sit further apart because pointing the acceptor back at
the donor without clashing needs room.

**If you think the acceptor's orientation matters for your chemistry, filter
on it.** `O_IC_rev` is a column in `poses.tsv`, so a one-line filter picks
the poses you want:

```bash
awk -F'\t' 'NR==1 || ($3=="ok" && $8>150)' poses.tsv    # $8 = O_IC_rev
```

**If you want it in the score, supply your own.** `channelbridge.geometry`
takes the coefficients as data (`DOCoefficients`, and `--do-coef` on the
command line), so a two-term rescoring with different constants needs no
code. A genuine third term does need code -- and it would no longer be the
published iMARS score, so give it another name and say what it is. We have
deliberately left `DO` as iMARS defined it, and report the acceptor-side
angle beside it rather than folding it in.

## What each stage reads and writes

| stage | reads | writes |
|---|---|---|
| 1 | PDB/CIF, ligand choices | prepared `<name>.pdb` per enzyme |
| 2 | prepared PDBs | CAVER output, normalised tunnel chains |
| 3 | prepared PDBs, tunnels | `<NAME>_residue_classes.tsv`, `<NAME>_hotspot_patches.tsv` |
| 4 | the above | `<prefix>_hotspot-map.pse` |
| 5 | the above | `poses.tsv`, `patch_pairs.tsv`, `geom<NNN>_*_complex.pdb` |
| 6a | stage 5's output | `<prefix>_<pose>_render.pse` |
| 6b | **any** session, including one you rotated and saved | `<prefix>_<pose>_render_<pair\|assembly>.svg` |
| 7 | stage 5's output | generator input JSON/PDB |

Every stage takes a `--prefix` and puts it at the front of everything it
writes, so two runs over one folder never collide. Stage 6 is split in two
because choosing the camera is a human act: 6a builds the session, 6b renders
from whatever session you point it at. See [figures.md](figures.md).

Stage 2 is not optional because everything after it is defined in terms of a
channel, and stage 5 is not optional because it is the method. Each of the
other five can be handed an answer from outside.

Tunnel normalisation happens inside stage 2 and is on by default: CAVER's exit
position, span and O_IC arm depend on how finely that particular build sampled
spheres, not on the protein. See [caver-builds.md](caver-builds.md).
