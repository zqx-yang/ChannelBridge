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

## What each stage reads and writes

| stage | reads | writes |
|---|---|---|
| 1 | PDB/CIF, ligand choices | prepared `<name>.pdb` per enzyme |
| 2 | prepared PDBs | CAVER output, normalised tunnel chains |
| 3 | prepared PDBs, tunnels | `<NAME>_residue_classes.tsv`, `<NAME>_hotspot_patches.tsv` |
| 4 | the above | one `.pse` |
| 5 | the above | `poses.tsv`, `patch_pairs.tsv`, `geom<NNN>_*_complex.pdb` |
| 6 | stage 5's output | `.pse` and `.svg` per chosen pose |
| 7 | stage 5's output | generator input JSON/PDB |

Stage 2 is not optional because everything after it is defined in terms of a
channel, and stage 5 is not optional because it is the method. Each of the
other five can be handed an answer from outside.

Tunnel normalisation happens inside stage 2 and is on by default: CAVER's exit
position, span and O_IC arm depend on how finely that particular build sampled
spheres, not on the protein. See [caver-builds.md](caver-builds.md).
