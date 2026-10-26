# Two CAVER builds, one tunnel, different numbers

Measured 2026-10-23 on the GDH primary route (PDB 1RWB assembly 1, prepared
tetramer, `config_GDH_primary_sr14`: probe 1.4 A, shell_radius 14,
shell_depth 6, clustering 3.5, 12 balls, seed 1, excluding HOH WAT BGC NAP).
Same structure file, byte-identical config, same seed. Two CAVER builds.

## What differed

| property | build A | build B |
|---|---|---|
| spheres in cluster 1 | 16 | 145 |
| mean sphere spacing | 0.98 A | 0.17 A |
| span (first to last sphere) | 12.1 A | 19.1 A |
| burial at the outermost sphere | 21 | 0 |
| burial verdict | surface mouth | past the surface |
| bottleneck radius | 2.915 A | 2.910 A |
| O_IC arm (last sphere step) | 0.91 A | 0.26 A |
| clusters kept | 3 | 9 |
| `FNR_primary` negative control | no tunnels | 2 clusters |

Build B is **CAVER 3.03 BETA** from [caver.cz](https://caver.cz/) (jar built
2016-04-24), run under Java 21.

Build A is the build that produced this project's earlier reference
measurements, on Windows. **Its exact version was not recorded** — which is
itself the first lesson of this document: record the CAVER version beside
every tunnel measurement. It is **older** than build B and behaves like a
stable release: SBGrid's CAVER listing offers 3.0.1 and 3.0.2.

Two independent signs that build B is the newer:

* Its jar contains the parameters `correct_voronoi_diagram` and
  `long_residue_names`, both introduced in version 3.01 per the `readme.txt`
  shipped beside it.
* `FNR_primary` is a deliberate negative control whose start point sits
  inside an FAD atom. Version 3.01's note says the starting-point
  optimisation "was corrected... the starting vertex was sometimes identified
  in too narrow cavity due to error in implementation". Build A failed, as
  documented; build B relocates the start and finds two clusters.

Note that `caver/license/license.txt` says "CAVER 3.0" in both. It is the
unchanged 2012 GPL notice and does not identify the point release, so it
cannot be used to tell builds apart.

### That build B is a BETA matters

The ~0.17 A sphere spacing and the chain running ~9 A past the surface look
like deliberate behaviour rather than breakage — the route is identical and
the bottleneck agrees to 0.005 A — but they are behaviour of a pre-release,
and no changelog for 3.03 is published alongside the download.

This is the strongest argument for normalising rather than depending on the
build. A pipeline that reads "the outermost sphere" is pinned to whichever
CAVER it was developed against; one that truncates by burial and resamples to
a fixed spacing gets the same answer from 3.0.x, from 3.03 BETA, and from
whatever follows. A methods section should still state the CAVER version and
that it was a beta.

## It is the same tunnel

The disagreement is not about where the channel runs.

* Build B's chain passes **0.078 A** from build A's exit.
* The burial there is **21**, build A's recorded value exactly.
* Build B's chain then continues a further **~9 A** into bulk solvent, where
  burial falls to 0.
* The bottleneck agrees to 0.005 A, because it is found in the interior where
  the shell parameters have no influence.

Build A stopped at the surface. Build B kept going. Everything else follows
from that plus the ten-fold difference in sphere spacing.

## Why it matters

Three numbers the pipeline depends on are therefore **not** properties of the
protein:

1. **The exit position**, and so **D_IC** directly. An exit 7 A further out
   changes D_IC by up to 7 A, which at the published coefficients is 0.055 in
   DO -- larger than the entire spread between the top ten poses.
2. **`burial_exit`** and the band verdict derived from it, which is how a
   tunnel end is judged to be a real mouth.
3. **The O_IC arm.** The published construction takes P1 as the
   second-to-last sphere. On build A that is a 0.91 A baseline; on build B,
   0.26 A. Sphere-placement jitter is amplified roughly four-fold, on an
   angle whose grid step is 15 degrees.

A pipeline that hard-codes "the outermost sphere" and "the second-to-last
sphere" silently inherits all three.

## What ChannelBridge does about it

`channelbridge.tunnel.Tunnel.normalised()` applies three steps, all on by
default and all recorded in the tunnel's `history`:

1. **`oriented(protein)`** -- order the chain inner-to-outer by **burial**,
   not by distance from the centroid, which fails on a multimer whose
   centroid is not inside the protomer.
2. **`truncate_at_surface(protein)`** -- take the exit where the chain first
   leaves the solvent-excluded surface, scanning **outward from the buried
   inner end** and stopping where burial drops below 15 and stays there for
   2 A. Scanning outward matters: a neighbouring subunit or an additive out
   in the solvent also raises burial, and scanning inward from the far end
   would accept it and place the exit past the real mouth.
3. **`resample(1.0)`** -- re-space the chain at a fixed arc length, both
   endpoints preserved exactly. 1.0 A is the density the published iMARS
   construction was calibrated on.

and `axis_point(arm)` takes P1 a fixed arc length inward from the exit,
interpolated along the path, so the O_IC arm is a physical length rather than
one sphere step.

### Result

| | raw | normalised |
|---|---|---|
| exit separation between builds | **8.94 A** | **1.44 A** |
| build A verdict | surface mouth | surface mouth |
| build B verdict | past the surface | surface mouth |
| O_IC arm, build A / build B | 0.91 / 0.26 A | 0.99 / 1.00 A |

Reproduce it with:

```bash
chbridge tunnel \
  /path/to/buildA/out_GDH_primary_sr14 \
  /path/to/buildB/out_GDH_primary_sr14 \
  --structure gdh.pdb --exclude HOH,WAT,BGC,NAP --cluster 1 --compare
```

`tests/test_crossbuild_regression.py` asserts all of the above against both
real outputs, and is skipped when they are absent.

## What is *not* fixed

D_IC is reconciled. **O_IC is not made build-exact**, and should not be
claimed to be.

After identical normalisation, the two builds' mouth directions still differ
by **15.5 degrees** at a 1 A arm. The cause is physical: this tunnel curves as
it opens, and the two builds trace slightly different paths through that
curve. Sweeping the arm shows the disagreement is not even monotone:

| O_IC arm | angle between the two builds' mouth directions |
|---|---|
| 0.5 A | 13.9 deg |
| 1.0 A | 15.5 deg |
| 1.5 A | 28.7 deg |
| 2.0 A | 35.6 deg |
| 3.0 A | 26.9 deg |
| 5.0 A | 12.0 deg |
| 8.0 A | 5.3 deg |

So "use a longer arm" is not a rule to apply blindly. The mid-range is the
worst, because that is the scale on which build A's chain curves and build
B's does not: build A's own axis rotates 20.7 degrees between a 1 A and a 2 A
arm, build B's rotates 0.07 degrees.

`Tunnel.axis_spread()` reports exactly that, as a quality flag in every
tunnel report (`oic_arm_spread_deg`). A large spread means the tunnel's mouth
direction is scale-dependent and the O_IC should be quoted with that
uncertainty. At the published coefficients 15.5 degrees of O_IC is worth
2.3 A of D_IC, about 0.019 in DO.

### The 15 degrees is systematic, not noise

CAVER's `seed` changes the Voronoi construction enough to alter the cluster
count, so the obvious worry is that 15 degrees is simply stochastic scatter.
It is not. Running build B over seeds 1-8 on the same structure and config,
matching clusters on the bottleneck coordinate, and normalising each:

| quantity | spread across 8 seeds |
|---|---|
| exit position | **0.007 A** (max pairwise) |
| mouth axis at a 1 A arm | **0.06 deg** (max pairwise) |
| bottleneck radius | 0.0000 A |
| spheres before / after normalisation | 145 / 17 in every run |

Within one build the geometry is effectively deterministic. The 15.5 degrees
between builds is therefore a reproducible, version-specific difference in the
path traced through the curved mouth region -- a systematic offset, not an
error bar that shrinks with repetition.

Reproduce the seed check for another pair with:

```bash
for s in $(seq 1 8); do
  sed "s/^seed .*/seed $s/" config_X.txt > cfg_seed$s.txt
  java -Xmx8000m -cp "$CAVER/lib" -jar "$CAVER/caver.jar" -home "$CAVER" \
       -pdb ./in -conf cfg_seed$s.txt -out ./out_seed$s > log_seed$s.txt 2>&1
done
chbridge tunnel out_seed*/ --structure in/target.pdb --exclude HOH,WAT --compare
```

Match clusters on the bottleneck coordinate, not the cluster number: the
numbering is by priority and is not stable across seeds (7 to 12 clusters
appeared over the eight runs above, all containing the same primary route).

Two practical consequences:

* **Do not sweep the seed.** It buys nothing; `seed 1` is as good as any.
  Sweep `shell_radius` instead, which is what actually moves the exit.
* **A published O_IC is version-dependent at the ~15 degree level** on a
  curved mouth, which is 2.3 A of D_IC or 0.019 in DO. The CAVER version
  belongs in the methods next to the number.

The honest statement for a methods section is therefore: D_IC is reproducible
across CAVER builds once the chain is normalised, and reproducible to 0.007 A
within a build; O_IC is deterministic within a build but carries a systematic
tunnel-specific offset of order 15 degrees between builds on a curved mouth,
which should be propagated into DO rather than ignored.

## Consequence for the pose search

Run end to end on build B, the whole pose search picks a different answer from
the raw chain and build A's answer from the normalised chain.

| | build A (reference) | build B, normalised | build B, raw |
|---|---|---|---|
| best pose DO | 0.8401 | **0.8502** | 0.8695 |
| D_IC | 16.0 A | **17.0 A** | 10.0 A |
| O_IC | 150 deg | **165 deg** | 135 deg |
| bridged epitope pair | GDH_p4 + FNR_p1 | **GDH_p4 + FNR_p1** | GDH_p3 + FNR_p1 |

The normalised run reproduces build A's **hotspot patch table exactly** -- all
four upstream patches identical in score and geometry factor (4.54 / 4.39 /
4.04 / 3.41, every geometry factor 1.00) with identical residue lists, and
five of six downstream patches likewise -- and bridges the same epitope pair.

The raw chain distorts that table to 4.21 / 3.90 / 3.41 / 3.04 with geometry
factors 0.88 / 0.91 / 1.00 / **0.00**, the zero falling on the very patch
build A bridges. With the mouth ~7 A out in bulk solvent, the correct epitope
looks geometrically inaccessible, so a different pair wins.

The residual difference against build A is 1.0 A of D_IC and 15 deg of O_IC
(one grid step), +0.0101 in DO. Both match the independent tunnel
measurements: the normalised exits sit 1.44 A apart, and the normalised mouth
directions differ by 15.5 deg. Nothing unexplained remains.

### A higher DO from a raw chain is an artefact

Note that the raw run scores **highest**. It measures D_IC from an exit in
bulk solvent, which lets the downstream enzyme sit closer in the scored
geometry while the real channel mouths are further apart. DO is therefore not
comparable between a normalised and an un-normalised run, and a higher DO from
an untreated chain means the exit escaped the protein, not that the
channelling geometry improved.

## References

CAVER: Chovancova E, Pavelka A, Benes P, Strnad O, Brezovsky J, Kozlikova B,
Gora A, Sustr V, Klvana M, Medek P, Biedermannova L, Sochor J, Damborsky J
(2012). CAVER 3.0: A Tool for the Analysis of Transport Pathways in Dynamic
Protein Structures. *PLoS Comput Biol* 8: e1002708.
doi:10.1371/journal.pcbi.1002708
