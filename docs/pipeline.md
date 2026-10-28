# The pipeline

Seven stages. **Five of the seven are optional**, and each optional stage can
be replaced by the user's own answer rather than skipped — the point is that
nothing downstream has to care which route the input took.

(There is also a **stage 0**, the ChannelsDB lookup. It is numbered 0 rather
than counted among the seven because it is a place to *get* a tunnel, not a
step in the chain — see [Stage 0](#stage-0-look-before-you-compute).)

| # | stage | optional | you can instead |
|---|---|---|---|
| 1 | structure selection and preparation | yes | hand in prepared PDBs |
| 2 | CAVER: tunnels, channel residues, shadow | no | — |
| 3 | hotspot patch selection | yes | name your own patches, or edit the suggested ones |
| 4 | hotspot map session | yes | — (a figure, nothing depends on it) |
| 5 | DO sampling | no | — |
| 6 | complex session and SVG | yes | — (figures) |
| 7 | generator input prep (RFdiffusion3, RFD4-Proteina, BindCraft, ...) | yes | write the tool's spec by hand |

```
 1 prep ──► 2 CAVER ──► 3 patches ──► 4 hotspot map (look, decide)
                             │
                             ▼
                         5 DO sampling ──► 6 session + SVG
                             │
                             └──────────► 7 generator input
```

## What stage 7 must not loosen

The generated spec holds **both** target chains fixed *relative to each
other*, per atom. That is not a stylistic choice about rigidity: the relative
placement of the two enzymes is the entire output of stages 2 to 5, so a
condition that keeps each chain rigid on its own while leaving its global pose
to the model discards D_IC and O_IC and produces designs against a geometry
nobody searched for — without failing, and with plausible-looking output.

For RFD4-Proteina that means `C_CRD` on both chains and never `C_CRDGRP` per
chain, since `C_CRDGRP` preserves a group's internal geometry and explicitly
frees its global pose. `channelbridge.generators` writes the reason into every
file it produces, and refuses tools that have no way to express a fixed
two-chain frame (protpardelle-1c, RFDpoly) rather than emitting a spec they
cannot honour.

Hotspot residue numbers and the contig are renumbered **together** for any
tool that restarts numbering per chain (RFD4-Proteina does: "Output res_id
restarts at 1 per chain"). Renumbering one and not the other is the failure
mode that matters here — on another campaign a hotspot set 184 residues off
still produced designs, against the wrong epitope.

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

On the worked GDH–FNR example, run on the CAVER build this repository is
measured against, that is **D_IC 6.0 Å, O_IC 180°, DO 0.9551** (tilt 30,
azimuth 60, spin 120; closest approach 3.35 Å). 6.0 Å is the scan's own
`d_min`, so the real floor is lower still and this placement is bounded by
the grid rather than by the proteins.

An earlier prototype reported 9.0 Å / DO 0.9314 for the same pair. The difference
is the tunnels, not the placement: that run used a different CAVER build, and
[caver-builds.md](caver-builds.md) measures what that does to an exit. Treat
the absolute number as build-specific and the construction as the thing being
documented.

Two things this placement deliberately does **not** do:

* it applies no bridgeable-pair, binder-gap or occupancy test. Those are what
  the search is for, and a map that refused to be drawn because no binder
  fits would hide exactly the fact worth seeing. On the worked example the ideal geometry is
  one where no patch pair has room for a binder — which is why the search's
  best real pose is at D_IC 15 Å, not 6 Å.
* it is not a pose. It is a frame to look at epitopes in.

## Input preparation: which structure, before any of this

Every stage below measures distances to atoms, so the atoms decide the answer.
Getting them is a stage with its own decisions, and it is the one place where
a wrong choice is invisible afterwards -- a tunnel computed from the wrong
protein looks exactly like a tunnel.

**Preference order, and it is about identity, not resolution:**

1. **An experimental structure of this protein.** Confirm it holds *this*
   protein rather than resembling it: UniProt's `xref_pdb` and RCSB's
   `/core/uniprot/{entry}/1` both say which protein an entry contains.
2. **An AlphaFill model.** An AlphaFold DB model with cofactors transplanted
   from a homologous holo structure, reported with donor, local identity,
   RMSD and clash count. The cofactor matters here: the start point of a
   channel is derived from where the channelled species sits.
3. **An AlphaFold DB model**, apo.
4. **An AlphaFold3 prediction with MSAs.** For anything the databases do not
   hold. Single-sequence mode is the weakest mode for exactly the plant and
   bacterial enzymes that tend to be missing.

**Holding the protein is not the same as holding all of it.** Check how much
of the reference sequence actually has coordinates, over the entry's aligned
span, and reject a large unobserved run: a missing domain measures as a wide
channel. Among several entries of the same protein, prefer a bound cargo, then
wild type, then resolution -- an active-site mutant is a different channel.

**A homolog is not a fallback.** It is a different protein, and the channel is
the property being measured. In the iMARS example set the structure actually
used for most enzymes was a homolog, because the search walked *down* the
sequence-search ranks until it found an entry with ChannelsDB annotation --
trading protein identity for the convenience of a precomputed tunnel. That
is what the first version of `examples/iMARS` did: 22 of its 24 enzymes were
represented by an entry that does not contain them, because the search walked
down the ranks looking for one ChannelsDB had annotated. Every enzyme there
now holds its own protein, and that cost the ChannelsDB branch: it applies to
2 of the 24. Record for every input which of
the four cases it is, and the evidence; `examples/iMARS/structures.tsv` is
that record for the example set.

**Watch the sequence-search score.** RCSB's `score` is this hit's bitscore
divided by the best hit's bitscore in the same search. It is not sequence
identity, and reading it as identity makes a distant homolog look like a
near-identical match.

**Oligomeric state is a decision.** A monomer is the cheap answer and is wrong
wherever the active site is completed across the interface: a type III
polyketide synthase builds its initiation cavity partly from the *partner*
subunit, so a monomeric model has an incomplete site and therefore an
incomplete tunnel. The opposite case is just as real -- the
*trans*-prenyltransferase fold carries a complete site in each subunit, so a
monomer is enough even though the protein dimerises. Record which and why;
once the model exists the choice is invisible.

**Then convert.** Everything above arrives as mmCIF; CAVER and everything
downstream of it read PDB. `chbridge prep` is the one stage that reads both,
and the conversion has traps of its own -- two numbering schemes, and ligands
that are not in their protein's chain. See
[pitfalls](pitfalls.md#input-preparation-mmcif-predicted-and-transplanted-models).

## Stage 0: look before you compute

The first question is not "what does this structure's channel look like" but
"has someone already computed it". ChannelsDB 2.0 holds tunnels for deposited
entries, annotated by what they lead to and computed by both MOLE and CAVER.

It is stage 0 rather than a branch of stage 2 because **it takes no structure
as input**. You search by enzyme name or PDB id, and the answer decides which
entry you go on to prepare -- so preparing first would be preparing the wrong
thing, or preparing something you did not need.

Two consequences worth stating before anyone relies on it:

* **Coverage is per ENTRY.** 4-coumarate:CoA ligase has thirteen PDB entries
  and two carry tunnels; 5BST has none, 5BSR has 25. Search the enzyme, not
  the structure you happen to hold -- and search past the first hit. Over the
  24 enzymes of [the iMARS example set](../examples/iMARS/README.md), checking
  198 candidate entries, only **4** of the 24 top-ranked hits were annotated,
  while **22** had an annotated entry somewhere in their top twelve. Going
  down the ranks is the difference between the branch being available and
  not. It is also a weaker claim: a rank-11 hit at 0.52 sequence-search score
  is evidence about the family, not about that enzyme, so record the rank and
  the score with the choice.
* **Taking their tunnel means taking their structure.** The sphere chain is
  in that entry's frame, and nothing re-fits it. Prep may extract chains and
  delete waters; it must not move coordinates. Expanding a biological
  assembly moves them, and the result looks like a 19-40 A disagreement
  rather than like the mistake it is.

If nothing is annotated, run CAVER (stage 2) -- which is the normal case for
a model, a mutant, or anything unpublished.

**And check what you fetched before measuring from it.** The database
truncates a chain by its own surface definition, which is not this package's:
measured over the first, homolog-based version of the iMARS set -- which had
21 annotated enzymes, against 2 today -- 13 of those 21 had no chain whose
ends reached bulk solvent as `truncate_at_surface` defines it. D_IC and O_IC are
both measured at the mouth, so the branch is sound for finding and
cross-checking a route and needs care before it is used to measure an exit.
The normalise step's verdict column says which case you are in.

## Stage 2b: normalisation, and why it is its own step

CAVER gives you a chain of probe spheres. Three things about that chain are
properties of the RUN rather than of the protein, and all three feed straight
into D_IC and O_IC:

* **Direction.** The chain's order is not guaranteed to run inner to outer,
  and everything downstream takes the last sphere as the exit. `oriented()`
  fixes the direction from burial.
* **Where it stops.** CAVER runs the chain out to `shell_radius`, which is
  deliberately set large enough that it certainly leaves the protein -- and
  therefore ends some arbitrary distance out in solvent.
  `truncate_at_surface()` cuts it back to where it actually crosses the
  surface. Measured over a `shell_radius` sweep: raw D_IC drifts by up to
  **19 A** across shell radii 3-16 (34.7 / 28.7 / 21.8 / 15.6 A on one
  construct) while the
  normalised value stops moving once the chain clears the surface (25.76 vs
  26.25 A between shell 10 and 16, where the raw value still moved 6.2 A).
* **Sphere spacing.** The O_IC arm is defined as an arc length, so its
  direction depends on how densely that particular build sampled.
  `resample()` puts the spheres at a uniform 1 A spacing.

**Why it is separate from the run.** Two reasons, and the second is the one
that shapes the pipeline:

1. CAVER does not do it. It is our post-processing, and keeping it separate
   means the raw output stays on disk exactly as the tool wrote it.
2. **It is not CAVER-specific.** A tunnel from ChannelsDB needs it just as
   much: MOLE's endpoint convention is not ours either, and the same chain
   from the two sources differs in length (our GDH cluster 1 is 24.2 A raw,
   the matching ChannelsDB tunnel 15.9 A). So normalisation sits *after* the
   two sources merge, not inside the CAVER branch.

## Stage 3: which surface a binder may use

Nothing reaches the pose search without this stage. The scan places a binder
on a (patch, patch) pair and marks any pose with none as
`no_bridgeable_pair`, so there is no route from a tunnel to a pose that skips
patch selection.

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

## Stage 5: the pose search

The method. Everything before it describes two enzymes; this searches for a
relative placement a bispecific binder could hold, and scores each candidate
in iMARS terms. It is the one stage that cannot be handed an answer from
outside, and it needs no PyMOL -- it reads what stages 1-3 wrote and is pure
numpy.

```bash
chbridge pose --pair-file examples/gdh_fnr/pair.yaml \
    --up-pdb GDH.pdb --dn-pdb FNR.pdb \
    --tunnel-up GDH_normalised.pdb --tunnel-dn FNR_normalised.pdb \
    --patches-up GDH_hotspot_patches.tsv --patches-dn FNR_hotspot_patches.tsv \
    --out-dir out-complex_GDH_FNR
```

The pair file supplies the grid, the DO coefficients and the O_IC arm; any
flag overrides it, so asking "what if" costs no edit to the record of what a
run was.

### Five parameters, and what each one decides

The downstream enzyme moves and the upstream one never does. The
parameterisation is the metric's, not a docking program's:

| | |
|---|---|
| `d_ic` | separation of the two exits, held exactly |
| `o_ic` | the published angle at the **upstream** exit |
| `azimuth` | where on that cone the downstream enzyme sits; degenerate at `o_ic = 180`, and collapsed there |
| `tilt` | how far the downstream channel turns from facing back |
| `tilt_azimuth` | which way it tilts -- see below |
| `spin` | roll about the downstream channel: changes no scored quantity, decides the downstream epitope |

### The axis a pose is built along

A pose is built along `Tunnel.axis` -- the exit minus the point one O_IC arm
back -- which is the same pair of points `geometry.o_ic` measures the result
with. `Tunnel.fitted_axis` is steadier and sits 9-12 degrees away from it on
real tunnels (see [caver-builds.md](caver-builds.md)), so building on one and
scoring on the other leaves the grid asking for one O_IC and recording
another. Building on the scored direction makes "requested" and "realised"
the same number. On the iMARS pairs they agree to 1e-14 angstrom in D_IC and
1e-6 degree in O_IC, and `tests/test_pose.py` asserts it.

### Four filters, cheapest first

A pose is feasible only if it survives all of them, and `poses.tsv` records
which one stopped it:

1. `clash` -- the two enzymes interpenetrate.
2. `no_binder_gap` -- they are closer than `binder_gap`, so nothing can run
   between them.
3. `no_bridgeable_pair` -- no epitope on one is the right distance from an
   epitope on the other *with a clear line between them*. This is the filter
   the whole stage exists for: a geometry can score beautifully and have
   nowhere to put a binder.
4. `occupancy_clash` -- where the upstream enzyme is a homomultimer, the pose
   must work on every protomer at once. One binder per protomer means that
   many copies of the downstream enzyme around one particle, and nothing in
   D_IC, O_IC or the single-pose clash test sees whether they fit. Run last,
   because it is much the dearest.

   This filter is not a formality. Measured on the two equivalent protomers
   of GPPS (3KRF chains A and D) against LIC, **24 of 25 otherwise-feasible
   poses fail it** -- every one because the two downstream copies would
   interpenetrate, at 0.1-0.2 A of separation. The survivor sits at D_IC
   36 A with 4.3 A between copies: requiring the pose to work on both
   protomers at once pushes the pair from DO 0.84 to 0.68. A search without
   this test would have reported the 0.84 geometry.

The grid is ordered **best-DO-first**, so the first feasible pose is the best
one and `--max-ok` truncates the scan without biasing the answer. Only
feasible poses are issued a `geomNNN` id: numbering rejections alongside
successes once made `geom050145` the best design in a 100,000-point grid,
which tells you nothing.

### How the proximity test stays exact

The scan asks "is anything closer than *x*" tens of thousands of times.
Exactly, that is 4,000 x 3,000 distances per pose -- about 0.3 s, so hours of
grid. A voxel lookup is one array index per atom, and the usual pattern of
conservative-prefilter-then-exact-check is the wrong way round here: nearly
every pose in a best-DO-first grid interpenetrates, so "conservative" would
send nearly every pose to the exact check and save nothing.

So each threshold gets two grids, eroded and dilated by one voxel
half-diagonal. A hit in the eroded grid means *definitely* closer; a miss in
the dilated grid means *definitely* not; between them the exact distance is
measured. **The verdict is therefore exact everywhere**, and the exact
measurement only runs near the boundary.

Where it does run, it runs on the atoms that can matter: the downstream atoms
within 12 A of the upstream enzyme, against the upstream atoms in a box
around those. That is still exact for every distance a verdict turns on, and
on the first iMARS pair it took the scan from 95.5 s to 8.6 s — **11 times
faster, with byte-identical results**: same status counts, same best pose,
same DO. When nothing at all is within the radius the narrowing has nothing
to offer and the full computation runs, rather than a bound being reported as
a measurement.

One consequence shows in the output: `min_interprotein_dist` is `-` for a
pose rejected without measurement. The verdict is exact and the number was
never needed; writing a voxel bound into a column headed "distance" would be
a worse answer than leaving it out.

### The orientation the default grid does not sample

O_IC constrains where the downstream *exit* goes and says nothing at all
about the downstream enzyme's own orientation, which is three free
parameters. `tilt` and `spin` are two. The third is `tilt_azimuth`, which way
the channel tilts -- and its default is a **single value**, so the default
grid reaches the downstream channel directions on one arc rather than over
the sphere.

That default is for grid size, not because the axis does not matter: sweeping
it costs twelve times the poses. It is the cheapest dimension to add when a
pair comes back with nothing feasible, because no scored quantity changes
along it -- only the clash and epitope tests. `examples/iMARS/README.md`
records what it bought there.

## What the score does not see, and what to do about it

`DO = -0.0079 D_IC + 0.0012 O_IC + 0.7865` has two terms, and **O_IC is
one-sided**: it is the angle P1-P2-P3 measured at the *upstream* exit, so it
asks only whether the donor channel is aimed at the acceptor's exit. Which
way the **acceptor** channel opens is not in the score at all. A receiving
mouth pointing 90 degrees sideways, or straight back into its own protein,
scores exactly like one facing the donor squarely.

The pose search samples that orientation anyway -- it is the `tilt`
parameter, gridded `0,30,60,90` by default -- and records the resulting
angle in `poses.tsv` as **`O_IC_rev`**.

| tilt | `O_IC_rev` | |
|---|---|---|
| 0 | 180 | the two mouths face each other |
| 30 | 150 | |
| 60 | 120 | |
| 90 | 90 | the acceptor opens sideways to the line joining the exits |

Those are exact, not approximate, and that is a property worth having: it
holds because the pose is *built* along the same O_IC arm the angle is
*measured* on. Build it along the fitted axis instead and the column drifts
9-12 degrees away from `180 - tilt` -- see
[The axis a pose is built along](#the-axis-a-pose-is-built-along).

So every pose carries the third quantity and nothing scores it. How much
that costs is a per-pair measurement, not a constant, and on the worked
GDH-FNR example it is **small**: of 58,130 grid points scanned before the
50-pose budget ran out, 19 of the 50 feasible poses have the acceptor facing
the donor (`O_IC_rev > 150`), and the best of them is **rank 5, at a cost of
0.0136 DO** -- 1.7 A of equivalent exit separation. Six of the top 14 point
the acceptor 90-120 degrees away; the rest face.

An earlier prototype reported this much more strongly on the same pair: all 14
best by DO pointing 90-120 degrees away, and the best facing pose at rank 20
costing 0.0597 DO. That does **not** reproduce here, and the likely reason is
the one running through this whole project -- a different CAVER build gives a
different exit, so the set of reachable geometries differs
([caver-builds.md](caver-builds.md)). Two things follow. The qualitative
point stands on its own: the metric is one-sided, so filter on `O_IC_rev` if
the acceptor's orientation matters to your chemistry. The quantitative claim
does not travel: measure it on your own pair and build.

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
| input prep | an enzyme name or accession | a structure of *that* protein, and the record of which of the four cases it is |
| 0 | an enzyme name or a PDB id | whether tunnels already exist, and for which entry |
| 1 | PDB **or mmCIF**, chain and ligand choices | prepared `<name>.pdb` per enzyme |
| 2a | prepared PDBs, a starting point | CAVER output (`summary.txt`, `data/clusters_timeless/`) |
| 2b | CAVER output or a fetched chain | normalised tunnel chains |
| 3 | prepared PDBs, tunnels | `<NAME>_residue_classes.tsv`, `<NAME>_hotspot_patches.tsv` |
| 4 | the above | `<prefix>_hotspot-map.pse` |
| 5 | the above | `poses.tsv`, `patch_pairs.tsv`, `geom<NNN>_*_complex.pdb` |
| 6a | stage 5's output | `<prefix>_<pose>_render.pse` |
| 6b | **any** session, including one you rotated and saved | `<prefix>_<pose>_render_<pair\|assembly>.svg` |
| 7 | stage 5's output, plus stage 3's patch tables | one spec per pose per tool (`p_<pose>.json` / `.yaml` / `.cif`) |

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
