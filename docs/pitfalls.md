# Pitfalls

Failure modes this pipeline has hit, each one either guarded by code in this
package or -- where the cause lies in a tool outside it -- stated here so it is
not rediscovered. They are rules, not anecdotes: everything below changes what
a user should do.

## CAVER

**Launch it with `-jar`, and pass `-home`.** The entry point is named in
the jar's own manifest (`caver.ui.Launcher`); there is no `CaverRunner` class
to put on a classpath. `-home <the caver directory>` is not optional either,
because CAVER reads `bin/atom_radii.csv` relative to it. The working command
is `java -Xmx6g -jar caver.jar -home <dir> -pdb <dir> -conf <f> -out <dir>`.

**A single structure writes `data/clusters_timeless/`, not `data/clusters/`.**
The name follows whether the run had a time dimension, and the files carry a
trailing snapshot index (`tun_cl_001_1.pdb`, not `tun_cl_001.pdb`). Code that
looks only in `clusters/` reports every static run as having found nothing --
which is indistinguishable from a bad starting point.

**There is no default starting point and no complaint about missing one.**
At least one of `starting_point_coordinates`, `starting_point_residue` and
`starting_point_atom` is required; they COMBINE as an equally-weighted centre
of gravity rather than overriding each other. `write_config` refuses to write
a config without one.

**`starting_point_residue` averages whole residues, which can land inside an
atom.** A catalytic triad's residues face each other, so their centre of
gravity is in the middle of the side chains: 0.61 A from the nearest atom
centre for stilbene synthase (3A5S) and 0.49 A for curcuminoid synthase
(3ALE). CAVER looks for a Voronoi vertex within `max_distance` (3 A), fails,
and reports **zero tunnels** -- at probe 0.9 and 1.4 and at shell radius 3, 5
and 10, with nothing in the log about the starting point. Use
`chbridge caver config --start-cavity` to move the point into the pocket
first; 3.7 A and 4.0 A respectively, and both then give four clusters.

**A start point from lining residues is not interchangeable with one from a
ligand.** On GDH the two sit 4.05 A apart and that swapped clusters 1 and 2,
which were tied at priority 0.773. The route survived -- cluster 2 of the
second run exits 0.00 A from cluster 1 of the first -- but the "top tunnel"
exit moved 10.02 A, worth 0.08 in DO. Record which you used.

**`-pdb` must be a directory, never a file.** Handing it a file makes Java's
`File.listFiles()` return null and CAVER dies with a bare
`NullPointerException at chemistry.pdb.PdbUtil.getTrajectoryFiles`, with no
explanation.

**One PDB per `-pdb` directory.** Every PDB in that folder is treated as a
frame of the same trajectory, so a folder holding two structures silently
analyses the wrong thing.

**Every run needs its own `-out`.** CAVER overwrites without refusing and
without warning. Two runs pointed at one `-out` leave a full-looking folder
of results for the wrong structure.

**Line endings must be LF.** A PDB with CRLF gives "No PDB files found" and
**exit code 0**, so a runner believes it succeeded. `channelbridge.io.pdb`
writes LF unconditionally and `has_crlf()` checks an input.

**CAVER reads PDB only.** No mmCIF. `PairFile.validate_against_disk()`
rejects a `.cif` target.

**Always read `warnings.txt`.** CAVER exits 0 either way.

**A start point inside an atom yields no tunnels** on older builds, and a
relocated start point on 3.01 and later. Check clearance first; see
[caver-builds.md](caver-builds.md).

**A zero-tunnel run is a legitimate outcome**, not a failure: a probe of
2.5 A on a route that needs the protein to move finds nothing, and that is a
result. A runner must not abort the remaining runs when `data/clusters/` is
absent -- a sweep runner written under `set -euo pipefail` will, and loses
the remaining structures.

**Burial must be counted against the atoms CAVER actually searched** -- the
input minus `exclude_residue_names`. Counting a cofactor CAVER never saw
inflates it. `AtomArray.exclude_resnames()` exists for this; the `--exclude`
flag on `chbridge tunnel` must match the config.

**A prosthetic group is never excluded.** FAD, a haem, a metal: these are part
of the protein. Strip FAD from FNR and CAVER finds a large false cavity
exactly where the chemistry happens. The pair file has a `prosthetic` field
that is protected in every config.

**The bottleneck is quotable across builds; the span and the exit are not.**
The bottleneck sits in the interior where the shell parameters have no
influence. See [caver-builds.md](caver-builds.md).

**Match tunnels across a `shell_radius` sweep on the bottleneck coordinate,**
not on the cluster number. Cluster numbering is by priority and is not stable
when the shell parameters change.

**Which cluster to use is a decision, not a default.** Keep clusters whose
spheres lie within ~4 A of the excluded cofactor (the genuine-channel test),
then choose between survivors empirically.

## ChannelsDB

**`Distance` is not the same quantity in the two methods.** In MOLE entries
it is cumulative arc length, as you would expect. In the CAVER-derived
entries it is non-monotonic and saturates near 1 -- a per-step or normalised
value. Reading `Distance[-1]` as a length therefore reported every
CAVER-derived tunnel as about 1 A long, which is plausible enough that
nothing contradicted it; the symptom was a selection rule preferring
two-sphere "tunnels" over real ones. `DbTunnel.length` measures the
coordinates, and `DbTunnel.distance_is_arc_length` says whether the reported
series can be used.

**A fetched chain often stops before the surface.** The database truncates by
its own surface definition, not this package's. Measured over the first,
homolog-based version of `examples/iMARS`, which put 31 enzyme sides on that
branch: only 10 had a candidate chain whose ends fell below burial 15, and on
5ICE every Procognate chain has a minimum burial of 29-71 atoms along its
whole length. D_IC and O_IC are both measured at the
mouth, so check the verdict (`mouth` usable; `buried`, `recessed`, `surface`
not) before quoting geometry from a fetched tunnel.

**What to do about it: seed CAVER from it, do not measure it.** The database
identifies the route (its category says what the route is FOR, which CAVER
cannot) and CAVER produces a chain that reaches solvent (which the database
does not). Run CAVER from the database chain's buried end, then check the
result back against that chain. Two cautions, both measured over the 13
seeded enzymes of that first version: seed from the chain END, not from a
bound ligand -- seeding from the ligand first dropped annotation agreement
from 6 enzymes to 2, because CAVER then finds whatever leaves the ligand
site; and judge the check with
`Match.recovery()`, not `Match.verdict()`, because your chain is longer by
design and the exit-led verdict calls a recovered route "different tunnel".

**A fetched tunnel may belong to a different protomer.** A report annotates
the whole deposited assembly and never says which chain each tunnel is in.
Prepping chain A of 3KRF and taking its best Procognate tunnel gave a chain
32 A from the protein at BOTH ends: correct coordinates, in a frame whose
protein had been deleted. It passes a "does it leave the protein" test
trivially. Check that a fetched chain comes within a few angstroms of the
chain you kept, and if it does not, prepare the chain the annotation is for.

**Coverage is per entry, not per protein.** Only 4 of 24 rank-1 sequence-
search hits in the iMARS set are annotated; going down to rank 12 found 22.
Search the enzyme, then pick an annotated entry -- but note that a rank-11
hit is a weaker claim than a rank-1 hit at the same score.

**Two tunnels can share an entry across "different" enzymes.** 4'OMT and
7OMT in the iMARS set have the same single annotated candidate. Putting both
on it would have modelled two enzymes on one set of coordinates, which looks
like two measurements and is one.

## Geometry

**O_IC is one-sided.** It is measured at the upstream exit only, so
`upstream:` in the pair file changes the score. A short *downstream* tunnel is
not a defect: only the upstream tunnel's outer direction orients anything.

**Do not take P1 as "the second-to-last sphere".** That is a sphere step, not
a length, and it varies ten-fold between CAVER builds. Use
`Tunnel.axis_point(arm)`.

**A curved mouth has no single direction.** Check
`oic_arm_spread_deg` before quoting an O_IC to better than ~15 degrees.

**The two published DO constant sets disagree by up to 0.032** over the grid
this pipeline scans, which is larger than the spread between top-ranked
poses. `do_coefficients:` records which was used.

**DO is a linear fit and is not clamped.** A value outside `[0, 1]` means the
pose is outside the range the fit was made on.

## The pose search

**Build a pose along the axis you will score it with.** `Tunnel.axis` (the
O_IC arm) and `Tunnel.fitted_axis` differ by 9-12 degrees on real tunnels. If
the grid places the downstream exit along one and `geometry.o_ic` then
measures against the other, every row of `poses.tsv` records an O_IC that was
never requested, and nothing in the output looks wrong. Assert the two agree;
`tests/test_pose.py` does.

**O_IC constrains the downstream exit, not the downstream enzyme.** Its
orientation is three free parameters and the grid samples `tilt` and `spin`,
which are two. The third -- which way the channel tilts -- is
`tilt_azimuth`, and it defaults to a single value, so the default grid
reaches the downstream channel directions on an arc rather than over the
sphere. It is the first thing to sweep when a pair returns nothing feasible.

**A pose that scores well may have nowhere to put a binder.** DO sees two
exits. It does not see whether any pair of epitopes is both the right
distance apart and has a clear line between them. On the worked example the
best-DO geometry is at D_IC 9 A and the best *buildable* pose is at 17 A --
so a search without the bridgeable-pair test returns the wrong answer
confidently.

**A homomultimer presents the binder once per protomer.** One good pose on
one protomer can be impossible on all four at once, and neither D_IC, O_IC
nor the single-pose clash test sees it: that test already covers one
downstream enzyme against the whole upstream assembly, but never one
downstream copy against another.

**Judge a voxel clash test at its boundary, not on average.** A grid built by
growing whole voxels around snapped atom positions carries about 1.7 A of
slop at 1 A spacing, which is half a clash cut-off. Erode and dilate by one
voxel half-diagonal and measure exactly between the two, or the verdict is
wrong by a few tenths of an angstrom in a way no output reveals.

**Number only the poses you could build.** Issuing ids across the whole grid
made `geom050145` the best design in a 100,000-point scan: a number that
encodes where in a grid a pose sat, not how good it is.

**Do not edit the code a running scan is calling.** Each pair here is a fresh
subprocess, so a mid-run edit is picked up by every pair after it -- while
the ones before it ran the old code. Fifteen of nineteen pairs crashed this
way, on an argument added to a function before it was added to the parser.

## Re-running a stage after its inputs change

**Key a skip on timestamps, not on existence.** A stage that reuses its
output because the file is there cannot see that its input changed, so a
corrected run is answered out of the cache it was meant to replace. This has
bitten this project three times -- the tunnel/patch/pose stages, the CAVER
run directory, and the SVG's helper renders -- and the fix is one line each
time: compare mtimes against every input, including the config.

**Never judge success by the existence of an output file.** Judge it by what
the step reported. A stage that printed `ok` for 38 runs that had all exited
127 is the worst case, and the tell was that the `ok` lines carried no
counts.

## Input preparation (mmCIF, predicted and transplanted models)

**An entry can hold the right protein and not hold all of it.** Identity is
necessary and not sufficient. 4CL's only experimental structure is a
4CL::stilbene-synthase *fusion* at 3.1 A in which 4CL's last 102 residues were
never modelled -- and 4CL's C-terminal domain is the one that closes over its
active site. A tunnel computed there runs out through a hole where a domain
should be, and reports as a wide channel. Require a minimum fraction of the
aligned reference span observed, and a maximum unobserved run. Measure
coverage over the **aligned span**, not the whole UniProt length: a cleaved
transit peptide is correct biology, not 10 % missing.

**The first cross-referenced entry is not the right entry.** Coverage is a
veto, not a ranking. Among entries that pass it, prefer a bound cargo (the
start point comes from where the channelled species sits), then wild type (an
active-site mutant is a different channel -- which is why someone made it),
then resolution. Taking the first cross-reference for BAS gives the apo 3A5Q
over the holo wild-type 3A5R at 1.60 A, and can give 3A5S, the I207L/L208F
mutant.

**AlphaFill only outranks a plain model if its transplant is usable.** It is
preferred *because* a ligand-derived start point survives, so a transplant
that clashes into the model's own side chains cannot carry that preference:
4CL's ATP, AMP and CoA arrive from a 0.69-identity homolog at 2.6, 4.1 and 4.4
clashes per heavy atom. Normalise the clash count by atom count before
comparing ligands, and count only the **best-identity tier** of donors --
accepting any usable transplant at any identity is picking by clash count,
which lands on buffer.

**Keep overrides in the code, never as an edit to a generated table.** A
hand-edit to a file produced by a script is lost the next time anyone runs
the script, and lost silently. An override belongs where the script will read
it, with the reason and what the generic rule would have chosen.

**A domain-split score that ignores intra-block confidence finds floppy
termini.** Maximising inter-block minus mean intra-block PAE looks like it
identifies two separately-placed domains and does not: a disordered
80-residue terminus has a high PAE to everything *including itself*, so it
wins. On four models here that version put one dimer's boundary 322 residues
from the truth with a higher score than another dimer's correct one, and
invented a split in a single-chain protein. Require both blocks to be
confidently placed within themselves and score by `inter - max(intra_a,
intra_b)`; then both dimers land on their chain boundary exactly. The score
still separates two-body from one-body cases by only 1.2 against 0.9, so
report it against known controls rather than as a verdict.

**Read `auth_*`, not `label_*`.** An mmCIF carries two numbering schemes.
`label_asym_id` is the structure-model one: each distinct chemical entity gets
its own asym id, so a deposited protein is `A` and *each* of its cofactors,
metals and waters becomes `B`, `C`, `D`... `auth_*` carries the depositor's
chain letters and residue numbers -- the ones in the paper, in ChannelsDB and
in the `--chain A` a user will type. Reading `label_*` therefore breaks the
pipeline in a way that looks like biology: `prep --chain A --keep HEM` keeps
the protein, silently drops the haem because it is `label_asym` B, reports
"ligands kept: none", and the channel then measured runs straight through the
empty cofactor site. Residue numbers diverge too -- `label_seq_id` counts
from 1 along the entity and is null for anything non-polymer, so a residue
named in a paper is a different number here, or none at all.

**A `label_seq_id` of `.` is normal, not corrupt.** It is null for every
HETATM. A reader that calls `int()` on it crashes on the first ligand, which
is also the first file where it matters.

**A quoted value in `atom_site` shifts every later column.** Atom names like
`"O5'"` are quoted, and mmCIF ends a quoted value at quote-then-whitespace,
not at the first quote. A `str.split()` reader works until it meets one, and
then reads coordinates out of the wrong fields -- a wrong shape, not an error.

**A ligand in a modelled structure is not in its protein's chain.** AlphaFill
gives every transplant a chain of its own: one model in the iMARS set has the
protein in `A` and ligands in `B`..`Z` and `AA`..`AW`. Selecting by chain
alone keeps the protein and loses the ligand whose site is the whole point, so
`chbridge prep` adopts a `--keep` ligand from any chain when it comes within
`--ligand-contact` of the kept one, and reports every copy it saw with its
contact count.

**Several copies of one ligand are the normal case there, and each spare one
is a wall.** AlphaFill transplants from every donor it finds: four copies of
GPP in one model, with 40, 35, 18 and 3 protein contacts -- the last sitting
half on the surface. `--ligand-copies best` keeps the most-contacted one.
The default keeps them all, because a protein with two real cofactor sites is
also the normal case, and which one you have is not something the file says.

**Renaming the chain collides the ligand with residue 1.** AlphaFill numbers
every transplant `1`, so moving it into the protein's chain with `--relabel`
gives it the address of the protein's first residue. Downstream a residue is
selected by chain and number, so this is not cosmetic; prep renumbers past the
last protein residue and says so in both the report and the PDB header.

**mmCIF has no column widths and PDB does.** A two-character chain id, a
residue number past 9999, a five-character atom name or a serial past 99999
does not fail to write -- it writes into the neighbouring field, and is then
parsed without complaint as something else. Checked before writing.

## PyMOL (the prep and patch stages)

**`iterate` has no builtins** in its namespace, so an expression that works in
a normal interpreter can fail there. Supply the names you use explicitly (a
`_ns`/`_iterate`/`_iterate_state`/`_alter` shim) and check an expression before
relying on it.

**`cmd.centerofmass()` fails on pseudoatoms** (element PS). Use
`cmd.get_coords().mean(axis=0)`.

**`index` is renumbered after every `remove`.** Build the whole selection
first.

**PyMOL segfaults on interpreter shutdown** in headless mode (seen on 3.1),
*after* a script has finished and printed its result. A driver that
reads the exit code will see 139 and call a passing self-test a failure. Judge
by the script's own output marker, and use line-buffered output
(`stdbuf -oL`) or the marker is lost in the buffer when the process dies.

**Self-tests that write scratch folders into the current directory** leave
litter in the work folder. Give them a temp directory.

## RFdiffusion3

**Ligands must be off polymer chains,** one ligand residue per chain, named by
index token (e.g. `R315`).

**A declared ligand is something to pack the binder against.** Declaring a
ligand that merely happens to be nearby distorts the design.

**Occupancy-0.00 atoms crash it** (`could not broadcast (23,3) into (27,3)`).
Note that occupancy 0.00 can be genuine: 22 of A2P's 27 atoms are 0.00 in
1FNC itself, so the fix is to drop the ligand from the declaration, not to
raise its occupancy.

**Do not set `allow_ligand_on_existing_chain: true`.** It disables both
ligand checks and renumbers globally.

**Context is quadratic.** Trimming context chains to residues within 20 A of
the target costs about 0.25x.
