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
its own surface definition, not this package's. Of the 31 database-branch
enzyme sides in `examples/iMARS`, only 10 have a candidate chain whose ends
fall below burial 15; on 5ICE every Procognate chain has a minimum burial of
29-71 atoms along its whole length. D_IC and O_IC are both measured at the
mouth, so check the verdict (`mouth` usable; `buried`, `recessed`, `surface`
not) before quoting geometry from a fetched tunnel.

**What to do about it: seed CAVER from it, do not measure it.** The database
identifies the route (its category says what the route is FOR, which CAVER
cannot) and CAVER produces a chain that reaches solvent (which the database
does not). Run CAVER from the database chain's buried end, then check the
result back against that chain. Two cautions from doing it over 13 enzymes:
seed from the chain END, not from a bound ligand -- seeding from the ligand
first dropped annotation agreement from 6 enzymes to 2, because CAVER then
finds whatever leaves the ligand site; and judge the check with
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

## PyMOL (the prep and pose stages)

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

## A CAVER tunnel is too dense to hand straight to CaverDock

A current CAVER build emits sphere chains at **~0.17 Å spacing** (185 spheres
over 31.7 Å on a typical cluster). CaverDock's own example tunnel is 13 spheres
at **0.956 Å** — about six times coarser.

Feed the raw chain to CaverDock's `discretizer` and it produces **2 disks from
95 spheres**. The job then **completes normally and prints plausible energies**
(−8.2 / −7.5 kcal/mol in one observed case), because the job wrapper only checks
that the `.dsd` file is non-empty. Nothing in the exit code, the log, or the
presence of output files reveals that the trajectory is meaningless.

This is the same root cause as the finding this whole tool exists to encode —
**sphere density is a property of the CAVER build, not of the protein** — and it
propagates into every downstream tool that consumes a tunnel.

**Re-space to ~1 Å before the discretizer, and assert the disk count** against
the tunnel's length. A 2-disk `.dsd` from a 16 Å tunnel is a hard error, not a
result.

### The preparation is not the same as for D_IC

Counter-intuitive, and easy to get backwards:

| | truncate at the surface? | re-space? |
|---|---|---|
| **D_IC / O_IC** | **yes** — the exit *is* the measurement | yes, to a fixed arc length |
| **CaverDock** | **no** — the ligand has to reach bulk solvent to leave | yes, to ~1 Å |

The same chain needs two incompatible preparations, and feeding either tool the
other's version produces output that looks entirely normal. Keep the raw chain
alongside whatever you derive from it.

### Crystal ligands often do not match their PDB chemical component

Meeko needs bond orders, hydrogens and formal charges; a ligand copied out of a
crystal structure has none of them. Matching the component's ideal SDF onto the
coordinates with RDKit's `AssignBondOrdersFromTemplate` supplies them — but only
if the atom counts agree. In 3A5R the HC4 ligand has **11 of the component's 12
heavy atoms** (one carboxyl oxygen unmodelled), the match fails, and chemistry
falls back to geometry-based bond inference. Docking energies depend on the atom
typing that follows, so record which route each ligand took and discount the
fallbacks rather than treating all the numbers alike.

## Ranking tunnels with CaverDock: cost is driven by the search

Runtime is set by the upper-bound search, not by the size of the input.
Measured over six runs on one project, elapsed time spanned **19x** (under
6 min to nearly 2 h) and the slowest run was on the **shortest** tunnel in the
set. Neither tunnel length, disk count nor ligand size predicted it.

So **do not price a batch from a median**. Report the distribution and size the
batch from the tail. Cost rises steeply with ligand flexibility rather than
with ligand size, so price the most flexible ligand in the set first.

`--final_state LB` stops after the lower bound and is ~8x cheaper, but the
lower bound optimises each disk independently, so **its maximum is not a
barrier the ligand must continuously cross** -- which is the quantity one would
be ranking tunnels by. Reasonable as a cheap pre-filter over many tunnels; not
a substitute for the upper bound in a result.

**Output-file age is not a liveness signal.** CaverDock writes nothing for the
whole smoothed upper-bound phase, which is the great majority of a run, so a
healthy job can look stalled for hours. If it runs behind a wrapper that
mirrors results on a timer, the output directory also looks fresh when the tool
is producing nothing. Judge progress from the tool's own log
(`cScheduler::step(): conformations computed: ...`), never from file
timestamps.
