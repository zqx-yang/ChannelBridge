# Pitfalls

Things that cost days on the P5 prototype. Each one is either guarded by code
in this package or, where it is outside the package, stated here so it is not
rediscovered.

## CAVER

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
absent. The P5 `run_caver_all.sh` did, under `set -euo pipefail`.

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

**`iterate` has no builtins** in its namespace; the prototype scripts carry
`_ns`/`_iterate`/`_iterate_state`/`_alter` shims and a
`check_expression_portability()` check.

**`cmd.centerofmass()` fails on pseudoatoms** (element PS). Use
`cmd.get_coords().mean(axis=0)`.

**`index` is renumbered after every `remove`.** Build the whole selection
first.

**PyMOL 3.1 segfaults on interpreter shutdown** in headless mode on this
cluster, *after* a script has finished and printed its result. A driver that
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

Our cluster's CAVER build emits sphere chains at **~0.17 Å spacing** (185
spheres over 31.7 Å on a typical cluster). CaverDock's own example tunnel is
13 spheres at **0.956 Å** — about six times coarser.

Feed the raw chain to CaverDock's `discretizer` and it produces **2 disks from
95 spheres**. The job then **completes normally and prints plausible energies**
(−8.2 / −7.5 kcal/mol in the case we hit), because the job wrapper only checks
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

## Output-file age is not a liveness signal for a wrapped job

The house stall rule is "log idle 30 min **and** no new output → investigate".
Applied to CaverDock it is wrong in both directions, and following it would
have meant cancelling six healthy multi-hour jobs.

* **CaverDock writes nothing during the smoothed upper-bound phase**, which is
  ~86 % of a run (one measured run: LB 11 min, UB 1 h 37 min). Output can sit
  untouched for a long time while the job is perfectly healthy.
* **The job wrapper rsyncs results every 30 s**, so the output directory's
  mtime tracks the *syncer*, not the tool. It can look seconds-fresh while the
  tool is producing nothing at all.

So the directory's timestamp is decoupled from progress in both directions.
**Judge liveness from the mirrored log only**, which shows the search
advancing (`cScheduler::step(): conformations computed: …`).

This generalises to any tool run behind a wrapper that syncs on a timer: the
sync makes "recent output" meaningless, and a long internal phase makes
"no recent output" meaningless.

## Cost is driven by the search, not by the size of the input

Measured over six CaverDock runs: 5:41, 9:33, 10:57, 19:37, 47:50 and
**1:48:14**. The 1 h 48 outlier was the **shortest** tunnel in the set
(16.3 Å, 57 disks, a 20-heavy-atom ligand). Neither tunnel length, disk count
nor ligand size predicts runtime; the upper-bound search does.

Quoting a single "~20 min per run" from the median hides a 19× spread. Report
the distribution, and size a batch from the tail rather than the middle.

`--final_state LB` stops after the lower bound and is ~8× cheaper, but the
lower bound optimises each disk independently, so **its maximum is not a
barrier the ligand must continuously cross** — which is the quantity one would
be ranking tunnels by. Reasonable as a cheap pre-filter over many tunnels;
not a substitute for the upper bound in a result.
