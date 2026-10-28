# CaverDock ranking (stage 2c): measured, and not used

ChannelBridge selects a tunnel by CAVER's own geometry — priority, then
bottleneck. The obvious improvement is to rank the candidate tunnels by how
hard the substrate actually finds them: dock the ligand along each one with
[CaverDock](https://loschmidt.chemi.muni.cz/caverdock/) and prefer the
tunnel with the lowest energy barrier. That is stage 2c. It was **measured
and rejected**, and this page is why.

The short version: on the only channels the metric is measured on, CaverDock
did not return an answer at all, and on the channels where it did, the
answer moved further when we changed nothing than it did between the
tunnels it was supposed to tell apart.

Measured on three enzyme domains, CaverDock 1.2 with AutoDock Vina, over
about 3,700 CPU-hours. A small set — so read the findings as "this did not
work here", not as a general claim about CaverDock.

## What CaverDock computes

Two profiles per tunnel, both as energy against position along it:

| | how | cost |
|---|---|---|
| **lower bound (LB)** | each disc optimised independently | minutes to an hour |
| **upper bound (UB)** | one smoothed, continuous trajectory | hours to days |

The barrier is `E_a = max(E) − E(first disc)`, read from the
`REMARK CAVERDOCK RESULT:` lines of the output PDBQT. The LB is a bound, not
a trajectory: it can pass through states no single conformation reaches. The
UB is the physical one, and it is the expensive one.

## Finding 1 — no upper bound on the upstream channels, ever

The DO metric's O_IC term is measured on the **upstream** enzyme's channel.
That is the side that has to work.

Nine upstream runs were submitted (three domains × three CAVER clusters),
each with the 4-coumaroyl-CoA ligand, 39 atoms, `TORSDOF 12`. **Zero
produced an upper bound.** Every run directory holds a `-lb.pdbqt` and no
`-ub.pdbqt`, and every `-ranges.dat` is `3.40282e+38 / -3.40282e+38 / 0` —
the uninitialised sentinel, one row per disc.

Every one of the nine ran for days and then ended without an upper bound:
five were killed when they exhausted the memory they had been given, one
failed inside CaverDock itself (`CaverDock execution failed: trajectory was
not smoothed`) after spending 93 % of its wall clock in the UB phase, and the
rest were still running when this was written. More memory would plausibly
let some finish. That does not change the conclusion below, because of
finding 2.

## Finding 2 — the UB barrier is seed noise

Five runs of one tunnel, byte-identical receptor, ligand and tunnel, differing
only in `--seed`:

| run | UB `E_a` | LB `E_a` |
|---|---|---|
| seed 1 | 2.0 | 2.9 |
| seed 2 | 3.6 | 2.9 |
| seed 3 | 3.3 | 2.9 |
| seed 4 | 3.3 | 2.9 |
| seed 5 | 4.8 | 2.9 |

UB `E_a`: **sd 1.00, range 2.80 kcal/mol**. `E_max` ranges 4.10. The LB is
bit-identical across all five — it is seed-deterministic, which is worth
knowing and is not a substitute, because the LB is the bound and not the
trajectory.

Against that, the quantity the ranking has to resolve — the gap between the
best and second-best cluster of a domain — is **0.5 kcal/mol** for one domain
and **0.7** for the other. The noise is roughly 4× the signal, and about 6×
the 0.5 kcal/mol threshold the analysis had declared for calling two tunnels
tied.

So a UB ranking of these tunnels reports the seed. Five seeds of one tunnel
would have ranked it first, third and last.

## Finding 3 — LB and CAVER's own ranking disagree, with nothing to break the tie

On the two domains with a valid downstream set, the LB barrier picks CAVER's
cluster 1 in one of the two; in the other it prefers cluster 2 by
0.2 kcal/mol. With no reliable UB there is nothing to adjudicate that, and
0.2 is inside any honest error bar on the LB.

## Check `TORSDOF` before you read any barrier

A PDBQT with `TORSDOF 0` is a **rigid** ligand. It is well-formed, Vina
accepts it, the run completes in minutes and the profile looks smooth and
plausible. Nothing anywhere says the ligand could not rotate a single bond.
We lost six runs to this.

```sh
grep TORSDOF ligand.pdbqt      # 0 means rigid: the profile is meaningless
```

## What ChannelBridge does instead

Tunnel selection stays on CAVER's priority and bottleneck
(`channelbridge.caver.select`), and the tunnel geometry is normalised
(`chbridge tunnel`) so that the published D_IC and O_IC do not depend on the
CAVER build or on which database produced the chain. See
[caver-builds.md](caver-builds.md).

`channelbridge.caver.caverdock_note` still writes the per-selection note for
anyone who wants to run CaverDock themselves. If you do:

* budget the UB, not the LB — it is the whole cost, and it is days;
* ask for well over 8 GB per run, and measure rather than assume;
* run **several seeds of the same tunnel** before you compare two tunnels,
  and publish that spread beside any barrier you quote;
* check `TORSDOF`.

On this evidence CaverDock cannot rank these channels, and a number that
moves by 2.8 kcal/mol under a reseed should not be in a selection rule. It may well work on a smaller, more rigid ligand through a tighter tunnel.
It did not work here.

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
