# Changelog

All notable changes to this project are documented here.
This project follows [Semantic Versioning](https://semver.org/).

## [Unreleased]

Nothing is tagged yet, so everything below is the initial content of the
package and the changes made to it since the first commit.

### Added

Core:
- `channelbridge.geometry`: D_IC, O_IC and the DO score, with both published
  iMARS coefficient sets selectable and recorded per run.
- `channelbridge.tunnel`: CAVER sphere chains as geometry, with the
  normalisation that makes the exit, span and O_IC arm independent of which
  CAVER build produced them (`oriented`, `truncate_at_surface`, `resample`,
  `axis_point`), plus `axis_spread` as a quality flag for a curved mouth.
- `Tunnel.fitted_axis` / `channelbridge.tunnel.fitted_axis`: the mouth
  direction fitted over the outer spheres, which the residue classification
  is measured against. Deliberately NOT `Tunnel.axis`, the O_IC arm: the two
  differ by 9-12 degrees and only one of them is what O_IC is defined on.
- `channelbridge.config`: the pair file, replacing the prototype's ~40
  command-line flags and hard-coded chain ids. Unknown keys are errors.
- `channelbridge.io.pdb`: fixed-column PDB reading, CAVER's
  radius-in-B-factor convention, LF-only writing and CRLF detection.
- **Stage 3**: `channelbridge.patches` (residue classes and hotspot patch
  selection, pure and unit-tested) with `render/run_patches.py` as its PyMOL
  shim and `chbridge patches` as the command. Reproduces the reference
  implementation exactly on the worked example: geometry columns agree to
  0.00 over all 261 GDH residues, patch table identical on rank, score, seed,
  n_res and anchors.

Figures:
- `render/annotation.py`: the D_IC line, O_IC arm, arc and labels, shared by
  the pose session and the hotspot map instead of duplicated.
- `render/frame.arm_frame` and `clearest_roll`: roll the camera about the
  O_IC arm to clear the drawn channels off the D_IC line, against a stated
  budget on how far the drawn O_IC and D_IC may stray. Both are reported.
- `chbridge session` / `chbridge render` as separate steps, `--prefix` on
  every figure command, and `--axis-px` so a rolled frame's angle marker is
  measured rather than assumed.

CLI, docs and packaging:
- `chbridge` CLI: `validate`, `tunnel` (with `--compare`) and `do`.
- `docs/caver-builds.md`: the measurement showing two CAVER builds give the
  same tunnel with exits 8.9 A apart, and what normalisation does and does
  not fix. Includes the seed sweep establishing that the residual ~15 deg
  O_IC difference is systematic between builds (0.06 deg within a build over
  8 seeds), not stochastic.
- `docs/pitfalls.md`: the failure modes carried over from the prototype.
- `tests/test_crossbuild_regression.py`: asserts the cross-build claim
  against both real CAVER outputs; skipped when they are absent.
- `environment.yml`, a dev container, and `.claude/skills/dev_channelbridge`.
- **Stage 1**: `channelbridge.prep` and `chbridge prep` -- chains, ligands,
  waters and altlocs, with the guarantee that coordinates are never moved
  (verified to 0.000000 A), and `--relabel` for the common case of a pair
  built from two single-chain entries.
- **Stage 2a**: `caver.find_caver` / `caver.run_caver` and `chbridge caver
  run`, plus `starting_point_residue` / `starting_point_atom` support and
  `caver.cavity_start_point` (`--start-cavity`), which derives a startable
  point from residues that name a site.
- `examples/iMARS/`: 19 enzyme pairs from the iMARS paper over 24 enzymes,
  exercising both tunnel branches end to end (stages 0-4). Decisions are
  committed (`enzymes.tsv`, `pairs.tsv`, 19 pair files, `run_examples.py`);
  the structures are fetched, not bundled.
- `Match.db_coverage` and `Match.recovery()`: the other direction of the
  tunnel comparison, for asking whether OUR chain found THEIR route when
  ours is deliberately longer. `Match.verdict()` is exit-led and answers a
  different question; using it here called six recovered routes "different
  tunnel", one of them an exact match over 100 % of the database chain.

### Changed
- The hotspot map is built for toggling: a translucent per-enzyme base with
  every class and patch a solid object over it, and two scenes
  (`classes_pair`, `classes_assembly`) holding the same content at different
  multiplicity. Binder envelopes are placed on every protomer.
- Figure scene names are `pair` / `assembly`; `monomer` / `multimer` survive
  as CLI aliases only.
- `SvgStyle.mouth_ring_width` defaults to 0 (off).
- The published docs were split from the project's own post-mortems: rules a
  user acts on stayed in `README.md` and `docs/`, and the campaign numbers and
  debugging histories behind them moved out of the repository.

### Removed
- The cross-section scene and its options. It was not meaningfully different
  from `render_pair` and its cut copy mis-coloured a ligand.

### Fixed
- `run_caver` invoked a `CaverRunner` class that does not exist and passed no
  `-home`, so stage 2a could never have run. CAVER's entry point is its jar
  manifest (`caver.ui.Launcher`) and `-home` is required -- without it CAVER
  cannot read its own `bin/atom_radii.csv`. Found by running it.
- `_find_cluster` looked only in `data/clusters/`. CAVER writes
  `data/clusters_timeless/` for a single structure -- the only kind this
  pipeline runs -- and names the files `tun_cl_001_1.pdb`, so every static
  run looked like a run that had found nothing.
- `DbTunnel.length` read ChannelsDB's `Distance` series as cumulative arc
  length. It is, in the MOLE entries; in the CAVER-derived ones it is
  non-monotonic and saturates near 1, so every one of those tunnels was
  reported as about 1 A long -- a plausible number that nothing contradicted.
  Length is now measured from the coordinates, and
  `DbTunnel.distance_is_arc_length` reports whether the series is usable.
- `write_config` would write a config with no starting point at all. CAVER
  has no default and refuses to run, so the failure surfaced far from its
  cause; it is now refused where it is made.
- The O_IC angle arc is a slerp dividing by `sin(phi)` and was guarded only
  near 0. At exactly 180 degrees every arc point collapsed onto the vertex,
  drawing a spike that read as an arc. Both ends are guarded; at 180 the
  plane is chosen and the code says so.
- `ic_line_frame` put the upstream enzyme on the right. Both frames now test
  where P1 lands rather than applying a fixed flip.

### Known limits
See [Known limits](README.md#known-limits) in the README for the method's
limits. Not yet exposed through `chbridge`: the pose search itself, CaverDock
ranking and the generator-input stage; the table under
[Layout](README.md#layout-and-what-is-built) is the current state per stage.

On the ChannelsDB branch a fetched chain usually stops before it reaches the
surface as this package defines it: only 10 of the 31 database-branch enzyme
sides in `examples/iMARS` have one that does. D_IC and O_IC are both measured
at the mouth, so such a chain has no exit to measure. The resolution is to
use each source for what it is good at -- the database identifies the route,
CAVER produces a chain that reaches solvent -- by re-running CAVER seeded
from the database chain's buried end and checking the result back against it.
That takes the set from 15 usable mouths to 36 of 38. On 6 of the 13 seeded
enzymes CAVER recovers the annotated route; on the rest it does not, and the
tunnel is then CAVER's own choice with no category attached.
