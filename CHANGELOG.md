# Changelog

All notable changes to this project are documented here.
This project follows [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- **Stage 3**: `channelbridge.patches` (residue classes and hotspot patch
  selection, pure and unit-tested) with `render/run_patches.py` as its PyMOL
  shim and `chbridge patches` as the command. Reproduces the author's
  prototype exactly on the worked example: geometry columns agree to 0.00
  over all 261 GDH residues, patch table identical on rank, score, seed,
  n_res and anchors.
- `Tunnel.fitted_axis` / `channelbridge.tunnel.fitted_axis`: the mouth
  direction fitted over the outer spheres, which the residue classification
  is measured against. Deliberately NOT `Tunnel.axis`, the O_IC arm: the two
  differ by 9-12 degrees and only one of them is what O_IC is defined on.
- `render/annotation.py`: the D_IC line, O_IC arm, arc and labels, shared by
  the pose session and the hotspot map instead of duplicated.
- `render/frame.arm_frame` and `clearest_roll`: roll the camera about the
  O_IC arm to clear the drawn channels off the D_IC line, against a stated
  budget on how far the drawn O_IC and D_IC may stray. Both are reported.
- `chbridge session` / `chbridge render` as separate steps, `--prefix` on
  every figure command, and `--axis-px` so a rolled frame's angle marker is
  measured rather than assumed.
- `environment.yml`, a dev container, and `.claude/skills/dev_channelbridge`.

### Changed
- The hotspot map is built for toggling: a translucent per-enzyme base with
  every class and patch a solid object over it, and two scenes
  (`classes_pair`, `classes_assembly`) holding the same content at different
  multiplicity. Binder envelopes are placed on every protomer.
- Figure scene names are `pair` / `assembly`; `monomer` / `multimer` survive
  as CLI aliases only.
- `SvgStyle.mouth_ring_width` defaults to 0 (off).

### Removed
- The cross-section scene and its options. It was not meaningfully different
  from `render_pair` and its cut copy mis-coloured a ligand.

### Fixed
- The O_IC angle arc is a slerp dividing by `sin(phi)` and was guarded only
  near 0. At exactly 180 degrees every arc point collapsed onto the vertex,
  drawing a spike that read as an arc. Both ends are guarded; at 180 the
  plane is chosen and the code says so.
- `ic_line_frame` put the upstream enzyme on the right. Both frames now test
  where P1 lands rather than applying a fixed flip.
- `channelbridge.geometry`: D_IC, O_IC and the DO score, with both published
  iMARS coefficient sets selectable and recorded.
- `channelbridge.tunnel`: CAVER sphere chains as geometry, with the
  normalisation that makes the exit, span and O_IC arm independent of which
  CAVER build produced them (`oriented`, `truncate_at_surface`, `resample`,
  `axis_point`), plus `axis_spread` as a quality flag for a curved mouth.
- `channelbridge.config`: the pair file, replacing the prototype's ~40
  command-line flags and hard-coded chain ids. Unknown keys are errors.
- `channelbridge.io.pdb`: fixed-column PDB reading, CAVER's
  radius-in-B-factor convention, LF-only writing and CRLF detection.
- `chbridge` CLI: `validate`, `tunnel` (with `--compare`) and `do`.
- `docs/caver-builds.md`: the measurement showing two CAVER builds give the
  same tunnel with exits 8.9 A apart, and what normalisation does and does
  not fix. Includes the seed sweep establishing that the residual ~15 deg
  O_IC difference is systematic between builds (0.06 deg within a build over
  8 seeds), not stochastic.
- `docs/pitfalls.md`: the prototype's hard-won failure modes.
- `tests/test_crossbuild_regression.py`: asserts the cross-build claim
  against both real CAVER outputs; skipped when they are absent.

### Known limits
- The prep, pose-scan and RFdiffusion3 stages are not yet exposed through
  `chbridge`; they remain the prototype PyMOL scripts.
- Normalisation reconciles D_IC across CAVER builds but not O_IC, which
  retains a ~15 degree tunnel-specific uncertainty on a curved mouth.
