# Changelog

All notable changes to this project are documented here.
This project follows [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
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
