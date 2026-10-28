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
- `channelbridge.io.mmcif`: mmCIF reading, so the prep stage can convert.
  Reads the `auth_*` numbering rather than `label_*` -- the scheme a paper,
  ChannelsDB and a user's `--chain A` all mean. Reading `label_*` puts each
  cofactor in an asym id of its own, so a kept ligand is silently dropped and
  its site measures as empty space. Validated against the two example entries
  present in both formats: 7,129 and 26,081 atoms identical atom for atom,
  including altlocs, insertion codes and HETATM flags. `pdb_column_problems`
  checks what a fixed-column PDB cannot hold before writing it, because a
  two-letter chain id does not fail to write -- it corrupts the field beside
  it and is then parsed without complaint.
- `examples/iMARS/find_structures.py` checks two things beyond protein
  identity, each of which changed an answer: **coverage** (an entry can hold
  the protein and not all of it -- 4CL's only experimental structure is a
  fusion missing its last 102 residues, where the C-terminal domain closes
  the active site, so a tunnel there runs out through a hole), measured over
  the aligned span so a cleaved transit peptide does not read as missing; and
  **which entry** among several of the same protein, by bound cargo, then
  wild type, then resolution -- the first cross-reference would give BAS the
  apo 3A5Q over the holo wild-type 3A5R, or the I207L/L208F mutant 3A5S.
  AlphaFill now also has to have a *usable* best-identity-tier transplant
  (<= 2.0 clashes per heavy atom, one definition shared with
  `pick_ligands.py`) to outrank a plain model, since that preference exists
  only because a ligand-derived start point survives. `structures.tsv` gains
  a `reason` column, and overrides live in the script, never as an edit to
  the generated table.
- `examples/iMARS/verify_4cl.py`: checks the one predicted-by-override model
  against the only experimental structure of it, in two parts that are not
  merged -- a residue-for-residue superposition on the fusion's observed
  N-terminal domain, and the prediction's own inter-domain PAE for the part
  the fusion never modelled. The two-body split is taken from the PAE matrix
  and validated on the four models whose answer is known: both dimers land on
  their chain boundary exactly, both single-chain enzymes score below them.
- `chbridge prep` adopts a `--keep` ligand from another chain when it comes
  within `--ligand-contact` of the kept one, with `--ligand-copies best` for
  the transplanted case. A modelled structure does not put a ligand in its
  protein's chain: AlphaFill gives every transplant its own chain id (one
  model in the iMARS set has the protein in A and ligands in B..Z, AA..AW),
  and one has four copies of GPP with 40, 35, 18 and 3 protein contacts.
  Every copy and its fate is reported, and a ligand moved into the protein's
  chain by `--relabel` is renumbered past the last protein residue rather
  than left colliding with residue 1.
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
- `environment.yml` and a dev container.
- **Stage 1**: `channelbridge.prep` and `chbridge prep` -- chains, ligands,
  waters and altlocs, with the guarantee that coordinates are never moved
  (verified to 0.000000 A), and `--relabel` for the common case of a pair
  built from two single-chain entries.
- **Stage 2a**: `caver.find_caver` / `caver.run_caver` and `chbridge caver
  run`, plus `starting_point_residue` / `starting_point_atom` support and
  `caver.cavity_start_point` (`--start-cavity`), which derives a startable
  point from residues that name a site.
- `examples/iMARS/`: 19 enzyme pairs from the iMARS paper over 24 enzymes,
  exercising both tunnel branches end to end (stages 0-7). Decisions are
  committed (six decision tables, 19 pair files, nine scripts and the 19
  published figures); the structures are fetched, not bundled, and
  `report.tsv` is an output rather than a decision.
- **Stage 5**: `channelbridge.pose` and `chbridge pose` -- the pose search,
  which is the method. Scans the (D_IC, O_IC, azimuth, tilt, spin) grid
  best-DO-first and applies four filters in increasing cost: clash, binder
  gap, a bridgeable epitope pair with a clear line between the two patches,
  and -- where the upstream enzyme is a homomultimer -- whether the pose can
  be realised on every protomer at once. Writes `poses.tsv`,
  `patch_pairs.tsv`, `patch_pair_best.tsv` and a complex PDB per chosen pose.
  Pure numpy: unlike the prototype it needs no PyMOL, reading instead what
  stages 1-3 wrote.
  - Poses are **built along the O_IC arm**, the same two points the angle is
    measured on, so requested and realised agree to 1e-14 A in D_IC and 1e-6
    degrees in O_IC on real structures, and `O_IC_rev` is exactly
    `180 - tilt` at the grid's tilt values. Building along `fitted_axis`
    instead would leave the grid asking for one O_IC and recording another,
    9-12 degrees away.
  - `pose.PoseParams.tilt_azimuth`: the third, previously unsampled degree of
    freedom of the downstream enzyme's orientation. The prototype's grid
    tilted the downstream channel in one fixed plane, so its downstream
    channel directions lay on an arc rather than over the sphere. The axis
    now exists and is documented; it defaults to a single value, which keeps
    the default grid the size it was.
  - Only feasible poses are issued a `geomNNN` id. Numbering rejections
    alongside successes once made `geom050145` the best design in a
    100,000-point grid.
  - Protomer copies for the occupancy test are superposed here (Kabsch over
    matching CA atoms) rather than by `cmd.super`, and are taken from the
    structure's own protomers rather than an assumed symmetry operator.
  - Point-to-segment distances are solved analytically instead of by
    sampling the segment every angstrom, which is both faster and strictly
    more correct: a 1 A sample can step over a thin obstruction.
- `Match.db_coverage` and `Match.recovery()`: the other direction of the
  tunnel comparison, for asking whether OUR chain found THEIR route when
  ours is deliberately longer. `Match.verdict()` is exit-led and answers a
  different question; using it here called six recovered routes "different
  tunnel", one of them an exact match over 100 % of the database chain.
- `channelbridge.generators` (stage 7): generator input from a finished pose
  search. Writes RFdiffusion3 JSON, RFD4-Proteina YAML + mmCIF and BindCraft
  JSON; writes the job-script argument lines for flexcraft and RFdiffusion,
  which take no spec file; and REFUSES protpardelle-1c and RFDpoly with the
  reason, rather than emitting a spec neither tool can honour. Exposed as
  `chbridge generate`, with `--list-tools`.
  The two targets are conditioned with **`C_CRD`**, per atom, on both chains.
  `C_CRDGRP` is the obvious-looking choice and is wrong here: it preserves a
  group's internal geometry and leaves its global pose free, so with two
  targets it would release the inter-enzyme placement -- which IS the output
  of stages 2-5 -- and discard D_IC and O_IC without failing. Every generated
  file carries that sentence.
  Hotspot residue numbers and the contig are renumbered **together** for a
  tool that restarts numbering per chain, because they were not once on
  another campaign and a hotspot set 184 residues off still produced designs,
  against the wrong epitope.
- `channelbridge.caver.read_log` / `Loaded`: what CAVER actually put in its
  Voronoi diagram, from its own `log.txt`. The gap between the file and the
  measurement is silent in both directions -- `exclude_residue_names` removes
  a residue without a word, a starting point inside a residue removes that
  one too, and a clashing AlphaFill transplant is dropped by CAVER itself --
  so the count is now read and bounded rather than assumed.
- `docs/caverdock.md`: why stage 2c is not built. CaverDock produced no upper
  bound at all on nine upstream runs (five killed OUT_OF_MEMORY at 2-3 days,
  one tool failure, ~3,700 CPU-hours), and on the channels where it did, five
  reseeds of ONE tunnel spread the barrier by 2.8 kcal/mol against a 0.5-0.7
  margin between the tunnels it was meant to rank.
- `chbridge session` writes the pose's geometry beside the `.pse`
  (`<session>_geometry.json`: D_IC, O_IC, DO, the projected O_IC axis and the
  frame) and then CHECKS itself against the pose it was built from, refusing a
  session whose own D_IC is more than 1 A from the pose's. `chbridge render`
  reads that file when `--dic/--oic/--do` are not given.
- `run_examples.py figures`: stage 6 for the iMARS set, one session and one
  editable SVG per pair, with the per-pair arguments that make the figure
  describe its own pose.
- `chbridge generate --mode`: **monospecific** (one job per enzyme, the
  other absent) and **context** (one job per enzyme with the other present
  but not targeted, cropped to the residues within `--context-cut` of the
  epitope) beside the default **bispecific**. The three are a ladder a
  campaign can walk up. `generators.crop_runs` does the trimming BindCraft
  asks for as a measurement; `generators.memory_note` reports the job's size
  against BindCraft's own guidance (about 550 target+binder residues on a
  32 GB card, 950 on an 80 GB one). On the iMARS set this takes a pair from
  1,128 residues to 620.
- `patches.pick_hotspots`: the five residues of a patch a generator should
  be told about -- the exposed F/Y/W/I/L/M, spread along the sequence so they
  span the epitope rather than one turn of a helix. The worked example's
  rule, and BindCraft's advice for choosing a site.

### Changed
- Internal project codenames were scrubbed from `src/`, `tests/` and the
  published docs, and `CONTRIBUTING.md` no longer has to warn contributors
  that part of the source is private shorthand.
- `make_af3_inputs.py` no longer defaults to one cluster's paths: `--job`
  requires `$CHANNELBRIDGE_AF3_RUN`, the scheduler settings are
  overridable, and the emitted script says which lines are site-specific.
- **Figure colours (user's defaults, 2026-10-27; every role still
  overridable).** The D_IC/O_IC annotation — the line, the exit dots, the
  O_IC arm, arrowhead and arc, and the segment names — is **magenta**; the
  channel stays red. They were one red, so the measurement was drawn in the
  colour of the thing it measures. The hotspot patch the pose bridges is
  **red** and now has its own `patch` role instead of borrowing
  `annotation`, and it is drawn in the figure rather than only in the
  inspection scene, so a pose render shows the epitope pair it was built on.
  The cofactor role moves from magenta to purple: the annotation took pure
  magenta and the two would otherwise have been the same colour, which the
  palette's own contrast check reported.
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
- `channelbridge.rfd3`, an empty placeholder, superseded by
  `channelbridge.generators`.
- The cross-section scene and its options. It was not meaningfully different
  from `render_pair` and its cut copy mis-coloured a ligand.

### Fixed
- **A chain with unmodelled gaps got the wrong hotspot numbers.** The
  renumbering offset was computed from the span, so a tool that renumbers
  the PRESERVED residues 1..N counted the missing ones as present: measured
  on a 1-100 chain missing 50-60, a hotspot at 70 was emitted as 70 where
  the tool calls it 59. The same class of error as the 184-residue offset,
  reached from the other direction. `Target.output_index` is now the one
  place a residue number is converted, correct for a dense chain, a gapped
  chain and a cropped chain alike, and `offset` refuses the cases it cannot
  express.
- **A monospecific job's file still contained the enzyme the mode removes,
  and a cropped job's file was not cropped for any tool but one.** The crop
  reached the mmCIF and the contig and nothing else -- so BindCraft, the
  tool whose memory guidance motivated cropping, was handed the whole
  partner while `memory_note` described a file that was never written. One
  filter (`job_atoms`) now decides what a job's input contains, and every
  tool that takes a structure path is given a written-out `_target.pdb` that
  matches its job.
- **Every generated file carried the two-chain frame warning, including
  single-chain ones.** The suppression existed and no call site passed the
  spec, so a monospecific file read as deliberately two-chain.
- A hotspot on Ala, Gly or Pro is refused with the residue named, instead of
  failing the whole RFdiffusion3 batch at load time.
- `design_jobs` refuses to return two jobs with one name; it used to be
  possible for two labels sharing a first word, and `write_for` would
  overwrite the first without a word.
- `pick_ligands.py --no-evidence` wrote `start_mode` from a clash count it
  had not fetched, reading "no evidence" as "zero clashes" -- which silently
  turned the three worst transplants in the set (2.5-5.2 clashes per atom)
  from `cavity` into `ligand`. It now writes `unknown` and `build_inputs.py`
  refuses that row.
- `run_examples.py fetch` built its download list from `enzymes.tsv`'s
  superseded `entry` column -- a fourth reader of the table nothing should
  read. It downloaded 20 entries of which 12 were read by nothing, while the
  6 that mattered (the site-transfer donors) arrived by accident. It now
  names what it needs.
- `make_af3_inputs.py`: `--collect` iterated `COPIES` rather than the table's
  `af3` rows, so it hunted an output folder that was never made and exited 1
  on a run where everything succeeded; `stem_name` took the folder name
  rather than the input name, so a timestamped rerun copied no confidence
  file and `verify_4cl.py` silently printed "PART 2 skipped"; and a missing
  oligomeric-state decision defaulted to a monomer instead of refusing.
- `derive_start_sites.py` truncated a site to 40 residues while the evidence
  string on the same row counted all of them, and its `route` column was a
  second copy of `family`. It now writes the route that found the site
  (`pks-triad`, `terpene-motifs`, `homolog-<id>`).
- `pairs.tsv`'s `branches` column disagreed with the computed answer in 19
  of 19 rows and was read by nothing; it is gone.
- **Whole patches were being handed over as hotspots.** A patch is 7-14
  residues and most of them are not something a binder can grip; the list
  therefore included glycines, alanines and prolines, which are exactly the
  three residue types RFdiffusion3's `TIP` keyword is undefined for -- so
  the generator carried its own per-residue tip-atom table purely to emit
  `CA` and `CB` for residues that should never have been hotspots.
  `pick_hotspots` chooses anchors, and `select_hotspots` is now `TIP` for
  every entry, which is RFdiffusion3's own table rather than a second copy
  of it. A patch with fewer than `--n-hotspots` anchors yields fewer and
  says so instead of padding.
- **Two complexes for one pose id, silently.** The pose stage encodes the
  geometry in the file NAME, so re-running a search over a folder that
  already held one left both — `geom001_DO0.874_...` beside
  `geom001_DO0.856_...`, same id, different pose, different structure.
  Nothing overwrote anything and nothing complained, and since `--pse top1`
  selects by DO the **superseded** pose won wherever it scored higher. On
  the iMARS set that was 36 stale complexes across 17 pairs, several
  computed on structures the set no longer uses; a design spec built from
  one of them was caught only by the generator's span check.
  `figures.scan_output` now refuses the collision and names both files, and
  `run_examples.py poses` clears the previous search before running.
- The PNG preview drew the channel outline in the ANNOTATION colour where
  the SVG uses the channel colour. Invisible while the two roles were the
  same red; the moment the annotation went magenta the preview showed a
  magenta channel the SVG does not contain.
- **`--pse-arg` was dropped in silence.** `cmd_session` and `cmd_figures`
  handed an `action="append"` list to a splitter written for a single
  comma-separated string, whose `if not value` guard made the empty list
  falsy. Every `--pse-arg` vanished, and the session builder then failed with
  "no chains given for the upstream enzyme" -- naming neither the option nor
  the reason. Splitting on commas would have been wrong too: these are
  `key=value` tokens whose values contain commas (`up_chains=B,C,D,E`), so
  each occurrence is now passed on verbatim.
- **A figure could be labelled `D_IC = 0.0 A, O_IC = 0.0 deg, DO = 0.000`.**
  Those flags default to 0.0, which is a real-looking value, and the session
  -- the only place the numbers exist -- printed them and threw them away. It
  now writes them down; `chbridge render` picks them up, an explicit flag
  still wins, and with no numbers from either source it says so instead of
  labelling zero.
- **A session could draw both channels in the wrong frame and say nothing.**
  Each tunnel is in the frame of the PREPPED enzyme, not of the pose, so it
  has to be carried over -- which needs `tunnel_up_ref`/`tunnel_dn_ref`.
  Without them the tunnels are drawn where the unposed enzymes had them and
  the session measures between those: on one pair, **D_IC 38.0 A and DO 0.538
  against the pose's own 6.0 A and 0.937**, with every step exiting 0.
  `chbridge session` now compares its own geometry with the pose's and
  refuses the mismatch, naming the likely cause.
- **The SVG step cached its helper renders on existence.** The mask, mouth
  markers, footprint and transparent base were reused whenever the file was
  there, so a corrected session rendered with the previous session's helpers:
  corrected D_IC, O_IC and DO labels over exit markers at byte-identical
  pixel coordinates from the old frame -- a picture that disagrees with its
  own caption. They are keyed on the session's mtime.
- **`build_session` crashed on a channelled ligand that matches nothing.**
  `cmd.get_coords` returns `None` for an empty selection and
  `np.asarray(None, float)` is a 0-d array -- not `None`, and `len()` of it
  raises `TypeError: len() of unsized object`. The "no atoms matched the
  channelled ligand" branch written for exactly this case was unreachable.
- **`examples/iMARS` tunnel start points, the same error in a second place.**
  `_seeded_caver_side` averaged the cargo over every copy, which
  `_caver_side` had already been fixed not to do. Taxadiene synthase's 3P5P
  holds two A3C copies 8.20 A apart in one chain; the pooled centroid sat
  4.10 A from each, in neither site, and CAVER ran from it without complaint
  and returned tunnels of the space between them. Both branches now share
  `_cargo_seed`, which picks one copy and says in its own sentence when that
  choice is arbitrary.
- **The CAVER cache could not see a corrected starting point.** Staleness was
  keyed on the structure's mtime only, so changing the start point -- the
  whole content of the fix above -- left the old clusters in place and the
  report printed the new start point beside the old tunnel. The cache would
  have hidden its own fix. It is now keyed on the structure AND the config's
  CONTENT, with the config stamped into the run directory after CAVER
  succeeds, and the fallback seed gets its own run directory so a cached
  answer from either seed stays cached.
- **The cargo was excluded from CAVER only when it also supplied the start
  point.** Four enzymes keep a cargo and start from lining residues instead,
  so they were measured with nothing but waters named -- the cargo asked to
  sit in the channel whose bottleneck was being measured. On this set it
  changed no geometry, and that is the point: the GST transplants are dropped
  by CAVER anyway and the A3C was removed as a starting-point residue, both
  accidents of this data. The cargo is now excluded by name in both branches
  and in every generated `pair.yaml`, and `_check_loaded` bounds what CAVER
  loaded so an exclusion that does not take is an error instead of a silence.
- **`report.tsv` and all 19 `pair.yaml` credited the wrong coordinates.** Both
  printed `entry` from `enzymes.tsv`, the homolog of the first structure
  search, for all 24 enzymes -- so four results were attributed to PDB entries
  no stage in the run opened (4CL to 5BSM, really `af3_4CL.cif`; STS to 3A5S;
  LIC to 1N1B; MYS to 5UV1). They now carry `structure` and `source` from
  `inputs.tsv`, with `superseded_entry` keeping the history where it cannot be
  mistaken for provenance. `branch_changed` was computed and written to no
  column at all; it is now in the header. A pair file also no longer claims
  "ChannelsDB for both" for a pair whose upstream enzyme became a model with
  no annotation to fetch, and no longer names a `channel_ligand` that is not
  in the file.
- **`run_examples.py check` reported duplication that no longer existed.** It
  grouped on `enzymes.tsv`'s `(entry, chain)` and so answered for the first
  version of the set, printing four shared structures and "16 distinct
  computations" after `inputs.tsv` had given every one of those enzymes a
  model of its own. The answer is 19 of 19. It now keys on the whole
  `inputs.tsv` computation -- file, chains, cargo, kept ions, trim, start mode
  and start residues -- so two enzymes reading one file are only one
  computation when every decision about it also matches. This fed
  `duplicate_pairs`, which had been dropping three real pairs from the
  comparison against the paper.
- **`examples/iMARS/find_structures.py` crashed in its documented default
  mode.** `decide` grew a third return value, `rows` grew with it, and two
  unpacks did not: `./find_structures.py --only AS` died with
  `too many values to unpack (expected 8)` before printing anything, so the
  counts summary and the "need an AlphaFold3 run" list had never executed.
- The hotspot map was cached on existence, not timestamps, so a map built
  from a superseded tunnel survived a corrected `tunnels` run and was
  reported as current -- the cache `_is_stale` exists to abolish, at a call
  site it had missed.
- `pairs.tsv`'s `branches` column is no longer read: it was written when the
  set first ran and went stale when `inputs.tsv` replaced 19 of the 24
  structures. The branch each side takes is computed from the merged table.
- `figures.scan_output` anchored the pose id at the start of the file name, so
  a folder written with `chbridge pose --prefix` -- which puts the prefix in
  front of everything -- reported **zero poses**. A folder full of complexes
  that reads as an empty search is the worst possible failure mode for this,
  since it is indistinguishable from a pose stage that found nothing.
- The epitope-pair tag in a pose's file name is now built from enzyme names
  stripped to letters and digits. The figure stage parses
  `<up>p<N>-<dn>p<M>` out of the name, and the iMARS set's `enz_4CL` style
  split the field: both patch ids were silently lost from every name.
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
See [Known limits](README.md#known-limits) for the method's limits, and the
stage table under [Layout](README.md#layout-and-what-is-built) for the state
of each stage. CaverDock ranking (stage 2c) was measured and rejected; see
[docs/caverdock.md](docs/caverdock.md). Every other stage runs through
`chbridge`.

The ChannelsDB branch applies only where an enzyme's own coordinates happen
to be an annotated PDB entry, which in `examples/iMARS` is 2 of 24. Where it
does apply, a fetched chain usually stops before the surface as this package
defines it, and D_IC and O_IC are both measured at the mouth -- so each
source is used for what it is good at: the database identifies the route,
CAVER produces a chain that reaches solvent. Neither of the two seeded runs
in the example set recovered the annotated route, so both tunnels are CAVER's
own choice with no category attached.
