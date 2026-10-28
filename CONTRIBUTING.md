# Contributing to ChannelBridge

The tool is under construction and published early so the method can be read
and criticised. Issues pointing at a wrong number are more valuable than
issues pointing at a missing feature.

## Layout, and the one architectural rule

```
pure Python, pytest-testable          tunnel.py  geometry.py  patches.py
  (no PyMOL, no cluster, no I/O)      channelsdb.py  caver/  reference_pose.py
                                      pose.py
                                      render/{frame,cut,svg,style,select}.py
          ▲
          │ called by
          │
needs PyMOL                           render/{scenes,hotspots,annotation}.py
  (measures things, draws things)     render/run_*.py  ← entry points run
                                                          INSIDE PyMOL
```

Decisions live in the pure layer; measurement and drawing live in the PyMOL
layer. When adding to a stage, ask which side a piece belongs on: **if it is
arithmetic, it goes where a test can reach it without a renderer.**

PyMOL code cannot be reached by pytest, so each entry point carries its own
`--selftest`. All of these must pass before a change is done:

```bash
pytest -q
ruff check src tests
pymol -cq src/channelbridge/render/run_pse.py      -- --selftest
pymol -cq src/channelbridge/render/run_hotspots.py -- --selftest
pymol -cq src/channelbridge/render/run_patches.py  -- --selftest
pymol -cq src/channelbridge/render/run_validate.py -- --selftest
```

**Judge a selftest by the printed `PASSED`/`FAILED`, not by the exit code.**
PyMOL segfaults during interpreter teardown after the test has already
finished. Run with `PYTHONUNBUFFERED=1` or the marker can be lost in the
buffer when it does.

Three PyMOL facts that each cost a run to discover, and are worked around in
every `run_*.py`:

* `pymol -m` is not Python's `-m`: it runs a **file**, so an entry point has
  to be a file.
* A file PyMOL `execfile()`s gets `__file__ = pymol/__init__.py`, not the
  script. The package is located from `sys.argv[0]` instead.
* Such a file never gets `__name__ == "__main__"`, so an `if __name__` guard
  silently does nothing and exits 0.

## Conventions

* **Defaults are a contract.** Every threshold lives in one dataclass
  (`PatchParams`, `RenderStyle`, `SvgStyle`, `ClassPalette`, `CaverParams`)
  with its value and the reason beside it. Changing one silently makes old
  and new outputs incomparable; change it deliberately and record it in
  `CHANGELOG.md`.
* **An unknown configuration key is an error**, never a silently ignored
  line. A figure that quietly keeps a default you thought you changed is
  worse than one that refuses to build.
* **`--prefix` on every stage**, at the front of every file it writes, so two
  runs over one folder cannot collide. Anything that *parses* such a name
  must therefore find its field inside the name rather than anchor at the
  start -- `figures.scan_output` anchored, and read a prefixed pose folder as
  an empty search.
* **`pair` / `assembly`, never `monomer` / `multimer`** — the "monomer" scene
  holds two chains. The old spellings survive as CLI aliases only.
* PyMOL **scenes store object visibility and the camera, never global
  settings**. `ray_trace_mode`, `float_labels` and `orthoscopic` must be
  re-applied per scene at render time.

## Things this project has got wrong, as rules

Each of these cost real time here. They are stated as rules because the
general form is what transfers.

1. **Read the reference implementation before porting.** Not "reason about
   what it probably does". The residue-classification axis is an SVD fit over
   a tunnel's outer spheres, not its last two — consecutive CAVER spheres sit
   ~1 Å apart, so a pair is placement jitter, and the two differ by 9–12°.
   Not guessable.
2. **A label nothing asserts on is only documentation.** A column recording
   which chemistry route produced a ligand was correct the whole time and
   nobody checked it; five docking runs used a molecule with 11 heavy atoms
   instead of 12, no aromatic ring, and zero torsions. Validate at the
   boundary and **raise**.
3. **Measure a rate over the regime you will use it in.** Several cost
   estimates were wrong in the same direction with one cause: a quantity
   measured in one regime used to price another. Finished runs are a biased
   sample *because* they finished.
4. **Verify a fix on the artefact, not in the source.** A marker was written
   to an SVG correctly and drawn in the same colour as what sat under it —
   present, correct and invisible. Look at the picture; read the number back
   out of the session.
5. **Withhold a statistic the data cannot carry.** ρ = +0.655 at n = 6 is
   p ≈ 0.16: indistinguishable from no relationship *and* from a strong one.
   Report the denominator beside every verdict.
6. **Check the invariant; do not reason about the sign.** Two camera frames
   build their x axis from different vectors, so "apply the 180° flip" means
   opposite things in the two. Both now *test* where the upstream enzyme
   lands and turn only if needed.
7. **Handle degenerate cases; do not avoid them.** An angle arc drawn as a
   slerp divides by `sin(φ)` and was guarded only near 0. At exactly 180° —
   which is the normal case for one of the two sessions — every arc point
   collapsed onto the vertex.

## Before opening a pull request

- [ ] `pytest -q` and `ruff check src tests` clean
- [ ] every affected `--selftest` prints `PASSED`
- [ ] the change looked at in the rendered artefact, not just the source
- [ ] a new threshold has its value, its reason and its override path
- [ ] numbers quoted with their denominator and the regime they were
      measured in
