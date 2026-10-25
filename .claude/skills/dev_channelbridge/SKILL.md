---
name: dev_channelbridge
description: "How to work on ChannelBridge (the channel-aware bispecific pose-search tool at programs/ChannelBridge/production, with P6_CBD/p6a as its worked example and p6z as the iMARS reproduction). Covers the layout and the seven pipeline stages, what is built and what is not, the house conventions (naming, prefixes, defaults, provenance), how to test code that needs PyMOL or a cluster, and the ten failure modes this project has actually hit -- each with the evidence that produced it. Use when porting a stage, adding a figure or a metric, changing a default, writing or reviewing anything in src/channelbridge, sizing a compute campaign, or reporting a number that came out of one."
---

# dev_channelbridge — working on this tool

ChannelBridge places two enzymes so a channelled intermediate can pass between
them, scores the geometry with the iMARS DO metric, and draws the result.
Package at `programs/ChannelBridge/production`, worked example at
`projects/P6_CBD/p6a`, iMARS reproduction at `projects/P6_CBD/p6z_iMARS-repro`.

Related skills: `pd_render` (the user's house style for complex figures — this
tool implements it), `E14-cluster_slurm` (running jobs), `pd_pipeline_orga`
(campaign layout and naming).

---

## 1. What this is, structurally

Seven stages, five optional, in `docs/pipeline.md`. Built: tunnel
normalisation (`tunnel.py`), the DO metric (`geometry.py`), residue classes
and hotspot patches (`patches.py`), the reference placement
(`reference_pose.py`), the sessions and SVGs (`render/`). **Not built:**
structure prep (`prep/`), running CAVER (`caver/`), the pose search
(`pose/`), generator input (`rfd3/`) — those four package directories are
empty stage homes, not dead code.

**The layering is the one architectural rule that matters:**

```
pure Python, pytest-testable          tunnel.py geometry.py patches.py
  (no PyMOL, no cluster, no I/O)      reference_pose.py render/{frame,cut,svg,style,select}.py
          ▲
          │ called by
          │
thin shims that run INSIDE PyMOL      render/run_{pse,png,hotspots,patches}.py
  (measure things, draw things)       render/{scenes,hotspots,annotation}.py
```

Decisions live in the pure layer; measurements and drawing live in the shim.
When adding to a stage, ask which side a piece belongs on — if it is
arithmetic, it goes where a test can reach it without a renderer.

PyMOL code cannot be reached by pytest, so each shim carries its own
`--selftest`. **Both must pass before anything is called done:**

```bash
pytest -q && ruff check src tests
pymol -cq src/channelbridge/render/run_pse.py      -- --selftest
pymol -cq src/channelbridge/render/run_hotspots.py -- --selftest
pymol -cq src/channelbridge/render/run_patches.py  -- --selftest
```

Judge those by the printed `PASSED`/`FAILED`, not the exit code: PyMOL
segfaults during teardown after the test has already finished. Run with
`PYTHONUNBUFFERED=1` or the marker is lost when it does.

PyMOL gotchas that each cost a run: `-m` runs a FILE, not a module; a file
PyMOL `execfile()`s gets `__file__ = pymol/__init__.py`, so locate the
package from `sys.argv[0]`; and such a file never gets
`__name__ == "__main__"`, so an `if __name__` guard silently does nothing and
exits 0.

---

## 2. House conventions

* **`--prefix` on every stage**, at the front of every file it writes:
  `<prefix>_hotspot-map.pse`, `<prefix>_geom001_render.pse`,
  `<prefix>_geom001_render_{pair,assembly}.svg`. Pose ids keep their decimals
  (`geom001_DO0.870_d010.0_...`), so never split a name on `.`.
* **pair / assembly, never monomer / multimer** — the "monomer" scene holds
  two chains. Old spellings survive as CLI aliases only.
* **Upstream enzyme green, on the left**, in every frame and at every roll.
  Judged on P1 (which points into the upstream body), not on the exits,
  which sit together at the interface.
* **Defaults are a contract.** Every threshold lives in one dataclass
  (`PatchParams`, `RenderStyle`, `SvgStyle`, `ClassPalette`) with its value
  and the reason beside it. Changing one silently makes old and new outputs
  incomparable; change it deliberately and say so in `CHANGELOG.md`.
* **An unknown key is an error**, not a silently ignored line — a figure that
  quietly keeps a default you thought you changed is worse than one that
  refuses to build.
* **Scenes store object visibility and the CAMERA, never global settings.**
  `ray_trace_mode`, `float_labels`, `orthoscopic` must be re-applied per
  scene at render time. A scene that differs from another only by what is
  switched on is a toggle, not a scene — delete it.
* **Say which numbers a figure's reader can measure off it.** A rolled frame
  foreshortens O_IC; report `o_ic_true` and `o_ic_drawn` and pass the
  projected arm to the SVG as `axis_px`.

---

## 3. The ten failure modes this project has actually hit

Each of these cost real time or real CPU here. They are listed with their
evidence because the general form is only convincing with the instance.

**1. Read the reference implementation before porting.** Not "reason about
what it probably does". The classification axis is an SVD fit over the outer
spheres, not the last two spheres — the original says why (consecutive CAVER
spheres are ~1 Å apart, a pair is jitter; the two differ by 9–12°). SHADOW is
tested at *every* opening via `any()` while the table reports only the
primary. Neither is guessable. Four separate figure defects were found the
same way, after several rounds of adjusting colours against an imagined
target. Open the reference, read its settings back out, then write.

**2. A label nothing asserts on is only documentation.** The
`ligand_chemistry_route` column was correct the whole time and nobody checked
it; five CaverDock runs docked a molecule with 11 heavy atoms instead of 12,
no aromatic ring, no carboxylic acid and `TORSDOF 0`. Validate at the
boundary and **raise** — `validate_ligand()` compares canonical SMILES,
formula, heavy count and torsion count against the component reference, with
no warn-only path, and a regression test feeds the known-bad file back in and
requires the refusal.

**3. Measure a rate over the regime you will use it in.** Three budget
estimates were wrong in the same direction (970 → 3,400–5,300 → 12,365
CPU-hours) with one root cause: a quantity measured in one regime used to
price another. The mean of the finished runs (finished runs are a biased
sample *because* they finished); a per-disk rate taken over elapsed time that
includes a phase the disk counter does not cover; a rate from a rigid ligand
applied to a flexible one (6.7× out, measured on the same tunnel). State
which part of an estimate rests on measurement and which is extrapolation,
and give them separately.

**4. Verify a fix on the artefact, not in the source.** The mouth ring was
written to the SVG correctly and drawn in the tunnel colour on top of the
tunnel fill — present, correct and invisible. `float_labels 1` is the obvious
fix for labels hidden behind a surface and does nothing to a ray trace; the
re-render came back identical. Look at the picture, diff the bytes, read the
number back out of the session.

**5. Withhold a statistic the data cannot carry.** ρ = +0.655 at n = 6 is
p ≈ 0.16 — indistinguishable from no relationship *and* from a strong one.
A caveat still anchors a reader, so the summariser returns a refusal string
below n = 10 rather than a number. Report the denominator beside every
verdict, and add an **indistinguishable** category where a margin is inside
the method's own error (0.5 kcal/mol for CaverDock) — a tie is a result, and
collapsing ties into either camp overstates the evidence.

**6. Pre-register the interpretation.** Before running the seed-variation
test, write down what each outcome will mean — tighter threshold, looser
threshold, or "deterministic, which is an answer about determinism and not a
measurement of precision". Then the result cannot be read to suit.

**7. Fix the class, not the instance.** Three jobs died on a relative
`OUT_DIR`; that variable was fixed and two more died on a relative
`SCRIPT_DIR`. The generator now asserts *every* path-valued variable is
absolute. Same shape as #3.

**8. Check the invariant; do not reason about the sign.** `angle_frame`
builds x from P1−P2 (into the upstream protein) and needs a 180° flip;
`ic_line_frame` builds x from P3−P2 (away from it) and inherited the same
flip, shipping mirrored. Both now *test* where P1 lands and turn only if
needed. Likewise, measure clearance from the projected segment, not from a
camera axis — they coincide only when the line is on x, and the difference
flattered one frame family by 5×.

**9. Handle degenerate cases; do not avoid them.** The O_IC arc is a slerp
dividing by `sin(φ)`, guarded only at the low end. That was fine for years of
poses at 135–165° and collapsed the moment the hotspot map asked for its
reference placement at exactly 180°, where `sin(π) ≈ 1e-16` and every arc
point lands on the vertex. Guard both ends. Reusing one function in a second
place is what exposed it — a second copy would have hidden it.

**10. A known-bad input goes in an exclusion list, not in memory.** A void
run that had already been invalidated was included in a later spread
calculation and moved the headline from 19× to 66×. Put the id in the code.

---

## 4. Reporting numbers

This project's output is numbers other people will act on, so:

* **Reproduce a port against the reference before believing it.** The patches
  port was checked against the published p6a tables: geometry columns 0.00
  difference over 261 residues, patch table identical on rank/score/seed/
  n_res/anchors. Differences that remain must be *explained* — ours were
  exactly accounted for by inputs (no `--avoid` passed, CHANNEL from the
  geometric fallback rather than CAVER's list).
* **State what the metric cannot see.** O_IC is one-sided: it measures only
  whether the donor channel aims at the acceptor's exit, and is blind to
  which way the acceptor opens (90° on the p6a pose). The intermediate in
  transit is never modelled; the D_IC gap is empty. Both belong in any
  write-up that quotes DO.
* **Correct in public and carry the correction forward.** When a number
  changes, say which earlier statement it replaces. A quoted figure computed
  on the whole cluster does not survive a restriction to one node.

---

## 5. Cluster work

Follow `E14-cluster_slurm`. Specific to this tool:

* CaverDock needs a CAVER tunnel as **input** — it ranks tunnels, it cannot
  find them. "Use CaverDock instead of CAVER" is not a possible pipeline.
* Cost scales steeply with ligand size: ~0.0144·e^(0.196·T) h/disk in torsion
  count, so a 12-torsion ligand is ~100× a rigid one. Price the **expensive**
  class first; we deferred it and discovered its true cost 15 hours in.
* Pin nodes in the generator with a per-script assertion, never by hand.
* Submit **one** job and read its first real output lines before releasing
  the rest. This rule was skipped twice in one day and cost queue churn both
  times.
* A hold needs a documented owner and a documented release condition, and is
  released only on the user's word. Prefer replacing a bad job over releasing
  it: build the replacement, submit it, then cancel — never leave a window
  where the bad one could start.

---

## 6. Before saying something is done

- [ ] `pytest -q` and `ruff check src tests` clean
- [ ] every affected `--selftest` prints `PASSED`
- [ ] the change looked at in the rendered artefact, not just the source
- [ ] a new threshold has its value, its reason and its override path
- [ ] a ported stage reproduced against its reference, differences explained
- [ ] numbers quoted with their denominator and their regime
- [ ] **ask before committing** — the user does not want many commits
