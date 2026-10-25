# Figures

Two products from one pose search, and three commands: `chbridge session`
builds a session, `chbridge render` draws from one, and `chbridge figures`
does both in one go for the common case.

| | `chbridge session` + `render` | `chbridge hotspots` |
|---|---|---|
| asks | what does *this pose* look like | what did the search have to *choose from* |
| one per | pose | enzyme pair |
| enzymes | told apart by colour (upstream green, downstream orange) | painted **alike** |
| patches | only the one this pose bridges, in red | every one, as shades of a single hue |
| position | the pose | computed: smallest clash-free D_IC at O_IC 180 |
| runs | after the DO sampling | **before** it |

Figure generation is deliberately a **separate command from the pose search**.
Which poses deserve a figure is a question asked long after the search, often
several times with different answers, and re-running a search to change a
colour would be absurd.

Building and rendering are separate for the same reason one step down:

```
chbridge session <output-dir> --pse rank1 --prefix p6a      # -> p6a_geom001_render.pse
# open it, rotate to the frame you want, save it, then:
chbridge render p6a_geom001_render.pse --prefix p6a --variant both \
    --dic 17.0 --oic 165.0 --do 0.850
```

`render` takes **any** session — ours, or one you built yourself — so "turn it
until it looks right, then render that" is the normal path rather than
something to work around. A camera is a judgement, and no scan of clear-view
scores substitutes for someone looking at it.

---

## Naming

Every command takes `--prefix`, and it goes at the front of every file that
step writes. Two runs over the same poses then never overwrite each other, and
`ls` groups the output by run rather than by scene.

| file | what |
|---|---|
| `<prefix>_hotspot-map.pse` | the map. No pose in the name: it is built before any pose exists, and there is one per enzyme *pair* |
| `<prefix>_geom001_render.pse` | the pose session |
| `<prefix>_geom001_render_pair.svg` | the SVG of the two protomers |
| `<prefix>_geom001_render_assembly.svg` | the SVG of the full assemblies |
| `<prefix>_<full pose stem>_complex.pdb` | a structure copied out with `--pdb`. The whole stem, not `geom001`: this one leaves the project, so it carries its own geometry in its name |
| `intermediate/<prefix>_geom001_mask.png` | helper renders, kept so a re-render of one variant need not redo them |

Pose ids keep the decimals of the DO sampling
(`geom001_DO0.870_d010.0_o135_p3-p1_az000_ti030_sp270_ok`), so nothing in the
naming may split a name on `.`; the extension is taken from the end.

---

## Selecting poses

`--pse`, `--svg` and `--pdb` each take the same small language, and each is
independent, so "sessions for the top ten, an SVG for the best one, and the
structure of the one I want to show a colleague" is one command:

```
chbridge figures <output-dir> --pse top10 --svg rank1 --pdb geom007
```

| form | means |
|---|---|
| `none` | nothing |
| `all` | every pose |
| `5` or `top5` | the five best by DO |
| `rank3`, `#3` | the third best |
| `rank2-4` | the second to fourth |
| `geom007` | that pose, whatever its rank |
| `top3,geom019` | mixtures; the result comes back in rank order |

An SVG needs a session to render from, so asking for one adds the session
unless it is already on disk. The command says which it added and which it
reused. Ray-traces run **one at a time**: they are minutes to tens of minutes
each on a shared machine, and two at once is slower than two in a row.

---

## What is in a pose session

One `.pse` per camera frame, every variant an object *and* a scene. Scenes
share the frame's camera, so helper renders of one frame serve every SVG of
that frame.

| scene | shows | drawn |
|---|---|---|
| `occupancy` | the whole assembly, one binder envelope **per upstream protomer** | translucent (`inspect`) |
| `base` | proteins and ligands only | translucent (render) |
| `mask`, `mouth`, `footprint` | flat black helpers the SVG overlay is traced from | flat |
| `mask_assembly`, `mouth_assembly`, `footprint_assembly` | the same three stored with the **assembly** camera | flat |
| `render_pair` / `render_assembly` | the illustration, with D_IC, O_IC and DO drawn | translucent (render) |
| `base_pair` / `base_assembly` | the same frames with the 3D channel and annotation **off** | translucent (render) |

`pair` is one protomer of each enzyme; `assembly` is every chain of both. The
earlier names were `monomer` and `multimer`, which was wrong on the face of
it — the "monomer" scene holds two chains. Both old spellings are still
accepted on the command line.

The helper scenes exist **twice** because the two variants are framed
differently. The first assembly SVG traced its overlay from the pair-framed
helpers, which put the channel outline and the D_IC/O_IC lines off the
proteins entirely; the objects had not moved at all (chains identical to
0.000 Å) — only the camera the helpers were rendered with.

There is no cross-section scene. One was built and dropped: it was not
meaningfully different from `render_pair`, which already shows the tunnels and
the relative positions, and the cut copy coloured part of the downstream
ligand black. The cutting geometry stays in `render/cut.py`, tested on its own.

Ligands appear in **every** render scene. One bug here is worth naming because
it never showed up as an error: the assembly object only ever had
`show surface` for polymer, so the assembly render silently lost every
cofactor.

The `base_*` scenes exist for one reason: **an SVG is assembled on them, not
on the matching `render_*` scene.** The render scenes carry the channel as a
3D object and D_IC/O_IC/DO as PyMOL labels; the SVG adds a traced vector
channel and editable text. Build the SVG on a render raster and every number
appears twice, with the 3D channel's solvent-side smear showing under the
traced outline.

A pose session carries **no classification scene** by default. It answers
"what does this pose look like"; the residue classification belongs to the
hotspot map, which is built once per pair and shows every patch. Pass
`classes=1` to put it back.

Three things here were learned the hard way:

* **The drawn channel wraps the channelled ligand.** A CAVER sphere chain
  starts at the starting point and runs *outward*, so on its own it leaves
  the cofactor sitting outside the red volume — but the ligand is *in* the
  channel. The drawn volume is the chain unioned with that enzyme's
  channelled species wrapped at `ligand_wrap` (3 Å), as p5c does. Only the
  channelled species: wrapping every cofactor would paint red over prosthetic
  groups the channel never reaches (for GDH–FNR, NAP and A2P are wrapped;
  the prosthetic FAD is not).

  This is **visual only**. Two sets are kept internally — the CAVER chain
  alone feeds every number (P1, P2, P3, D_IC, O_IC, the mouth marker) and the
  wrapped set feeds only what is painted and the SVG mask. A selftest builds
  the session with and without the wrap and asserts D_IC and O_IC come out
  identical; if wrapping ever moved an exit, every measurement in the project
  would shift silently.

* **The inspection scenes are translucent, not solid.** They were solid at
  first, on the reasoning that an inspection scene should show the surface
  plainly. It hid the cofactors — and a hotspot patch only means anything
  next to the ligand the channel serves. They are now `transparency.inspect`
  (0.55), less see-through than the render scenes but not opaque.
* **`occupancy` draws one capsule per upstream protomer**, paired to the
  nearest downstream copy, because the only reason the scene exists is to
  show whether neighbouring binders collide. A single capsule cannot answer
  that. The scene reports the closest axis-to-axis distance; capsules of
  radius 11 Å touch below 22 Å. On the p6a GDH–FNR pose it reports 30.8 Å
  against the pose scan's own `binder_binder` of 30.9 Å — an independent
  cross-check of the search's clash test.

---

## What is in an SVG

The proteins stay a raster — there is no honest vector form of a ray-traced
molecular surface — and everything the figure *claims* is vector with an id:

| id | what |
|---|---|
| `background` | transparent rect; raise `fill-opacity` for a solid colour |
| `focus_outline` | the two protomers the geometry belongs to, picked out of an assembly; suppressed when they are the whole picture |
| `mouth_rings` | an exit whose outline has no gap to mark it, drawn as the disc really projects. **Off by default** (`mouth_ring_width`) |
| `proteins` | the embedded transparent render |
| `tunnels` | `tunnel_fill_<i>` (30 % red) and `tunnel_outline_<i>_<j>` |
| `channel` | `channel_line`, `exit_upstream`, `exit_downstream` |
| `angle` | `segment_X`, its arrowhead, `angle_arc` |
| `segment_names` | `label_X` (O_IC), `label_Y` (D_IC) |
| `labels` | `label_DIC`, `label_OIC`, `label_DO`, with real `tspan` subscripts |

Why the overlay is traced in 2D rather than rendered in 3D: `ray_trace_mode 1`
outlines the tunnel's bumpy surface as a field of dots and creases, and the
part of a channel that reaches into solvent ray-traces as a red smear over
nothing. Tracing the silhouette from the flat `mask` render and **clipping it
to the enzymes' own `footprint`** gives one clean curve per channel and
guarantees no red is drawn outside a protein. The outline is then cut open at
the `mouth` disc, so a channel reads as something that opens onto the surface
rather than as a closed cavity.

Label blocks are placed by scanning for the position that covers the least
protein, with every corner clear of the channel line. "Above" and "below" are
taken perpendicular to the line rather than in screen y, so a near-vertical
channel still has two sides.

The frame is cropped in the `viewBox`, not on the raster, so the overlay's
coordinates stay the render's coordinates and can be checked against the PNG.

### The one approximation

In the default `rotated0` frame the plane P1–P2–P3 is face-on, so the drawn
upstream-axis arm is derived from the reported O_IC rather than measured off
the render, and the assembler says so in its output. For any other camera,
pass the projected axis explicitly as `axis_px`.

### Two frames, and what each one costs

| | `angle_frame` (default) | `ic_line_frame` |
|---|---|---|
| on the x axis | the O_IC **arm** (P2→P1) | the **D_IC line** (P2→P3) |
| P1–P2–P3 plane | face-on | tilted by the roll |
| O_IC as drawn | exactly true | foreshortened, by a reported amount |
| free parameters | **none** | the roll about the line |

The angle frame is right for reading the angle, and it spends every degree of
freedom doing so. That is fine until a channel's drawn volume happens to lie
across the D_IC line — on the worked GDH–FNR pose it passes **0.97 Å** from
it, with two sphere centres inside 2 Å — because then there is nothing left
to adjust.

`ic_line_frame` keeps the D_IC line horizontal and rolls the camera about it.
The channels swing clear; the angle foreshortens. `clearest_roll` picks the
roll that maximises clearance **subject to a stated error budget**, so the
figure can never quietly trade its own measurement for a tidier picture:

| budget | roll | clearance | O_IC drawn (true 165.0°) |
|---|---|---|---|
| — (angle frame) | 0° | 0.97 Å | 165.0° |
| 1° | 20° | 2.25 Å | 165.9° |
| 2° | 210° | 2.79 Å | 166.9° |
| 5° | 225° | 3.44 Å | 169.3° |

Pass `oic_error_budget=1.0` to the session builder to use it. How much a roll
costs depends on the pose: near O_IC 180° the arm is almost antiparallel to
the line and hardly moves, near 90° it moves a great deal. So the roll is
chosen per pose against the budget, never fixed as a house constant — and
the session reports `roll_deg`, `clearance_A`, `o_ic_true` and `o_ic_drawn`
so a caption can say which frame it is.

Under a rolled frame the SVG **must** be given `axis_px`; the session reports
it. Letting the assembler fall back to the `rotated0` rule there would draw
the arm at an angle the picture does not contain.

---

## The hotspot map

Shown **before** the DO sampling, so the epitopes can be judged and changed
before committing to a search over them ([pipeline.md](pipeline.md)). It
therefore takes the two prepared enzymes and their tunnels, not a pose:

```
chbridge hotspots --prefix p6a \
    --up-pdb GDH.pdb --dn-pdb FNR.pdb \
    --tunnel-up GDH_cl1_norm.pdb --tunnel-dn FNR_cl2_norm.pdb \
    --up-name GDH --dn-name FNR --ligands NAP,BGC,FDA,A2P \
    --classes-up GDH_residue_classes.tsv --patches-up GDH_hotspot_patches.tsv \
    --classes-dn FNR_residue_classes.tsv --patches-dn FNR_hotspot_patches.tsv
```

With no `--out` that writes `p6a_hotspot-map.pse`.

The relative position is computed here — the smallest clash-free D_IC at
O_IC 180°, which is the highest-DO geometry the pair can reach. On the p6a
GDH–FNR pair: **D_IC 9.0 Å, O_IC 180°, DO 0.9314**, closest approach 3.38 Å.
Pass `--complex` instead to use a pose that already exists.

Four choices distinguish it from a pose session:

* **Both enzymes are painted alike.** The subject is the classification, and
  two surface colours for one class would read as two classes.
* **Every element is its own object**, so any of them can be switched off: one
  per residue class per enzyme, one per patch, plus ligands, mouths, tunnels
  and envelopes. The p5d original colours the classes in place on one surface
  and keeps 31 named *selections* beside it. A selection can be picked but not
  hidden, so with seven translucent classes on one object every colour is seen
  through every colour behind it and nothing can be isolated — which is what
  made the first versions of this session unreadable. The classes make the
  split clean because they **partition** the residues: the union of the
  per-class objects is the whole surface with nothing left over (checked on
  p6a — GDH 7 objects = 7,884 atoms, FNR 6 = 2,352, 0 unclassified). The
  patches are solid objects over the translucent base, never painted into
  it: a tinted base means switching a patch off changes nothing on screen.
* **The patches are shades of one hue.** Labels are built but start OFF:
  ten names on two crowded faces cannot all be placed beside the patch they
  name, and every patch is its own object (`<name>_patch_01`, and its
  envelope `<name>_patch_01_binder`), so the object list is the legend. A rank is an ordering,
  not a set of unrelated categories. `patch_01` is palest and the last rank
  deepest by default (`best_is_darkest: true` flips it). Shades are mixed with
  white rather than scaled to black: a dark red and a dark grey are hard to
  tell apart on a shaded surface, a pale red still reads as red.
* **Every patch carries a minimum binder, on every protomer.** The envelope
  object holds one blob chain per copy of the enzyme, so the assembly scene
  answers whether binders on neighbouring protomers would collide — which is
  what it is for. It used to place one, taken from a `resi -> coord` map
  where residue numbers collide across chains and one protomer silently won. Not a sphere: the same chain of
  spheres every 2 Å that the pose session draws, spanning the pose stage's own
  `binder_gap` (6 Å) — four blobs — pointing out along the patch's *local*
  surface normal. So the figure shows how far and which way a binder would
  reach from each epitope, before any pose exists to constrain it.

**The map is drawn with an orthographic camera**, unlike the pose session.
A map compares epitopes, so an equal patch must not look larger for being
nearer. It also makes the label handling work: see below.

**Labels collide in two different ways, and one fix cannot do both.**

* *In the camera plane.* Anchors are relaxed apart on the two screen axes,
  anisotropically -- a line of text is far wider than it is tall, so two
  labels at the same height need much more room than two stacked vertically.
  The relaxation is bounded by the drawn pair plus a margin, because without
  bounds it trades an overlap for a label out in the frame's margin, which is
  what happened to `FNR p05`.
* *In depth.* A label sits at its own anchor's depth, so one naming a patch on
  the far side of an enzyme renders **behind** the near surface and washes out
  to unreadable. The in-plane relaxation cannot fix that -- it deliberately
  never changes depth -- so the anchors are lifted along the view axis until
  they clear the frontmost atom. PyMOL's own `float_labels` looks like the
  answer and is not: it applies to the interactive viewport, not to
  ray-traced output. That was confirmed on a real render, where setting it
  changed nothing.

The orthographic camera is what makes the lift free: under it, moving a point
along the view axis does not move it on screen, so a label clears the surface
without sliding away from the patch it names. Under perspective it drifts,
and one label ended up over the *other* enzyme.

What remains: with six patches crowded on one face, separation can still put
a label in open space rather than beside its patch. A leader line is the
honest fix and is not implemented.

Transparency is **one** thing: a translucent whole-enzyme base per enzyme,
with every element — each residue class, each patch — a solid object on top
of it. That is what makes toggling work: switch a patch on and it reads as
"this, here, on this surface". The base is never tinted with the patch
colours; it was, so that the scene still read with the patches hidden, which
meant switching a patch off changed nothing on screen.

Patch surfaces use a slightly larger probe (1.7 Å against PyMOL's 1.4) so
they sit just outside the coincident class surface beneath them and always
win the depth test. A rendering radius, not a coordinate — nothing moves.

**Two scenes, and they hold the same things**: `classes_pair` — one
protomer of each enzyme — and `classes_assembly`, which is that pair
replicated onto every other protomer. Same classes, same patches, same
binder envelopes, same annotation; only the number of copies differs, which
is what makes them comparable. `classes_assembly` is stored last, so the
session opens on it.

The downstream copies have to be **generated**: a prepared downstream enzyme
is often a single chain while the upstream one is an oligomer. Each copy is
the same placement carried onto another protomer by that protomer's own
symmetry operator, fitted on matched CA atoms with the rms reported — the
same construction a pose search uses when it writes an occupancy model, and
a construction rather than a measurement.

There is no separate `patches` or `envelopes` scene. Those were this scene
with objects switched off, which is a toggle, not a scene; what a scene
stores that an object cannot is the camera.

Default class colours follow the p5d map script:

| class | colour | |
|---|---|---|
| `ALLOWED` | cyan | solvent-exposed and usable |
| `AVOID` | grey60 | excluded, and the footprint around it |
| `SHADOW` | pink | inside the cylinder a binder must not block |
| `CHANNEL` | magenta | the tunnel wall |
| `LIGAND_CORE` | purple | the residues holding the cofactor |
| `IFACE` | grey30 | buried by another chain of the assembly |
| `BURIED` | wheat | everything else |

Exposed hydrophobes are **not** split out of the surface colour here: the map
is about where a binder may sit, and that split belongs to the per-pose view.

## Changing any of it

Everything above is a default. The pair file's `render:` section carries
`colors`, `transparency`, `background`, `stick_radius`, an `svg:` block (fill
opacity, stroke widths, font size, arc radius, margins, which numbers to
print) and a `hotspots:` block (class colours, the patch hue and its gradient
ends, the envelope colour). An unknown key is an error rather than a silently
ignored line — a figure that quietly keeps a default you thought you changed
is worse than one that refuses to build. See `examples/gdh_fnr/pair.yaml`,
where every key is present and commented.
