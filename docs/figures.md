# Figures

Two products, built by two commands, from one pose search.

| | `chbridge figures` | `chbridge hotspots` |
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
| `render_monomer` / `render_multimer` | the illustration, with D_IC, O_IC and DO drawn | translucent (render) |
| `render_cs` | the same, cut open through the channel axis | cut copy solid, ghost at 85 % |
| `base_monomer` / `base_multimer` / `base_cs` | the same frames with the 3D channel and annotation **off** | translucent (render) |

Ligands appear in **every** render scene. Two bugs here are worth naming
because neither showed up as an error: the multimer object only ever had
`show surface` for polymer, so the assembly render silently lost every
cofactor; and the cross-section selected `and polymer`, excluding ligands
outright. The cut copies now come with an **uncut** `xsec_ligands` object —
the cofactor is what the section exists to reveal, and cutting it leaves a
meaningless fragment.

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

The cross-section cuts by **deleting atoms**, never with a clipping plane: a
clipping plane also cuts the tunnel spheres, so a thick slab leaves protein
over the channel and a thin one turns the spheres into hollow rings. A copy
with its near half removed is drawn solid while the original stays at high
transparency, which is what keeps the enzyme from looking broken into pieces.

---

## What is in an SVG

The proteins stay a raster — there is no honest vector form of a ray-traced
molecular surface — and everything the figure *claims* is vector with an id:

| id | what |
|---|---|
| `background` | transparent rect; raise `fill-opacity` for a solid colour |
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
pass the projected axis explicitly.

---

## The hotspot map

Shown **before** the DO sampling, so the epitopes can be judged and changed
before committing to a search over them ([pipeline.md](pipeline.md)). It
therefore takes the two prepared enzymes and their tunnels, not a pose:

```
chbridge hotspots --out map.pse \
    --up-pdb GDH.pdb --dn-pdb FNR.pdb \
    --tunnel-up GDH_cl1_norm.pdb --tunnel-dn FNR_cl2_norm.pdb \
    --up-name GDH --dn-name FNR --ligands NAP,BGC,FDA,A2P \
    --classes-up GDH_residue_classes.tsv --patches-up GDH_hotspot_patches.tsv \
    --classes-dn FNR_residue_classes.tsv --patches-dn FNR_hotspot_patches.tsv
```

The relative position is computed here — the smallest clash-free D_IC at
O_IC 180°, which is the highest-DO geometry the pair can reach. On the p6a
GDH–FNR pair: **D_IC 9.0 Å, O_IC 180°, DO 0.9314**, closest approach 3.38 Å.
Pass `--complex` instead to use a pose that already exists.

Four choices distinguish it from a pose session:

* **Both enzymes are painted alike.** The subject is the classification, and
  two surface colours for one class would read as two classes.
* **The patches are coloured in place, not built as separate surfaces.** A
  surface built from ten residues alone renders as a free-floating blob that
  z-fights with the parent surface — which is what made the first version of
  this session unreadable. Each patch is still a named selection, so it stays
  pickable.
* **The patches are shades of one hue, and labelled.** A rank is an ordering,
  not a set of unrelated categories. `patch_01` is palest and the last rank
  deepest by default (`best_is_darkest: true` flips it). Shades are mixed with
  white rather than scaled to black: a dark red and a dark grey are hard to
  tell apart on a shaded surface, a pale red still reads as red.
* **Every patch carries a minimum binder.** Not a sphere: the same chain of
  spheres every 2 Å that the pose session draws, spanning the pose stage's own
  `binder_gap` (6 Å) — four blobs — pointing out along the patch's *local*
  surface normal. So the figure shows how far and which way a binder would
  reach from each epitope, before any pose exists to constrain it.

Transparencies matter here and are not uniform: the surfaces are at
`inspect` (0.55), the tunnels at 0.10 because a translucent tunnel inside a
translucent surface disappears, and the envelopes at 0.70 because they sit in
front of things worth seeing.

Scenes: `map` (everything), `classes` (classification + tunnels), `patches`
(epitopes, labelled), `envelopes` (the binders).

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
