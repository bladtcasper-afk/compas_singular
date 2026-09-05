# What we changed in `compas_singular`

A complete inventory of the work done on this fork, measured against the upstream
repository it was forked from. Written as the basis for a presentation.

---

## 0. The baseline we are measuring against

| | |
|---|---|
| Upstream | `BlockResearchGroup` / `BRG-research` — `compas_singular` |
| Fork point | commit `66024561`, *"rename vertex_index as vertex_topo_index"*, Robin Oval, **13 April 2022** |
| Fork | `github.com/bladtcasper-afk/compas_singular` |
| Our first commit | `f102cd05`, *"Update of compas_singular to compas-v2"*, **5 August 2026** |
| Work period | 5 – 28 August 2026 |

Everything below is new since April 2022. Upstream had been dormant for four years and
**did not run at all** on a current COMPAS install — that is where the work started.

### Headline numbers

| | Upstream (2022) | This fork | Change |
|---|---|---|---|
| Library code (`src/`) | 9 229 lines | 26 877 lines | **×2.9** |
| Python modules in `src/` | 67 | 97 | +30 |
| Example scripts | 5 | 54 | **×11** |
| Tests | 1 (trivial import) | 3 files, 27 tests | — |
| Design / reference documents | 0 | 8 major `.md` | — |
| COMPAS version supported | 0.x / 1.x | **2.15.1** | port |

---

## 1. Layer one — making it run again: the COMPAS 2 port

*Commit `f102cd05`. This is the precondition for everything else.*

Upstream targets the COMPAS 0.x/1.x API. On COMPAS 2.x it fails at import time and, where
it does import, silently calls functions with the wrong signature.

**What was broken and what was done:**

| Problem | Fix |
|---|---|
| ~20 functions moved or were deleted from `compas.datastructures` / `compas.utilities` (`mesh_weld`, `trimesh_face_circle`, `network_polylines`, `mesh_explode`, `geometric_key`, …) | rewired to their COMPAS 2 homes (`compas.itertools`, `compas.tolerance.TOL`, `compas.datastructures.graph.operations`); the ones COMPAS 2 dropped outright were **re-implemented in `datastructures/mesh/operations.py`** (`trimesh_face_circle`, `mesh_weld`, `meshes_join`, `meshes_join_and_weld`, +171 lines) |
| COMPAS 2 changed every edge method from `f(u, v)` to `f((u, v))` — 12+ call sites across the library | a compatibility shim on our `Mesh` subclass (`_as_edge`) that **accepts both calling conventions**, so the rest of the library keeps the old signatures instead of being rewritten line by line |
| `Mesh.data` setter did hand-rolled version gating with `distutils.LooseVersion` (removed from Python 3.12) | deleted; COMPAS 2 serialisation used directly |
| `Polyline.point(t)` renamed to `point_at(t)` | updated |
| `__init__(self)` signatures blocked COMPAS 2's `from_*` constructors | every datastructure now takes `*args, **kwargs` and forwards them |
| `from .foo import *` re-exported the **module object** `foo`, shadowing same-named functions further up the chain | every package `__init__` now filters module objects out of `__all__` |
| `setup.py` shipped only the top-level package — subpackages were missing from an install | `find_packages(where="src")` |

An intermediate `_compat.py` shim (347 lines) carried the port while it was in progress;
it has since been **retired**, with its contents either resolved to real COMPAS 2 imports
or promoted into `datastructures/mesh/operations.py`. The library now imports COMPAS 2
directly, with no shim layer.

**Talking point:** this is not cosmetic. Upstream `compas_singular` is unusable today;
this fork is the only version that runs on the current COMPAS.

---

## 2. Layer two — a new front end: the frame field

*`src/compas_singular/framefield/`, 17 modules, **8 400 lines**. The single largest addition.*

Upstream has exactly one way to get from a boundary to a quad layout:
`SkeletonDecomposition` — the **medial axis** of a constrained Delaunay triangulation.
The medial axis is decided entirely by the outline. You cannot tell it where you want the
mesh lines to run.

We built a **second, drop-in front end** where the layout is generated *from a cross field*,
and cables / force lines are **inputs that shape that field**:

```python
d    = FieldDecomposition.from_boundary(outline, guides=cables)
mesh = d.quad_mesh(target_length=0.5)
```

`FieldDecomposition` deliberately mirrors `SkeletonDecomposition`'s method names, so an
existing script changes by two lines. Nothing downstream — `collect_strips`,
`densification`, `add_strip`, the grammar — was modified.

### The pipeline, module by module

| Module | Lines | Step |
|---|---|---|
| `background.py` | 314 | a triangulation a field can live on (upstream's puts vertices *only* on the boundary — useless for a field) |
| `field.py` | 422 | the cross field itself: one complex `u = exp(4iθ)` per vertex, so "smoothest cross field" becomes a plain sparse linear solve |
| `constraints.py` | 232 | pinning the field to a direction — boundaries and guide curves |
| `trace.py` | 754 | separatrix tracing, with guards for limit cycles, non-termination and near-singularity drift |
| `repair.py` | 1 239 | turning raw traces into a network `from_polylines` will accept |
| `arrangement.py` | 804 | the planar arrangement `from_polylines` *assumes* it is given |
| `decomposition.py` | 1 457 | `FieldDecomposition` — the public class |
| `densify.py` | 777 | **field-aware densification**: integrate the field inside a patch instead of blending its four edges |
| `relax.py` | 366 | global relaxation across patch seams, which per-patch relaxation cannot reach |
| `quality.py` | 281 | the two-tier element-quality gate |
| `symmetry.py` | 706 | keeping a symmetric input symmetric all the way to the layout |
| `edit.py` | 426 | the one seam: layout goes out to Rhino, comes back edited |
| `cache.py` | 748 | a solve cache that knows when it is stale |
| `guides.py` | 257 | directional metrics — does the mesh run *along* the guide, and mid-block? |
| `block.py`, `viz.py`, `__init__.py` | 332 | block dual, viewer scenes, exports |

### Four findings worth presenting

1. **`from_polylines` never computes an intersection.** Every curved domain was falling
   back to a non-field route. The stated reason ("the separatrix network does not close")
   was wrong: the network *does* close. `CoarsePseudoQuadMesh.from_polylines` hands raw
   segments to `Mesh.from_lines`, which keys vertices off segment endpoints only — so two
   separatrices crossing in their interiors got no node at the crossing and the recovered
   face visited the same vertex twice. The tell was there all along: the hexagon failed at
   background 0.5 and 0.4 and **succeeded at 0.35**, and a geometric limit does not come
   and go with resolution. Fixed in `arrangement.py`; **every curved domain tested now
   takes the field route** — disc, hexagon, ellipse, stadium.

2. **A 45° turn *is* a corner to a cross field.** `comb-plate` raised at every resolution
   and `rect with slot` produced a 0.00°/180.00° face while passing every structural check
   there was. One cause, and it was the fourth instance of the same species: *a tolerance
   scaled to a global quantity, applied to a domain whose local feature is smaller than
   that scale* — `boundary_corners`' spacing default of perimeter/25 swallowing the comb's
   2-unit teeth. The fix (`repair.SHARP_TURN`) is immune to that failure mode because
   **an angle has no length scale**: above 45° a vertex is a corner unconditionally.
   Measured margins: sampled curves reach 18°, real corners start at 60° — the threshold
   sits 2.5× above one and 1.33× below the other, and the harness fails if either side
   ever comes within 1.25×.

3. **"Manifold, all-quad, 100% coverage" was never the guarantee it sounded like.** Those
   three properties certified a mesh carrying a literal **180° angle** and an aspect ratio
   of 38. Robustness is now stated in **element quality** instead, in two tiers: a
   permanent hard floor (no degenerate angles, no non-finite aspect ratio) checked on every
   `quad_mesh()` call, and a committed 64-row regression baseline.

4. **`cable → field` always worked; `field → mesh` did not.** A square with a diagonal
   cable came out as a plain axis-aligned grid — *byte-identical to the unguided mesh* —
   because `densification` fills every patch with a bilinear Coons blend that never
   consults the field. `densify.py` closed this: 36.9° off the cable → **26.9°**, and
   25.2° over the middle of the cable where it is clear of a wall. What remains is a
   **ceiling, not a gap**: patch boundaries are held fixed so patches still weld, and on a
   plain square that boundary *is* the four walls.

### Verification

- **64-run locked baseline** — 16 domains × 4 background resolutions, every metric committed
  to `baseline.json` and diffed on every run.
- **One external check** — the disc's singularity radius against the published value of
  0.85. The default solver measures ~0.45 and *drifts under refinement*; `solve(relax=True)`
  measures 0.85. This is recorded as a known inaccuracy in the default solver, not hidden.
- **Head-to-head against the front end it replaces** — the same 16 shapes through both.
- **Solve cache correctness** — 44 assertions, one per way an entry can go stale.

### Honest state (as recorded in the suite)

| | Rows |
|---|---|
| route `field`, quality floor clear | 57 |
| route `triangulation`, floor clear | 4 |
| route `polygon`, **floor breached** | 2 |
| route `field`, **floor breached** | 1 |
| raises outright | **0** |

---

## 3. Layer three — the singularity as a *block*, not a joint

*`src/compas_singular/blocks.py`, **1 088 lines**, plus `framefield/block.py`.*

A conceptual addition rather than a port. In upstream, a point feature becomes a **pole**:
a vertex where mesh lines gather, its incident faces fanned into pseudo-quads. **A block
system cannot build that** — blocks meet along faces, never in a point.

`blocks.py` supplies the other half of the pair: a picked point that becomes the **centre
of an n-gon face with quads all round it** — the same singularity, carried by an element
instead of a joint.

The mechanism is one line of difference: pass the point to `boundary_triangulation` so the
layout wraps patches around it, but **withhold it from `decomposition_mesh(poles)`** so it
is never fanned. It then survives densification as a genuine valence-*n* joint with an
all-quad 1-ring, at the picked location, and `block_points` truncates it into the block.

Public API: `block_points`, `block_at`, `blocks_at`, `blocks_any_spin`, `blockable`,
`pole_at`, `pole_blocks`, `joints`, `joint_at`, `head_strips`, `relax`, `seam_edges`,
`index_sum`, `interior_valences`, `quality`.

Two properties the tests pin down: **topological index is conserved** through the
truncation, and a refusal always comes back **with a reason** rather than a silent failure.

---

## 4. Layer four — guide lines: making a curve steer the mesh

*`src/compas_singular/guide_lines.py`, **2 319 lines**.*

Two self-contained implementations, either liftable out of the file on its own:

| | What it does |
|---|---|
| `guide_line_mesh(coarse, guides)` | **GUIDE F** — snap existing polyedges onto the curves. No topology change. The guide becomes a chain of mesh edges. |
| `guide_band_mesh(coarse, guides)` | **GUIDE B** — inflate polyedges into bands with `add_strip`, snap each band's two new rails onto parallel offsets of the curve. The guide runs **down the middle of a row of faces**. |

The distinction that matters is not curved-versus-oblique — it is **where on the quads the
guide lands**. For a cable driving a block layout the answer is almost always "down the
middle", because a cable sitting on the shared edge between two rows is a cable in the joint.

Since the curve-feature restoration (§4a) there is a third:

| | What it does |
|---|---|
| `guide_feature_mesh(outer, guides)` | **GUIDE C** — build the layout *around* the guides by passing them as `polyline_features`. The guide is a course by construction, not by snapping. |

The first two are **post-processing**, and their docstring says so: *refinement subdivides
the lines that exist, it does not invent lines in new directions.* That limitation is what
motivated the frame-field work in §2. The third is not post-processing — it lets the guide
decide the layout, at the cost of no longer choosing that layout yourself, and of a
singularity wherever a guide ends in the interior.

---

## 4a. Layer four-a — curve features: a five-year regression, found and fixed

*`src/compas_singular/algorithms/decomposition.py`, `triangulation.py`.*

`polyline_features` — Oval's thesis §4.3.2, Figs 4.17–4.22 — raised `KeyError` in
`collect_strips` for any curve feature but a lucky single segment. It was not a missing
feature and it was **not a COMPAS 2 port regression**: it worked continuously from February
2019 to March 2021.

Commit `7d69ca0b` (16 March 2021, Robin Oval) commented out the call to
`quadrangulate_polygonal_faces` — the Fig 4.20 seam-propagation step. Its own message says
why: *"fix consequences of new method mesh.vertices_on_boundaries"*. COMPAS had renamed
`vertices_on_boundary()` **and changed its return from a flat vertex list to a list of
loops**. Of the three call sites in that file, two were migrated with the flattening idiom
and the third — inside `quadrangulate_polygonal_faces` — was missed; the call was disabled
rather than fixed. The 2026 port then renamed `mesh_explode` and `geometric_key` *inside the
dead method* without running it, so the un-migrated line survived that pass too.

Underneath sat a second, 2019-vintage bug: the method welded **inside** its
`for mesh in exploded()` loop, so `self.mesh` became the last component only. Fixing only
the migration returns a zero-face mesh for the one case that previously worked. That bug is
what the method's `# WIP` marker was about.

Two further defects, independent of those:

* `boundary_triangulation` never welded its feature input, which `discrete_mapping` has
  always done in Rhino. Now `weld_polyline_features`.
* `store_pole_data` recorded a pole only for a triangle with a **point feature** at a
  corner. A triangle is a legitimate pseudo-quad, but only with its pole recorded; curve
  extremities and medial-axis degeneracies produce triangles no point feature covers. That
  is the `"pole missing"` print, and the `KeyError` right after it.

### Results

Nine cases in `examples/10_curve_features.py`, previously five crashes, now all densify:
Fig 4.17 (extremities on the boundary), Fig 4.18 (extremities off it), Fig 4.19 (several
features), a mid-edge cross, two parallel features and a closed loop.

**Nothing that already worked moved.** `quadrangulate_polygonal_faces` calls `mesh_weld`,
which renumbers every key, so it returns immediately unless a face has more than four
vertices. No features, one point feature, two point features and the single diagonal are
identical face-for-face and key-for-key; that is a test, not an assertion.

### Still open

Chew's constrained Delaunay (Fig 4.20) is still not implemented — harmless at present, since
0% of feature segments were missing from the unconstrained triangulation across seven
configurations. Fig 4.22's "unwanted triangle" operation is not implemented either, so a free
extremity stays a pole instead of collapsing to a two-valent singularity: layouts carry more
singularities than the thesis figures, and a closed feature still needs the fallback repair.

---

## 5. Layer five — constrained smoothing without Rhino

*`src/compas_singular/datastructures/mesh/smoothing.py`, **589 lines**, 18 tests.*

Upstream's constrained smoothing lives in `compas_singular.rhino` and **needs Rhino open**
to project a vertex onto its constraint. Reimplemented against `compas.geometry` so it runs
headless, plus new capability:

- `boundary_constrained_smoothing` — boundary vertices that **slide along** their boundary
  curve instead of being pinned, with corners optionally held.
- `automated_boundary_constraints` — matches **every** boundary loop to its curve, not just
  the longest one (the trap in the original).
- `smoothing_region` — smooth a *region* with damping tapered over a blend band.
- A **windowed projection search**, because the naive full search was the dominant cost.

A finding recorded alongside it, because it is counter-intuitive and was measured:
**smoothing makes joint-perpendicularity worse, not better** — it degrades feature-line
perpendicularity from ~90° to ~45–65° while barely improving boundary regularity. Do not
smooth a mesh whose purpose is joint orthogonality to a force line.

---

## 6. Layer six — hand-editing a layout, with no CAD in it

*`src/compas_singular/editing/`, **1 977 lines**.*

Every operation a user performs on a coarse layout in Rhino — drag a corner, cut with a
drawn curve, delete a strip, commit back to the decomposition — is defined **on the mesh**
and runs from a plain script with no Rhino open. The Rhino side supplies picks, previews
and prompts, and nothing else.

The split is not tidiness: this is the half that can be wrong in ways a user cannot see.

| Module | Lines | What |
|---|---|---|
| `coarse_layout.py` | 1 122 | `CoarseLayoutEditor` — move corner, cut, delete strip, commit |
| `guide_chain.py` | 803 | choosing the dense-mesh vertices that follow a guide curve |

**Two results worth presenting:**

- **Edit the mesh, never the baked polylines.** `face_polylines` writes one closed polyline
  per patch, so an interior corner where four patches meet is **four coincident points**.
  Moving one and leaving three is a tear — and the weld rounds at 3 decimals, so the moved
  copy silently becomes a *new* vertex. The layout comes back with a slit in it and nothing
  reports an error.
- **The guide chain is a piece of a polyedge, and nothing else.** Three ways of building a
  chain vertex-by-vertex were tried and all three failed, each only visibly when drawn: a
  *radius* collects a ragged band from both sides and folds faces (max face angle 180.00);
  the *nearest vertex per station* hops between columns and comes out a zigzag; a *walk*
  that prefers the straight continuation still turns, and a turn puts it on a different
  polyedge. `collect_polyedges` already had the answer — the chain is not constructed, it
  is **chosen**, as the longest run of one polyedge that still follows the guide. Purity is
  then true by construction.

---

## 7. Layer seven — the Rhino plugin support layer

*`src/compas_singular/rhino/` — 1 535 new lines across four modules.*

Upstream's `rhino/` subpackage did not import at all. New modules, none of which duplicate
what upstream had:

| Module | Lines | What |
|---|---|---|
| `mesh_ui.py` | 400 | draw a mesh as pickable objects, pick → key, drag with live preview, ask a question. **Interaction only, decides nothing.** Shared by two commands. |
| `coarse_curves.py` | 594 | give a coarse layout back the curvature the Rhino round trip took off it |
| `dual_mesh.py` | 537 | the dual of a quad mesh, **closed on the boundary** |
| `helpers/helpers.py` | 404 | conversions, layer handling |

**Three things worth presenting:**

- **Rhino 8's own CPython runs the whole stack, unmodified.** That was the question with the
  power to kill the project. Python 3.9.10 + compas 2.15.1 + numpy 2.0.2 + scipy 1.13.1, all
  already present. Six domains solved in-app, every one on route `field`, 0.1–0.9 s. The
  same interpreter can be driven from the command line, so **Rhino compatibility is testable
  without opening Rhino**.
- **COMPAS's own dual is wrong for a plate.** `mesh_conway_dual` emits one dual face per
  *interior* vertex, so a plate shrinks inward by one full ring and the outline is gone —
  right for a closed polyhedron, wrong for a slab. `Mesh.dual(include_boundary=True)` closes
  the boundary but inserts the primal vertex into *every* closure face, so every regular
  boundary vertex becomes a pentagon. `dual_mesh.py` inserts it **only at a valence-2 corner**,
  which is the one place the cell would otherwise be a triangle — so corners and regular
  boundary vertices both come out as quads.
- **A single-precision round trip sets every tolerance in the file.** `rs.AddMesh` stores
  vertices as `Point3f`, so a corner read back from the document is ~2e-6 off on a 20-unit
  plate. Nothing in `coarse_curves.py` matches by vertex key — a bake renumbers them anyway.

The Rhino **commands themselves** (`CMD_*`, `ff_*`) live in the sibling `compas_topology`
repository; what is in *this* repo is the backend they call.

---

## 8. Layer eight — extensions to the core datastructures

Not new subsystems, but new capability on classes upstream already had.

**`QuadMesh` (`mesh_quad.py`, +259 lines)** — polyedges became first-class:

- `number_of_polyedges`, `polyedge_vertices`, `polyedge_edges`, `polyedge_midpoint`,
  `polyedge_length`, `polyline(pkey)` — all new.
- `collect_polyedge` / `collect_polyedges` / `vertex_opposite_vertex` /
  `singularity_polyedge_decomposition` gained a `strict=` mode; `collect_polyedge` also
  gained `both_sides=` and `oriented=`; `collect_strip` gained `both_sides=`.

**`datastructures/mesh/operations.py` (+171)** — `trimesh_face_circle`, `mesh_weld`,
`meshes_join`, `meshes_join_and_weld`, replacing what COMPAS 2 deleted.

**`geometry/array.py` (+60)** — `circle_evaluate`, `archimedean_spiral_evaluate`.

**`_viewer_patches.py` (65, new)** — a runtime fix for a real `compas_viewer` 2.0.2 bug:
`Renderer.init()` includes `Group` objects in the GPU settings buffer and
`rebuild_buffers()` excludes them, so **the first `scene.add()` after the viewer has ever
been shown permanently breaks group show/hide** for the rest of the process. The sidebar
checkbox keeps setting `group.show` and nothing reads it. Diagnosed to the vertex shader
(`getEffectiveShow` walks a parent-index chain that resolves to −1), fixed, and documented
with before/after screenshots.

---

## 9. Layer nine — tests and verification harnesses

Upstream shipped one trivial import test.

| | Tests | Covers |
|---|---|---|
| `tests/test_pipeline.py` | 8 | densification, pole survival through JSON, skeleton decomposition, block-from-point index conservation, refusal-with-reason, lizard grammar manifoldness |
| `tests/test_smoothing.py` | 18 | constraint dispatch, boundary loops/polylines/corners, sliding boundaries, multi-loop meshes, region damping taper, windowed-projection equivalence |
| `tests/test_import.py` | 1 | upstream's |

Plus, in `examples/New approach/`, harnesses that are tests in everything but name:
the 64-row `15_baseline.py` regression diff, `19_field_accuracy.py`'s external check,
`25_cache.py`'s 44 staleness assertions, `13_limits.py`'s reproducer-per-limitation,
`17_vs_skeleton.py`'s head-to-head, `24_symmetry.py`'s determinism proof.

**A deliberate rule worth stating:** the baseline and the symmetry check prove determinism
by **re-solving**. Putting the solve cache underneath either would make that a tautology,
so they are explicitly excluded from it.

---

## 10. Layer ten — documentation

Upstream had a README and a changelog. Eight substantial documents were written, all of
them recording *why* and *what was measured*, not just *what*:

| Document | Size | What |
|---|---|---|
| `HOW_IT_WORKS.md` | 256 lines | the data model and pipeline, why it is shaped this way, and the gotchas |
| `BASIC_COMMANDS.md` | 191 lines | command reference in the order you actually call things |
| `RHINO_PLUGIN.md` | 26 KB | plan, feasibility, measured environment, work list, progress log |
| `COMPAS_QUAD_REVIEW.md` | 37 KB | file-by-file review of Robin Oval's `compas_quad` |
| `New approach/README.md` | — | the frame-field front end, its state stated honestly |
| `cross-field-quad-meshing-evaluation.md` | — | why a field front end, which representation, academic sources |
| `frame-field-front-end-module-layout.md` | — | module design, 8-step plan, risk register |
| `paper-diffusion-generated-crossfields.md` | — | review of Dai/Qiao/Wang 2026, and where our field disagrees with it |

Plus a generated **visual overview page** (`New approach/site/`) built from the baseline
data — every domain drawn from its real outline with its route, patch count and quality.

### One finding from the `compas_quad` review worth its own slide

`compas_quad` (Robin Oval, Princeton) looks like the newer upstream. **It is not an
upgrade — it is a sibling fork of the same ancestor**, trimmed to the quad-grammar core,
and on COMPAS 2 compatibility this fork is strictly further ahead:

| | this fork | `compas_quad` |
|---|---|---|
| Library code | ~16 100 lines (at review time) | ~5 400 |
| COMPAS 2 imports | complete | partial — `compas_quad.coloring` fails outright |
| COMPAS 2 *call signatures* | complete | **not migrated** — 12 sites use the removed 2-arg edge API |
| `from_json` / `to_json` | works | **broken for every datastructure** |

Wholesale replacement is off the table. What it *does* contain is ~a dozen genuine
improvements to code that exists here, plus one entirely new feature (the ATP lizard
grammar), all ranked for cherry-picking.

---

## 11. Examples: 5 → 54 scripts

| Group | Scripts | What |
|---|---|---|
| Upstream, ported | `00`–`04` | densification, poles, planar/curved decomposition, lizard |
| Vault / block topology | `05` (852 lines) | four scenarios boundary → +opening → +supports → +force line, each with correctness checks *and* advisory quality metrics, plus a lift onto a doubly-curved shell |
| Bug reproducers with fixes | `06`, `07` | interior-triangle pole crash; square-with-centre-pole |
| Blocks | `08`, `09` | skeleton visualisation; a picked point that becomes a block |
| Workflow demos | `FloorExample*`, `MoreControl`, `Editing QuadMesh`, `StripDelete`, `Smoothing*`, `BasicCommands` | the user-facing workflows |
| Frame field | `New approach/01`–`25` (25 scripts, ~7 000 lines) | the front end, its checks, its limits, its baseline |

A convention runs through all of them: one `scene.add_group` per stage, laid out on a grid,
camera framed to the whole set, `--no-view` for numbers only.

---

## 12. Timeline

| Date | Milestone |
|---|---|
| **2026-08-05** | COMPAS 2 port lands — the library runs again |
| ~08-06 → 08-09 | vault/block topology example, pole crash fix, block-from-point |
| **2026-08-10** | Rhino feasibility established — the whole stack runs in Rhino 8's own CPython; six-step command series built |
| **2026-08-11** | Step 5's backend found missing and written (`framefield/edit.py`, `edit_coarse`) |
| **2026-08-14** | Symmetry work adopted and on by default; pre-change tree kept as a backup |
| **2026-08-20** | `compas_quad` review completed |
| ~08-21 → 08-27 | constraint-input study (10 force-line patterns), guide chain, coarse layout editor, solve cache |
| **2026-08-28** | `framefield` promoted from `examples/` into `src/` |
| **2026-09-01** | curve features restored — a five-year regression traced to `7d69ca0b` (§4a) |

---

## 13. What is still open

Stated as recorded, not softened.

1. **Multiply-connected domains with no corner anywhere fall back.** A disc with a round
   hole has no singularity and no reflex corner, hence no separatrix and no layout to
   repair. Their meshes are correct in shape but ignore the field.
2. **Patch count is inflated on the field route.** `solve_non_quad_faces`' guaranteed step
   is global, so one stray triangle quad-splits the whole layout — the hexagon's 11 raw
   patches become 43.
3. **One domain is still broken and the baseline records it rather than hiding it** —
   `square+ring cable`, a 180.00° face at three resolutions. It is in the matrix precisely
   because it is the suite's only source of poles, which is what a naive quality metric
   divides by zero on.
4. **The constraint input is not the bottleneck; the layout is.** The field takes all ten
   force-line patterns (0.0–3.8° off the curves fed in, against 36.9° unguided). The layout
   takes about half: patterns the walls can absorb are free, patterns forcing a singularity
   cost minimum angle 90° → 7.2°, and three of ten lose the field route. One measured
   mechanism accounts for it: a singularity nearer a wall than `0.45 × target_length`
   launches **none** of its separatrices.
5. **`relax=True` is not the default**, although it is the only setting under which real
   principal-stress load paths beat the walls alone — because it moves baseline rows
   unrelated to cables. **Open decision.**
6. **Manual frame rotation** has no backend — the one genuine gap in the Rhino workflow.
7. **Two settings blobs** exist across the `ff_*` and `CMD_*` command families, silently
   dropping the user's spacing and guide-alignment choices.

---

## 14. Repository / packaging state

| Branch | Contents |
|---|---|
| `main` | pristine upstream at the 2022 fork point |
| `feature-carbcomn` | the COMPAS 2 port |
| `dev` | **the working branch, and what the fork publishes** — the full library: the port plus the frame field, blocks, guide lines, the editing package and the research examples, all under `src/compas_singular/` |
| `dev-compas2` | the port **trimmed for a clean public release** — new-method examples and docs removed, `pyproject.toml` + `environment.yml` added, installable with a single `pip` command |

**New modules belong in `src/compas_singular/` and are committed to the fork** (policy
changed 2026-08-28). The earlier rule — that the fork carried the ported original library
only, and that the frame field, blocks, guide lines and editing package stayed off it —
is retired: it cost a second import root, a manual exclusion checklist before every push,
and extra `sys.path` entries for the Rhino commands. `dev-compas2` still exists, and what
a public release carries is decided when one is prepared, not when a module is written.

---

## 15. If you present one slide

> Upstream `compas_singular` stopped working four years ago and could only put mesh
> singularities where the *outline* dictated. This fork got it running on current COMPAS,
> then added a **second front end that generates the layout from a cross field**, so cables
> and force lines become **inputs** rather than something you snap the mesh to afterwards —
> plus the machinery to make that usable in practice: blocks instead of poles, headless
> smoothing and editing, a Rhino backend proven to run in-app, and a 64-row locked baseline
> that states what works and what does not.
>
> **9 200 → 26 900 lines. 5 → 54 examples. 1 → 27 tests. 0 → 8 design documents.**
