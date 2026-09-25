# How compas_singular works

Working notes on this package's data model and generation pipeline, written while
building `examples/05_vault_block_topology.py` (a 3DCP vault block-topology
generator). Not official docs — see `docsource/` for those. This file explains
*why* the pipeline is shaped the way it is, and records the gotchas we hit.

## 1. What it's for

`compas_singular` generates **quad meshes with controlled topology** — meshes
where you choose where the irregular vertices (*singularities*) sit, rather than
getting them wherever a generic quadrangulation algorithm happens to put them.
For a fabrication use case (3DCP block coursing, panelisation, structural grid
shells) that control is the whole point: singularities are where courses
converge/diverge, and every mesh edge becomes a physical interface.

The name "singular" refers to *singularity-driven* topology design: you feed in
a boundary + internal constraints, and the library builds a mesh whose
singularities are pinned exactly where those constraints demand.

## 2. Core vocabulary

| Term | Meaning |
|---|---|
| **Singularity** | A vertex that isn't 4-valent (interior) or 3-valent (boundary) — i.e. not part of a regular grid. Poles count as singularities too. |
| **Strip** | The chain of quad faces you get by crossing a quad edge-to-opposite-edge, repeatedly, until you hit a boundary or it closes into a loop. A strip has one **density** (its uniform subdivision count). Deleting/adding a whole strip is a topological edit. |
| **Polyedge** | The chain of *edges* running straight through a sequence of 4-valent vertices, from one singularity/boundary to another. Roughly the "dual" of a strip: a strip flows *across* a polyedge. |
| **Pole** | A singularity realised as a face with a repeated corner, `[a, b, c, c]` (a triangle stored as a degenerate quad). Used to fan an odd number of strips into one point — e.g. a column head where 3 or 5 courses meet. Requires `PseudoQuadMesh`. |
| **Coarse mesh** | The topology skeleton: one quad face per "patch", each patch edge carrying a *strip density*. Cheap to edit topologically. |
| **Dense mesh** | The coarse mesh after **densification**: each coarse face is filled with an actual grid of quads (a discrete Coons patch) according to its edges' strip densities. This is the fabrication-ready mesh. |

Mental model: **coarse mesh = topology graph, dense mesh = geometry**. You edit
topology cheaply on the coarse mesh (add/delete strips, set densities), then
densify once you're happy.

## 3. Class hierarchy

```
Mesh (compas_singular's thin subclass of compas.datastructures.Mesh)
 └─ QuadMesh                     strips, polyedges, singularities
     ├─ CoarseQuadMesh           + density attributes, densification()
     ├─ PseudoQuadMesh           + poles ([a,b,c,c] faces), pole-aware strip/opposite-edge logic
     └─ CoarsePseudoQuadMesh(PseudoQuadMesh, CoarseQuadMesh)   both: coarse AND poles
```

Source: `src/compas_singular/datastructures/mesh_quad*/`.

For the vault/block use case you almost always end up with
`CoarsePseudoQuadMesh` as the coarse mesh (because internal point supports need
poles), and its `densification()` output is a plain-ish `PseudoQuadMesh`
(dense, with `face_pole` data still attached for any surviving poles).

## 4. The generation pipeline

This is the one pipeline for going from "a boundary + some constraints" to a
controlled-topology quad mesh (`src/compas_singular/algorithms/`):

```python
trimesh   = boundary_triangulation(outer, inners, features, points)
decomp    = SkeletonDecomposition.from_mesh(trimesh)
coarse    = decomp.decomposition_mesh(points)          # CoarsePseudoQuadMesh
coarse.collect_strips()
coarse.set_strips_density(...)  # or set_strips_density_target(t)
coarse.densification()
dense     = coarse.get_quad_mesh()
```

### 4.1 `boundary_triangulation(outer, inners, features, points)`

Constrained Delaunay triangulation of the 2D (XY-planar) design domain.

- `outer` — outer boundary polyline (closed, list of XYZ points).
- `inners` — list of inner boundary polylines (holes/openings). Faces whose
  circumcentre falls inside one are deleted.
- `features` — polylines the triangulation is **cut along** (topologically
  unwelded): the mesh boundary is split open along the curve, so it becomes an
  edge path in the mesh rather than passing invisibly through triangle
  interiors. This is how you force a specific line to survive as mesh edges.
- `points` — extra Delaunay vertices (this is also where **internal supports /
  poles** enter the pipeline — see §6).

All four are literally point lists in the XY plane at this stage; nothing here
knows about 3D yet (see §7 for how we later lift the result onto a shell).

### 4.2 `SkeletonDecomposition.from_mesh(trimesh)` + `.decomposition_mesh(points)`

This is the core idea of the library. It computes the **topological skeleton**
(medial axis) of the triangulated domain — the set of triangle-circumcentre
polylines connecting the domain's singular triangles (corners, feature-forced
splits, curvature kinks) — and uses that skeleton to partition the domain into
quadrilateral **patches**. Every patch becomes one coarse-mesh face.

Concretely (`src/compas_singular/algorithms/decomposition.py`):

1. `branches_singularity_to_singularity()` / `..._to_boundary()` / boundary
   splits — build the skeleton branch polylines.
2. `decomposition_polylines()` — resolve intersections/splits between branches
   into a clean network of polylines.
3. `CoarsePseudoQuadMesh.from_polylines(...)` — build a coarse mesh whose faces
   are the regions enclosed by those polylines.
4. `solve_triangular_faces()` — patch up any degenerate/triangular faces this
   process leaves *if they touch the boundary* (see §8 for the gap here).
5. `split_quads_with_poles(points)` / `store_pole_data(points)` — for every
   input point that landed as a triangle corner, mark that face+vertex as a
   pole (`face_pole` attribute) instead of merging it away.

Net effect: **the medial axis IS the strip structure**. Every skeleton branch
becomes a polyedge; every gap between branches becomes a strip. This is why the
topology respects the boundary automatically — the skeleton is defined *by* the
boundary geometry.

### 4.3 Density + `densification()`

The coarse mesh only has topology. `collect_strips()` finds all strips;
`set_strips_density(d)` / `set_strip_density(skey, d)` /
`set_strips_density_target(t)` (target edge length → density from average strip
edge length) assign an integer subdivision count per strip.

`densification()` (`CoarsePseudoQuadMesh.densification`,
`src/compas_singular/datastructures/mesh_quad_pseudo_coarse/mesh_quad_pseudo_coarse.py`)
then fills every coarse face with a **discrete Coons patch**
(`compas.geometry.discrete_coons_patch`) built from its four (subdivided) side
polylines, and welds all the per-face patches together. Pole faces get one side
collapsed to `None` so the Coons patch degenerates into a triangular fan instead
of a quad grid.

**This is the mechanism behind interface perpendicularity** — see §7.

## 5. Curve ("line") features — how they work, and how they were broken

Curve features are Oval's thesis §4.3.2, Figs 4.17–4.22: hand a curve to
`boundary_triangulation` as a `polyline_features` entry and the pattern is laid
out **around** it. The curve becomes a *topological cut* in the Delaunay mesh
(`mesh_unweld_edges`), so the medial axis cannot cross it, and the decomposition
puts patch boundaries along it. The curve ends up as a continuous course of mesh
edges by construction — not by snapping.

Extremities carry the cost, as the thesis sets out (p. 100):

| where the curve ends | what you get |
|---|---|
| on the boundary | a three-valent boundary vertex, **no singularity** |
| off the boundary | a singularity — a pole, or a two-valent vertex |

### What the pipeline accepts

`boundary_triangulation` coerces its geometry, so all of these are equivalent:

```python
from compas.geometry import Point, Polyline

boundary_triangulation(outer, [], [[[1.4, 1.4, 0], [8.6, 8.6, 0]]], [])   # as always
boundary_triangulation(outer, [], Polyline([Point(1.4, 1.4, 0), Point(8.6, 8.6, 0)]), [])
boundary_triangulation(outer, [], [polyline_a, polyline_b], [Point(5, 5, 0)])
```

`as_points` / `as_curves` do the coercion and are exported, so a caller can use
them directly. A list of `[x, y, z]` is returned **as the same object**, so
existing callers pay nothing and their Delaunay tie-breaks cannot move. A single
curve may be passed without wrapping it in a list.

The closed-loop convention differs by role and is handled for you: an
`outer_boundary` or hole **drops** a final point coincident with the first, a
closed `polyline_features` entry **keeps** it, because the cut has to come back
round to where it started.

`point_features` is a flat list of `[x, y, z]` or `Point` — not a list of lists.

### What was broken, 2021 → 2026

Curve features raised `KeyError` in `collect_strips` for five years. Four
separate defects, worth keeping straight because only one of them was ever a
design limitation.

**1. The seam-propagation step was switched off.** Fig 4.20: the cut leaves
pentagonal and higher-valency faces, which "become quad faces by propagating the
seams of the discrepancies on the curve features". That is
`SkeletonDecomposition.quadrangulate_polygonal_faces`, and it had been called
since Robin Oval wrote it in February 2019.

Commit `7d69ca0b` (16 March 2021) commented the call out. Its own message says
why: *"fix consequences of new method mesh.vertices_on_boundaries"*. COMPAS had
renamed `vertices_on_boundary()` **and changed its return from a flat list of
vertices to a list of boundary loops**. There are three call sites in
`algorithms/decomposition.py`; two were migrated with the same flattening idiom
(`for bdry in ...: for vkey in bdry:`), the third — inside
`quadrangulate_polygonal_faces` — was missed, and the call was disabled instead
of fixed. The 2026 COMPAS 2 port then renamed `mesh_explode` → `.exploded()` and
`geometric_key` → `TOL.geometric_key` *inside the dead method* without ever
running it, so the un-migrated line survived that pass too.

**This failure mode is worth recognising**: a rename-only migration of an API
whose *return type* also changed, in code that no test reaches because its
caller is commented out.

**2. The method also welded inside its own loop.** `self.mesh = mesh_weld(mesh)`
sat inside `for mesh in supermesh.exploded()`, so `self.mesh` ended up as the
last component only. Fixing just defect 1 makes the single-diagonal case return
a mesh with **zero faces**. This is the bug behind the method's `# WIP` marker
and it dates from 2019 — the method was never trusted, which is probably why
disabling it looked cheap.

**3. Feature input was not welded.** `RhinoSurface.discrete_mapping` has always
pushed selected curves through `Network.from_lines` + `graph_polylines` before
triangulating; `boundary_triangulation` did not. Two polylines meeting at a
junction with no shared point leave the cut leaking through it, and the
decomposition returns an 11-gon. Now done inside `boundary_triangulation` by
`weld_polyline_features`. Note the limit, in Rhino too: curves that *cross*
without sharing a point cannot be welded, because there is no node. Supply a
crossing as the chains that meet at it — an X as four arms from the centre.

**4. Poles were sourced only from point features.** A triangular coarse face is
not a defect: it is the pseudo-quad `(p, a, b, p)` with one side collapsed at
the pole `p`, and `face_opposite_edge`, `collect_strips` and `densification` all
handle it. What breaks is a triangle with **no pole recorded** —
`face_opposite_edge` then raises `KeyError`. `store_pole_data` only matched a
triangle corner against the `poles` argument, so triangles from a curve
extremity, or from a medial-axis degeneracy, got none: that is the
`"pole missing"` print. It now assigns a pole to every triangle, preferring a
point feature and otherwise picking the corner whose two edges are closest to
equal, so the fan densifies evenly.

### Guarding the restoration

`quadrangulate_polygonal_faces` calls `mesh_weld`, which **rebuilds the mesh and
renumbers every key** — and the editing and agent layers address faces by key.
It therefore returns immediately unless some face has more than four vertices.
`tests/test_pipeline.py::test_curve_feature_repair_leaves_working_inputs_alone`
is the tripwire: no features, one point feature, two point features and the
single diagonal must come out with identical face counts, vertex counts and face
sizes.

### Two latent crashes in the corrections (2026-09-03)

Both are original defects, reachable without any curve feature, found from a
crash report on a real plate.

**`splits[0]` on a set.** `branches_splitting_collapsed_boundaries` built
`splits` as a `set` and then indexed it, so whenever a boundary loop has exactly
ONE split it raised `TypeError: 'set' object is not subscriptable`. Reached by a
square with two symmetric curve features -- **10 of 48 configurations crashed**.
`splits` is now a list in polyedge order with duplicates dropped, which also
removes an iteration-order dependency from the two-split branch. All 48 pass.

**`Polyline(None)`.** `solve_triangular_faces` case 2 passed
`decomposition_polyline(...)` straight to `Polyline`, and there is not always a
branch joining the two boundary vertices -- `TypeError: 'NoneType' object is not
iterable`. The merge is the point of that case and the polyline is only consulted
for WHERE to put the merged vertex, so it now falls back to the centroid of the
pair, which are coincident anyway.

Case 2 is hard to reach deliberately: instrumented across 25 featureless
configurations -- polygons, holes, discs, ellipses, thin strips -- it fired
**zero** times, with the concavity gate both on and off. So neither crash could
be attributed to the 2026-09-02 gate, and neither was.

### The collapsed edge, and the corner that is still cut off (2026-09-03)

Both found by pulling a real coarse mesh out of Rhino (`DemoAle.3dm`, 10 faces,
2 guides) over the MCP link.

**FIXED -- a singularity reading as a tiny edge.** `solve_triangular_faces` case 1
duplicates a vertex, leaving a zero-length edge that densification divides by, and
the tail of the method opens the pair again by moving each copy toward its own
neighbour centroid. The fraction was **0.1, too small to make an edge**: measured
0.019 on the Rhino plate -- aspect ratio **376**, a **178.6 degree** corner, and
the two irregular vertices `v:489.686,-34.028` / `v:489.701,-34.040` sitting
0.019 apart. Removing the nudge is not an option: the free-guide case then raises
`ZeroDivisionError`. At **0.5** the same case opens to 0.70 with the layout
unchanged. Exposed as `collapsed_edge_opening` and guarded by
`test_a_collapsed_edge_is_opened_enough_to_be_an_edge`.

**NOT FIXED -- a corner a feature lands on is cut off.** `corner_vertices` tests
`len(vertex_neighbors(v)) == 2`; at such a corner the wall gives two neighbours
and the feature a third, so CLOSING never splits the boundary there and the patch
cuts the corner off. Measured: a square with a diagonal guide keeps **2 of its 4
corners**.

Three attempts, all measured, none kept:

| attempt | result |
|---|---|
| discount feature edges | flags **142 vertices instead of 2** -- after the unweld the cut is a boundary too, so every vertex along it reads two-valent. 6 -> 164 patches |
| ...and require wall membership (`wall_gkeys` recorded at triangulation) | corners 4/4, but the split creates pentagons: diagonal 6 -> 8 patches with 2 triangles, both 20 -> **79** through the global repair |
| extended set for the boundary CUT only, plain set for the re-join | layout stays good, but corners back to **2/4** |

The last one localises it: **the corner only becomes a layout node when it is a
`splits=` for `graph_polylines`, and that is exactly what creates the pentagon.**
So this is the same blocker as everything else in R1/R2 -- a pentagon nothing can
turn into quads locally. Fix that first.

### Curve features are discretised like the walls (2026-09-03)

Thesis eq. 4.1 is stated for **every curve**, not only the walls.
`from_boundary` discretised the loops and handed `polyline_features` through
untouched -- and `curve_points` takes a polyline **at its own vertices**, so a
straight guide drawn in Rhino arrives as **two points**.

A feature is cut into the Delaunay along its OWN segments. One longer than the
surrounding sampling is not a Delaunay edge at all, so the cut never happens.
Measured on the Fig 4.17 diagonal against a wall sampled at 0.2:

| guide | segments missing | coarse layout |
|---|---|---|
| 2 points | **1 of 1** | 4 END faces, 4 patches, 0 singularities — *the no-feature layout* |
| 3 points | 2 | 4 END, 12 patches, 4 singularities |
| 26 points | 2 (the ones at the corners) | 4 END, 10 patches, 4 |
| 72 points | 0 | **6 END, 6 patches, 2** |

`geometry.discretise_line(line, spacing)` is a **separate function** from
`discretise_boundary`, not a flag on it: a boundary is a closed loop whose
closing segment is subdivided like any other, while a line is open and its two
**extremities must survive untouched** — whether a feature's end lands on a wall
decides whether the layout gets a node there. The two conventions do not belong
in one function.

`from_boundary` now discretises the features at the same spacing as the walls,
resolved against the same bounding-box diagonal (`_wall_spacing`), so a short
feature is not sampled far denser than the wall beside it. Every input point
survives, so a kink stays a kink. Guarded by
`test_a_two_point_guide_is_discretised_like_the_walls`.

**This is what stopped the fixes reaching Rhino.** `CMD_coarse_mesh` wires the
guides through correctly, but they arrived as two points and the cut did nothing.
(`CMD_boundary_selection` has a separate problem: line 280 sets
`polyline_features = [[2.5,3.5,0]]` and line 283 never passes it, so its guides
are dropped entirely.)

### The crossing fix (2026-09-03) -- all three benchmarks all-quad

Two curve features crossing in their interiors share no point, so the crossing is
no vertex, and **no triangulation can hold both crossing segments**. One is simply
absent, the cut leaks through the crossing, and the branches that pruning should
have removed run through it instead.

Chew (1989) cannot help here: the CDT is defined for "a set of vertices together
with a set of NONCROSSING edges", so a crossing breaks its precondition too. The
fix is a planar arrangement, not a stronger triangulator.

Measured on Fig 4.19, whose features cross at (4.7, 4.7): exactly **one** feature
segment is missing from the Delaunay -- the one containing the crossing -- and the
two surviving branches pass within **0.026** of it.

`arrange_polyline_features` puts a vertex at every crossing, before the weld.

| | before | after |
|---|---|---|
| 4.17 | 6 patches `{4:6}`, 2 singularities | unchanged |
| 4.18 | 10 patches `{4:10}`, 4 | unchanged |
| 4.19 | 20 patches `{4:18, 3:2}`, **10** | **20 patches `{4:20}`, 4** |

**This same change was tried on 2026-09-02 and made 4.19 far worse** (28 -> 115
patches, global quad-split). The difference is that pruning was broken then: an
adjacency across a feature still counted, so completing the cut at the crossing
only added structure the pruning could not remove. Order mattered -- the
arrangement is only useful once `real_neighbors` exists.

All three thesis benchmarks are now all-quad:

    4.17   6 patches, 2 interior singularities
    4.18  10 patches, 4
    4.19  20 patches, 4

### The grafting fix (2026-09-03) -- shared nodes on a feature

Thesis S4.2.2: "grafting adds branches between the singular face circumcentres
and their three vertices". Where two singular faces sit either side of a curve
feature, each grafts to its own nearest SAMPLE of that feature, and the Delaunay's
sampling puts those two samples one step apart rather than at the same place. The
patch spanning them comes out a triangle.

Measured on Fig 4.18: the singular faces at (2.75, 2.75) and (6.90, 6.90) grafted
to (4.77, 4.63) and (4.63, 4.77) -- 0.2 apart -- when both project to the
segment's midpoint (4.7, 4.7).

`SkeletonDecomposition.merge_graft_targets` collapses grafts that land on the same
feature at **adjacent samples** onto one shared sample.

| | before | after |
|---|---|---|
| 4.17 | 6 patches, 2 singularities | unchanged |
| 4.18 | 10 patches, **6 singularities** | 10 patches, **4 singularities** |
| 4.19 | 20 patches `{4:18, 3:2}`, 10 | unchanged |

**Adjacency along the chain is the test, not a distance.** The discretisation
decides what counts as the same place, so there is nothing to tune. Targets two
samples apart -- 4.19's (3.36, 6.04) and (3.08, 6.32) -- are deliberately left
alone: widening the rule to reach them gives 4 triangles at radius 0.5 and 82
patches at 0.8. Guarded by `test_graft_merge_only_joins_ADJACENT_samples`.

**Tested and NOT adopted: projecting a graft target onto the feature.** The
proposal was to move a target to the point of the feature closest to its
singularity. Implemented and measured: **0 of 12 targets move on 4.18 and 0 of 18
on 4.19** -- they are already the closest sample -- and on 4.17 it shifts two
targets 0.199 along the diagonal, doubling drift 0.0996 -> 0.1999 for no
topological gain. Reverted. The visual record is
`examples/images/12_graft_projection_test.png`; the adopted rule is
`13_graft_merge_proposal.png`.

### The pruning fix (2026-09-03) -- Fig 4.17 now matches the thesis

**A skeleton branch was crossing the curve feature.** `Skeleton.lines()` connects
the circumcentres of *adjacent* faces. `mesh_unweld_edges` never cuts the **first
and last segment** of a feature chain -- their end vertices are not split -- so
the two faces either side of that segment stay adjacent and the skeleton runs
straight across the feature. That is exactly Fig 4.20a, the thing the topological
cut exists to prevent.

Everything about Fig 4.17 followed from it:

- the corner where the feature lands never became an **END face** -- 2 were found
  where two triangles need 6 -- so **pruning could not remove the branches running
  to it**;
- those branches survived by wrapping *around* the corner (midpoints measured at
  (9.96, 9.90) and (0.10, 0.04)) and became extra divisions in the coarse mesh.

**The fix is classification-only.** `Skeleton.real_neighbors(fkey)` returns the
adjacent faces excluding any across a curve-feature edge, and `singular_faces`,
`corner_faces` and `lines` use it. The mesh is untouched: no vertex splitting, no
boundary splitting. With no features the neighbour set is unchanged, so a
featureless domain is bit-identical -- guarded by
`test_real_neighbors_is_inert_without_curve_features`.

| | before | after |
|---|---|---|
| 4.17 | 12 patches `{4:9, 3:3}`, 3 singularities | **6 patches `{4:6}`, 2 singularities** |
| 4.18 | 10 patches, 6 | unchanged -- no extremity lands on a wall |
| 4.19 | 26 patches `{4:21, 3:5}`, 11 | **20 patches `{4:18, 3:2}`, 10** |

**4.17 is now the thesis figure**: two Y-nodes, each with two arms to the walls
and one meeting the diagonal, six all-quad patches, two singularities. Guarded by
`test_no_skeleton_branch_crosses_a_curve_feature`.

The concavity gate below is still load-bearing -- with it disabled, 4.17 goes to
10 patches with 4 triangles and 4.19 to 82 through the global quad-split. The two
fixes are needed together.

### Following the thesis more closely (2026-09-02)

§4.2 of the thesis specifies the decomposition as **three operations** plus
**four ordered corrections**. The three operations map to the code 1:1, which is
why the skeletons are right. The corrections did not. What changed:

**1. The concavity correction now discriminates convex from concave.**
§4.2.3.1: *"The skeleton marks convex but not concave kinks ... branches are
added to include [the concavities]."* Only a concavity should get a branch.
`branches_splitting_boundary_kinks` used unsigned `angle_vectors`, so it could
not tell them apart, and corrected convex corners too. That is invisible on a
plain polygon — a convex corner is two-valent and the guard above catches it —
but a curve feature landing on a wall makes the corner **three-valent**, it slips
the guard, and gets corrected as though it were a concavity.

New `SkeletonDecomposition.boundary_interior_angle(vkey)` measures the interior
angle as the sum of incident face angles: **convex below π, concave above**. A
square corner reads 90°, an L-plate's reentrant corner 270°, and a vertex the
domain wraps entirely around — the end of a curve-feature cut — 360°. Measured
this way rather than from the turn of the boundary walk, whose sign flips between
the outer loop and a hole and is undefined on the loop around a cut, which
encloses no area.

Measured on the benchmarks:

| | before | after |
|---|---|---|
| 4.17 | 14 patches, 4 interior singularities | **12 patches, 3** |
| 4.18 | 10 patches, 6 | unchanged |
| 4.19 | 28 patches, 12 | **26 patches, 11** |

No fallback fires, and all of `examples/10_curve_features.py` still densifies.
Guarded by `test_boundary_interior_angle_separates_convex_from_concave` and
`test_convex_corner_carrying_a_curve_feature_is_not_corrected`.

**2. The corrections now run in the thesis's order.** §4.2.3 orders them
concavities → unwanted triangles → flipped patches → collapsed, and says why:
*"correcting collapsed boundaries occurs last"*, and flipped patches *"requires a
quad decomposition, so it occurs after correcting unwanted triangles"*. The code
ran collapsed second. Reordered. It changes no result today — the four generators
each read the Delaunay, not each other's output, so concatenating them in a
different order is inert — but the code no longer contradicts the specification
it implements, and a future correction that *is* state-dependent will be in the
right place.

**3. Thesis references in the functions.** Every method that implements a named
thesis operation now cites it — section number and the thesis's own word for the
operation — so a reader can open the right paragraph. The class docstring of
`SkeletonDecomposition` carries the full map, and `boundary_triangulation` cites
§4.2.1's discretisation rule and §4.3.2's topological cut.

**What was implemented, measured, and deliberately NOT kept:**

- **§4.2.3.2's rule for which corner gets the branch.** The thesis says: *"If two
  of the three patch corners are on the boundary, the branch is inserted at the
  other corner"* — split the interior corner, the mirror of case 1. The code
  merges the two boundary corners instead. Implemented behind a flag, then
  removed: across a square, an L-plate, a hole, a point feature and all three
  curve-feature benchmarks, **case 2 never fires once** (case 0 fires 12 times,
  case 1 four times). The two rules cannot be told apart, so shipping the branch
  would have been unexercised code. The finding is recorded at the case-2 comment.
- **Registering a wall extremity as a corner** (§4.3.2's "three-valent boundary
  vertex"). Gives 4.17 both extremities as proper nodes and halves its drift to
  0.037 — but takes 4.19 from 26 patches to 101 through the global quad-split.
  Reverted; see the attempts table below.

**The gap that matters now is case 0** — a triangle with all three corners
interior. It fires 12 times where case 2 fires none, it is what a free curve
extremity produces, and **neither the code nor §4.2.3.2 specifies it**. That is
the thing to solve before any further work on R1/R2.

### The specification, and what is present

Oval states the contract outright in §4.3.2 and on p. 100, so it can be tested
rather than eyeballed. `examples/11_thesis_curve_features.py` runs the thesis's
own three benchmarks (Figs 4.17–4.19) in its own four panels and reports each
rule.

| rule | requirement | status |
|---|---|---|
| R1 | extremity **on** the boundary → a three-valent boundary vertex, **no singularity** | partial — often not a layout vertex at all |
| R2 | extremity **off** the boundary → exactly **one** singularity, two-valent (or a pole, corrected to two-valent by Fig 4.22) | **not met** — a pole, or the tip is missed entirely; **zero** two-valent vertices are produced |
| R3 | pattern aligned along / around the curve | met — exact along the curve, drifting ≤0.14 within one patch of a free tip |
| R4 | coarse decomposition all quads | met |
| R5 | no skeleton branch crosses a feature | met — the unweld |
| R6 | features combined and superimposed | met |

And the four pipeline steps the text names:

| step | present? |
|---|---|
| points subdivide the curves, their edges **constrain** the Delaunay (Chew 1989) | **never, in any revision** |
| topological cut (unweld) along the curves | yes |
| seam propagation of the discrepancies → quads (Fig 4.20d) | yes, restored 2026-09-01 |
| correction: missed corners | yes, and it fires at the tips |
| correction: **unwanted triangles** (Fig 4.22) | exists, but **never fires at a cut** — it tests `is_vertex_on_boundary` on the coarse mesh, and `from_polylines` welds the cut shut, so every triangle at a feature reads case 0 and is skipped |

### Four attempts at R1/R2, measured and reverted

All four were implemented, measured against the three benchmarks, and **backed
out** on 2026-09-02. Recorded so they are not retried blind.

| attempt | result |
|---|---|
| **planar arrangement** — split features at their crossings, since Chew's CDT is defined only for *noncrossing* edges and a crossing costs exactly one missing Delaunay edge | 4.17/4.18 unchanged; **4.19 went 28 → 115 patches** and triggered the global quad-split. Helps a symmetric X of two full diagonals (96 → 32); wrecks the thesis benchmark |
| **pin free extremities** as decomposition splits | 4.19 finally produced two-valent tips and met 3 of 4 tip rules — at 28 → 123 patches. 4.18 got worse on every count |
| **cut-as-boundary** — let the unwanted-triangle correction see the cut | the most promising: 4.17 4 → 2 singularities, 4.18 10 → 8 patches and 6 → 2 singularities. But **4.19 went from 28 clean patches to 103 with the global fallback** |
| feeding the pinned tips into seam propagation as sources | no effect whatsoever |

**They all fail the same way.** Any extra node creates pentagons; seam
propagation cannot take them (it needs exactly four non-source corners) and the
*local* fan repair refuses them too (non-convex), so `solve_non_quad_faces`
falls through to the **global** topological quad split, which roughly quadruples
the patch count.

So the blocker is not the constraint and not the pinning: **it is that nothing
can turn those pentagons into quads locally.** Until that is solved, any work on
R1/R2 trades singularities for a 4× patch explosion. That is the thing to attack
next, or the current baseline is the sensible place to stop.

### What still is not implemented

* **The constrained Delaunay** of Fig 4.20 (Chew 1989). Currently harmless:
  across seven feature configurations, 0% of feature segments were missing from
  the unconstrained triangulation, because a feature is sampled at the same
  density as everything else. It would matter for a sparsely sampled feature.
* **Fig 4.22's "unwanted triangle" operation**, which collapses the pole at a
  free extremity into a two-valent singularity. Extremities here stay poles, so
  layouts carry more singularities than the thesis figures, and a guide can
  drift from the mesh edges within one patch of a free tip (measured: exact
  along the guide, ≤0.03 within one patch of the tip on a 10-unit plate).
* **A closed curve feature** goes through the fallback repair
  (`repair_polygonal_faces` → `framefield.repair.solve_non_quad_faces`) rather
  than seam propagation, because a face carrying more than one seam on the same
  side has fewer than four real corners. It densifies, but the layout around it
  is visibly poorer than an open feature's.

## 6. Internal supports → poles

A support point becomes a pole by being passed as one of the `points` /
`poles` arguments through the whole pipeline
(`boundary_triangulation(..., points)` → `decomposition.decomposition_mesh(points)`
→ `split_quads_with_poles` / `store_pole_data`). Concretely it needs to land
*exactly* on a Delaunay vertex, which is why it's passed into the triangulation
step, not added afterwards.

The result: the support point becomes a genuine mesh singularity — an odd
number of courses (typically 3 or 5) converge on it, exactly like a real
compression-only support needs (no ambiguous 4-way tension/compression split at
a point load).

## 7. Why interfaces come out perpendicular — and why NOT to smooth

This is the property that matters most for a compression-block assembly: mesh
edges (= block interfaces / joints) should cross a boundary, a force line, or a
support radially — perpendicular to it — so the thrust crosses every joint in
pure compression.

**It's already true in the raw densified mesh**, and here's the mechanism: a
discrete Coons patch fills a quad region with grid lines that run
**perpendicular to straight patch sides**. Every coarse-mesh edge is (locally)
straight, and every domain boundary / feature curve *is* a patch side by
construction (§4.1–4.2). So:

- at every regular boundary vertex, the inward course meets the boundary at ~90°;
- every interface edge crossing an embedded feature curve meets it at ~90°.

We verified this numerically on `05_vault_block_topology.py`'s scenarios:
boundary-perpendicular deviation median ≈0–2°, and interfaces crossing an
embedded force line: **median 90.0°, 100% within 15° of perpendicular** — with
**zero smoothing**.

**Counter-intuitive but measured: applying a global boundary-constrained
Laplacian/centroid smoothing pass makes this WORSE**, not better — it degrades
feature-line perpendicularity from ~90° down to ~45–65° while barely improving
boundary regularity, because it relaxes the free interior courses toward an
even-but-non-orthogonal field. **Don't smooth when perpendicularity matters.**
Smoothing is only useful for evening out block sizes, and should never be
applied to a mesh whose whole purpose is joint orthogonality to a force line.

The trade-off you keep instead: some sliver quads near poles and near where a
feature curve meets the boundary (a singularity fan has to be somewhere). Those
are a downstream block-quality problem (chamfer/merge/planarize), not a
topology problem.

## 8. Correctness vs. quality — how to think about testing this

Two different kinds of properties come out of this pipeline, and it's worth not
conflating them when testing (see `05_vault_block_topology.py`'s `Report`
class, which keeps `check()` — hard, must-pass — separate from `target()` —
advisory quality gates):

- **Topological correctness** (hard): quad-or-pole faces only, manifold, Euler
  characteristic matches boundary-loop count, every input boundary/opening/
  support/feature survives as the corresponding mesh feature, interfaces
  perpendicular to force lines. These either hold structurally or they don't —
  a bug here means the geometry is wrong, not just ugly.
- **Fabrication quality** (soft/advisory): planarity of individual blocks
  (inherently nonzero on a doubly-curved shell — fixed by a planarisation pass
  later, not by this topology step), minimum corner angle (avoid slivers),
  aspect ratio. These are targets for the *next* pipeline stage, not
  assertions this stage should be held to.

## 9. Environment

Examples run in the conda env **`singular312`**
(`C:\Users\Casper\anaconda3\envs\singular312\python.exe`): Python 3.12, compas
2.15.1, compas_singular installed **editable**, plus compas_viewer 2.0.2.
`carbcomn-core`'s compas_viewer is currently broken there (missing
freetype.dll) — don't use it for viewer scripts. To screenshot a
`viewer.show()` call non-interactively, use the `compas-viewer-screenshot`
Claude skill with this same interpreter.

## 10. Source map

```
src/compas_singular/
  algorithms/
    triangulation.py     boundary_triangulation()
    decomposition.py     Skeleton, SkeletonDecomposition (the medial-axis engine)
    layout.py, propagation.py, twocoloring.py, isomorphism.py   supporting/advanced
  datastructures/
    mesh/                 Mesh (compas 2.x compat shims)
    mesh_quad/             QuadMesh: strips, polyedges, singularities
    mesh_quad_coarse/       CoarseQuadMesh: density + densification (Coons patch)
    mesh_quad_pseudo/       PseudoQuadMesh: poles ([a,b,c,c] faces)
    mesh_quad_pseudo_coarse/ CoarsePseudoQuadMesh: coarse + poles (the usual coarse mesh type)
    skeleton/               Skeleton (base class SkeletonDecomposition builds on)
    lizard/                 turtle-graphics-style topological editing DSL (see examples/04_lizard.py)
  rhino/                  Rhino-only IO/constraints/artists (not used in the headless pipeline)
```

## 11. Worked example

`examples/05_vault_block_topology.py` puts all of the above together for the
carbcomn vault/block use case: four scenarios (boundary only → +opening →
+internal supports → +force line), each generated with this pipeline and run
through a battery of correctness checks + advisory quality metrics, plus a
lift onto an illustrative doubly-curved shell to sanity-check block
printability. Run it with `--no-view` for the numeric report or without flags
to open the compas_viewer scene (screenshot it via the skill above to actually
look at the mesh, since the viewer blocks headless execution otherwise).
