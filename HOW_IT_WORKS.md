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

## 5. Feature curves in the headless (non-Rhino) pipeline: the gotcha

`polyline_features` are only lightly exercised outside the Rhino workflow (the
Rhino example warps the mesh onto a real surface afterwards and never stresses
this path in isolation). Under COMPAS 2.x, headless, we found: **any
non-crossing single feature works, but the medial-axis cut it introduces
sometimes leaves a triangular coarse face with all three corners interior**
(not touching any boundary). `solve_triangular_faces()` only fixes
boundary-touching triangles, so that face is silently left triangular,
`store_pole_data` prints `"pole missing"`, and `collect_strips` /
`densification` crash later with a `KeyError` / `NoneType` unpack.

Fix used in `examples/05_vault_block_topology.py`: a `RobustSkeletonDecomposition`
subclass overrides `store_pole_data` to turn **any** leftover triangular face
into a pole (assign `face_pole[fkey] = face_vertices[0]`), whether or not it
matches an input point. This is a legitimate valence-3 singularity, not a hack —
it just needed the base class to also handle the "orphan interior triangle"
case, not only the "triangle produced by an intentional point pole" case. Needs
**uniform strip density** so the pole's fan of quads stays balanced.

Multiple *crossing* features still produce unresolved polygonal (>4-sided)
faces that this shim doesn't fix — keep features few and non-crossing.

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
  _compat.py              COMPAS 0.x/1.x -> 2.x shims compas_singular depends on
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
