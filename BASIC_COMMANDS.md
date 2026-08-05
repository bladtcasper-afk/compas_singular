# Basic commands in compas_singular

A command reference for `compas_singular`, organised in the order you'd actually
call things: build a mesh, densify it, optionally edit or analyse it, then view
it. For the *why* behind this pipeline (medial axis, poles, perpendicularity),
see [HOW_IT_WORKS.md](HOW_IT_WORKS.md) — this file is the short "what do I call
and what do I get back" version.

Run examples with the `singular312` conda env (Python 3.12, compas 2.15.1,
compas_singular installed editable, compas_viewer 2.0.2) — see
[HOW_IT_WORKS.md §9](HOW_IT_WORKS.md#9-environment).

## 0. What compas_singular does

`compas_singular` turns a planar boundary description (outer/inner curves,
optional feature curves and point supports) into a **coarse quad mesh** via a
topological skeleton (medial axis), then **densifies** that coarse mesh into a
fine quad mesh. Separately, it offers **grammar operations** to edit an
existing quad mesh strip-by-strip, and **analysis algorithms** (isomorphism,
two-coloring) that work on any quad mesh.

## 1. Data structures (what the commands operate on)

Inheritance chain, each level adding capability:

| Class | Adds | Source |
|---|---|---|
| `Mesh` | base compas mesh + strip/coloring helpers | `src/compas_singular/datastructures/mesh/mesh.py` |
| `QuadMesh` | strips, polyedges, grammar ops (add/delete/split strip) | `src/compas_singular/datastructures/mesh_quad/mesh_quad.py` |
| `PseudoQuadMesh` | quads collapsed to triangles = **poles** (singularities) | `src/compas_singular/datastructures/mesh_quad_pseudo/mesh_quad_pseudo.py` |
| `CoarseQuadMesh` | strip **densities** + `densification()` -> dense `QuadMesh` | `src/compas_singular/datastructures/mesh_quad_coarse/mesh_quad_coarse.py` |
| `CoarsePseudoQuadMesh` | coarse + poles combined (`from_polylines`, used by decomposition) | `src/compas_singular/datastructures/mesh_quad_pseudo_coarse/mesh_quad_pseudo_coarse.py` |
| `Skeleton` | topological medial-axis skeleton of a triangle mesh | `src/compas_singular/datastructures/skeleton/skeleton.py` |

Plus `Network` (edge graph, used internally for the skeleton) and `Lizard` (a
"turtle" for editing quad meshes by walking rules).

**Input/Output for all:** every class supports `.from_json(path)` /
`.to_json(path)` and `.from_vertices_and_faces(vertices, faces)` /
`.to_vertices_and_faces()`, since they're compas mesh subclasses.

```python
mesh = CoarseQuadMesh.from_json("coarse_quad_mesh.json")   # in: JSON, out: CoarseQuadMesh
```

## 2. Main pipeline: boundary curves -> dense quad mesh

This is the generative pipeline (see `examples/02_decomposition_discrete_planar.py`
and the annotated `examples/05_vault_block_topology.py`):

### Step A -- `boundary_triangulation(...)`

Constrained Delaunay triangulation of the design domain.

- **Input:** `outer_boundary` (closed polyline, list of XYZ points),
  `inner_boundaries` (list of closed polylines = holes), `polyline_features`
  (open/closed polylines the mesh must be cut/aligned along, e.g. force
  lines), `point_features` (points that will become mesh singularities/poles),
  optional `delaunay` callable.
- **Output:** a triangle `Mesh` (`trimesh`), boundaries respected, feature
  polylines topologically unwelded.

```python
trimesh = boundary_triangulation(outer_boundary, inner_boundaries, polyline_features, point_features)
```

Source: `src/compas_singular/algorithms/triangulation.py`

### Step B -- `SkeletonDecomposition.from_mesh(trimesh)` -> `.decomposition_mesh(poles)`

Computes the medial axis of `trimesh` and turns it into a coarse quad-patch
layout.

- **Input:** `trimesh` (from Step A); `poles` = the same `point_features` list
  (points to force as singularities).
- **Output:** `CoarsePseudoQuadMesh` -- a coarse quad mesh where each face is a
  Coons-patch placeholder, poles are pseudo-quads (triangles).

```python
decomposition = SkeletonDecomposition.from_mesh(trimesh)
coarsemesh = decomposition.decomposition_mesh(point_features)
```

Source: `src/compas_singular/algorithms/decomposition.py`

### Step C -- density + `densification()` -> `get_quad_mesh()`

Fills each coarse patch with a discrete Coons patch at the requested edge
count.

- **Input:** `coarsemesh.collect_strips()` first (groups parallel edges into
  "strips"); then one of:
  - `set_strips_density(d)` -- same integer subdivision count `d` on every strip
  - `set_strips_density_target(t)` -- pick per-strip density so edge length ~= `t`
  - `set_strips_density_func(func, func_args, skeys=None)` -- custom per-strip rule: calls
    `func(skey, func_args)` for each strip and stores the result via `int(...)` (truncates,
    doesn't round/ceil -- add `ceil()` inside `func` yourself if you want that). `func_args`
    is passed through untouched, so bundle whatever the rule needs into it (a dict keyed by
    `skey`, a target length, a tuple including the mesh itself if `func` needs geometry via
    closure isn't convenient)
  - `set_strip_density(skey, d)` / `set_strip_density_func(skey, func, func_args)` -- override
    a single strip (key from `coarse.strips()`)
  - then `densification()` (no args -- reads the density attributes just set)
- **Output:** `densification()` mutates the coarse mesh in place; `get_quad_mesh()`
  returns the resulting dense `QuadMesh`.

```python
coarsemesh.collect_strips()
coarsemesh.set_strips_density_target(0.5)   # or set_strips_density(3)
coarsemesh.densification()
densemesh = coarsemesh.get_quad_mesh()      # -> QuadMesh
```

Source: `src/compas_singular/datastructures/mesh_quad_coarse/mesh_quad_coarse.py`

> **Fragility note:** non-empty `polyline_features` in Step A can leave an
> unresolved interior triangle that Step B doesn't fix, crashing Step C. A
> single non-crossing feature works if you subclass `SkeletonDecomposition` to
> assign leftover interior triangles as pseudo-quad poles, and use uniform (not
> target-length) density. See `HOW_IT_WORKS.md §5` and the
> `RobustSkeletonDecomposition` shim in `examples/05_vault_block_topology.py`.

## 3. Editing an existing quad mesh (grammar operations)

Once you have any `QuadMesh` (dense or coarse), you can edit it strip-by-strip:

| Command | Input | Output |
|---|---|---|
| `add_strip(mesh, polyedge)` | a `QuadMesh` + a vertex polyedge (path) to duplicate as a new strip | mutates `mesh`; returns updated polyedge data |
| `delete_strip(mesh, skey)` | a `QuadMesh` + a strip key (from `mesh.strips()`) | mutates `mesh` (collapses the strip); returns `skey_to_skeys` remap dict or `None` |
| `split_strip` / `split_strips` | mesh + strip key(s) | subdivides a strip into two |

Source: `src/compas_singular/datastructures/mesh_quad/grammar/add_strip.py`,
`.../grammar/delete_strip.py`

### `Lizard` -- turtle-style rule editing

A "turtle" that walks the mesh edges and applies grammar rules from a string,
see `examples/04_lizard.py`:

- **Input:** a `QuadMesh` (with `collect_strips()` already called), a rule
  string of characters `t`/`p`/`a`/`d`.
- **Output:** mutates the mesh in place.

```python
lizard = Lizard(mesh)
lizard.initiate()                     # picks a starting directed edge (tail, head)
lizard.from_string_to_rules('atta')   # t=turn, p=pivot, a=add strip, d=delete strip
```

Source: `src/compas_singular/datastructures/lizard/lizard.py`

## 4. Analysis algorithms (on any quad mesh)

| Command | Input | Output |
|---|---|---|
| `are_meshes_isomorphic(mesh_i, mesh_j)` | two `QuadMesh` objects | `bool` |
| `are_strips_isomorphic(mesh_i, mesh_j)` | two `QuadMesh` objects | `bool` -- compares strip graphs only |
| `matches_between_ismorphic_meshes(mesh_i, mesh_j)` | two isomorphic meshes | vertex/face correspondence mapping |
| `TwoColourableProjection` | a `QuadMesh` | two-coloring of strips/polyedges (checks/produces bipartite coloring) |
| `quadrangulate_mesh(mesh, sources)` | a mesh with non-quad (polygonal) faces + source vertices | mutates mesh, converting polygons into quads |

Source: `src/compas_singular/algorithms/isomorphism.py`, `.../twocoloring.py`,
`.../propagation.py`

## 5. Visualization

Every data structure above is a compas mesh, so it drops straight into
`compas_viewer`:

```python
from compas_viewer import Viewer
viewer = Viewer()
viewer.scene.add(mesh, show_points=True, show_lines=True, show_faces=True)
viewer.show()
```

## Putting it together (the whole pipeline in one line each)

```
boundary_triangulation(...)                 curves        -> trimesh (Mesh)
SkeletonDecomposition.from_mesh(trimesh)    trimesh       -> decomposition (skeleton)
decomposition.decomposition_mesh(poles)     poles         -> coarse (CoarsePseudoQuadMesh)
coarse.collect_strips()                     -             -> strip data stored on coarse
coarse.set_strips_density_target(t)         target length -> per-strip density stored
coarse.densification()                      -             -> coarse mutated with dense sub-faces
coarse.get_quad_mesh()                       -             -> dense (QuadMesh)  <- final output
```

Everything downstream of `dense` (grammar edits, isomorphism checks, viewing)
is optional post-processing.
