# `pole missing` / `KeyError` on `face_pole` when using `polyline_features`

See it visually: `06_interior_triangle_pole_fix.py` builds a minimal domain that
reproduces this exact crash and shows the fix in `compas_viewer`, laid out as a
2x2 grid -- (1) the trimesh with the floating feature line, (2) the naive
decomposition with the unresolved interior triangle in red (and the crash
printed to console), (3) the same triangle registered as a valid pole after the
fix, and (4) the final densified quad mesh with that pole's singularity
highlighted in red among the mesh's other (pre-existing) singularities in
orange.

## Symptom

Running a headless decomposition with a non-empty `polyline_features` (a
feature/force curve, as opposed to `point_features`) prints:

```
pole missing
pole missing
```

and later crashes with:

```
File ".../mesh_quad_pseudo.py", line 95, in face_opposite_edge
    pole = self.attributes['face_pole'][fkey]
KeyError: 22
```

typically from `coarse_mesh.collect_strips()` (or later, `densification()`).

## Cause

`boundary_triangulation` ([triangulation.py](../src/compas_singular/algorithms/triangulation.py))
adds every `polyline_features` point as a Delaunay vertex and then topologically
cuts (unwelds) the mesh along that polyline. The medial-axis decomposition of
this slit domain can leave a triangular patch whose three corners are **all
interior vertices** (no vertex touches the outer/inner boundary).

`SkeletonDecomposition.solve_triangular_faces`
([decomposition.py:343](../src/compas_singular/algorithms/decomposition.py#L343))
only repairs triangles that touch the mesh boundary (1 or 2 boundary
vertices). A fully-interior triangle is left untouched.

`SkeletonDecomposition.store_pole_data`
([decomposition.py:454](../src/compas_singular/algorithms/decomposition.py#L454))
then tags every triangular face with a "pole" vertex, but only recognizes
poles that coincide with a `point_features` location. The leftover interior
triangle has none, so it prints `pole missing` and is never added to
`mesh.attributes['face_pole']`.

Later, `collect_strips` → `collect_strip` → `face_opposite_edge`
([mesh_quad_pseudo.py:95](../src/compas_singular/datastructures/mesh_quad_pseudo/mesh_quad_pseudo.py#L95))
looks up `face_pole[fkey]` unconditionally for every triangular face crossed
while walking a strip, so the untagged face raises `KeyError`.

This reproduces with any single, non-crossing `polyline_features` entry —
it is not specific to one script or geometry.

## Solution (used in `05_vault_block_topology.py`)

Subclass `SkeletonDecomposition` and override `store_pole_data` so that any
triangular face left without a `point_features` pole is registered as a
**pseudo-quad pole** instead (`face_poles[fkey] = mesh.face_vertices(fkey)[0]`).
Geometrically this is a legitimate valence-3 quad-mesh singularity — exactly
the kind of irregular vertex a feature line is meant to introduce — and it
densifies cleanly as a fan of quads *provided* the strips meeting at the pole
share the same density (i.e. use a uniform strip density, not per-strip
targets, when a feature line is present).

```python
class RobustSkeletonDecomposition(SkeletonDecomposition):
    """SkeletonDecomposition that also resolves fully-interior triangular faces
    left behind by a polyline_features cut."""

    def store_pole_data(self, poles):
        mesh = self.mesh
        pole_map = tuple(geometric_key(pole) for pole in poles)
        face_poles = {}
        for fkey in mesh.faces():
            if len(mesh.face_vertices(fkey)) == 3:
                chosen = None
                for vkey in mesh.face_vertices(fkey):
                    if geometric_key(mesh.vertex_coordinates(vkey)) in pole_map:
                        chosen = vkey
                        break
                if chosen is None:
                    # interior triangle from a feature cut -> make it a pole
                    chosen = mesh.face_vertices(fkey)[0]
                face_poles[fkey] = chosen
        mesh.attributes['face_pole'] = face_poles
```

Use `RobustSkeletonDecomposition.from_mesh(trimesh)` in place of
`SkeletonDecomposition.from_mesh(trimesh)`, and keep the strip density
uniform (e.g. `coarse_mesh.set_strips_density(d=...)`).

## Caveats

- Multiple/crossing `polyline_features` entries can still produce unresolved
  *polygonal* (>4-sided) faces that this shim does not address.
- Perpendicularity of interfaces to the feature line is only guaranteed "for
  free" when the feature runs along a full straight side of a coarse patch.
  An interior feature (embedded mid-patch, as here) does not get that
  guarantee — treat any perpendicularity requirement as advisory in that case.
