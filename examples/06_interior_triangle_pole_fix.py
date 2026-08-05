"""
06 - Interior-triangle pole crash: cause and fix
=================================================

Companion script to ``INTERIOR_TRIANGLE_POLE_FIX.md``. Reproduces -- then
fixes -- the ``pole missing`` / ``KeyError`` crash that happens when a
``polyline_features`` feature line floats free inside the domain (neither end
touching a boundary).

What the viewer shows
----------------------
    1. trimesh              -- the constrained Delaunay triangulation, with
                                the floating feature line in blue.
    2. naive decomposition  -- SkeletonDecomposition's coarse quad-patch mesh.
                                The feature cut leaves a fully-interior
                                triangular face (highlighted red) that
                                solve_triangular_faces cannot resolve (it only
                                fixes triangles that touch the boundary).
                                collect_strips() on this mesh raises the
                                KeyError reproduced in the console output.
    3. fixed decomposition   -- the same trimesh, decomposed with
                                RobustSkeletonDecomposition instead. The former
                                red triangle is now a valid pseudo-quad pole
                                (orange point) -- a valence-3 singularity --
                                and the mesh no longer crashes.
    4. densified quad mesh   -- the final quad mesh, with the singularity
                                coming from the fixed pole highlighted.

Run
---
    python "06_interior_triangle_pole_fix.py"            # opens the viewer
    python "06_interior_triangle_pole_fix.py" --no-view  # console only
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import sys

from compas.geometry import centroid_points

from compas_singular.algorithms import boundary_triangulation
from compas_singular.algorithms import SkeletonDecomposition
from compas_singular._compat import geometric_key


# =============================================================================
# The fix (identical to the one in 05_vault_block_topology.py)
# =============================================================================

class RobustSkeletonDecomposition(SkeletonDecomposition):
    """SkeletonDecomposition that also resolves fully-interior triangular faces
    left behind by a polyline_features cut. See INTERIOR_TRIANGLE_POLE_FIX.md."""

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


# =============================================================================
# Minimal domain that reproduces the problem
# =============================================================================

def segment(a, b, n):
    """Sample the open segment [a, b) with ``n`` points (endpoint b excluded)."""
    return [[a[0] + (b[0] - a[0]) * i / n,
             a[1] + (b[1] - a[1]) * i / n,
             a[2] + (b[2] - a[2]) * i / n] for i in range(n)]


def open_segment(a, b, n):
    """Sample the closed segment [a, b] with ``n`` + 1 points."""
    return segment(a, b, n) + [list(b)]


def sample_rectangle(w, l, spacing):
    """A closed rectangular boundary polyline sampled at ~``spacing``."""
    corners = [[0.0, 0.0, 0.0], [w, 0.0, 0.0], [w, l, 0.0], [0.0, l, 0.0]]
    pts = []
    for a, b in zip(corners, corners[1:] + corners[:1]):
        length = ((b[0] - a[0]) ** 2 + (b[1] - a[1]) ** 2) ** 0.5
        n = max(1, int(round(length / spacing)))
        pts += segment(a, b, n)
    return pts


def build_domain():
    outer = sample_rectangle(12.0, 8.0, 0.5)
    # both ends float free inside the domain -- neither touches the outer
    # boundary. A free-ended interior feature is exactly what leaves a
    # fully-interior triangular face behind in the medial-axis decomposition.
    # (An axis-aligned/centred segment here also produces an extra >4-sided
    # coarse face -- a separate, harder problem the shim does not address --
    # so this uses the same oblique placement as 05_vault_block_topology.py's
    # S5 scenario, which cleanly isolates the single-triangle case.)
    feature = open_segment([6.0, 2.0, 0.0], [10.0, 6.0, 0.0], 24)
    return outer, feature


# =============================================================================
# Pipeline: build both the naive and the fixed decomposition from one trimesh
# =============================================================================

def interior_triangles(mesh):
    """Triangular faces still present right after decomposition_mesh().

    solve_triangular_faces() already ran inside decomposition_mesh() and fixes
    every triangle that touches the boundary, so anything still triangular
    here is fully interior -- the case neither solve_triangular_faces() nor
    the base store_pole_data() (without a matching point_features pole) can
    resolve.
    """
    return [fkey for fkey in mesh.faces() if len(mesh.face_vertices(fkey)) == 3]


def run():
    outer, feature = build_domain()
    trimesh = boundary_triangulation(outer, [], [feature], [])

    # --- naive path: reproduce the crash ------------------------------------
    naive = SkeletonDecomposition.from_mesh(trimesh)
    naive_coarse = naive.decomposition_mesh([])  # prints "pole missing" below
    problem_faces = interior_triangles(naive_coarse)
    print("\nNaive SkeletonDecomposition: {} unresolved interior triangle(s): {}"
          .format(len(problem_faces), problem_faces))

    crash = None
    try:
        naive_coarse.collect_strips()
    except KeyError as e:
        crash = e
        print("Reproduced the crash: collect_strips() -> KeyError({})  <-- expected".format(e))
    else:
        print("collect_strips() did not crash here -- widen/move the feature "
              "in build_domain() so both ends stay clear of the boundary.")

    # --- fixed path ----------------------------------------------------------
    fixed = RobustSkeletonDecomposition.from_mesh(trimesh)
    fixed_coarse = fixed.decomposition_mesh([])
    fixed_coarse.collect_strips()
    fixed_coarse.set_strips_density(4)
    fixed_coarse.densification()
    dense = fixed_coarse.get_quad_mesh()
    pole_faces = interior_triangles(fixed_coarse)
    print("\nRobustSkeletonDecomposition: the same {} triangle(s) are now "
          "registered as pseudo-quad poles at vertices {} -- "
          "collect_strips()/densification() succeed."
          .format(len(pole_faces), [fixed_coarse.attributes['face_pole'][f] for f in pole_faces]))

    # trace the fixed pole(s) forward into the densified mesh (geometrically,
    # by coordinate match) so the viewer can single out the singularity that
    # came from the fix, versus the other singularities the layout has anyway
    pole_gkeys = {geometric_key(fixed_coarse.vertex_coordinates(fixed_coarse.attributes['face_pole'][f]))
                  for f in pole_faces}
    dense_pole_vertices = [v for v in dense.vertices() if geometric_key(dense.vertex_coordinates(v)) in pole_gkeys]

    return dict(outer=outer, feature=feature, trimesh=trimesh,
                naive_coarse=naive_coarse, problem_faces=problem_faces, crash=crash,
                fixed_coarse=fixed_coarse, dense=dense, pole_faces=pole_faces,
                dense_pole_vertices=dense_pole_vertices)


# =============================================================================
# Visualisation
# =============================================================================

def translate(mesh, dx, dy, dz=0.0):
    out = mesh.copy()
    for v in out.vertices():
        x, y, z = out.vertex_coordinates(v)
        out.vertex_attributes(v, 'xyz', [x + dx, y + dy, z + dz])
    return out


def face_centroid(mesh, fkey):
    return centroid_points([mesh.vertex_coordinates(v) for v in mesh.face_vertices(fkey)])


def view(result):
    from compas_viewer import Viewer
    from compas.geometry import Point, Polyline, Polygon
    from compas.colors import Color

    viewer = Viewer()
    pitch = 14.0

    # 2x2 grid: [1 trimesh | 2 naive] over [3 fixed coarse | 4 fixed dense],
    # each panel keeping the domain's own footprint (x in [0,12], y in [0,8]).
    # Point the camera at the middle of the grid and pull it back far enough
    # to see all 4 panels at once (the default camera is framed for one object).
    cx, cy = 0.5 * pitch + 6.0, -0.5 * pitch + 4.0
    dist = 2.0 * pitch
    viewer.renderer.camera.target = [cx, cy, 0.0]
    viewer.renderer.camera.position = [cx - dist, cy - dist, dist]

    # 1. trimesh + feature line -------------------------------------------------
    g1 = viewer.scene.add_group(name="1. trimesh + floating feature line")
    g1.add(result['trimesh'], show_points=False, show_lines=True, show_faces=True, name="trimesh")
    g1.add(Polyline([[p[0], p[1], p[2] + 0.05] for p in result['feature']]),
           linecolor=Color.blue(), linewidth=5, name="feature line (both ends free)")

    # 2. naive decomposition: interior triangle highlighted, in red ------------
    dx, dy = pitch, 0.0
    g2 = viewer.scene.add_group(name="2. naive decomposition -- pole missing, crashes")
    g2.add(translate(result['naive_coarse'], dx, dy), show_points=False, show_lines=True,
           show_faces=True, name="coarse mesh (naive)")
    for fkey in result['problem_faces']:
        pts = [Point(*result['naive_coarse'].vertex_coordinates(v)) for v in result['naive_coarse'].face_vertices(fkey)]
        pts = [Point(p.x + dx, p.y + dy, p.z + 0.03) for p in pts]
        g2.add(Polygon(pts), facecolor=Color.red(),
               name="unresolved interior triangle (fkey {})".format(fkey))
        outline = pts + pts[:1]
        g2.add(Polyline(outline), linecolor=Color.red(), linewidth=6,
               name="interior triangle outline")
    if result['crash'] is not None:
        print("Group 2 ('{}'): collect_strips() -> KeyError({})".format(g2.name, result['crash']))

    # 3. fixed decomposition: same triangle is now a pseudo-quad pole -----------
    dx2, dy2 = 0.0, -pitch
    g3 = viewer.scene.add_group(name="3. RobustSkeletonDecomposition -- fixed")
    g3.add(translate(result['fixed_coarse'], dx2, dy2), show_points=False, show_lines=True,
           show_faces=True, name="coarse mesh (fixed)")
    for fkey in result['pole_faces']:
        c = face_centroid(result['fixed_coarse'], fkey)
        g3.add(Point(c[0] + dx2, c[1] + dy2, c[2] + 0.05), pointcolor=Color.orange(), pointsize=20,
               name="pseudo-quad pole (was fkey {})".format(fkey))

    # 4. densified quad mesh ------------------------------------------------------
    dx3, dy3 = pitch, -pitch
    g4 = viewer.scene.add_group(name="4. densified quad mesh (fixed)")
    dense_t = translate(result['dense'], dx3, dy3)
    g4.add(dense_t, show_points=False, show_lines=True, show_faces=True, name="dense quad mesh")
    for vkey in dense_t.singularities():
        x, y, z = dense_t.vertex_coordinates(vkey)
        is_the_fix = vkey in result['dense_pole_vertices']
        g4.add(Point(x, y, z + 0.05),
               pointcolor=Color.red() if is_the_fix else Color.orange(),
               pointsize=22 if is_the_fix else 12,
               name="singularity from the fix" if is_the_fix else "other singularity (pre-existing)")

    viewer.show()


if __name__ == '__main__':
    result = run()
    if "--no-view" not in sys.argv:
        view(result)
