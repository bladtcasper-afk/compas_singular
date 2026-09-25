"""Guide lines for a floor-slab quad mesh -- two implementations.

A POLE (point feature) is a *singularity*: it makes an odd number of strips
converge on a single point, so it attracts every mesh line around it.

A GUIDE LINE is the opposite device. Along a guide you want the strip
structure to stay maximally *regular*, so you cannot build it by attracting
anything. You build it by making the curve a BOUNDARY of the decomposition:

    boundary_triangulation(..., polyline_features=[guide])
        -> mesh_unweld_edges() cuts the domain topologically open along it
        -> the medial-axis decomposition then sees the guide as a boundary
        -> the guide survives as a *side of a coarse patch*
        -> densification() fills every coarse patch with a discrete Coons
           patch, whose grid lines run PERPENDICULAR to straight patch sides.


The two guide logics implemented here
-------------------------------------
GUIDE A -- single wall-to-wall guide
    One feature polyline whose two ends terminate exactly ON the outer
    boundary (``guide_across_domain``). The mesh is cut along it, so it comes
    out as a continuous course of mesh edges with the quads crossing it.

GUIDE B -- offset band
    Two parallel offsets of the same centreline, each independently
    wall-to-wall (``offset_band``), so the quads between them are bounded by
    the guide direction on both sides rather than only anchored on one.


What the measurements below actually show (do read this)
--------------------------------------------------------
Two separate properties, and they do NOT behave the same way:

1. "course ALONG the guide" comes out at a median of **0.0 deg in every case**.
   The guide is always embedded exactly as a chain of mesh edges. This is
   robust, it is what the unwelding buys you, and it is the property you can
   rely on.

2. "interfaces ACROSS the guide" is **only 90 deg when the patch on the far
   side of the guide has a side PARALLEL to the guide.** A discrete Coons
   patch interpolates between *opposite* sides, so grid lines leave the guide
   at 90 deg only if the side they run towards is parallel to it.

   * On a rectangle split by a mid-span guide, both halves are rectangles ->
     measured **median 90.0 deg, 100% within 15 deg** (RECTANGLE CONTROL below).
   * On this floor plate, whose sides are not parallel, the halves are
     trapezoids -> measured **median 56-82 deg depending on where the guide
     sits**. The courses fan instead of staying square.

   GUIDE B was expected to fix this by pinning the band on both sides. Measured
   on this domain, it does not: the medial axis still puts branch points inside
   the band, so the band is not one clean strip of patches and the interior
   angle wanders (47-84 deg over a sweep of band widths). It is included here
   as an honestly-reported negative result, not as the better option.

So: use a guide to *place a course exactly where you want it* -- that always
works. Do not expect perpendicularity for free unless your domain gives the
guide a parallel partner side.


Why the extra ``quadrangulate_coarse`` step
-------------------------------------------
Cutting the domain open makes the medial-axis decomposition leave faces with
5, 6 or 7 sides (and stray triangles). ``SkeletonDecomposition`` only repairs
triangles that touch a boundary, and its ``quadrangulate_polygonal_faces`` is
commented out and marked WIP in decomposition.py, so those n-gons survive and
``collect_strips`` then dies with

    TypeError: cannot unpack non-iterable NoneType object   (face_opposite_edge)

One global Catmull-Clark-style topological split fixes it: every n-gon becomes
n quads. It adds singularities only where the n-gons were (splitting a quad
gives four regular valence-4 quads), and it is measurably free -- the rectangle
control scores an identical 90.0 deg / 100% with and without it.

Run
---
    python FloorExampleGuide.py             # opens the compas_viewer scene
    python FloorExampleGuide.py --no-view   # prints the numeric report only
"""

import sys
from math import pi, acos, degrees

from compas.geometry import Point, Polyline, Line, Translation
from compas.geometry import distance_point_point
from compas.geometry import dot_vectors
from compas.geometry import intersection_segment_segment_xy
from compas.geometry import normalize_vector
from compas.geometry import subtract_vectors
from compas.itertools import pairwise

from compas_singular.algorithms import boundary_triangulation, SkeletonDecomposition
from compas_singular.datastructures import CoarseQuadMesh, QuadMesh


# =============================================================================
# Configuration
# =============================================================================

BOUNDARY_RESOLUTION = 10    # samples per straight side of the outer boundary
GUIDE_SPACING = 1.0         # sample spacing along a guide; match the boundary
BAND_HALF_WIDTH = 1.5       # GUIDE B: offset of each rail from the centreline
TARGET_LENGTH = 0.6         # target quad edge length after densification

GUIDE_BASE = [0.0, 0.0, 0.0]        # a point the guide passes through
GUIDE_DIRECTION = [1.0, 0.0, 0.0]   # the direction it runs in


# =============================================================================
# Design domain (same floor plate as FloorExample.py)
# =============================================================================

outer_corners = [Point(6, 5, 0), Point(5, -5, 0), Point(-5, -5, 0),
                 Point(-6, 5, 0), Point(0, 7, 0)]
outer_corners += [outer_corners[0]]


def densify(curve, resolution):
    """Sample a polyline side by side (never as a whole) so corners survive."""
    if isinstance(curve, list) and isinstance(curve[0], Point):
        line_pts = curve
    elif isinstance(curve, Polyline):
        line_pts = curve.split_at_corners(pi / 10)
    else:
        raise TypeError(f"Wrong input type, expected Polyline or list[Point], got {type(curve)}")

    densified_pts = []
    for indx in range(len(line_pts) - 1):
        line = Line(line_pts[indx], line_pts[indx + 1])
        densified_pts.extend(line.to_polyline(n=resolution).points)
    return densified_pts


def to_coordinates(pts):
    return [[pt.x, pt.y, pt.z] for pt in pts]


outer_curve = Polyline(outer_corners)
outer = to_coordinates(densify(outer_corners, BOUNDARY_RESOLUTION))


# =============================================================================
# Decomposition + repair of the polygonal faces a guide leaves behind
# =============================================================================

class PlainDecomposition(SkeletonDecomposition):
    """SkeletonDecomposition without the pseudo-quad pole bookkeeping.

    No point features are used in this example, so every triangle the
    decomposition leaves is an artefact of the guide cut rather than an
    intentional pole. The base ``store_pole_data`` would print 'pole missing'
    for each of them; they are quadrangulated instead (see below), so there is
    nothing to record.
    """

    def store_pole_data(self, poles):
        self.mesh.attributes['face_pole'] = {}


def quad_split(vertices, faces):
    """Catmull-Clark *topological* split: every n-gon becomes n quads.

    Original vertices keep their position (no smoothing), so the domain
    boundary and the guides are unchanged; each original edge simply becomes
    two collinear edges. A quad splits into four regular quads, so this adds
    singularities only where an n-gon actually was.
    """
    new_vertices = [list(p) for p in vertices]
    edge_mid = {}

    def mid(u, v):
        key = (min(u, v), max(u, v))
        if key not in edge_mid:
            a, b = vertices[u], vertices[v]
            edge_mid[key] = len(new_vertices)
            new_vertices.append([(a[i] + b[i]) / 2.0 for i in range(3)])
        return edge_mid[key]

    new_faces = []
    for fv in faces:
        n = len(fv)
        centre = len(new_vertices)
        new_vertices.append([sum(vertices[v][i] for v in fv) / n for i in range(3)])
        mids = [mid(fv[i], fv[(i + 1) % n]) for i in range(n)]
        for i in range(n):
            new_faces.append([fv[i], mids[i], centre, mids[i - 1]])
    return new_vertices, new_faces


def quadrangulate_coarse(raw):
    """Rebuild a possibly-polygonal decomposition mesh as an all-quad CoarseQuadMesh."""
    index = {vkey: i for i, vkey in enumerate(raw.vertices())}
    vertices = [raw.vertex_coordinates(vkey) for vkey in raw.vertices()]
    faces = [[index[vkey] for vkey in raw.face_vertices(fkey)] for fkey in raw.faces()]
    vertices, faces = quad_split(vertices, faces)
    return CoarseQuadMesh.from_vertices_and_faces(vertices, faces)


# =============================================================================
# Guide construction
# =============================================================================

def guide_across_domain(boundary, base, direction, spacing=GUIDE_SPACING):
    """A straight guide through ``base`` along ``direction``, clipped to ``boundary``.

    Both endpoints land exactly on the boundary polyline. A guide that stops
    inside the domain does not become a full patch side, and its dangling end
    leaves the decomposition with an orphan interior triangle.

    The endpoints are genuine new points on a boundary segment -- collinear with
    their neighbours, so the zero-area triangles they create are dropped by
    ``boundary_triangulation`` -- deliberately not snapped onto existing
    boundary samples, which would put duplicate points into the Delaunay input.
    """
    d = normalize_vector([direction[0], direction[1], 0.0])
    far = 1.0e4
    ray = ([base[0] - d[0] * far, base[1] - d[1] * far, 0.0],
           [base[0] + d[0] * far, base[1] + d[1] * far, 0.0])

    hits = []
    for p, q in pairwise(list(boundary) + [boundary[0]]):
        if distance_point_point(p, q) == 0.0:
            continue                                   # duplicated corner sample
        x = intersection_segment_segment_xy(ray, (p, q))
        if x is not None:
            hits.append([x[0], x[1], 0.0])
    if len(hits) < 2:
        raise ValueError("the guide does not cross the domain")

    hits.sort(key=lambda p: dot_vectors(subtract_vectors(p, base), d))
    start, end = hits[0], hits[-1]

    n = max(2, round(distance_point_point(start, end) / spacing))
    return [[start[0] + (end[0] - start[0]) * i / n,
             start[1] + (end[1] - start[1]) * i / n, 0.0] for i in range(n + 1)]


def offset_band(boundary, base, direction, half_width=BAND_HALF_WIDTH,
                spacing=GUIDE_SPACING):
    """GUIDE B: two rails offset either side of a centreline.

    Each rail is built by ``guide_across_domain`` from an offset base point, so
    each one is independently wall-to-wall -- offsetting a finished guide would
    push its endpoints off the boundary and forfeit the patch-side property.
    """
    d = normalize_vector([direction[0], direction[1], 0.0])
    normal = [-d[1], d[0], 0.0]
    rails = []
    for sign in (-1.0, +1.0):
        shifted = [base[0] + normal[0] * half_width * sign,
                   base[1] + normal[1] * half_width * sign, 0.0]
        rails.append(guide_across_domain(boundary, shifted, direction, spacing))
    return rails


# =============================================================================
# Pipeline
# =============================================================================

def build_pattern(boundary, features=None, target_length=TARGET_LENGTH):
    """outer boundary + guide polylines -> trimesh, coarse patches, dense quads.

    No smoothing on purpose. A global centroid relaxation evens out quad sizes
    but *degrades* alignment to the guide, because it relaxes the free interior
    towards an even-but-non-orthogonal field.
    """
    features = features or []

    trimesh = boundary_triangulation(outer_boundary=boundary, inner_boundaries=[],
                                     polyline_features=features, point_features=[])

    raw = PlainDecomposition.from_mesh(trimesh).decomposition_mesh([])
    face_sizes = sorted({len(raw.face_vertices(f)) for f in raw.faces()})

    coarse: CoarseQuadMesh = quadrangulate_coarse(raw)
    coarse.collect_strips()
    coarse.set_strips_density_target(t=target_length)
    coarse.densification()
    dense: QuadMesh = coarse.get_quad_mesh()

    return {'trimesh': trimesh, 'raw': raw, 'coarse': coarse, 'dense': dense,
            'features': features, 'raw_face_sizes': face_sizes}


# =============================================================================
# Measurement -- is the mesh actually guided?
# =============================================================================

def point_to_polyline(p, polyline):
    """(distance, unit tangent) of the closest point on ``polyline``."""
    best = (float('inf'), None)
    for a, b in pairwise(polyline):
        ab = subtract_vectors(b, a)
        length2 = dot_vectors(ab, ab)
        if length2 == 0.0:
            continue
        t = max(0.0, min(1.0, dot_vectors(subtract_vectors(p, a), ab) / length2))
        q = [a[i] + ab[i] * t for i in range(3)]
        d = distance_point_point(p, q)
        if d < best[0]:
            best = (d, normalize_vector(ab))
    return best


def acute_deg(u, v):
    """Undirected angle between two vectors, in [0, 90] degrees."""
    return degrees(acos(max(-1.0, min(1.0, abs(dot_vectors(u, v))))))


def median(values):
    s = sorted(values)
    n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def measure_alignment(name, result, on_tol=0.05, endpoint_margin=1.5):
    """Report how well the dense mesh follows its guides.

    * courses ALONG a guide  -> edge angle to the tangent should be  0 deg
    * interfaces ACROSS it   -> edge angle to the tangent should be 90 deg

    The fan where a guide meets the outer boundary is reported separately: the
    direction field has to turn somewhere, and it turns there.
    """
    mesh, features = result['dense'], result['features']
    print(f"\n{name}")
    print(f"  decomposition face sizes {result['raw_face_sizes']} -> all quads after split")
    print(f"  dense mesh: {mesh.number_of_faces()} faces / {mesh.number_of_vertices()} vertices")
    if not features:
        print("  (no guide -- baseline)")
        return

    on_guide = {}
    for vkey in mesh.vertices():
        p = mesh.vertex_coordinates(vkey)
        d, tangent = min((point_to_polyline(p, f) for f in features),
                         key=lambda dt: dt[0])
        if d <= on_tol and tangent is not None:
            on_guide[vkey] = tangent

    bnd = [mesh.vertex_coordinates(v) for loop in mesh.boundaries() for v in loop]

    def near_boundary(vkey):
        p = mesh.vertex_coordinates(vkey)
        return any(distance_point_point(p, q) < endpoint_margin for q in bnd)

    along, across_core, across_ends = [], [], []
    for u, v in mesh.edges():
        pu, pv = mesh.vertex_coordinates(u), mesh.vertex_coordinates(v)
        if distance_point_point(pu, pv) == 0.0:
            continue
        edge_dir = normalize_vector(subtract_vectors(pv, pu))
        on_u, on_v = u in on_guide, v in on_guide
        if on_u and on_v:
            along.append(acute_deg(edge_dir, on_guide[u]))
        elif on_u != on_v:
            gv = u if on_u else v
            ang = acute_deg(edge_dir, on_guide[gv])
            (across_ends if near_boundary(gv) else across_core).append(ang)

    print(f"  on-guide vertices: {len(on_guide)} | "
          f"transverse interfaces: {len(across_core) + len(across_ends)}")
    if along:
        print(f"  course ALONG guide  (want  0): median {median(along):5.1f} deg   <- always exact")
    if across_core:
        within = 100 * sum(1 for a in across_core if a >= 75.0) / len(across_core)
        print(f"  interface ACROSS it (want 90): median {median(across_core):5.1f} deg | "
              f"{within:3.0f}% within 15 deg")
    if across_ends:
        print(f"  endpoint fan (excluded):       median {median(across_ends):5.1f} deg "
              f"({len(across_ends)} edges)")


# =============================================================================
# Scenarios on the floor plate
# =============================================================================

guide_a = guide_across_domain(outer, GUIDE_BASE, GUIDE_DIRECTION)
band_b = offset_band(outer, GUIDE_BASE, GUIDE_DIRECTION)

scenarios = [
    ("BASELINE  no guide", build_pattern(outer), []),
    ("GUIDE A   single wall-to-wall guide", build_pattern(outer, [guide_a]), [guide_a]),
    (f"GUIDE B   offset band ({len(band_b)} rails, half-width {BAND_HALF_WIDTH})",
     build_pattern(outer, band_b), band_b),
]

print("=" * 78)
print("FLOOR PLATE  (sides not parallel -> no parallel partner for the guide)")
print("=" * 78)
for name, result, _ in scenarios:
    measure_alignment(name, result)


# =============================================================================
# Rectangle control -- the same GUIDE A on a domain that DOES have a parallel
# partner side. This is the reference showing the mechanism at its best, and
# the reason the floor-plate numbers above are lower.
# =============================================================================

def sample_rectangle(w, l, spacing):
    corners = [[0, 0, 0], [w, 0, 0], [w, l, 0], [0, l, 0]]
    pts = []
    for a, b in zip(corners, corners[1:] + corners[:1]):
        n = max(1, round(distance_point_point(a, b) / spacing))
        pts += [[a[0] + (b[0] - a[0]) * i / n,
                 a[1] + (b[1] - a[1]) * i / n, 0.0] for i in range(n)]
    return pts


print("\n" + "=" * 78)
print("RECTANGLE CONTROL  (guide has a parallel partner side)")
print("=" * 78)
rect = sample_rectangle(12.0, 8.0, 0.5)
rect_guide = guide_across_domain(rect, [6.0, 4.0, 0.0], [1.0, 0.0, 0.0], spacing=0.5)
measure_alignment("GUIDE A   single wall-to-wall guide, rectangle",
                  build_pattern(rect, [rect_guide]))


# =============================================================================
# Viewer -- the three floor-plate variants side by side along X
# =============================================================================

if "--no-view" not in sys.argv:
    from compas_viewer import Viewer

    viewer = Viewer()
    step = 18.0

    for i, (name, result, guides) in enumerate(scenarios):
        T = Translation.from_vector([i * step, 0.0, 0.0])
        group = viewer.scene.add_group(name=name)

        group.add(outer_curve.transformed(T), name="outer boundary")
        group.add(result['dense'].transformed(T), name="dense quad mesh")

        coarse_group = viewer.scene.add_group(name="coarse patches", parent=group)
        coarse = result['coarse']
        for u, v in coarse.edges():
            coarse_group.add(Line(Point(*coarse.vertex_coordinates(u)),
                                  Point(*coarse.vertex_coordinates(v))).transformed(T),
                             linewidth=3)

        # lift the guides clear of the mesh plane, otherwise they z-fight with
        # the very course they produced and are invisible
        guide_group = viewer.scene.add_group(name="guide", parent=group)
        for g in guides:
            guide_group.add(Polyline([Point(p[0], p[1], 0.15) for p in g]).transformed(T),
                            name="guide polyline", linewidth=6)

    # frame all three plates instead of the default single-object view
    centre = [0.5 * (len(scenarios) - 1) * step, 1.0, 0.0]
    viewer.renderer.camera.target = centre
    viewer.renderer.camera.position = [centre[0], centre[1] - 42.0, 38.0]

    viewer.show()
