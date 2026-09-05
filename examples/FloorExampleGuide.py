"""Guide lines for a floor-slab quad mesh -- two implementations.

A POLE (point feature) is a *singularity*: it makes an odd number of strips
converge on a point, so it attracts every mesh line around it. A GUIDE LINE is
the opposite device -- along it the strip structure should stay maximally
*regular* -- so it is not built by attracting anything.


Why not `polyline_features` (the obvious answer, and the wrong one)
-------------------------------------------------------------------
``boundary_triangulation(..., polyline_features=[guide])`` unwelds the domain
along the curve, so the guide does survive as a chain of mesh edges. But the
medial axis is computed *from the domain*, so slitting the domain rewrites the
ENTIRE patch layout -- not just near the guide. On this floor plate it turns a
clean 7-quad decomposition into 20-odd patches including 5-, 6- and 7-sided
faces, and the mesh comes out visibly worse than with no guide at all:

    slit domain, guide at y =  0.0   interfaces across guide: median 80.5 deg, 72% within 15
    slit domain, guide at y = -1.5   interfaces across guide: median 50.1 deg, 32% within 15

(it is also erratic: how badly it does depends on where the slit happens to cut
the medial axis, which is not something you can design against)

The scenario is kept at the bottom of this file, measured, as GUIDE C.


The two guide logics implemented here
-------------------------------------
Both keep the CLEAN decomposition and edit the coarse mesh afterwards, so the
patch layout stays intact and only the guide region changes.

GUIDE A -- snap an existing polyedge  (``guide_by_snapping``)
    The decomposition usually already contains a polyedge running roughly where
    you want the guide. Project its vertices onto the guide line. No topology
    change at all.

GUIDE B -- add a strip           (``guide_by_strip``)
    ``add_strip`` inflates one polyedge into a band of new faces, giving two
    fresh parallel polyedges. Snap those onto the two rails of the guide band.
    This is the general answer: it *creates* the guide rather than reusing one,
    so you are not limited to the polyedges the decomposition happened to give
    you, and the band's two rails are parallel by construction.

Measured on this floor plate (see the report the script prints):

    GUIDE A  snap        interfaces across guide: median 88.1 deg, 100% within 15
    GUIDE B  add_strip   interfaces across guide: median 88.1 deg, 100% within 15

versus 50.1 deg / 32% for the slit at the same position. The guide is embedded
as a course at median 0.0 deg in every case -- that part was never the problem.


Oblique guides -- where the two logics stop being equivalent (GUIDE D)
----------------------------------------------------------------------
Run axis-aligned, A and B score the same, and B's only advantage is that it
*creates* a slot. Tilt the guide and they separate, because B's two rails stay
parallel TO EACH OTHER however the guide is oriented, so the patches inside the
band keep their parallel partner side; A's single line has the (non-parallel)
domain boundary on both sides. Measured over a tilt sweep:

    tilt      A: snap                    B: add_strip
     0.0 deg  88.1 deg, 100% within 15   88.1 deg, 100% within 15
     5.7 deg  84.8 deg, 100%             86.8 deg, 100%
    11.3 deg  79.6 deg,  65%             84.6 deg,  75%
    16.7 deg  75.4 deg,  50%             82.2 deg,  74%
    24.2 deg  67.5 deg,  46%             rejected by assert_not_folded

So: **prefer add_strip for anything off-axis.** Past ~24 deg the guide's ends
can no longer reach the boundary without the outline folding, and the guard
refuses it rather than returning a slab that silently overlaps itself -- on this
plate that is the honest limit, not a bug to tune away.

GUIDE D and GUIDE E below are the same 8.5 deg guide as a band and as a single
line: 85.8 deg / 81% within 15 versus 82.2 deg / 100% within 15. The band holds
a better median, the single line a tighter spread.


Curved guides (GUIDE F)
-----------------------
A guide can only bend where the coarse rail carrying it HAS vertices --
``densification`` fills each patch with a discrete Coons patch, whose sides are
straight between coarse vertices. The raw slot on this plate has **three**
vertices, so snapping it to a curve gives a chevron, not a curve.

``refine_coarse`` fixes that: a topological quad split of an all-quad mesh
leaves every vertex valence 4, so it buys rail resolution without introducing a
single new singularity. Two splits take the rail from 3 vertices to 9, which
carries a gentle bow cleanly.

``bowed_guide`` parametrises the guide as ``point_at(t)`` with a half-sine
sagitta, so it leaves and meets the boundary flat instead of kinking into it,
and ``bulge=0`` reduces exactly to the straight case -- straight, oblique and
curved guides all run through one code path. Measured for GUIDE F: the course
follows the curve at median 1.6 deg (not 0.0 -- the rail is an 8-segment
approximation of a smooth arc) and interfaces cross at 82.9 deg, 100% within 15.


Scenarios
---------
    BASELINE  no guide
    GUIDE A   snap an existing polyedge      straight, axis-aligned
    GUIDE B   add_strip band                 straight, axis-aligned
    GUIDE C   slit the domain                [REJECTED -- for comparison]
    GUIDE D   add_strip band                 oblique ~8.5 deg
    GUIDE E   snap an existing polyedge      oblique ~8.5 deg
    GUIDE F   snap, on a refined coarse mesh curved, side to side


The constraint you have to respect
----------------------------------
1. Only an INTERIOR WALL-TO-WALL polyedge can carry a guide -- both ends on the
   boundary, but not itself part of the boundary (``guide_slots``). This floor
   plate's clean decomposition has exactly ONE: ``[4, 13, 8]``, running
   (5.5, 0) -> (-5.5, 0). That is the real reason GUIDE B matters: adding a
   strip is how you get another slot.

2. Snapping slides the slot's two end vertices ALONG the boundary. Push the
   guide past a neighbouring boundary vertex and the coarse outline folds back
   on itself. On this plate the right-hand end vertex sits between neighbours at
   y = 1.0 and y = -5.0, so the guide is only valid for roughly -5 < y < 1; at
   y = 1.5 the mesh silently develops 2 inverted faces. ``assert_not_folded``
   catches this instead of letting it through.

3. Still no perpendicularity *guarantee*: a discrete Coons patch interpolates
   between OPPOSITE sides, so courses leave the guide at exactly 90 deg only
   where the patch across the guide has a parallel side. 88 deg here is the
   residual of a domain whose sides are not parallel, not a bug.

Run
---
    python FloorExampleGuide.py             # opens the compas_viewer scene
    python FloorExampleGuide.py --no-view   # prints the numeric report only
"""

import sys
from math import pi, sin, acos, degrees

from compas.geometry import Point, Polyline, Line, Translation
from compas.geometry import distance_point_point
from compas.geometry import dot_vectors
from compas.geometry import intersection_segment_segment_xy
from compas.geometry import normalize_vector
from compas.geometry import subtract_vectors
from compas.itertools import pairwise

from compas_singular.algorithms import boundary_triangulation, SkeletonDecomposition
from compas_singular.datastructures import CoarseQuadMesh
from compas_singular.datastructures.mesh_quad.grammar.add_strip import (
    add_strip,
    is_polyedge_valid_for_strip_addition,
)


# =============================================================================
# Configuration
# =============================================================================

BOUNDARY_RESOLUTION = 10        # samples per straight side of the outer boundary
TARGET_LENGTH = 0.6             # target quad edge length after densification
BAND_HALF_WIDTH = 1.2           # GUIDE B: offset of each rail from the centreline

GUIDE_BASE = [0.0, -1.5, 0.0]       # a point the guide passes through
GUIDE_DIRECTION = [1.0, 0.0, 0.0]   # the direction it runs in
OBLIQUE_DIRECTION = [1.0, 0.15, 0.0]    # GUIDE D/E: ~8.5 deg off, see docstring
CURVE_BULGE = 1.1               # GUIDE F: sagitta of the bowed guide
CURVE_REFINEMENT = 2            # GUIDE F: quad splits, to give the rail vertices


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
# Geometry helpers
# =============================================================================

def guide_ends(boundary, base, direction):
    """Where the infinite guide line through ``base`` crosses the boundary."""
    d = normalize_vector([direction[0], direction[1], 0.0])
    far = 1.0e4
    ray = ([base[0] - d[0] * far, base[1] - d[1] * far, 0.0],
           [base[0] + d[0] * far, base[1] + d[1] * far, 0.0])
    hits = []
    for p, q in pairwise(list(boundary) + [boundary[0]]):
        if distance_point_point(p, q) == 0.0:
            continue                               # duplicated corner sample
        x = intersection_segment_segment_xy(ray, (p, q))
        if x is not None:
            hits.append([x[0], x[1], 0.0])
    if len(hits) < 2:
        raise ValueError("the guide does not cross the domain")
    hits.sort(key=lambda p: dot_vectors(subtract_vectors(p, base), d))
    return hits[0], hits[-1]


def bowed_guide(start, end, bulge=0.0):
    """Parametrise a guide from ``start`` to ``end`` as ``point_at(t)``, t in [0, 1].

    ``bulge`` is the sagitta: the guide bows that far off the straight chord at
    mid-span and comes back to the boundary at both ends (a half sine, so it is
    flat where it meets the boundary rather than kinking into it).
    ``bulge = 0`` gives the straight chord, so straight and curved guides go
    through exactly the same code path.
    """
    chord = subtract_vectors(end, start)
    normal = normalize_vector([-chord[1], chord[0], 0.0])

    def point_at(t):
        offset = bulge * sin(pi * t)
        return [start[i] + chord[i] * t + normal[i] * offset for i in range(3)]

    return point_at


def sample_guide(point_at, samples=48):
    """Polyline through a parametrised guide, for measuring and drawing."""
    return [point_at(i / samples) for i in range(samples + 1)]


def offset_base(base, direction, distance):
    d = normalize_vector(direction)
    return [base[0] - d[1] * distance, base[1] + d[0] * distance, 0.0]


# =============================================================================
# Coarse-mesh editing
# =============================================================================

class PlainDecomposition(SkeletonDecomposition):
    """SkeletonDecomposition without the pseudo-quad pole bookkeeping.

    No point features are used here, so there are no poles to record; the base
    class would only print 'pole missing' for stray triangles.
    """

    def store_pole_data(self, poles):
        self.mesh.attributes['face_pole'] = {}


def build_baseline():
    """The clean decomposition of the bare domain -- no features, no slit."""
    trimesh = boundary_triangulation(outer_boundary=outer, inner_boundaries=[],
                                     polyline_features=[], point_features=[])
    coarse = PlainDecomposition.from_mesh(trimesh).decomposition_mesh([])
    coarse.collect_strips()
    return coarse


def guide_slots(mesh):
    """Polyedges that can carry a guide: interior, and spanning wall to wall.

    ``is_polyedge_valid_for_strip_addition`` accepts anything with both ends on
    the boundary -- including the boundary chains themselves, which must be
    excluded or snapping would drag the outline of the slab around.
    """
    slots = {}
    for pkey, polyedge in mesh.collect_polyedges():
        if not is_polyedge_valid_for_strip_addition(mesh, list(polyedge)):
            continue
        if all(mesh.is_vertex_on_boundary(v) for v in polyedge):
            continue                                        # the boundary itself
        slots[pkey] = list(polyedge)
    return slots


def pick_slot(mesh, base, direction):
    """The slot most aligned with the guide direction and closest to it."""
    slots = guide_slots(mesh)
    if not slots:
        raise ValueError("no interior wall-to-wall polyedge to carry a guide")
    d = normalize_vector(direction)
    normal = [-d[1], d[0], 0.0]

    def score(polyedge):
        pts = [mesh.vertex_coordinates(v) for v in polyedge]
        chord = normalize_vector(subtract_vectors(pts[-1], pts[0]))
        alignment = abs(dot_vectors(chord, d))                      # 1 == parallel
        offset = sum(abs(dot_vectors(subtract_vectors(p, base), normal))
                     for p in pts) / len(pts)
        return alignment - 0.15 * offset

    pkey = max(slots, key=lambda k: score(slots[k]))
    return pkey, slots[pkey]


def snap_rail(mesh, rail, start, end, point_at):
    """Move a chain of coarse vertices onto the guide ``point_at`` spans.

    Each vertex keeps its position ALONG the guide -- its projection parameter
    on the straight chord -- and only moves across onto the guide itself. So a
    straight guide moves every vertex the minimum distance, and a curved one
    just adds the bow rather than sliding vertices along the span.

    The two ends are boundary vertices and must stay on the boundary, so they go
    to the guide's endpoints (t = 0 / 1), which are the boundary crossings.
    """
    chord = subtract_vectors(end, start)
    chord2 = dot_vectors(chord, chord)
    last = len(rail) - 1
    for i, vkey in enumerate(rail):
        p = mesh.vertex_coordinates(vkey)
        if i == 0 or i == last:
            t = 0.0 if distance_point_point(p, start) <= distance_point_point(p, end) else 1.0
        else:
            t = min(1.0, max(0.0, dot_vectors(subtract_vectors(p, start), chord) / chord2))
        q = point_at(t)
        mesh.vertex[vkey]['x'], mesh.vertex[vkey]['y'], mesh.vertex[vkey]['z'] = q


def assert_not_folded(mesh, label):
    """Guard: snapping must not push a boundary vertex past its neighbour.

    A folded coarse outline shows up as a face with non-positive signed area in
    XY. Densifying one produces a self-overlapping slab, silently.
    """
    folded = []
    for fkey in mesh.faces():
        pts = [mesh.vertex_coordinates(v) for v in mesh.face_vertices(fkey)]
        area2 = sum(p[0] * q[1] - q[0] * p[1] for p, q in zip(pts, pts[1:] + pts[:1]))
        if area2 <= 1e-9:
            folded.append(fkey)
    if folded:
        raise ValueError(
            f"{label}: the guide is outside the range its slot can reach -- "
            f"{len(folded)} coarse face(s) inverted. Move GUIDE_BASE back "
            f"towards the middle of the plate (see the module docstring)."
        )


def rail_side(mesh, rail, base, direction, exclude):
    """Which side of the guide the given rail's faces are on (sign of the offset)."""
    d = normalize_vector(direction)
    normal = [-d[1], d[0], 0.0]
    values = []
    for vkey in rail:
        for nbr in mesh.vertex_neighbors(vkey):
            if nbr in exclude:
                continue
            values.append(dot_vectors(subtract_vectors(mesh.vertex_coordinates(nbr), base),
                                      normal))
    return sum(values) / len(values) if values else 1.0


def quad_split(vertices, faces):
    """Catmull-Clark *topological* split: every n-gon becomes n quads.

    Original vertices keep their position (no smoothing). Used for two things:
    refining a coarse mesh so a guide can curve (``refine_coarse``), and making
    the slit decomposition's 5/6/7-gons walkable at all (``guide_by_slit``).
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


def refine_coarse(coarse, times=1):
    """Refine the coarse mesh by topological quad split, adding no singularities.

    A guide can only bend where the coarse rail carrying it HAS vertices. The
    raw slot on this plate has three, so snapping it to a curve would give a
    chevron, not a curve -- and ``densification`` would honour the chevron,
    because a discrete Coons patch interpolates its sides as straight segments
    between coarse vertices.

    Splitting an all-quad mesh leaves every vertex valence 4 (face centres and
    edge midpoints alike), so this buys rail resolution without introducing a
    single new singularity.
    """
    for _ in range(times):
        index = {vkey: i for i, vkey in enumerate(coarse.vertices())}
        vertices = [coarse.vertex_coordinates(v) for v in coarse.vertices()]
        faces = [[index[v] for v in coarse.face_vertices(f)] for f in coarse.faces()]
        vertices, faces = quad_split(vertices, faces)
        coarse = CoarseQuadMesh.from_vertices_and_faces(vertices, faces)
    coarse.collect_strips()
    return coarse


# --- the two guide logics ----------------------------------------------------

def guide_by_snapping(coarse, base, direction, bulge=0.0):
    """GUIDE A -- reuse an existing polyedge; no topology change."""
    pkey, rail = pick_slot(coarse, base, direction)
    start, end = guide_ends(outer, base, direction)
    point_at = bowed_guide(start, end, bulge)
    snap_rail(coarse, rail, start, end, point_at)
    assert_not_folded(coarse, "GUIDE (snap)")
    return [sample_guide(point_at)], f"reused polyedge {pkey} ({len(rail)} vertices) = {rail}"


def guide_by_strip(coarse, base, direction, half_width=BAND_HALF_WIDTH, bulge=0.0):
    """GUIDE B -- add a strip, then snap its two new rails onto the band."""
    pkey, rail = pick_slot(coarse, base, direction)

    # add_strip pops from the list it is given, so hand it a copy
    skey, old_to_new = add_strip(coarse, list(rail))
    left = [old_to_new[v][0] for v in rail]
    right = [old_to_new[v][1] for v in rail]

    # the new vertices are created ON TOP of the vertex they replace, so the
    # strip has zero width until it is positioned -- densifying it in that
    # state raises ZeroDivisionError in Polyline.point_at
    new_vertices = set(left) | set(right)
    sign = 1.0 if rail_side(coarse, left, base, direction, new_vertices) >= 0 else -1.0

    guides = []
    for chain, offset in ((left, sign * half_width), (right, -sign * half_width)):
        start, end = guide_ends(outer, offset_base(base, direction, offset), direction)
        point_at = bowed_guide(start, end, bulge)
        snap_rail(coarse, chain, start, end, point_at)
        guides.append(sample_guide(point_at))

    assert_not_folded(coarse, "GUIDE (add_strip)")
    coarse.collect_strips()          # refresh strip data after the topology edit
    return guides, f"added strip {skey} along polyedge {pkey} = {rail}"


# =============================================================================
# GUIDE C -- the rejected approach, kept for the measured comparison
# =============================================================================

def guide_by_slit(base, direction):
    """Cut the DOMAIN along the guide -- rewrites the whole patch layout.

    Needs ``quad_split`` afterwards because the slit decomposition leaves 5-, 6-
    and 7-sided faces that ``collect_strips`` cannot walk (it dies with
    ``TypeError: cannot unpack non-iterable NoneType`` in ``face_opposite_edge``);
    ``SkeletonDecomposition.quadrangulate_polygonal_faces`` is commented out and
    marked WIP in decomposition.py.
    """
    start, end = guide_ends(outer, base, direction)
    guide = sample_guide(bowed_guide(start, end), samples=22)

    trimesh = boundary_triangulation(outer_boundary=outer, inner_boundaries=[],
                                     polyline_features=[guide], point_features=[])
    raw = PlainDecomposition.from_mesh(trimesh).decomposition_mesh([])
    sizes = sorted({len(raw.face_vertices(f)) for f in raw.faces()})

    index = {vkey: i for i, vkey in enumerate(raw.vertices())}
    vertices = [raw.vertex_coordinates(v) for v in raw.vertices()]
    faces = [[index[v] for v in raw.face_vertices(f)] for f in raw.faces()]
    vertices, faces = quad_split(vertices, faces)

    coarse = CoarseQuadMesh.from_vertices_and_faces(vertices, faces)
    coarse.collect_strips()
    return coarse, [guide], f"slit domain; decomposition face sizes {sizes} -> quad split"


# =============================================================================
# Densification
# =============================================================================

def densify_coarse(coarse, target_length=TARGET_LENGTH):
    """Densify a coarse mesh into blocks. No smoothing: a global centroid
    relaxation evens out quad sizes but relaxes the interior towards an
    even-but-non-orthogonal field, which is exactly what the guide is for."""
    coarse.set_strips_density_target(t=target_length)
    coarse.densification()
    return coarse.get_quad_mesh()


# =============================================================================
# Measurement
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


def measure_alignment(name, mesh, guides, note, on_tol=0.05, endpoint_margin=1.5):
    """* courses ALONG a guide -> angle to the tangent should be  0 deg
       * interfaces ACROSS it  -> angle to the tangent should be 90 deg

    The fan where a guide meets the boundary is excluded: the direction field
    has to turn somewhere, and it turns there.
    """
    print(f"\n{name}")
    print(f"  {note}")
    print(f"  dense mesh: {mesh.number_of_faces()} faces / {mesh.number_of_vertices()} vertices")
    if not guides:
        print("  (no guide -- baseline)")
        return

    on_guide = {}
    for vkey in mesh.vertices():
        p = mesh.vertex_coordinates(vkey)
        d, tangent = min((point_to_polyline(p, g) for g in guides), key=lambda dt: dt[0])
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
        print(f"  course ALONG guide  (want  0): median {median(along):5.1f} deg")
    if across_core:
        within = 100 * sum(1 for a in across_core if a >= 75.0) / len(across_core)
        print(f"  interface ACROSS it (want 90): median {median(across_core):5.1f} deg | "
              f"{within:3.0f}% within 15 deg")
    if across_ends:
        print(f"  endpoint fan (excluded):       median {median(across_ends):5.1f} deg "
              f"({len(across_ends)} edges)")


# =============================================================================
# Scenarios
# =============================================================================

def scenario_baseline():
    coarse = build_baseline()
    slots = guide_slots(coarse)
    note = (f"clean decomposition: {coarse.number_of_faces()} patches | "
            f"guide slots available: {slots}")
    return coarse, [], note


def scenario_snap():
    coarse = build_baseline()
    guides, note = guide_by_snapping(coarse, GUIDE_BASE, GUIDE_DIRECTION)
    return coarse, guides, note


def scenario_strip():
    coarse = build_baseline()
    guides, note = guide_by_strip(coarse, GUIDE_BASE, GUIDE_DIRECTION)
    return coarse, guides, note


def scenario_slit():
    return guide_by_slit(GUIDE_BASE, GUIDE_DIRECTION)


def scenario_oblique_band():
    coarse = build_baseline()
    guides, note = guide_by_strip(coarse, GUIDE_BASE, OBLIQUE_DIRECTION)
    return coarse, guides, note


def scenario_oblique_line():
    coarse = build_baseline()
    guides, note = guide_by_snapping(coarse, GUIDE_BASE, OBLIQUE_DIRECTION)
    return coarse, guides, note


def scenario_curved():
    # refine first: a 3-vertex rail can only bend into a chevron
    coarse = refine_coarse(build_baseline(), CURVE_REFINEMENT)
    guides, note = guide_by_snapping(coarse, GUIDE_BASE, GUIDE_DIRECTION,
                                     bulge=CURVE_BULGE)
    return coarse, guides, f"{note}; refined x{CURVE_REFINEMENT}, bulge {CURVE_BULGE}"


scenarios = []
for title, factory in (
    ("BASELINE  no guide", scenario_baseline),
    ("GUIDE A   snap an existing polyedge (no topology change)", scenario_snap),
    ("GUIDE B   add_strip band (creates the guide)", scenario_strip),
    ("GUIDE C   slit the domain with polyline_features  [REJECTED]", scenario_slit),
    ("GUIDE D   oblique add_strip band (~8.5 deg off-axis)", scenario_oblique_band),
    ("GUIDE E   oblique single line, snapped (~8.5 deg off-axis)", scenario_oblique_line),
    ("GUIDE F   curved single line, side to side", scenario_curved),
):
    coarse, guides, note = factory()
    dense = densify_coarse(coarse)
    scenarios.append((title, coarse, dense, guides, note))

print("=" * 78)
print(f"FLOOR PLATE   guide through {GUIDE_BASE} along {GUIDE_DIRECTION}")
print("=" * 78)
for title, _, dense, guides, note in scenarios:
    measure_alignment(title, dense, guides, note)


# =============================================================================
# Viewer -- the four plates side by side along X
# =============================================================================

if "--no-view" not in sys.argv:
    from compas_viewer import Viewer

    viewer = Viewer()
    step = 16.0

    for i, (title, coarse, dense, guides, _) in enumerate(scenarios):
        T = Translation.from_vector([i * step, 0.0, 0.0])
        group = viewer.scene.add_group(name=title)

        group.add(coarse.transformed(T), name="course mesh")
        group.add(outer_curve.transformed(T), name="outer boundary")
        group.add(dense.transformed(T), name="dense quad mesh")

        coarse_group = viewer.scene.add_group(name="coarse patches", parent=group)
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

    span = (len(scenarios) - 1) * step + 14.0
    centre = [0.5 * (len(scenarios) - 1) * step, 1.0, 0.0]
    viewer.renderer.camera.target = centre
    viewer.renderer.camera.position = [centre[0], centre[1] - 0.92 * span, 0.58 * span]

    viewer.show()
