#! python3

# r: compas
# r: pydantic

"""**Step: read a coarse layout DRAWN over the domain boundaries.**

    reads   Outer, Inner                             the domain walls (step 1)
            TopologyProblem::Skeleton::EdgeCurves    the division lines
            TopologyProblem::Skeleton::Poles         collapsed corners, if any
    writes  TopologyProblem::Skeleton::{Mesh, Poles, Polylines}
            the session's layout
    never   touches Outer, Inner or EdgeCurves -- they are the input

The mirror of ``CMD_coarse_mesh``: that one GENERATES a layout, this one reads
one you drew. Draw the domain the way every other step expects it -- the outer
boundary on ``Outer``, holes on ``Inner`` -- and the coarse division on
``EdgeCurves``.

WALLS AND DIVISIONS
-------------------
All three layers are read as coarse EDGES of one network. A wall is an edge
like any other, split wherever a division line lands on it. What the layers
add is which edges are WALLS, and that decides the one thing a drawing cannot
say by itself: which closed region is a hole. A hexagon drawn on ``EdgeCurves``
is a hexagonal patch; the same hexagon on ``Inner`` is a hole. Both are right --
they are different layouts, and only the layer says which one was meant.

The walls also define the domain, so a division line must stay inside it: one
that overshoots the outer boundary, or runs into a hole, is refused and
selected rather than read as a patch nobody meant.

HOW CURVES BECOME EDGES
-----------------------
* a POLYLINE is split at its own corners -- a rectangle drawn as one closed
  polyline is its four sides, what ``_Explode`` gives. A vertex that does not
  turn is not a corner;
* a SMOOTH curve (arc, circle, interpolated curve) is one edge, because its
  points are samples, not corners;
* every curve is then split wherever another curve MEETS it: where a division
  lands on a wall, and where two divisions cross.

What is NOT a meeting: an end that stops short of the curve it was meant to
reach. Snap ends onto the curves they land on (End, Near, Int), or they are
refused as dangling.

WALLS ALREADY ON ``EdgeCurves``
-------------------------------
``CMD_coarse_mesh`` and ``CMD_edit_coarse_mesh`` bake the WHOLE layout onto
``EdgeCurves``, walls included. A piece of ``EdgeCurves`` lying along a wall is
therefore not a second edge there -- two edges between the same corners would
be refused -- but it does say where that wall has corners. So its two ends are
kept as corners on the wall and the piece itself is dropped. That is also how a
corner on a smooth wall with no division reaching it is drawn: a disc as ONE
patch is its circle on ``Outer`` and four arcs of it on ``EdgeCurves``.

IT REFUSES RATHER THAN REPAIRS
------------------------------
A dangling curve, a division outside the domain, a piece of the drawing that
touches nothing else (a hole no division reaches), a patch that is not a quad
or a triangle, or two curves between the same two corners all stop the command
with the offending curve selected. Each is seconds to fix on the screen, and a
command that repaired them would hand back a layout nobody drew. The one thing
it does silently is treat a triangle as a pseudo-quad, which is what the rest
of the workflow already does with them.

WITH NO ``Outer`` CURVE
-----------------------
The divisions then have to carry the outer boundary themselves, as they did
before the walls were read separately: the layout's boundary is whatever ends
up with one patch beside it. ``Inner`` curves are still walls and holes.
"""

import re

import rhinoscriptsyntax as rs

import compas_rhino as cr
from compas.tolerance import TOL
from compas_rhino.conversions import point_to_compas

from compas_singular.datastructures import CoarsePseudoQuadMesh
from compas_singular.datastructures import split_at_corners
from compas_singular.datastructures import split_at_junctions
from compas_singular.rhino.helpers import curve_points

from compas_singular.rhino.project import get_settings, resolve_spacing
from compas_singular.rhino.session import RhinoSession

#How finely a CURVED input curve is sampled, as a multiple of the background
#spacing. The same factor CMD_coarse_mesh samples its walls at, and for the same
#reason: these points ARE the layout's geometry from here on, so the final mesh
#follows the drawn curve only as closely as they do. A polyline is taken at its
#own vertices instead.
CURVE_SAMPLING_FACTOR = 0.25

#Most points a SMOOTH curve is sampled into, whatever the spacing says. The
#spacing is set for the domain in the settings, and a drawing at a different
#scale -- a 400-unit rectangle at a 0.5 spacing -- would otherwise put thousands
#of points on every arc, and the junction and validity checks compare segments
#pairwise. 256 holds a radius-150 arc to about 0.01 of its true shape.
MAX_CURVE_POINTS = 256

#How close a curve end must be to another curve to have LANDED on it. The weld
#resolution, so the same distance the constructor refuses a T-junction at. A
#snapped end (End, Near, Int) is exact; anything further is a near miss and is
#refused as dangling rather than guessed at.
ON_CURVE = 1e-3


# ----------------------------------------------------------------------------
# small geometry
# ----------------------------------------------------------------------------

def point_in_loop(point, loop):
    """Ray casting in XY. ``loop`` is closed, its last point not repeated."""
    x, y = point[0], point[1]
    hit = False
    n = len(loop)
    for i in range(n):
        ax, ay = loop[i][0], loop[i][1]
        bx, by = loop[(i + 1) % n][0], loop[(i + 1) % n][1]
        if (ay > y) != (by > y):
            t = (y - ay) / (by - ay)
            if x < ax + t * (bx - ax):
                hit = not hit
    return hit


def interior_point(loop):
    """A point strictly inside a closed loop, or ``None``.

    The centroid is right for anything convex and wrong for a crescent, whose
    centroid can sit outside it -- so fall back to midpoints between pairs of
    corners and take the first that is inside.
    """
    n = len(loop)
    centre = [sum(p[0] for p in loop) / n, sum(p[1] for p in loop) / n, 0.0]
    if point_in_loop(centre, loop):
        return centre
    for i in range(n):
        for j in range(i + 2, n):
            mid = [(loop[i][0] + loop[j][0]) / 2.0,
                   (loop[i][1] + loop[j][1]) / 2.0, 0.0]
            if point_in_loop(mid, loop):
                return mid
    return None


def midpoint(points):
    """The point halfway along a polyline, by arc length."""
    segments = list(zip(points, points[1:]))
    lengths = [((b[0] - a[0]) ** 2 + (b[1] - a[1]) ** 2) ** 0.5 for a, b in segments]
    half = sum(lengths) / 2.0
    for (a, b), length in zip(segments, lengths):
        if half <= length and length > 0.0:
            f = half / length
            return [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, 0.0]
        half -= length
    return list(points[-1])


def on_curve(guid, curve, point):
    """Is ``point`` within :data:`ON_CURVE` of the TRUE Rhino curve?"""
    t = rs.CurveClosestPoint(guid, point)
    if t is None:
        return False
    p = curve.PointAt(t)
    return ((p.X - point[0]) ** 2 + (p.Y - point[1]) ** 2) ** 0.5 <= ON_CURVE


def lies_on_wall(curve, walls):
    """Does this WHOLE curve lie along one of the walls?

    Tested before the curve is split at its corners, and that order is the whole
    point. ``CMD_coarse_mesh`` bakes a curved wall piece as a POLYLINE sampled
    from the wall; split at its corners first, it becomes a run of chords, each
    chord's middle a sagitta off the wall, and every one of them read as a
    division duplicating the wall -- measured on a disc, 48 "divisions" and a
    refusal. So the test uses only the points that are exact: a polyline's
    VERTICES, a smooth curve's samples.

    A two-point line is the one case where the vertices say too little -- a
    straight division with both ends on a round wall has both vertices on it --
    so its middle is tested as well. A straight wall piece passes that; a chord
    across a round wall does not.
    """
    ok, polyline = curve.TryGetPolyline()
    if ok:
        points = [[p.X, p.Y, 0.0] for p in polyline]
    else:
        points = [[q.X, q.Y, 0.0] for q in (curve.PointAt(t) for t in curve.DivideByCount(16, True))]
    if len(points) == 2:
        points.append(midpoint(points))
    return any(all(on_curve(guid, wall, q) for q in points) for guid, wall in walls)


# ----------------------------------------------------------------------------
# reading
# ----------------------------------------------------------------------------

def sample_smooth(guid, curve, sampling, ends):
    """A SMOOTH curve as one point list, with every end landing on it as a sample.

    Sampling puts points ON the curve, but the polyline between two samples is a
    chord, and a division snapped onto the real curve between them sits a
    sagitta off that chord -- on a radius-60 circle at 256 samples, 0.0045, four
    times :data:`ON_CURVE`. It would read as a near miss and be refused as
    dangling. So each end is located on the TRUE curve and, if it is on it,
    inserted at its own coordinates -- the polyline then passes exactly through
    it, and ``split_at_junctions`` cuts there.

    Crossings need no such help: two chord polylines crossing meet at one point
    that lies on both, whatever the sagitta.

    A closed curve comes back closed (first point repeated), for
    ``split_at_junctions`` to cut at its landings.
    """
    count = int(round(curve.GetLength() / max(sampling, 1e-6)))
    count = min(max(8, count), MAX_CURVE_POINTS)
    entries = []
    for t in curve.DivideByCount(count, True):
        p = curve.PointAt(t)
        entries.append((t, [p.X, p.Y, 0.0]))
    for end in ends:
        t = rs.CurveClosestPoint(guid, end)
        if t is None:
            continue
        p = curve.PointAt(t)
        if ((p.X - end[0]) ** 2 + (p.Y - end[1]) ** 2) ** 0.5 <= ON_CURVE:
            entries.append((t, [end[0], end[1], 0.0]))
    entries.sort(key=lambda entry: entry[0])
    points = [point for _t, point in entries]
    if curve.IsClosed and points:
        points.append(list(points[0]))
    return points


def read_curves(guids):
    """``(pieces, owner, smooth)``: polylines split at their corners, smooth curves kept.

    Polylines are exact and are read at once. A smooth curve waits: it has to be
    sampled with every curve END in the drawing to hand -- see
    :func:`sample_smooth` -- so it is returned for a second pass.

    **Not** ``curve_points``. That is a LOOP reader and it drops a closed curve's
    closing point: a rectangle drawn as one closed polyline came through as three
    sides, and the diagonals landing on its corners were reported as T-junctions.
    """
    pieces, owner, smooth = [], [], []
    for guid in guids:
        curve = rs.coercecurve(guid)
        if curve is None:
            continue
        ok, polyline = curve.TryGetPolyline()
        if ok:
            for piece in split_at_corners([[p.X, p.Y, 0.0] for p in polyline]):
                pieces.append(piece)
                owner.append(guid)
        else:
            smooth.append((guid, curve))
    return pieces, owner, smooth


def read_network(sampling):
    """Everything drawn, as one-curve-per-edge. ``(pieces, source, is_wall, marks)``.

    ``source[i]`` is the Rhino curve piece ``i`` came from -- the constructor
    numbers the pieces it is handed, and after the splits that is no longer the
    numbering on screen, so a refusal is mapped back through it and the curve
    selected. ``is_wall[i]`` says whether it came from ``Outer`` or ``Inner``.
    ``marks`` counts the wall pieces found on ``EdgeCurves`` and used as corners.
    """
    outer_guids = list(rs.ObjectsByLayer("Outer") or [])
    wall_guids = outer_guids + list(rs.ObjectsByLayer("Inner") or [])
    layer = rs.AddLayer(name="EdgeCurves", parent="Skeleton")
    division_guids = list(rs.ObjectsByLayer(layer) or [])
    if not division_guids and not outer_guids:
        print("nothing on 'Outer' or '{}' -- select the curves of the layout "
              "instead.".format(layer))
        division_guids = rs.GetObjects(message="Select the coarse layout's curves",
                                       filter=4, preselect=True) or []

    # A curve on EdgeCurves lying ALONG a wall is a wall piece baked by
    # CMD_coarse_mesh or CMD_edit_coarse_mesh. Kept, it would put two edges
    # between the same two corners; dropped outright, it would lose the corners
    # it marks on a SMOOTH wall. So its two ends become cut points on the wall and
    # the curve goes -- decided for the WHOLE curve, see :func:`lies_on_wall`.
    walls = [(guid, rs.coercecurve(guid)) for guid in wall_guids]
    walls = [(guid, curve) for guid, curve in walls if curve is not None]
    marks = []
    kept = []
    for guid in division_guids:
        curve = rs.coercecurve(guid)
        if curve is None:
            continue
        if walls and lies_on_wall(curve, walls):
            for p in (curve.PointAtStart, curve.PointAtEnd):
                marks.append([p.X, p.Y, 0.0])
        else:
            kept.append(guid)

    wall_pieces, wall_owner, wall_smooth = read_curves(wall_guids)
    divisions, owners, division_smooth = read_curves(kept)

    # Every point a smooth curve may be cut at has to be ON its sampled polyline
    # -- the marks included, or a corner on a round wall is missed by a sagitta.
    ends = [end for piece in wall_pieces + divisions for end in (piece[0], piece[-1])]
    for _guid, curve in wall_smooth + division_smooth:
        if not curve.IsClosed:
            for p in (curve.PointAtStart, curve.PointAtEnd):
                ends.append([p.X, p.Y, 0.0])
    ends.extend(marks)

    for smooth, pieces, owner in ((wall_smooth, wall_pieces, wall_owner),
                                  (division_smooth, divisions, owners)):
        for guid, curve in smooth:
            points = sample_smooth(guid, curve, sampling, ends)
            if len(points) >= 2:
                pieces.append(points)
                owner.append(guid)

    network = wall_pieces + divisions
    owner = wall_owner + owners
    flags = [True] * len(wall_pieces) + [False] * len(divisions)
    pieces, origin = split_at_junctions(network, tol=ON_CURVE, points=marks)
    return (pieces, [owner[i] for i in origin], [flags[i] for i in origin],
            len(marks) // 2)


def read_loops(layer, sampling):
    """The closed curves on a boundary layer, as point loops."""
    loops = []
    for guid in rs.ObjectsByLayer(layer) or []:
        curve = rs.coercecurve(guid)
        if curve is None or not curve.IsClosed:
            continue
        loop = curve_points(guid, sampling)
        if len(loop) >= 3:
            loops.append(loop)
    return loops


def read_poles():
    """The pole points already on the layer, if the layout had any."""
    layer = rs.AddLayer(name="Poles", parent="Skeleton")
    return [list(point_to_compas(cr.objects.find_object(guid).Geometry))
            for guid in rs.ObjectsByLayer(layer) or []]


# ----------------------------------------------------------------------------
# what the walls know
# ----------------------------------------------------------------------------

def outside_domain(pieces, is_wall, outer_loops, inner_loops):
    """``(index, why)`` for the first DIVISION outside the domain, or ``None``.

    Only with exactly one ``Outer`` loop -- with none the divisions are the
    boundary, and with several there is no single domain to test against. The
    test point is halfway along the piece: after splitting, a division that
    overshoots a wall is cut there, and the overshoot is a piece of its own
    lying wholly outside.
    """
    if len(outer_loops) != 1:
        return None
    for i, (piece, wall) in enumerate(zip(pieces, is_wall)):
        if wall:
            continue
        point = midpoint(piece)
        if not point_in_loop(point, outer_loops[0]):
            return i, "a division line runs OUTSIDE the Outer boundary -- trim it to the wall"
        if any(point_in_loop(point, loop) for loop in inner_loops):
            return i, "a division line runs INSIDE a hole -- a hole is not meshed"
    return None


def off_boundary(coarse, pieces, is_wall, has_outer):
    """``(index, why)`` for a piece on the wrong side of the layout's boundary, or ``None``.

    The layers say which edges are walls, so the finished layout can be held to
    it: its boundary must be exactly the walls. A wall INSIDE the layout means a
    hole was not recognised; a division ON its boundary means the patch beside it
    is missing. Neither should survive the checks before this one -- this is the
    statement of the rule, kept so that a case they miss is named, not baked.
    """
    index = {TOL.geometric_key(coarse.vertex_coordinates(v)): v for v in coarse.vertices()}
    for i, (piece, wall) in enumerate(zip(pieces, is_wall)):
        u = index.get(TOL.geometric_key(piece[0]))
        v = index.get(TOL.geometric_key(piece[-1]))
        if u is None or v is None or v not in coarse.halfedge[u]:
            continue
        boundary = coarse.is_edge_on_boundary((u, v))
        if wall and not boundary:
            return i, "this Outer/Inner curve ends up INSIDE the layout -- a hole was not recognised"
        if not wall and has_outer and boundary:
            return i, "this division line ends up on the EDGE of the layout -- the patch beside it is missing"
    return None


def select(guids):
    rs.UnselectAllObjects()
    rs.SelectObjects(list(dict.fromkeys(guids)))


def show_refusal(error, source):
    """Print a refusal, and select the Rhino curves it is about.

    The constructor's messages name edges by index ("polyline 3", "polylines 2
    and 5", "polyline(s) [4, 6]"). Every number in such a phrase is looked up in
    ``source`` and the curve it came from is selected, so the fix starts from the
    right object rather than from a coordinate. A message that names only a
    corner -- a dangling end -- carries its coordinates and selects nothing.
    """
    message = str(error)
    print("the drawing is not a coarse layout yet:")
    print("  {}".format(message))

    indices = set()
    for phrase in re.findall(r"polylines?(?:\(s\))?\s+([\d\s,\[\]and]+)", message):
        indices.update(int(n) for n in re.findall(r"\d+", phrase))
    guids = [source[i] for i in sorted(indices) if 0 <= i < len(source)]
    if guids:
        select(guids)
        print("  selected the curve(s) involved.")
    if "same corner" in message or "both run between" in message:
        print("  a closed curve is split into edges only where other curves land "
              "on it -- a hole or a round wall needs at least three.")
    if "separate pieces" in message:
        print("  if the selected curve is a hole, draw at least three division "
              "lines from it to the rest of the layout.")
    if "one curve only" in message:
        print("  if that end was meant to reach a curve, snap it onto the curve "
              "(End, Near or Int) -- a near miss is not a landing.")
    print("nothing was changed.")


# ----------------------------------------------------------------------------
# the command
# ----------------------------------------------------------------------------

def main():
    settings = get_settings()
    sampling = resolve_spacing(settings) * CURVE_SAMPLING_FACTOR

    # Everything is read BEFORE anything is written: the poles are read from the
    # layer this command rewrites.
    pieces, source, is_wall, marks = read_network(sampling)
    if not pieces:
        print("no curves to read.")
        return
    outer_loops = read_loops("Outer", sampling)
    inner_loops = read_loops("Inner", sampling)
    holes = [point for point in (interior_point(loop) for loop in inner_loops)
             if point is not None]
    poles = read_poles()

    print("read {} curve(s) as {} edge(s): {} wall, {} division; {} hole(s), "
          "{} pole point(s)".format(len(set(source)), len(pieces), sum(is_wall),
                                    len(pieces) - sum(is_wall), len(holes), len(poles)))
    if marks:
        print("  {} piece(s) of EdgeCurves lie along a wall -- read as corners on "
              "it, not as edges".format(marks))

    # ------------------------------------------------------------------
    # the layout
    # ------------------------------------------------------------------
    # A refusal is the useful output here, not a failure to be worked around --
    # it names and selects the curve, which is what makes the fix a few seconds
    # of drawing. Printed rather than raised so the command ends cleanly with
    # the document untouched.
    stray = outside_domain(pieces, is_wall, outer_loops, inner_loops)
    if stray is not None:
        print("the drawing is not a coarse layout yet:")
        print("  {}.".format(stray[1]))
        select([source[stray[0]]])
        print("  selected the curve involved.")
        print("nothing was changed.")
        return

    try:
        coarse = CoarsePseudoQuadMesh.from_coarse_polylines(
            pieces, poles=poles, holes=holes)
    except ValueError as error:
        show_refusal(error, source)
        return

    wrong = off_boundary(coarse, pieces, is_wall, bool(outer_loops))
    if wrong is not None:
        print("the drawing is not a coarse layout yet:")
        print("  {}.".format(wrong[1]))
        select([source[wrong[0]]])
        print("  selected the curve involved.")
        print("nothing was changed.")
        return

    sides = {}
    for fkey in coarse.faces():
        n = len(coarse.face_vertices(fkey))
        sides[n] = sides.get(n, 0) + 1
    print("layout: {} patch(es), {} corner(s), sides {}, {} strip(s)".format(
        coarse.number_of_faces(), coarse.number_of_vertices(), sides,
        coarse.number_of_strips()))
    if sides.get(3):
        print("  {} triangular patch(es), kept as pseudo-quads with a collapsed "
              "corner".format(sides[3]))

    # ``route`` rides in ``attributes`` and so travels with the session. It is
    # what a later step reads to know the layout was not generated: there is no
    # field behind it and no traced network to warp.
    coarse.attributes["route"] = "drawn"

    # ------------------------------------------------------------------
    # into the session -- and never onto Outer, Inner or EdgeCurves, the input
    # ------------------------------------------------------------------
    # ``shape_polylines`` is what the edge curves are rebuilt from, matching by
    # the geometric key of each edge's two ends. Every piece IS one coarse edge
    # with exactly those ends, so every edge finds its curve and no interior
    # edge densifies as a chord. The walls are included: a boundary edge's shape
    # comes from the Outer and Inner curves first, and these are its fallback.
    coarse.set_shape_polylines(pieces)

    # Recording draws the layout on Skeleton::Mesh, its poles on ::Poles and the
    # pieces on ::Polylines. EdgeCurves is left alone for a DRAWN layout: the
    # curves there are the ones just read.
    session = RhinoSession.current()
    session.coarse = coarse
    session.record("Read coarse mesh")
    print("layout stored in the session, and drawn")
    print("note: this replaces any previous layout -- densities (CMD_densities) and "
          "patterns (CMD_dense_pattern) set on it do not carry over.")
    print("next: CMD_densities, then CMD_quad_mesh.")


if __name__ == "__main__":
    main()
