"""**Give a coarse layout back the curvature the Rhino round trip took off it.**

A coarse edge is always a straight chord. That is not a defect of the layout --
the layout is a topological quad graph and has to stay one, because strips,
densities, poles and ``add_strip`` are all defined on it -- but it does mean the
SHAPE of each edge lives somewhere else. ``densification(edges_to_curves=...)``
is where it is handed back in: one polyline per coarse edge, sampled by
``Polyline.point_at`` at the strip's density. Without that argument every edge
densifies as the chord between its two corners, and on a curved domain the mesh
loses the area between chord and wall outright -- measured on the field route at
65% coverage on an ellipse and 74% on a disc.

**This module rebuilds that mapping from what is on the Rhino document alone.**
No cached decomposition, no session state, so it survives a save, a reopen and
running the steps out of order -- which matters because the coarse layout is
baked as a Rhino mesh, and a Rhino mesh is vertices and faces. It cannot carry a
per-edge polyline, and ``rs.AddMesh`` stores its vertices as ``Point3f``:
``compas_rhino``'s ``mesh_to_compas`` widens them back with ``float(vertex.X)``,
so a corner read back from the document is a SINGLE-PRECISION copy of the one
that was baked -- about 2e-6 off on a 20-unit plate. Every tolerance here is
sized against that, and nothing here matches by vertex key, which a bake
renumbers anyway.

Three ways an edge finds its curve, best first:

1. **wall arc** -- the edge is on the layout's boundary, so the piece of the
   domain wall between its two corners IS its shape. Re-derived from the input
   curve rather than looked up, which is the lesson of the ``guide_lines`` work:
   a corner that snapping, refining or editing moved is still ON the wall, so
   re-deriving works for any boundary edge in the final layout while matching a
   remembered curve only works for the ones nothing touched. It is also the more
   accurate of the two here -- the wall comes from the Rhino curve at whatever
   resolution is asked for, while a traced polyline is the wall as the
   background triangulation sampled it;
2. **traced polyline** -- the decomposition branch whose ends are this edge's
   ends, by geometric key. This is what carries interior separatrices, and it is
   why ``Skeleton::Polylines`` is worth baking;
3. **chord** -- a straight line, and the count is reported so a layout that
   quietly lost its curvature says so.

Which branch each edge took comes back as a tally, to be printed. ``chord`` is
the one that costs area -- but note that ``boundary`` falling to ``traced`` is
not free either, and the tally cannot show it: a traced polyline is the wall as
the background triangulation sampled it, so an edge that loses branch 1 to
branch 2 on a HOLE follows the triangulation's eight-or-so-sided sampling of the
circle instead of the wall. Measured on a radius-0.7 hole at 0.5 background
spacing, that is 0.043 units inside the true circle against 0.003 for the wall
arc -- a visible facet reported as ``chord: 0``.

**Holes are where branch 1 is hard, and the reason is scale.** A hole's loop is
short: one coarse edge can be half of it or more, its corners can be a tenth of
a mean coarse edge apart, and both of the things branch 1 has to decide -- which
way round the loop the edge goes, and whether another corner lies on the arc --
were being decided with the layout's GLOBAL mean edge length as the yardstick.
Both are handled at :meth:`BoundaryLoop.arcs` and :func:`_arc_is_one_edge`, and
the symptom they produced is worth recognising: several identical circles in one
domain, one of them meshed round and the others as polygons, because the arc
each hole needed happened to survive or not.

**Not covered: the warp.** ``FieldDecomposition.edges_to_curves`` has a fourth
branch that takes the separatrix an edge came FROM and warps it onto moved
corners, so a hand edit costs the nudge and not the curvature. That needs the
traced network and the live decomposition, so an interior edge whose corner was
dragged still falls to the chord here. Boundary edges do not -- branch 1 does
not care whether a corner moved, only that it is on a wall.

Nothing in this module touches ``rhinoscriptsyntax`` or ``Rhino``. That is
deliberate and is the same rule the lower half of ``edit_coarse`` keeps: this is
the half that can be wrong in ways a user cannot see, so it has to be runnable,
and testable, without Rhino open.

``compas_singular.rhino.coarse_curves`` re-exports this module's public names
so the ``CMD_`` commands keep importing from where they always have -- if a
search for ``coarse_curves.py`` lands you there instead, that file is the shim.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import annotations

from typing import TYPE_CHECKING

from compas_singular.geometry.polyline import project_on_polyline
from compas.geometry import distance_point_point
from compas.itertools import pairwise

from compas.tolerance import TOL

if TYPE_CHECKING:
    from compas_singular.datastructures import CoarseQuadMesh
    from compas_singular.datastructures import CoarsePseudoQuadMesh


__all__ = [
    'coarse_edges_to_curves',
    'snap_corners_to_walls',
    'mean_edge_length',
    'BoundaryLoop',
]


def mean_edge_length(coarse: "CoarseQuadMesh | CoarsePseudoQuadMesh") -> float:
    """Average coarse edge length -- the scale every tolerance here is in.

    Tolerances are relative to the layout rather than absolute because the same
    workflow runs on a 2 m detail and a 200 m plate, and an absolute tolerance
    that is right for one silently does nothing on the other.
    """
    lengths = [coarse.edge_length(edge) for edge in coarse.edges()]
    lengths = [length for length in lengths if length > 0.0]
    return sum(lengths) / len(lengths) if lengths else 1.0


def _clean(points: list[list[float]], tol: float = 1e-9) -> list[list[float]]:
    """Drop consecutive duplicates -- ``Polyline.point_at`` divides by them."""
    out = [list(points[0])]
    for point in points[1:]:
        if distance_point_point(out[-1], point) > tol:
            out.append(list(point))
    return out


class BoundaryLoop(object):
    """A closed domain wall, arc-length parametrised, with its seam unrolled.

    The ring is stored DOUBLED rather than indexed modulo. An arc that crosses
    the point where the loop closes is then a plain slice of a longer list, and
    there is no modular index arithmetic to get wrong -- the same trick
    ``guide_lines.boundary_arc_between`` uses.

    Parameters
    ----------
    points : list[[x, y, z]]
        The loop, closed. A repeated closing point is dropped: the loop is
        stored open and closed by the doubling, so keeping it would put a
        zero-length segment in the middle of every arc that crosses the seam.
    """

    def __init__(self, points: list[list[float]]) -> None:
        pts = [[float(p[0]), float(p[1]), 0.0] for p in points]
        while len(pts) > 1 and distance_point_point(pts[0], pts[-1]) < 1e-9:
            pts = pts[:-1]
        self.points = pts
        # two full turns, so a half-turn arc starting anywhere is contiguous
        self.ring = pts + pts + pts[:1]
        cum = [0.0]
        for a, b in pairwise(self.ring):
            cum.append(cum[-1] + distance_point_point(a, b))
        self.cum = cum
        self.total = cum[len(pts)]

    def project(self, point: list[float]) -> tuple[float, float, list[float] | None]:
        """``(distance, arclength, closest point)`` for a point near the loop.

        Searched over ONE turn -- the second turn is the same loop and would
        only ever tie.
        """
        best = (float('inf'), 0.0, None)
        for i in range(len(self.points)):
            a, b = self.ring[i], self.ring[i + 1]
            abx, aby = b[0] - a[0], b[1] - a[1]
            length2 = abx * abx + aby * aby
            if length2 == 0.0:
                continue
            t = ((point[0] - a[0]) * abx + (point[1] - a[1]) * aby) / length2
            t = max(0.0, min(1.0, t))
            q = [a[0] + abx * t, a[1] + aby * t, 0.0]
            d = distance_point_point(point, q)
            if d < best[0]:
                best = (d, self.cum[i] + t * (self.cum[i + 1] - self.cum[i]), q)
        return best

    def _point_at(self, s: float) -> list[float]:
        """The loop point at arclength ``s`` along the doubled ring."""
        for i in range(len(self.cum) - 1):
            if self.cum[i] <= s <= self.cum[i + 1]:
                seg = self.cum[i + 1] - self.cum[i]
                t = 0.0 if seg == 0.0 else (s - self.cum[i]) / seg
                a, b = self.ring[i], self.ring[i + 1]
                return [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, 0.0]
        return list(self.ring[-1])

    def _span(self, s0: float, s1: float, flip: bool, pa: list[float], pb: list[float]) -> list[list[float]] | None:
        """One way round, as a point list running ``pa`` -> ``pb``.

        **The loop's own points are kept, not resampled.** They are the input
        curve at the resolution it was read at, so passing them straight through
        is both the most faithful thing to do and gives each arc a point count
        proportional to its length for free.

        The two ends are then overwritten with ``pa`` and ``pb`` exactly. They
        are near the wall but not necessarily on it -- a corner sits whereever
        the background triangulation put it, which for a curved wall is on a
        chord of it -- and ``densification`` builds a Coons patch per face, so
        two patches sharing an edge must be handed curves with IDENTICAL ends or
        ``meshes_join_and_weld`` fails to weld there and the dense mesh comes
        back as disconnected fragments rather than an error.
        """
        inner = [list(self.ring[i]) for i in range(len(self.ring))
                 if s0 + 1e-9 < self.cum[i] < s1 - 1e-9]

        # An arc SHORTER than the loop's own sampling step straddles no loop
        # point, so it would come back as two points -- the chord -- and be
        # rejected for having fewer than three. Sample it instead. The shape
        # barely matters at that length (0.002 sagitta on a 0.1 arc of a
        # radius-0.7 hole), but "no arc found" is a lie the tally then repeats.
        if len(inner) < 3:
            inner = [self._point_at(s0 + (s1 - s0) * k / 4.0) for k in (1, 2, 3)]

        head, tail = (pb, pa) if flip else (pa, pb)
        arc = _clean([list(head)] + inner + [list(tail)])
        if flip:
            arc.reverse()
        if len(arc) < 2:
            return None
        arc[0], arc[-1] = [pa[0], pa[1], 0.0], [pb[0], pb[1], 0.0]
        return arc

    def arcs(self, pa: list[float], pb: list[float]) -> list[list[list[float]]]:
        """BOTH ways round the loop from ``pa`` to ``pb``, shorter first.

        Two, not one, because **length does not decide which way round is the
        edge**. The usual case says it does -- a coarse edge spans a fraction of
        the wall, so the shorter arc is the edge -- and that holds for the outer
        boundary of a plate, where one edge is a small part of a long loop. It
        fails on a HOLE, whose loop is short and whose ring may be as few as
        three corners: an edge that spans more than half the circle then has the
        shorter arc running the wrong way, through the ring's other corners.
        Measured on a hole with corners at 0, 200 and 280 degrees: the 0 -> 200
        edge was handed the 160-degree complement, :func:`_arc_is_one_edge`
        rejected it -- correctly, it runs through the corner at 280 -- and the
        edge densified as a chord while its two neighbours came out round. That
        is the whole "one circle is meshed correctly, the others are polygons"
        symptom, and there is nothing wrong with the geometry: the right arc was
        the one never offered.

        Returning both and letting :func:`_wall_arc` keep the first that no
        other corner lies on puts the decision where the answer actually is.
        Shorter stays first, so any edge that was already resolved resolves the
        same way.
        """
        _da, sa, _qa = self.project(pa)
        _db, sb, _qb = self.project(pb)

        forward = (sb - sa) % self.total
        backward = self.total - forward
        if forward <= 1e-12 or backward <= 1e-12:
            return []

        ahead = (sa, sa + forward, False)
        behind = (sb, sb + backward, True)
        order = (ahead, behind) if forward <= backward else (behind, ahead)
        return [arc for arc in (self._span(s0, s1, flip, pa, pb)
                                for s0, s1, flip in order) if arc is not None]

    def arc(self, pa: list[float], pb: list[float]) -> list[list[float]] | None:
        """The shorter way round the loop from ``pa`` to ``pb``, or ``None``."""
        found = self.arcs(pa, pb)
        return found[0] if found else None


def _arc_is_one_edge(arc: list[list[float]] | None, pa: list[float], pb: list[float], corners: list[list[float]], tol: float) -> bool:
    """Is this arc ONE edge of the layout, or several?

    Branch 1 assumes an edge with both ends on a wall IS the piece of wall
    between them. That is true of an ordinary boundary edge and false of a
    **chord across a corner of the domain** -- an edge whose two ends happen to
    sit on the same loop but whose patch lies inside it.

    The two are told apart by what the arc runs THROUGH. A genuine piece of wall
    runs between two adjacent corners of the layout, so no other corner lies on
    it. An arc that passes through one has gone round a corner of the domain and
    come back, and taking it makes the patch retrace its own other sides.

    ``corners`` must be the layout's BOUNDARY corners only. What the guard
    detects is an arc spanning several boundary EDGES, and only a boundary
    vertex can be one of their ends -- so an interior vertex is not evidence of
    anything, and letting one veto is a plain category error. It is not a rare
    one either: a hole near another feature puts interior corners just outside
    its wall, and one of those, 0.076 from a radius-0.7 hole against a ``tol``
    of 0.119, was what left that hole with two chord edges slicing 21% of the
    way across it while its two identical siblings came out round.

    Ported from ``FieldDecomposition._arc_is_one_edge``, where it is not
    hypothetical: on ``22_force_lines``' arch it was the difference between a
    195-face field mesh and two zero-area quads that got the whole run rejected.

    **THROUGH is measured along the arc, not across it.** A corner that is near
    the arc because it is near one of its ENDS is not one the arc runs through,
    and vetoing on distance alone made that mistake constantly on small holes:
    ``tol`` is a fraction of the mean coarse edge length of the WHOLE layout,
    which the big patches set, while the thing it has to resolve is the corner
    spacing on one small circle. Measured, a radius-0.7 hole in a layout of mean
    edge 2.4 (so ``tol`` = 0.12 against a 4.4-long loop): any two corners of the
    ring closer than 12 degrees vetoed BOTH of their outward neighbours, two of
    the ring's four edges, and the hole densified as a polygon while an
    identical hole a few units away came out round. Requiring the offending
    corner to sit a margin clear of both ends restores those without weakening
    the guard -- a chord across a domain corner has that corner in its middle.
    """
    if not arc or len(arc) < 3:
        return bool(arc)

    xs = [p[0] for p in arc]
    ys = [p[1] for p in arc]
    lo_x, hi_x = min(xs) - tol, max(xs) + tol
    lo_y, hi_y = min(ys) - tol, max(ys) + tol

    for point in corners:
        # the arc's own endpoints are not "another corner"
        if (distance_point_point(point, pa) < 1e-6
                or distance_point_point(point, pb) < 1e-6):
            continue
        if not (lo_x <= point[0] <= hi_x and lo_y <= point[1] <= hi_y):
            continue
        d, s, total = project_on_polyline(point, arc)
        if d >= tol:
            continue
        # capped so the margin can never swallow the arc it is protecting
        margin = min(tol, 0.25 * total)
        if s <= margin or s >= total - margin:
            continue
        return False
    return True


def _wall_arc(loops: list[BoundaryLoop], pa: list[float], pb: list[float], corners: list[list[float]], wall_tol: float, arc_tol: float) -> list[list[float]] | None:
    """The piece of wall between two corners, if they really are on one.

    Every loop the two corners could be on is tried, nearest first, and both
    ways round each -- see :meth:`BoundaryLoop.arcs`. The first arc that no
    other corner of the layout lies on wins.

    Rejecting used to end the search: the nearest loop's shorter arc was the
    only candidate, so an edge that spanned more than half of a hole, or one
    whose neighbour on the ring sat inside ``arc_tol``, fell straight to a
    chord. Both are hole-shaped problems -- a hole's loop is short enough for
    one edge to be most of it -- which is why they showed up as some circles
    round and others faceted in the same mesh.
    """
    candidates = []
    for loop in loops:
        d = max(loop.project(pa)[0], loop.project(pb)[0])
        if d <= wall_tol:
            candidates.append((d, loop))
    candidates.sort(key=lambda pair: pair[0])

    for _d, loop in candidates:
        for arc in loop.arcs(pa, pb):
            # Two points IS the chord. Falling through lets the traced branch
            # have a go.
            if len(arc) < 3:
                continue
            if _arc_is_one_edge(arc, pa, pb, corners, arc_tol):
                return arc
    return None


def snap_corners_to_walls(coarse: "CoarseQuadMesh | CoarsePseudoQuadMesh", loops: "list[list[list[float]]] | tuple[list[list[float]], ...]" = (), wall_tol: float | None = None) -> tuple[int, float]:
    """Put the layout's boundary corners ON the wall. Returns ``(moved, worst)``.

    **This MUTATES ``coarse``**, which is why it is a separate call and not part
    of :func:`coarse_edges_to_curves` -- that one is read-only, and a function
    that silently moved the layout it was asked to describe would be a trap.
    Call it BEFORE building the mapping; the arcs are anchored on corner
    positions, so afterwards is too late.

    A boundary corner of the layout is, by definition, a point of the domain
    boundary -- but nothing in the pipeline puts it there. It is placed by the
    background triangulation, and on a curved wall that means on a CHORD of it.
    The docstring of this module used to record the consequence as a known loss
    ("the dense boundary passes through them and bows out to the curve
    between"), on the estimate that a corner is a sagitta off. That estimate is
    right for a well-sampled wall and wrong exactly where it matters: measured
    on a deltoid plate with three radius-0.7 holes near the wall, the worst
    corner was **0.106 off its own hole, 15.1% of the radius**, and it was not
    a sampling artefact -- it held at 16, 32, 48 and 64 points per circle. With
    every edge correctly given its wall arc (``wall_missed`` 0), those corners
    were then the ONLY points of the dense hole ring not on the circle, because
    an arc is pinned to its endpoints.

    Only vertices on the layout boundary move, only onto the loop they are
    already nearest to, and only if the move is under ``wall_tol`` -- a corner
    further off than that is not a corner that lost its wall, it is a corner
    somewhere else, and moving it would be a guess.

    Parameters
    ----------
    coarse : :class:`CoarseQuadMesh` or :class:`CoarsePseudoQuadMesh`
        The layout. **Modified in place.**
    loops : list[list[[x, y, z]]], optional
        The domain walls, outer first, each closed -- the same ones handed to
        :func:`coarse_edges_to_curves`.
    wall_tol : float, optional
        Largest move allowed. Defaults to a quarter of the mean coarse edge
        length, matching :func:`coarse_edges_to_curves`.

    Returns
    -------
    tuple[int, float]
        How many corners moved, and the largest distance any of them moved.
    """
    walls = [BoundaryLoop(loop) for loop in (loops or []) if len(loop) >= 3]
    if not walls:
        return 0, 0.0
    if wall_tol is None:
        wall_tol = 0.25 * mean_edge_length(coarse)

    moved, worst = 0, 0.0
    for vertex in coarse.vertices():
        if not coarse.is_vertex_on_boundary(vertex):
            continue
        point = coarse.vertex_coordinates(vertex)
        best = None
        for wall in walls:
            d, _s, q = wall.project(point)
            if best is None or d < best[0]:
                best = (d, q)
        if best is None or best[1] is None:
            continue
        if best[0] <= 1e-12 or best[0] > wall_tol:
            continue
        coarse.vertex_attributes(vertex, 'xyz', best[1])
        moved += 1
        worst = max(worst, best[0])
    return moved, worst


def coarse_edges_to_curves(coarse: "CoarseQuadMesh | CoarsePseudoQuadMesh", loops: "list[list[list[float]]] | tuple[list[list[float]], ...]" = (),
                           polylines: "list[list[list[float]]] | tuple[list[list[float]], ...]" = (), wall_tol: float | None = None,
                           precision: int | None = None) -> tuple[dict[tuple[int, int], list[list[float]]], dict[str, int]]:
    """``({(u, v): polyline}, tally)`` -- the shape of every coarse edge.

    Hand the dict straight to ``densification(edges_to_curves=...)``.

    Parameters
    ----------
    coarse : :class:`CoarseQuadMesh` or :class:`CoarsePseudoQuadMesh`
        The layout to densify. Read only.
    loops : list[list[[x, y, z]]], optional
        The domain walls, outer first, each closed. Sample them from the Rhino
        curves at whatever resolution the final mesh should follow -- these
        points ARE the boundary from here on, so an arc handed over as two
        points is a chord.
    polylines : list[list[[x, y, z]]], optional
        The decomposition's branches, as baked on ``Skeleton::Polylines``.
        Matched by the geometric key of their two ends.
    wall_tol : float, optional
        How far from a wall a boundary corner may sit and still be taken as on
        it. Defaults to a quarter of the mean coarse edge length, which is
        generous on purpose: a corner placed by the background triangulation
        lies on a CHORD of the wall, so on a tightly curved boundary it is a
        sagitta off it -- 0.03 for a 0.5 background on a 1-unit radius -- and
        the branch is already gated on the edge being topologically on the
        layout boundary and on :func:`_arc_is_one_edge`.
    precision : int, optional
        Decimals for the geometric key of the traced branch. The default is
        COMPAS's own (3), which is far coarser than the single-precision error a
        Rhino mesh round trip introduces and so is safe.

    Returns
    -------
    tuple[dict, dict]
        The mapping, complete for every edge of ``coarse``, and
        ``{'boundary': n, 'traced': n, 'chord': n, 'wall_missed': n}``.

        ``wall_missed`` counts the edges that are ON the layout boundary and did
        NOT get branch 1 -- they fell to a traced polyline or to a chord. It is a
        subset of the other counts, not a fourth category, and it exists because
        ``chord`` alone cannot see the quiet half of this failure: a boundary
        edge that loses its wall arc to a traced branch still gets a curve, so
        the tally reads clean while the edge follows the background
        triangulation's sampling of the wall rather than the wall. On a hole that
        is the difference between round and faintly faceted. Expect 0 on a
        generated layout; anything else is a corner sitting off its wall.

    Notes
    -----
    The mapping is COMPLETE by construction, and has to be: once
    ``densification`` is given a mapping at all it looks every edge up in it and
    has no per-edge straight-chord branch, so a missing key raises ``KeyError``
    several frames away rather than falling back.
    """
    mean = mean_edge_length(coarse)
    if wall_tol is None:
        wall_tol = 0.25 * mean
    arc_tol = 0.05 * mean

    walls = [BoundaryLoop(loop) for loop in (loops or []) if len(loop) >= 3]

    lookup = {}
    traced = []
    for polyline in (polylines or []):
        points = [[float(p[0]), float(p[1]), float(p[2])] for p in polyline]
        if len(points) < 2:
            continue
        ka = TOL.geometric_key(points[0], precision)
        kb = TOL.geometric_key(points[-1], precision)
        # A two-point branch is kept even though it IS the chord, because the
        # tally is the point: 'chord' has to mean "no shape was found for this
        # edge", not "the shape found happens to be straight". Most interior
        # branches of a skeleton decomposition are straight, and counting them
        # as losses buries the boundary edge that is a real one.
        for key, curve in (((ka, kb), points), ((kb, ka), list(reversed(points)))):
            if len(lookup.get(key, ())) < len(curve):
                lookup[key] = curve
        traced.append(points)

    def _nearest_branch(pa: list[float], pb: list[float]) -> list[list[float]] | None:
        """The branch whose ends are these ends, allowing for a MOVED corner.

        The geometric key above is exact, and exact stops being right the moment
        a corner moves: :func:`snap_corners_to_walls` shifts boundary corners
        onto the wall, and a Rhino round trip shifts every corner by the
        single-precision error. Both leave the branch in place and unfindable.
        Measured on a deltoid plate with three holes, snapping alone moved 13
        interior edges from ``traced`` to ``chord`` -- each of them straight, so
        the mesh did not change, but ``chord`` then meant "the shape found was
        straight" instead of "no shape was found", which is the one thing this
        module's tally is for.

        Scored on the SUM of the two end distances so a branch has to match at
        both ends, and bounded by the same tolerance the snap respects.
        """
        best, found = wall_tol, None
        for points in traced:
            for curve in (points, list(reversed(points))):
                score = (distance_point_point(curve[0], pa)
                         + distance_point_point(curve[-1], pb))
                if score < best:
                    best, found = score, curve
        return found

    # BOUNDARY corners only -- see _arc_is_one_edge. An interior corner cannot
    # end a boundary edge, so it is not evidence that an arc spans several.
    corners = [coarse.vertex_coordinates(w) for w in coarse.vertices()
               if coarse.is_vertex_on_boundary(w)]
    tally = {'boundary': 0, 'traced': 0, 'chord': 0, 'wall_missed': 0}
    out = {}

    for u, v in coarse.edges():
        pa = coarse.vertex_coordinates(u)
        pb = coarse.vertex_coordinates(v)

        curve, how = None, 'chord'

        on_wall = bool(walls) and coarse.is_edge_on_boundary((u, v))
        if on_wall:
            curve = _wall_arc(walls, pa, pb, corners, wall_tol, arc_tol)
            if curve is not None:
                how = 'boundary'
            else:
                tally['wall_missed'] += 1

        if curve is None:
            curve = lookup.get((TOL.geometric_key(pa, precision),
                                TOL.geometric_key(pb, precision)))
            if curve is None:
                curve = _nearest_branch(pa, pb)
            if curve is not None:
                curve = [list(point) for point in curve]
                how = 'traced'

        if curve is None or len(curve) < 2:
            curve = [list(pa), list(pb)]
            how = 'chord'
        else:
            # Both branches must end EXACTLY on the corners or the Coons patches
            # of the two faces sharing this edge disagree and will not weld.
            curve[0], curve[-1] = list(pa), list(pb)

        tally[how] += 1
        out[u, v] = curve

    return out, tally
