from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from __future__ import annotations

from math import ceil

from compas.geometry import Polyline
from compas.geometry import length_vector
from compas.geometry import cross_vectors
from compas.geometry import distance_point_point
from compas.geometry import subtract_vectors
from compas.itertools import pairwise


__all__ = [
    'Polyline',
    'bounding_box_diagonal',
    'closest_on_polyline',
    'discretise_boundary',
    'discretise_line',
    'distance_to_loop',
    'distance_to_polyline',
    'project_on_polyline',
    'resample_loop',
]


def closest_on_polyline(point: list[float], points: list[list[float]], closed: bool = False) -> tuple[int, float, list[float], float]:
    """``(segment index, t, closest point, distance)`` on a polyline.

    The richest of this module's projection results: distance, the point itself,
    and WHICH segment it landed on. A caller wanting only one of those -- a
    distance, a tangent, an arc length -- derives it from here rather than
    running the segment loop again.

    ``t`` is clamped to ``[0, 1]``, so the result is the closest point on the
    polyline itself and not merely its nearest vertex.

    The projection is in 3D, which for a polyline lying in ``z = 0`` reduces
    EXACTLY to the planar one -- ``abz`` is then zero, so it contributes nothing
    to ``length2``, nothing to ``t``, and leaves the closest point at ``z = 0``
    (verified bit-identical over 30 000 random cases, open and closed). That is
    why there is no planar flag: the flat case is not a special case. A polyline
    at a constant ``z != 0`` does differ -- the result lands on the polyline
    rather than on its shadow, which is the intended reading.

    Parameters
    ----------
    point : [x, y, z]
    points : list[[x, y, z]]
        Three components each; unlike the rest of this module, ``z`` is read.
    closed : bool, optional
        Add the segment from the last point back to the first.

    Returns
    -------
    tuple
        ``(index, t, [x, y, z], distance)``. Distance is ``inf``, and the other
        three meaningless, when no segment has non-zero length.
    """
    # compas has closest_point_on_polyline_xy, but it returns the point ALONE --
    # no index and no t, so a tangent or an arc length cannot be recovered from
    # it -- and it is ~2x slower, sorting a per-segment cloud instead of keeping
    # a running min. Its distance is bit-identical to this one where it applies.
    # Polyline.tangent_at_point and Polyline.parameter_at look like the missing
    # half, but both search back for the point against INFINITE segment lines:
    # the first raises when it misses, the second returns a silently wrong
    # parameter on any outline with two non-contiguous collinear runs -- an
    # ordinary T- or L-plate.
    chain = list(points) + list(points[:1]) if closed else points
    best = (0, 0.0, list(points[0]), float('inf'))
    for i, (a, b) in enumerate(pairwise(chain)):
        abx, aby, abz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
        length2 = abx * abx + aby * aby + abz * abz
        if length2 == 0.0:
            continue
        t = ((point[0] - a[0]) * abx + (point[1] - a[1]) * aby
             + (point[2] - a[2]) * abz) / length2
        t = max(0.0, min(1.0, t))
        q = [a[0] + abx * t, a[1] + aby * t, a[2] + abz * t]
        d = distance_point_point(point, q)
        if d < best[3]:
            best = (i, t, q, d)
    return best


def project_on_polyline(point: list[float], points: list[list[float]]) -> tuple[float, float, float]:
    """``(distance, arclength, total length)`` for a point near an OPEN polyline.

    Open, not closed: joining last to first would invent a segment an arc does
    not have. :func:`distance_to_loop` is the closed version.

    The arclength is what tells "this corner sits at the END of the arc" apart
    from "this arc runs THROUGH that corner" -- a distinction the distance alone
    cannot make, and the reason this returns more than a distance.

    Parameters
    ----------
    point : [x, y, z]
    points : list[[x, y, z]]
        The open polyline.

    Returns
    -------
    tuple
        ``(distance, arclength of the closest point, total length)``.
    """
    best = (float('inf'), 0.0)
    travelled = 0.0
    for a, b in pairwise(points):
        abx, aby = b[0] - a[0], b[1] - a[1]
        length = (abx * abx + aby * aby) ** 0.5
        if length == 0.0:
            continue
        t = ((point[0] - a[0]) * abx + (point[1] - a[1]) * aby) / (length * length)
        t = max(0.0, min(1.0, t))
        q = [a[0] + abx * t, a[1] + aby * t, 0.0]
        d = distance_point_point(point, q)
        if d < best[0]:
            best = (d, travelled + t * length)
        travelled += length
    return best[0], best[1], travelled


def distance_to_polyline(point: list[float], points: list[list[float]]) -> float:
    """Shortest distance from a point to an OPEN polyline."""
    return project_on_polyline(point, points)[0]


def distance_to_loop(p: list[float], loop: list[list[float]]) -> float:
    """Shortest distance from a point to a CLOSED polyline, in XY.

    Distance to the polyline itself, not to its nearest vertex: ``t`` is clamped
    to ``[0, 1]`` so the foot of the perpendicular is used whenever it falls
    inside a segment, and only otherwise an endpoint.

    The loop is stored open -- the closing segment is added here.

    Parameters
    ----------
    p : [x, y, z]
    loop : list[[x, y, z]]
        Open, the last point is not the first.

    Returns
    -------
    float
        ``inf`` for a loop with no segment of non-zero length.
    """
    # Same answer from compas (import ``closest_point_on_polyline_xy``), but ~3.3x slower -- it sorts a per-segment cloud instead of keeping a running min:
    # return distance_point_point(p, closest_point_on_polyline_xy(p, list(loop) + list(loop[:1])))
    best = float('inf')
    for a, b in pairwise(list(loop) + list(loop[:1])):
        abx, aby = b[0] - a[0], b[1] - a[1]
        length2 = abx * abx + aby * aby
        if length2 == 0.0:
            continue
        t = ((p[0] - a[0]) * abx + (p[1] - a[1]) * aby) / length2
        t = max(0.0, min(1.0, t))
        q = [a[0] + abx * t, a[1] + aby * t, 0.0]
        best = min(best, distance_point_point(p, q))
    return best


def bounding_box_diagonal(*loops: list[list[float]]) -> float:
    """The diagonal ``D`` of the TOTAL bounding box of every loop given, in XY.

    The scale the thesis measures a discretisation against -- see
    :func:`discretise_boundary`. Total, so a hole that pokes outside its outer
    boundary still counts; for a valid domain the inners are inside the outer
    and this is the outer's own bounding box.

    Parameters
    ----------
    *loops : list[[x, y, z]]
        Any number of point lists. Empty ones are ignored.

    Returns
    -------
    float
        ``0.0`` if no points were given.
    """
    points = [p for loop in loops for p in (loop or [])]
    if not points:
        return 0.0
    xs = [float(p[0]) for p in points]
    ys = [float(p[1]) for p in points]
    return ((max(xs) - min(xs)) ** 2 + (max(ys) - min(ys)) ** 2) ** 0.5


def _clean_loop(points: list[list[float]]) -> list[list[float]]:
    """A loop as XY floats, with consecutive duplicates and a closing point dropped."""
    pts = [[float(p[0]), float(p[1]), float(p[2]) if len(p) > 2 else 0.0]
           for p in points]
    if not pts:
        return []
    out = [pts[0]]
    for p in pts[1:]:
        if distance_point_point(p, out[-1]) > 1e-9:
            out.append(p)
    if len(out) > 1 and distance_point_point(out[0], out[-1]) < 1e-9:
        out.pop()
    return out


def _subdivide(loop: list[list[float]], spacing: float) -> list[list[float]]:
    """Subdivide EACH segment of a closed loop so that none exceeds ``spacing``."""
    out = []
    for a, b in pairwise(loop + loop[:1]):
        n = max(1, int(ceil(distance_point_point(a, b) / spacing)))
        # the end point is dropped -- the next segment starts there, and the
        # last segment's end is the loop's first point
        for i in range(n):
            t = float(i) / n
            out.append([a[k] + (b[k] - a[k]) * t for k in range(3)])
    return out


def _subdivide_chain(chain: list[list[float]], spacing: float) -> list[list[float]]:
    """Subdivide EACH segment of an OPEN chain so that none exceeds ``spacing``.

    The loop version closes the chain and drops every segment's end point. An
    open chain must keep its last point: a curve feature's extremity decides
    whether it lands on a wall, and moving it is exactly what must not happen.
    """
    out = []
    for a, b in pairwise(chain):
        n = max(1, int(ceil(distance_point_point(a, b) / spacing)))
        for i in range(n):
            t = float(i) / n
            out.append([a[k] + (b[k] - a[k]) * t for k in range(3)])
    out.append([float(c) for c in chain[-1]])
    return out


def _clean_chain(points: list[list[float]]) -> list[list[float]]:
    """An open chain as XY floats, consecutive duplicates dropped, ends kept."""
    pts = [[float(p[0]), float(p[1]), float(p[2]) if len(p) > 2 else 0.0]
           for p in points]
    if not pts:
        return []
    out = [pts[0]]
    for p in pts[1:]:
        if distance_point_point(p, out[-1]) > 1e-9:
            out.append(p)
    return out


def discretise_line(line: list[list[float]], spacing: float | None) -> list[list[float]]:
    """Discretise ONE open curve at ``spacing``, per thesis eq. 4.1.

    Deliberately separate from :func:`discretise_boundary` rather than a flag on
    it. A boundary is a closed loop whose closing segment is subdivided like any
    other and whose first point is also its last; a line is open, and its two
    EXTREMITIES must survive untouched -- whether a curve feature's end lands on
    a wall decides whether the layout gets a node there, so moving it is exactly
    what must not happen. The two conventions do not belong in one function.

    Pass the spacing the walls were discretised at. A feature is cut into the
    Delaunay along its own segments, so a segment longer than the surrounding
    sampling is not a Delaunay edge at all and the cut does not happen there.

    Measured on the Fig 4.17 diagonal against a wall sampled at 0.2: a 2-point
    guide -- what a straight polyline drawn in Rhino gives, since
    ``curve_points`` takes a polyline at its own vertices -- leaves its single
    segment missing from the triangulation, and the layout comes back as though
    there were no feature at all.

    Parameters
    ----------
    line : list[[x, y, z]]
        The open curve.
    spacing : float
        Maximum segment length. ``None`` returns the curve unchanged apart from
        dropping consecutive duplicates.

    Returns
    -------
    list[[x, y, z]]

    """
    points = _clean_chain(line)
    if spacing is None or len(points) < 2:
        return points
    spacing = float(spacing)
    if spacing <= 0.0:
        raise ValueError('spacing must be positive, got {}'.format(spacing))
    return _subdivide_chain(points, spacing)


def discretise_boundary(
    outer: list[list[float]],
    inners: list[list[list[float]]] | None = None,
    alpha: float | None = 0.04,
    d_min: int | None = 5,
    spacing: float | None = None,
) -> tuple[list[list[float]], list[list[list[float]]]]:
    """Discretise closed boundary loops, per Oval's thesis eq. 4.1.

    The number of points of each curve ``i`` is

        d_i = max(d_scale, d_min),  d_scale = ceil(l_i / (alpha * D))

    with ``l_i`` the length of the curve, ``D`` a scale of the surface -- here
    the diagonal of the total bounding box -- and ``alpha`` a percentage of it.
    The thesis reports good results for ``d_min`` between 5 and 10 and ``alpha``
    between 0.01 and 0.05.

    **Every input point survives.** Each SEGMENT is subdivided on its own, which
    is what keeps a corner a corner; this is deliberately not compas's
    ``Polyline.divide_by_length``, which divides by arclength from one end and
    walks straight past a corner, so a resampled square comes back with its
    corners rounded off the point list.

    Both front ends discretise their walls with this. The skeleton route reads
    the medial axis off a Delaunay triangulation of the boundary POINTS, so a
    four-point square gives a two-triangle mesh with no interior structure and a
    caricature decomposition. The field route wants the same sampling for a
    different reason: a sampled arc turning more than 45 degrees between two
    points reads as a CORNER to a cross field.

    ``ceil``, not ``round``: ``d_scale`` is an upper integer value, so
    ``target_length`` is a bound and not an average. Rounding to nearest left
    segments up to 1.5x the target -- measured 1.12x on a 28-point disc at 0.5.

    Parameters
    ----------
    outer : list[[x, y, z]]
        The outer loop, stored OPEN -- the last point is not the first. A
        repeated closing point is dropped, so a closed list is accepted too.
    inners : list[list[[x, y, z]]], optional
        The holes, same convention. They are discretised against the SAME ``D``
        as the outer loop, which is the point of taking them here rather than
        loop by loop: a small hole measured against its own bounding box would
        come back far denser than the wall beside it.
    alpha : float, optional
        Fraction of ``D`` to use as the target length when ``target_length`` is
        ``None``. ``None`` for both returns the loops unchanged -- the explicit
        opt-out, for a caller who has already sampled its walls.
    d_min : int, optional
        Fewest points per loop, whatever the target length says. ``None`` or
        ``0`` disables the floor.
    spacing : float, optional
        Maximum segment length, given directly. Overrides ``alpha`` when not
        ``None``.

    Returns
    -------
    tuple
        ``(outer, inners)``, both open, with the closing segment subdivided like
        any other.
    """
    outer = _clean_loop(outer)
    inners = [_clean_loop(loop) for loop in (inners or [])]

    if spacing is None and alpha is None:
        return outer, inners

    if spacing is None:
        diagonal = bounding_box_diagonal(outer, *inners)
        if diagonal <= 0.0:
            return outer, inners
        spacing = alpha * diagonal

    spacing = float(spacing)
    if spacing <= 0.0:
        raise ValueError(
            'spacing must be positive, got {}'.format(spacing))

    def one(loop: list[list[float]]) -> list[list[float]]:
        if len(loop) < 2:
            return loop
        # ``d_min`` as a spacing rather than a count: the subdivision is
        # per-segment, so the floor has to reach the segments to be felt.
        # ``sum(ceil(l_j / s)) >= ceil(perimeter / s) = d_min`` at this spacing.
        perimeter = sum(distance_point_point(a, b)
                        for a, b in pairwise(loop + loop[:1]))
        # Its own name: assigning to ``spacing`` in here would make it local to
        # ``one`` and unbound on this very line.
        step = spacing
        if d_min and perimeter > 0.0:
            step = min(spacing, perimeter / float(d_min))
        return _subdivide(loop, step)

    return one(outer), [one(loop) for loop in inners]


def resample_loop(points: list[list[float]], spacing: float | None = None) -> list[list[float]]:
    """One loop through :func:`discretise_boundary`, with no scale rule.

    Kept for callers outside this repository. ``spacing`` is
    ``target_length``, and ``None`` returns the loop unchanged. New code should
    call :func:`discretise_boundary`, which also applies eq. 4.1.
    """
    loop, _ = discretise_boundary(points, spacing=spacing,
                                  alpha=None, d_min=None)
    return loop


class Polyline(Polyline):

    def __init__(self, points: list[list[float]]) -> None:
        super(Polyline, self).__init__(points)

    def vertex_curvature(self, i: int) -> float | None:
        """Discrete polyline curvature.

        Parameters
        ----------
        i : int
            Vertex index.

        Returns
        -------
        curvature : float, None
            Curvature at the vertex.
            None if index out of range

        References
        ----------
        .. [1] Lionel Du Peloux. Modeling of bending-torsion couplings in active-bending structures.
               Application to the design of elastic gridshells. PhD thesis, Universite Paris Est, Ecole des Ponts ParisTech. 2017.
               Available at: https://tel.archives-ouvertes.fr/tel-01757782/document.

        """

        n = len(self.points)

        if i < 0 or i > n - 1:
            return None

        if i == 0 or i == n - 1:
            return 0.0

        a, b, c = self.points[i - 1: i + 2]
        ab = subtract_vectors(b, a)
        bc = subtract_vectors(c, b)
        ac = subtract_vectors(c, a)

        return 2 * length_vector(cross_vectors(ab, bc)) / (length_vector(ac) * length_vector(ab) * length_vector(bc))


# ==============================================================================
# Main
# ==============================================================================

if __name__ == '__main__':
    pass

    # import compas

    # points = [[0.0, 0.0, 0.0], [1.0, 5.0, 0.0], [2.0, 0.0, 0.0], [3.0, -5.0, 0.0], [4.0, 0.0, 0.0]]
    # polyline = Polyline(points)

    # for i in range(len(points)):
    #     print(polyline.vertex_curvature(i))
