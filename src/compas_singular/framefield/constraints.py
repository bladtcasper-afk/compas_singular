"""Step 4 -- directional constraints on the field.

A constraint pins the field at one vertex to one direction. Because a cross is
4-fold symmetric, the pinned quantity is ``exp(i * 4 * theta)``: all four arms of
the cross map to the same complex number, which is exactly what makes the solve
in ``field.py`` a plain linear system.

Constraints live on VERTICES, not faces -- ``field.py`` carries one unknown per
vertex.
"""
from __future__ import annotations

from cmath import exp as cexp
from cmath import phase
from collections import namedtuple
from math import atan2
from math import cos
from math import pi
from math import sin
from numbers import Number
from typing import Any
from typing import TYPE_CHECKING

from compas.geometry import distance_point_point
from compas.geometry import normalize_vector
from compas.geometry import subtract_vectors
from compas.itertools import pairwise

if TYPE_CHECKING:
    from compas_singular.framefield.background import BackgroundMesh


__all__ = ['Constraint', 'as_curve_list', 'from_boundary', 'from_curves',
           'representation']


#: ``weight`` of ``None`` means hard (eliminated from the solve); a float means a
#: least-squares penalty, so several constraints can disagree gracefully.
Constraint = namedtuple('Constraint', 'vkey direction weight')


def representation(direction: list[float]) -> complex:
    """The 4-fold-symmetric complex representation of a direction.

    ``exp(i * 4 * theta)``. Note ``representation(d) == representation(rotate90(d))``
    -- that identity IS the cross field, and it is what makes ``mode='tangent'``
    and ``mode='perpendicular'`` below the same constraint.
    """
    theta = atan2(direction[1], direction[0])
    return cexp(4j * theta)


def as_curve_list(curves: Any) -> list[Any]:
    """A list of curves, whether one curve or several were passed.

    The test is on the DEPTH, not on the leaf type: ``curves[0][0]`` is a number
    only when ``curves[0]`` is a point, i.e. when ``curves`` is itself a single
    curve. The version that tested ``isinstance(curves[0][0], (list, tuple))``
    instead misread a list of compas ``Polyline``s -- or of lists of ``Point``s
    -- as ONE curve, wrapped it, and then failed in ``from_curves`` on
    ``float(Point)``, because ``p[0]`` of a Polyline is a Point rather than a
    coordinate. Rhino input arrives in exactly that shape.

    Points may be lists, tuples or compas ``Point``s; curves lists or compas
    ``Polyline``s. Nothing is copied or coerced here -- callers do that.
    """
    if curves and isinstance(curves[0][0], Number):
        return [curves]
    return list(curves or [])


def from_boundary(
    background: BackgroundMesh,
    weight: float | None = None,
    corner_tolerance: float = 0.1,
) -> list[Constraint]:
    """Pin the field tangent to every boundary vertex.

    This is the constraint that makes quads sit flush against the wall, and it is
    the only one needed to reproduce a skeleton-quality layout (milestone 1).

    A boundary vertex has TWO adjacent edge directions, and they are combined
    here in the 4th-power representation rather than as a geometric average --
    see ``BackgroundMesh.boundary_tangents`` for why the geometric one is wrong.
    In representation space a right-angle corner's two edges are literally the
    same number, so the average is exact and the corner is pinned to the cross
    both walls agree on.

    Where the two disagree the average shrinks, and at an interior angle near 45
    or 135 degrees it cancels to nothing: the two walls want crosses 45 degrees
    apart and no single cross satisfies both. That corner is left FREE rather
    than pinned to an arbitrary compromise -- the field then decides for itself,
    which is the honest outcome and usually means putting a singularity nearby.

    Parameters
    ----------
    background : BackgroundMesh
    weight : float, optional
        ``None`` (default) for hard constraints.
    corner_tolerance : float, optional
        Averaged representations shorter than this leave the vertex free.

    Returns
    -------
    list[Constraint]
    """
    out = []
    for vkey, tangents in background.boundary_tangents().items():
        acc = sum((representation(t) for t in tangents), 0j) / len(tangents)
        if abs(acc) < corner_tolerance:
            continue                       # ~45 degree corner: no cross satisfies both walls
        theta = phase(acc) / 4.0
        out.append(Constraint(vkey, [cos(theta), sin(theta), 0.0], weight))
    return out


def from_curves(
    background: BackgroundMesh,
    curves: Any,
    mode: str = 'perpendicular',
    band: float | None = None,
    weight: float | None = 1.0,
    skip_boundary: bool = True,
) -> list[Constraint]:
    """Pin the field along guide curves -- cables, force lines.

    Parameters
    ----------
    background : BackgroundMesh
    curves : list[list[[x, y, z]]]
        The guide curves. Points may be lists, tuples or compas ``Point``s, and
        a curve may be a list or a compas ``Polyline`` -- everything is coerced
        to plain XY floats here, exactly as ``background._as_open_loop`` does
        for the boundaries. A single curve may also be passed directly.
    mode : {'perpendicular', 'tangent'}, optional
        Whether mesh edges should run ACROSS the curve or ALONG it.

        For a cross field these are the SAME CONSTRAINT: rotating a direction by
        90 degrees leaves ``exp(i*4*theta)`` unchanged. The parameter is kept
        because it stops being an identity the moment ``warp.py`` exists and the
        frame is no longer orthogonal -- at which point "perpendicular to a
        cable" and "along a cable" pick different arms of the frame. Removing it
        as dead code now would be a real regression later.
    band : float, optional
        Only vertices within this distance of a curve are constrained. Defaults
        to ``background.target_length``, i.e. the curve steers its own
        neighbourhood rather than fighting the boundary across the domain.
    weight : float, optional
        Least-squares weight. ``None`` for hard. Soft by default: a guide that
        disagrees with the boundary should bend, not break the solve.
    skip_boundary : bool, optional
        Leave boundary vertices to ``from_boundary``. A guide reaching the wall
        would otherwise fight the wall-tangency constraint at the one vertex
        where the wall must win.

    Returns
    -------
    list[Constraint]

    Notes
    -----
    **The result does not depend on the order of ``curves``, and it used to.**
    Every curve within ``band`` of a vertex contributes its own constraint and
    the solve accumulates them; there is no first-come claim. The version that
    had one gave each vertex to the FIRST curve reaching it, so two guides 45
    degrees apart produced a different mesh depending on list order, silently --
    6 of 91 constraints dropped on the measured case. On four arcs that are
    exact mirror images of each other it was worse than silent: with a perfectly
    symmetric background the constraint set still mismatched by 1.145 under
    every non-identity symmetry, which is the full magnitude of a
    representation, i.e. a different constraint entirely.
    """
    curves = as_curve_list(curves)

    radius = background.target_length if band is None else band
    boundary = background.boundary_vertices() if skip_boundary else set()

    out = []
    for curve in curves:
        pts = [[float(p[0]), float(p[1]), 0.0] for p in curve]
        for vkey in background.mesh.vertices():
            if vkey in boundary:
                continue
            p = background.mesh.vertex_coordinates(vkey)
            best, tangent = float('inf'), None
            for a, b in pairwise(pts):
                ab = subtract_vectors(b, a)
                length2 = ab[0] * ab[0] + ab[1] * ab[1]
                if length2 == 0.0:
                    continue
                t = ((p[0] - a[0]) * ab[0] + (p[1] - a[1]) * ab[1]) / length2
                t = max(0.0, min(1.0, t))
                q = [a[0] + ab[0] * t, a[1] + ab[1] * t, 0.0]
                d = distance_point_point(p, q)
                if d < best:
                    best, tangent = d, normalize_vector(ab)
            if tangent is None or best > radius:
                continue
            direction = tangent
            if mode == 'perpendicular':
                direction = [-tangent[1], tangent[0], 0.0]
            elif mode != 'tangent':
                raise ValueError("mode must be 'perpendicular' or 'tangent'")
            out.append(Constraint(vkey, direction, weight))
    return _merge_hard(out) if weight is None else out


def _merge_hard(constraints: list[Constraint]) -> list[Constraint]:
    """One hard constraint per vertex, averaged in the 4th-power representation.

    Soft constraints need no help: ``field.solve`` accumulates them per vertex,
    and its own comment calls that "the right behaviour where two guides
    overlap". Hard ones it puts in a dict, so the LAST one written wins and the
    answer depends on curve order -- the same defect the ``claimed`` set used to
    have, just moved one function along.

    Averaging in the representation is the same averaging ``from_boundary`` does
    at a corner, and it cancels for the same reason: two guides crossing at 45
    degrees want crosses no single cross satisfies, and that vertex is left free
    rather than pinned to an arbitrary one of them.
    """
    grouped = {}
    for c in constraints:
        grouped.setdefault(c.vkey, []).append(c)

    out = []
    for group in grouped.values():
        if len(group) == 1:
            out.append(group[0])
            continue
        acc = sum((representation(c.direction) for c in group), 0j) / len(group)
        if abs(acc) < 1e-9:
            continue
        theta = phase(acc) / 4.0
        out.append(Constraint(group[0].vkey, [cos(theta), sin(theta), 0.0], None))
    return out


def _assert_mode_identity() -> bool:
    """The 4-fold identity the ``mode`` docstring claims. Called by the tests."""
    d = normalize_vector([0.37, 0.93, 0.0])
    perp = [-d[1], d[0], 0.0]
    assert abs(representation(d) - representation(perp)) < 1e-12
    assert abs(representation(d) - representation([-d[0], -d[1], 0.0])) < 1e-12
    return True


# a cross's period: the angle you can rotate by without changing the object
PERIOD = pi / 2.0
