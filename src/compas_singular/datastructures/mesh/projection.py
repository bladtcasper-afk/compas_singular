"""**Smoothing a mesh ON a surface** -- the 3D half of constrained smoothing.

:mod:`.smoothing` projects onto points and curves, which is all a planar domain
needs. A mesh draped over a curved surface needs two more things: every interior
vertex held ON the surface, and every boundary vertex held on the surface's
BORDERS rather than on a polyline of the mesh's own boundary. This module adds
those two, on top of :func:`.smoothing.constrained_smoothing`.

Ported from ``compas_singular.rhino.constraints``, which took Rhino GUIDs and
stopped importing under COMPAS 2. **Nothing here imports Rhino.** A surface or a
curve is anything :func:`.smoothing.closest_point_on_constraint` can project onto
-- ``compas_rhino``'s ``RhinoSurface`` / ``RhinoCurve`` inside Rhino, an
OCC-backed one outside it, a polyline -- and the borders and kinks the old Rhino
surface wrapper found for itself are passed in.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import annotations

from typing import Any
from typing import Sequence
from typing import TYPE_CHECKING

from compas.geometry import Point
from compas.geometry import closest_point_in_cloud
from compas.geometry import distance_point_point

from compas_singular.datastructures.mesh.smoothing import _closest_curve
from compas_singular.datastructures.mesh.smoothing import _split_loop_at_corners
from compas_singular.datastructures.mesh.smoothing import closest_point_on_constraint
from compas_singular.datastructures.mesh.smoothing import constrained_smoothing
from compas_singular.datastructures.mesh.smoothing import mesh_boundary_loops

if TYPE_CHECKING:
    from compas_singular.datastructures import Mesh


__all__ = [
    'automated_smoothing_surface_constraints',
    'automated_smoothing_constraints',
    'surface_constrained_smoothing',
]


class _NearestOf(object):
    """Several curves read as ONE constraint: a point goes to whichever is nearest.

    Re-chosen on every projection, so a boundary vertex that slides past the end
    of one border is picked up by the next.
    """

    def __init__(self, curves: Sequence[Any]) -> None:
        self.curves = list(curves)

    def closest_point(self, point: Any) -> list[float] | None:
        xyz = list(point)[:3]
        best, minimum = None, None
        for curve in self.curves:
            projected = closest_point_on_constraint(curve, xyz)
            distance = distance_point_point(xyz, projected)
            if minimum is None or distance < minimum:
                best, minimum = projected, distance
        return best


def automated_smoothing_surface_constraints(mesh: Mesh, surface: Any, borders: Sequence[Any], kinks: Sequence[list[float]] | None = None) -> dict[int, Any]:
    """Constrain every vertex of a mesh to a surface, its borders and its kinks.

    Interior vertices go to the surface. Each boundary vertex goes to the border
    curve it is nearest to NOW, and keeps that curve for the whole smoothing.
    Each kink pins the boundary vertex nearest to it, so the corners of the
    surface stay corners of the mesh.

    Parameters
    ----------
    mesh : :class:`compas.datastructures.Mesh`
        The mesh to constrain.
    surface : Any
        The surface, as anything with a ``closest_point`` method.
    borders : sequence
        The border curves of the surface, as polylines or curves.
    kinks : sequence[[float, float, float]], optional
        The points along the borders where the tangent jumps.

    Returns
    -------
    dict
        Vertex keys pointing to the surface, a border curve or a :class:`compas.geometry.Point`.

    """
    constraints = {vertex: surface for vertex in mesh.vertices()}

    boundary = [vertex for loop in mesh_boundary_loops(mesh) for vertex in loop]
    for vertex in boundary:
        xyz = mesh.vertex_coordinates(vertex)
        constraints[vertex] = min(
            borders, key=lambda curve: distance_point_point(xyz, closest_point_on_constraint(curve, xyz)))

    if kinks and boundary:
        cloud = [mesh.vertex_coordinates(vertex) for vertex in boundary]
        for xyz in kinks:
            index = closest_point_in_cloud(xyz, cloud)[2]
            constraints[boundary[index]] = Point(*list(xyz)[:3])

    return constraints


def automated_smoothing_constraints(mesh: Mesh, points: Sequence[list[float]] | None = None, curves: Sequence[Any] | None = None, surface: Any | None = None) -> dict[int, Any]:
    """Constrain the vertices of a mesh to points, curves and a surface.

    Every vertex goes to ``surface``. The boundary is then split at the vertices
    the ``points`` pin, and each run between two of them goes to the curve it is
    closest to on average -- :func:`.smoothing.automated_boundary_constraints`'s
    rule, with the pinned vertices as the corners. Pins win over everything.

    Parameters
    ----------
    mesh : :class:`compas.datastructures.Mesh`
        The mesh to constrain.
    points : sequence[[float, float, float]], optional
        Each pins the mesh vertex nearest to it, anywhere in the mesh.
    curves : sequence, optional
        Curves for the boundary vertices, as polylines or curves.
    surface : Any, optional
        A surface -- or a mesh -- for every vertex, as anything with a
        ``closest_point`` method.

    Returns
    -------
    dict
        Vertex keys pointing to a point, curve or surface.

    """
    constraints = {}

    pinned = {}
    if points:
        vertices = list(mesh.vertices())
        cloud = [mesh.vertex_coordinates(vertex) for vertex in vertices]
        for xyz in points:
            pinned[vertices[closest_point_in_cloud(xyz, cloud)[2]]] = Point(*list(xyz)[:3])

    if surface is not None:
        constraints.update({vertex: surface for vertex in mesh.vertices()})

    if curves:
        for loop in mesh_boundary_loops(mesh):
            for segment in _split_loop_at_corners(loop, pinned):
                curve = _closest_curve(mesh, segment, curves)
                for vertex in segment:
                    constraints[vertex] = curve

    constraints.update(pinned)
    return constraints


def surface_constrained_smoothing(mesh: Mesh, surface: Any, borders: Sequence[Any], kmax: int = 100, damping: float = 0.5, algorithm: str = 'centroid') -> dict[int, Any]:
    """Smooth a mesh while it stays on a surface and its boundary on the surface's borders.

    Unlike :func:`automated_smoothing_surface_constraints`, a boundary vertex is
    not tied to one border: every iteration projects it onto the nearest one.
    Boundary vertices with two neighbours -- the corners of a quad mesh -- are
    fixed.

    Parameters
    ----------
    mesh : :class:`compas.datastructures.Mesh`
        A mesh to smooth, modified in place.
    surface : Any
        The surface, as anything with a ``closest_point`` method.
    borders : sequence
        The border curves of the surface, as polylines or curves.
    kmax : int, optional
        Number of iterations for smoothing. Default is ``100``.
    damping : float, optional
        Damping value for smoothing between 0 and 1. Default is ``0.5``.
    algorithm : {'centroid', 'area', 'centerofmass'}, optional
        Type of smoothing algorithm to apply. Classic centroid by default.

    Returns
    -------
    dict
        The constraints that were applied. The mesh itself is modified in place.

    """
    edge = _NearestOf(borders)
    constraints = {}
    for vertex in mesh.vertices():
        constraints[vertex] = edge if mesh.is_vertex_on_boundary(vertex) else surface

    fixed = [vertex for loop in mesh_boundary_loops(mesh) for vertex in loop
             if mesh.vertex_degree(vertex) == 2]

    constrained_smoothing(
        mesh, kmax=kmax, damping=damping, constraints=constraints, algorithm=algorithm, fixed=fixed)

    return constraints
