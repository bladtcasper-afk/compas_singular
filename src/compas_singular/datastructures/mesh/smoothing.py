from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

from math import pi

from compas.datastructures import mesh_smooth_centerofmass
from compas.geometry import Curve
from compas.geometry import Point
from compas.geometry import Polyline
from compas.geometry import angle_vectors
from compas.geometry import closest_point_on_polyline
from compas.geometry import closest_point_on_segment
from compas.geometry import distance_point_point
from compas.geometry import length_vector
from compas.geometry import subtract_vectors
from compas.itertools import flatten
from compas.itertools import pairwise


__all__ = [
    'closest_point_on_constraint',
    'mesh_boundary_loops',
    'mesh_boundary_polylines',
    'mesh_boundary_corners',
    'constrained_smoothing',
    'automated_boundary_constraints',
    'boundary_constrained_smoothing',
    'smoothing_region',
]


# ==============================================================================
# Closest point
# ==============================================================================

def _is_xyz(item):
    """Return whether an item is a bare sequence of three numbers."""
    return (
        isinstance(item, (list, tuple))
        and len(item) == 3
        and all(isinstance(i, (int, float)) for i in item)
    )


def closest_point_on_constraint(constraint, xyz, discretisation=128):
    """Project a point onto a constraint object, using COMPAS geometry only.

    Parameters
    ----------
    constraint : None | :class:`compas.geometry.Point` | [float, float, float] | :class:`compas.geometry.Polyline` | sequence[point] | :class:`compas.geometry.Curve` | Any
        The geometry to project onto:

        * ``None`` -- no constraint, None is returned;
        * a :class:`compas.geometry.Point` or a bare ``[x, y, z]`` -- the point itself, i.e. a pin;
        * a :class:`compas.geometry.Polyline` or any sequence of points -- the closest point
          on that polyline, via :func:`compas.geometry.closest_point_on_polyline`;
        * any other :class:`compas.geometry.Curve` -- a :class:`~compas.geometry.Line`, a
          :class:`~compas.geometry.Circle`, a curve from ``compas_rhino``'s
          ``curve_to_compas`` or from ``compas_occ`` -- the closest point ON THE CURVE, from
          its own ``closest_point``. See :func:`_closest_point_on_curve` for the curves that
          have none and fall back on their ``to_polyline`` discretisation;
        * any other object with a working ``closest_point`` method (a
          :class:`~compas_singular.editing.GuideCurve`, a Rhino or OCC surface, ...) -- the
          result of that method;
        * anything else with a ``to_polyline`` -- the closest point on that discretisation.
    xyz : [float, float, float]
        The coordinates of the point to project.
    discretisation : int, optional
        Number of segments used to discretise a curve that has no analytical closest point.
        Default is ``128``.

    Returns
    -------
    list[float] | None
        The XYZ coordinates of the closest point on the constraint,
        or None if the constraint is None.

    Raises
    ------
    TypeError
        If the closest point on the constraint cannot be computed.

    """
    if constraint is None:
        return None

    # a point constraint pins the vertex
    if isinstance(constraint, Point) or _is_xyz(constraint):
        return [float(constraint[0]), float(constraint[1]), float(constraint[2])]

    # a polyline, or any bare sequence of points read as one. BEFORE the curve branch: a
    # Polyline is a Curve too, and its own closest_point raises NotImplementedError
    if isinstance(constraint, Polyline):
        return closest_point_on_polyline(xyz, constraint)
    if isinstance(constraint, (list, tuple)):
        return closest_point_on_polyline(xyz, [list(point) for point in constraint])

    # a parametric curve: onto the curve itself, not onto a sampling of it
    if isinstance(constraint, Curve):
        return _closest_point_on_curve(constraint, xyz, discretisation)

    # anything else that knows how to project a point: a GuideCurve, Rhino or OCC surfaces, ...
    method = getattr(constraint, 'closest_point', None)
    if callable(method):
        try:
            return list(method(Point(*xyz)))[:3]
        except NotImplementedError:
            pass

    # curves without an analytical projection: discretise and project on the polyline
    method = getattr(constraint, 'to_polyline', None)
    if callable(method):
        return closest_point_on_polyline(xyz, method(n=discretisation))

    raise TypeError('Cannot compute the closest point on a constraint of type {}.'.format(type(constraint)))


#: Curve types whose missing ``closest_point`` has already been explained by
#: :func:`_closest_point_on_curve`. Once per TYPE, not per call: smoothing projects every
#: constrained vertex at every iteration, and a print per call would put thousands of
#: identical lines on Rhino's command line.
_EXPLAINED_FALLBACKS = set()


def _closest_point_on_curve(curve, xyz, discretisation):
    """The closest point on a compas curve, from the curve itself wherever it can say.

    ``closest_point(point=...)`` is exact on a ``Line``, a ``Circle``, a ``compas_rhino``
    curve (Rhino's ``ClosestPoint``) and a ``compas_occ`` curve
    (``GeomAPI_ProjectPointOnCurve``). Two ways it gives no answer, and both fall back on
    the curve's ``to_polyline`` -- a sagitta off the curve rather than no point at all:

    * ``Arc``, ``Ellipse`` and ``Bezier`` raise ``NotImplementedError`` (compas 2.15). The
      fallback is printed, once per curve type;
    * a Rhino curve returns ``None`` when ``ClosestPoint`` fails.
    """
    try:
        closest = curve.closest_point(point=Point(*xyz))
    except NotImplementedError:
        if type(curve) not in _EXPLAINED_FALLBACKS:
            _EXPLAINED_FALLBACKS.add(type(curve))
            print('note: {} does not implement closest_point, so points are projected onto a '
                  '{}-segment polyline of it instead -- up to a sagitta off the curve itself. '
                  'Shown once per curve type.'.format(type(curve).__name__, discretisation))
        closest = None
    if closest is None:
        return closest_point_on_polyline(xyz, curve.to_polyline(n=discretisation))
    return [float(closest[0]), float(closest[1]), float(closest[2])]


def _is_polyline_like(constraint):
    """Return whether a constraint is a polyline, or a sequence of points read as one."""
    if isinstance(constraint, Polyline):
        return True
    return isinstance(constraint, (list, tuple)) and not _is_xyz(constraint)


class _PolylineProjector(object):
    """Project points onto a polyline, searching around the previous result.

    Projecting on a polyline costs one segment projection per segment, which the
    smoothing loop would pay again for every vertex at every iteration. A vertex barely
    moves between two iterations, so the search is limited to a window of segments
    around the one it was projected on last time. The window is only a shortcut: the
    first projection, and any projection that lands on the edge of the window, searches
    the whole polyline.
    """

    def __init__(self, polyline, window=10):
        points = [list(point) for point in polyline]
        if len(points) < 2:
            raise ValueError('A polyline constraint needs at least two points.')
        self.segments = list(pairwise(points))
        self.closed = distance_point_point(points[0], points[-1]) < 1e-9
        self.window = window

    def _search(self, xyz, indices):
        closest, minimum, index = None, None, None
        for i in indices:
            point = closest_point_on_segment(xyz, self.segments[i])
            distance = distance_point_point(xyz, point)
            if minimum is None or distance < minimum:
                closest, minimum, index = point, distance, i
        return closest, index

    def closest_point(self, xyz, hint=None):
        """Return the closest point on the polyline and the index of its segment."""
        count = len(self.segments)

        if hint is None or 2 * self.window + 1 >= count:
            return self._search(xyz, range(count))

        if self.closed:
            offsets = range(-self.window, self.window + 1)
            indices = [(hint + offset) % count for offset in offsets]
        else:
            indices = list(range(max(0, hint - self.window), min(count, hint + self.window + 1)))

        closest, index = self._search(xyz, indices)

        # the window was too small to contain the closest segment: search all of it
        if index == indices[0] or index == indices[-1]:
            return self._search(xyz, range(count))

        return closest, index


# ==============================================================================
# Boundary
# ==============================================================================

def mesh_boundary_loops(mesh):
    """Collect the boundaries of a mesh as ordered, non-repeating loops of vertices.

    Parameters
    ----------
    mesh : :class:`compas.datastructures.Mesh`
        A mesh.

    Returns
    -------
    list[list[int]]
        A list of boundaries as ordered lists of vertex keys.
        The first vertex is not repeated at the end, the loops are implicitly closed.

    """
    loops = []
    for boundary in mesh.vertices_on_boundaries():
        loop = list(boundary)
        while len(loop) > 1 and loop[-1] == loop[0]:
            loop.pop()
        if len(loop) > 2:
            loops.append(loop)
    return loops


def mesh_boundary_polylines(mesh):
    """Collect the boundaries of a mesh as closed polylines.

    Parameters
    ----------
    mesh : :class:`compas.datastructures.Mesh`
        A mesh.

    Returns
    -------
    list[:class:`compas.geometry.Polyline`]
        A closed polyline per mesh boundary.

    """
    polylines = []
    for loop in mesh_boundary_loops(mesh):
        points = [mesh.vertex_coordinates(vertex) for vertex in loop]
        polylines.append(Polyline(points + points[:1]))
    return polylines


def mesh_boundary_corners(mesh, corner_angle=pi / 6):
    """Collect the boundary vertices of a mesh that sit on a kink.

    The kink is measured as the deviation from straight between the two boundary
    edges at the vertex: zero along a straight boundary, pi / 2 at a square corner.

    Parameters
    ----------
    mesh : :class:`compas.datastructures.Mesh`
        A mesh.
    corner_angle : float, optional
        Threshold deviation angle, in radians. Default is ``pi / 6`` (30 degrees).

    Returns
    -------
    list[int]
        The keys of the boundary vertices at a kink.

    """
    corners = []
    for loop in mesh_boundary_loops(mesh):
        count = len(loop)
        for index, vertex in enumerate(loop):
            a = mesh.vertex_coordinates(loop[index - 1])
            b = mesh.vertex_coordinates(vertex)
            c = mesh.vertex_coordinates(loop[(index + 1) % count])
            u = subtract_vectors(b, a)
            v = subtract_vectors(c, b)
            if length_vector(u) == 0 or length_vector(v) == 0:
                continue
            if angle_vectors(u, v) > corner_angle:
                corners.append(vertex)
    return corners


def _split_loop_at_corners(loop, corners):
    """Split a closed loop of vertices into the segments delimited by its corner vertices.

    The corner vertices are shared by the two segments they delimit.
    """
    indices = [index for index, vertex in enumerate(loop) if vertex in corners]
    if len(indices) < 2:
        return [loop]

    count = len(loop)
    segments = []
    for i, start in enumerate(indices):
        end = indices[(i + 1) % len(indices)]
        segment = [loop[start]]
        index = start
        while index != end:
            index = (index + 1) % count
            segment.append(loop[index])
        segments.append(segment)
    return segments


def _closest_curve(mesh, vertices, curves):
    """Return the curve that is closest, on average, to a run of mesh vertices."""
    sample = vertices[1:-1] or vertices
    points = [mesh.vertex_coordinates(vertex) for vertex in sample]

    closest, minimum = None, None
    for curve in curves:
        distance = sum(
            distance_point_point(point, closest_point_on_constraint(curve, point))
            for point in points) / len(points)
        if minimum is None or distance < minimum:
            closest, minimum = curve, distance
    return closest


# ==============================================================================
# Smoothing
# ==============================================================================

def constrained_smoothing(mesh, kmax=100, damping=0.5, constraints=None, algorithm='centroid', fixed=None,
                          symmetric=None):
    """Constrained smoothing of a mesh. Constraints can be points, curves or surfaces.

    The projections go through :func:`closest_point_on_constraint`, which builds on
    :mod:`compas.geometry`, so this runs headless, without Rhino. For a mesh on a
    surface, :mod:`.projection` builds the constraints.

    Parameters
    ----------
    mesh : :class:`compas.datastructures.Mesh`
        A mesh to smooth, modified in place.
    kmax : int, optional
        Number of iterations for smoothing. Default is ``100``.
    damping : float, optional
        Damping value for smoothing between 0 and 1. Default is ``0.5``.
    constraints : dict, optional
        Dictionary of constraints as vertex keys pointing to COMPAS geometry
        (points, polylines, curves or surfaces). Empty by default.
        See :func:`closest_point_on_constraint` for the accepted geometry.
    algorithm : {'centroid', 'area', 'centerofmass'}, optional
        Type of smoothing algorithm to apply. Classic centroid by default.
    fixed : sequence[int], optional
        Vertices that the smoothing algorithm must not move at all. Default is None.
    symmetric : bool, optional
        Keep a mesh made by ``expand_symmetrically`` exactly symmetric: every
        iteration averages each vertex over its orbit (read from
        ``attributes['orbits']``, not searched for) before projecting onto the
        constraints. ``None`` (default) means: when the mesh has orbits. Without
        it, each iteration drifts by floating point and by the constraint
        projections, and the symmetry the expansion built slowly erodes.

    Returns
    -------
    None
        The mesh is modified in place.

    """
    constraints = constraints or {}
    if symmetric is None:
        symmetric = bool(mesh.attributes.get('orbits'))
    orbit_maps = None
    if symmetric:
        from compas_singular.symmetry.replicate import orbit_maps as _orbit_maps
        from compas_singular.symmetry.replicate import symmetrise_positions
        orbit_maps = _orbit_maps(mesh)

    # polylines get a projector that only searches around the previous result
    projectors = {}
    for constraint in constraints.values():
        if _is_polyline_like(constraint) and id(constraint) not in projectors:
            projectors[id(constraint)] = _PolylineProjector(constraint)
    hints = {}

    def project(vertex, constraint, xyz):
        projector = projectors.get(id(constraint))
        if projector is None:
            return closest_point_on_constraint(constraint, xyz)
        point, hints[vertex] = projector.closest_point(xyz, hints.get(vertex))
        return point

    def callback(k, args):
        mesh, constraints = args
        if orbit_maps is not None:
            symmetrise_positions(mesh, orbit_maps)
        for vertex, constraint in constraints.items():
            if constraint is None:
                continue
            xyz = project(vertex, constraint, mesh.vertex_coordinates(vertex))
            if xyz is None:
                continue
            mesh.vertex_attributes(vertex, 'xyz', xyz)

    if algorithm == 'area':
        mesh.smooth_area(fixed=fixed, kmax=kmax, damping=damping, callback=callback, callback_args=[mesh, constraints])
    elif algorithm == 'centerofmass':
        mesh_smooth_centerofmass(mesh, fixed=fixed, kmax=kmax, damping=damping, callback=callback, callback_args=[mesh, constraints])
    elif algorithm == 'centroid':
        mesh.smooth_centroid(fixed=fixed, kmax=kmax, damping=damping, callback=callback, callback_args=[mesh, constraints])
    else:
        raise ValueError(f'{algorithm} is not a recognised smoothing algorithm. Pick area, centerofmass or centroid instead')
        # mesh.smooth_centroid(fixed=fixed, kmax=kmax, damping=damping, callback=callback, callback_args=[mesh, constraints])


def automated_boundary_constraints(mesh, curves=None, corner_angle=pi / 6, fix_corners=True):
    """Automatically constrain the boundary vertices of a mesh to the boundary.

    Every boundary vertex is constrained to the curve it belongs to, so that it slides
    along the boundary during smoothing instead of being pinned to its start position.
    Without input curves, the boundary of the mesh itself -- in its current state -- is
    the constraint. Vertices at a kink are pinned, so that corners survive the smoothing.

    Parameters
    ----------
    mesh : :class:`compas.datastructures.Mesh`
        The mesh to constrain.
    curves : sequence, optional
        The boundary geometry to snap the boundary vertices to, as COMPAS polylines or
        curves. Each run of boundary vertices between two corners is matched to the curve
        it is closest to. Default is None, in which case the mesh boundary polylines are used.
    corner_angle : float, optional
        Threshold deviation angle for a boundary kink, in radians.
        Default is ``pi / 6`` (30 degrees).
    fix_corners : bool, optional
        If True, pin the boundary vertices at a kink to their current position.
        Default is True.

    Returns
    -------
    dict
        A dictionary of mesh constraints for smoothing, as vertex keys pointing to
        point, polyline or curve objects.

    """
    constraints = {}
    corners = set(mesh_boundary_corners(mesh, corner_angle)) if fix_corners else set()

    for loop in mesh_boundary_loops(mesh):
        if not curves:
            points = [mesh.vertex_coordinates(vertex) for vertex in loop]
            polyline = Polyline(points + points[:1])
            for vertex in loop:
                constraints[vertex] = polyline
        else:
            for segment in _split_loop_at_corners(loop, corners):
                curve = _closest_curve(mesh, segment, curves)
                for vertex in segment:
                    constraints[vertex] = curve

    for vertex in corners:
        constraints[vertex] = Point(*mesh.vertex_coordinates(vertex))

    return constraints


def boundary_constrained_smoothing(mesh, curves=None, kmax=100, damping=0.5, algorithm='centroid',
                                   corner_angle=pi / 6, fix_corners=True, constraints=None):
    """Smooth a mesh with its boundary vertices automatically constrained to the boundary.

    The boundary vertices are identified and constrained by
    :func:`automated_boundary_constraints`, then :func:`constrained_smoothing` snaps them
    back onto that boundary after every iteration. The interior vertices relax freely, the
    boundary vertices slide along the boundary, and the corners stay put.

    Parameters
    ----------
    mesh : :class:`compas.datastructures.Mesh`
        A mesh to smooth, modified in place.
    curves : sequence, optional
        The boundary geometry to snap the boundary vertices to, as COMPAS polylines or
        curves. Default is None, in which case the mesh boundary polylines are used.
    kmax : int, optional
        Number of iterations for smoothing. Default is ``100``.
    damping : float, optional
        Damping value for smoothing between 0 and 1. Default is ``0.5``.
    algorithm : {'centroid', 'area', 'centerofmass'}, optional
        Type of smoothing algorithm to apply. Classic centroid by default.
    corner_angle : float, optional
        Threshold deviation angle for a boundary kink, in radians.
        Default is ``pi / 6`` (30 degrees).
    fix_corners : bool, optional
        If True, pin the boundary vertices at a kink to their current position.
        Default is True.
    constraints : dict, optional
        Additional constraints, as vertex keys pointing to COMPAS geometry, applied on
        top of the automated ones. Default is None.

    Returns
    -------
    dict
        The constraints that were applied, as vertex keys pointing to point, polyline or
        curve objects. The mesh itself is modified in place.

    """
    boundary_constraints = automated_boundary_constraints(
        mesh, curves=curves, corner_angle=corner_angle, fix_corners=fix_corners)

    if constraints:
        boundary_constraints.update(constraints)

    constrained_smoothing(
        mesh, kmax=kmax, damping=damping, constraints=boundary_constraints, algorithm=algorithm)

    return boundary_constraints

def boundary_smoothing(mesh, kmax=100, damping=0.5):
    """Only the boundary vertices are smoothened by sliding them along the boundary curves."""

    vertices = list(mesh.vertices())
    boundary_vertices = list(flatten(mesh.vertices_on_boundaries()))

    fixed = [vkey for vkey in vertices if vkey not in boundary_vertices]

    boundary_constraints = automated_boundary_constraints(mesh)

    constrained_smoothing(mesh, kmax=kmax, damping=damping, constraints=boundary_constraints, algorithm='area', fixed=fixed)


def smoothing_region(mesh, vertices, kmax=50, damping=0.5, blend=3, constraints=None):
    """Area-weighted smoothing of a region of a mesh, with per-vertex damping.

    Only the selected region moves: every vertex outside it is fixed. The damping is
    tapered from full strength on the selected core down to zero over ``blend`` rings of
    neighbours around it, so the region blends into the rest of the mesh. A region
    smoothed at full strength up to a hard edge leaves a CREASE at that edge -- a fixed
    ring next to a fully relaxed one is a kink, which is the sort of defect this is meant
    to remove. The taper is what makes it blend.

    The smoothing is area-weighted, not centroid: centroid equalises edge lengths, which
    fights the grading a frame-field mesh is supposed to have.

    This is the region counterpart of :func:`boundary_constrained_smoothing`, and the
    headless form of the region branch of the Rhino ``CMD_smoothen`` command.

    Parameters
    ----------
    mesh : :class:`compas.datastructures.Mesh`
        A mesh to smooth, modified in place.
    vertices : sequence[int] | dict[int, float]
        The region to smooth, either as the vertex keys of its core -- which are then
        tapered over ``blend`` rings -- or as a ready-made dictionary of vertex keys
        pointing to damping weights between 0 and 1, in which case ``blend`` is ignored.
        A vertex that is absent, or whose weight is not positive, is fixed.
    kmax : int, optional
        Number of iterations for smoothing. Default is ``50``.
    damping : float, optional
        Damping value for smoothing between 0 and 1, scaled per vertex by its weight.
        Default is ``0.5``.
    blend : int, optional
        Rings of vertices outside the core over which the damping falls to zero.
        Default is ``3``. Ignored if ``vertices`` is a dictionary of weights.
    constraints : dict, optional
        Dictionary of constraints as vertex keys pointing to COMPAS geometry, applied to
        the vertices of the region after every iteration. Default is None, in which case
        :func:`automated_boundary_constraints` is used: a region touching the boundary
        then slides along it instead of being dragged inward.
        See :func:`closest_point_on_constraint` for the accepted geometry.

    Returns
    -------
    dict[int, float]
        The damping weights that were applied, as vertex keys pointing to a weight
        between 0 and 1. The mesh itself is modified in place.

    """
    def taper(core, rings):
        """Per-vertex damping weight: 1 in the core, falling to 0 outside the blend.

        A vertex with no weight is fixed. The ring just beyond the blend is therefore
        held and the ones just inside it barely move, so there is no step in the result.
        """
        weights = {vertex: 1.0 for vertex in core}
        frontier = set(core)
        for ring in range(1, rings + 1):
            frontier = {n for vertex in frontier for n in mesh.vertex_neighbors(vertex)} - set(weights)
            if not frontier:
                break
            weight = 1.0 - ring / float(rings + 1)
            for vertex in frontier:
                weights[vertex] = weight
        return weights

    def area_target(vertex):
        """The area-weighted centroid of the faces around a vertex, or None."""
        areas, centroids = [], []
        for face in mesh.vertex_faces(vertex, ordered=True):
            if face is None:
                continue
            areas.append(mesh.face_area(face))
            centroids.append(mesh.face_centroid(face))
        total = sum(areas)
        if not areas or total <= 0.0:
            return None
        return [sum(a * c[i] for a, c in zip(areas, centroids)) / total for i in range(3)]

    weights = dict(vertices) if isinstance(vertices, dict) else taper(set(vertices), blend)

    if constraints is None:
        constraints = automated_boundary_constraints(mesh)
    constraints = {vertex: constraint for vertex, constraint in constraints.items()
                   if weights.get(vertex, 0.0) > 0.0}

    for _ in range(kmax):
        # every target is collected before any vertex moves, so the result does not
        # depend on the order the region is walked in
        targets = {}
        for vertex, weight in weights.items():
            if weight <= 0.0:
                continue
            target = area_target(vertex)
            if target is not None:
                targets[vertex] = target

        for vertex, target in targets.items():
            xyz = mesh.vertex_coordinates(vertex)
            step = damping * weights[vertex]
            mesh.vertex_attributes(vertex, 'xyz', [
                xyz[i] + step * (target[i] - xyz[i]) for i in range(3)])

        # a region touching the outline must slide along it, not be dragged inward
        for vertex, constraint in constraints.items():
            moved = closest_point_on_constraint(constraint, mesh.vertex_coordinates(vertex))
            if moved is not None:
                mesh.vertex_attributes(vertex, 'xyz', moved)

    return weights


# ==============================================================================
# Relaxation
# ==============================================================================


def relaxation(mesh, fixed="corners", fixed_vertices = [], constraints=None, q_factor=100, algorithm='forcedensity'):
    from compas_fd.solvers import fd_constrained_numpy

    if fixed == "corners":
        fixed = mesh_boundary_corners(mesh)
    elif fixed == "boundary":
        boundaries = mesh_boundary_loops(mesh)
        fixed = list(flatten(boundaries))
    elif fixed=="manual":
        fixed = fixed_vertices
    else:
        raise ValueError("This type of fixing is not recognised. Choose corners, boundary or manual with fixed_vertices.")

    vertices = mesh.vertices_attributes("xyz")
    edges = list(mesh.edges())
    constraints = []
    q = []
    q_baseline = 1.0
    for edge in edges:
        if mesh.is_edge_on_boundary(edge):
            q.append(q_factor*q_baseline)
        else:
            q.append(q_baseline)
    loads = [[0, 0, 0] for _ in range(mesh.number_of_vertices())]

    result = fd_constrained_numpy(
        vertices=vertices,
        fixed=fixed,
        edges=edges,
        forcedensities=q,
        loads=loads,
        constraints=constraints,
        )

    for vertex, attr in mesh.vertices(data=True):
        attr["x"] = result.vertices[vertex, 0]
        attr["y"] = result.vertices[vertex, 1]
        attr["z"] = result.vertices[vertex, 2]

    return mesh


# ==============================================================================
# Main
# ==============================================================================

if __name__ == '__main__':
    pass
