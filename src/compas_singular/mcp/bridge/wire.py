"""**Meshes and curves as plain, index-based JSON.**

The encoding on the spool::

    {"vertices": [[x, y, z], ...],
     "faces":    [[i, j, k, l], ...],     indices into vertices
     "poles":    [[x, y, z], ...]}        POINTS, never keys

**Poles travel as points because keys renumber.** Every topological edit
renumbers vertices, and a pole recorded as key 37 is a pole somewhere else after
the next weld. This is the same decision ``helpers.read_coarse`` already made
when it rebuilds a layout with
``CoarsePseudoQuadMesh.from_vertices_and_faces_with_poles(vertices, faces,
poles)`` -- the points survive the trip, the keys do not.

**Nothing compas-typed crosses the wire, and that is load-bearing.**
``CMD_start.import_compas_singular`` purges every ``compas_singular*`` module
from ``sys.modules``, so running any other command in the family while the link
is attached leaves two copies of the package alive at once. Anything that
pickled, or that checked ``isinstance`` across the boundary, would fail with the
``not the same object as compas_singular.datastructures.mesh.mesh.Mesh`` error
recorded in ``RHINO_PLUGIN.md``. Plain lists of floats have no identity to get
wrong, so the link survives a purge it never notices.

**Do not reach for** ``mesh.to_vertices_and_faces()`` **here.** This library
overrides it with ``keep_keys=True`` as the DEFAULT, which returns dicts keyed by
vertex and face key -- while upstream compas returns lists. A mesh pulled out of
Rhino is a plain compas ``Mesh`` and a mesh built here is this library's, so the
same call means two different things depending on which side made the object.
Both are JSON-serialisable, so the mistake would not raise; it would just put
face-key-indexed nonsense on the wire. The index map below is built explicitly
for that reason.

**Coordinates are not rounded.** The server holds the authoritative mesh in
double precision; Rhino's copy is already single precision by the time it is a
``Point3f`` inside ``bake_mesh``. Rounding here would lose precision a second
time, for nothing.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


__all__ = [
    'mesh_to_wire',
    'mesh_from_wire',
    'points_of',
    'curve_to_wire',
    'curves_to_wire',
]


class WireError(ValueError):
    """The wire payload is not a mesh. Raised where the fault is, not later."""


def points_of(curve):
    """Any curve-ish thing as a plain list of ``[x, y, z]``.

    Accepts what the rest of the library hands around: a compas ``Polyline``, a
    sequence of compas ``Point`` objects, or a sequence of bare triples.

    Parameters
    ----------
    curve : Polyline | sequence

    Returns
    -------
    list[list[float]]
    """
    points = getattr(curve, 'points', None)
    if points is None:
        points = curve
    return [[float(p[0]), float(p[1]), float(p[2])] for p in points]


def curve_to_wire(curve):
    """One curve as a point list. Alias of :func:`points_of`, named for the wire."""
    return points_of(curve)


def curves_to_wire(curves):
    """A sequence of curves as a list of point lists."""
    return [points_of(curve) for curve in curves or []]


def _pole_points(mesh):
    """The mesh's poles as points, or ``[]`` if it has no notion of one.

    A plain compas ``Mesh`` -- what comes back out of Rhino -- has no ``poles``,
    and that is not an error: it means the mesh carries no pole data, so there is
    nothing to preserve.
    """
    poles = getattr(mesh, 'poles', None)
    if poles is None:
        return []
    try:
        keys = list(poles())
    except Exception:
        return []
    return [list(mesh.vertex_coordinates(vkey)) for vkey in keys]


def mesh_to_wire(mesh):
    """Encode a mesh for the spool.

    Parameters
    ----------
    mesh : Mesh
        Any compas mesh, including this library's pseudo-quad subclasses.

    Returns
    -------
    dict
        ``vertices``, ``faces`` and ``poles``.
    """
    keys = list(mesh.vertices())
    index = dict((key, i) for i, key in enumerate(keys))
    vertices = [list(mesh.vertex_coordinates(key)) for key in keys]
    faces = [[index[key] for key in mesh.face_vertices(fkey)]
             for fkey in mesh.faces()]
    return {'vertices': vertices, 'faces': faces, 'poles': _pole_points(mesh)}


def mesh_from_wire(data, cls=None):
    """Rebuild a mesh from the spool.

    Parameters
    ----------
    data : dict
        As produced by :func:`mesh_to_wire`.
    cls : type, optional
        The class to build. Defaults to
        :class:`~compas_singular.datastructures.QuadMesh`, or to
        :class:`~compas_singular.datastructures.CoarsePseudoQuadMesh` when the
        payload carries poles -- a mesh with poles is not a plain quad mesh and
        building it as one silently discards the collapsed sides.

    Returns
    -------
    Mesh
    """
    if not isinstance(data, dict):
        raise WireError('mesh payload is {}, not a dict'.format(type(data).__name__))
    vertices = data.get('vertices')
    faces = data.get('faces')
    if not isinstance(vertices, list) or not isinstance(faces, list):
        raise WireError("mesh payload needs 'vertices' and 'faces' lists; got "
                        "keys {}".format(sorted(data)))
    if not vertices:
        raise WireError('mesh payload has no vertices')

    count = len(vertices)
    for i, face in enumerate(faces):
        if not isinstance(face, list) or len(face) < 3:
            raise WireError('face {} is not a list of at least 3 indices'.format(i))
        for key in face:
            if not isinstance(key, int) or key < 0 or key >= count:
                raise WireError('face {} refers to vertex {}, but the payload '
                                'has {} vertices'.format(i, key, count))

    poles = data.get('poles') or []

    # Imported here, not at module scope: the bridge is imported by the Rhino
    # side too, and keeping it free of package-level imports is what lets it be
    # loaded after ``ensure_paths()`` without dragging the datastructures in.
    if cls is None:
        if poles:
            from ...datastructures import CoarsePseudoQuadMesh as cls
        else:
            from ...datastructures import QuadMesh as cls

    builder = getattr(cls, 'from_vertices_and_faces_with_poles', None)
    if poles and builder is not None:
        return builder(vertices, faces, poles)
    return cls.from_vertices_and_faces(vertices, faces)
