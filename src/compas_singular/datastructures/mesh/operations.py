from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

from compas.geometry import cross_vectors
from compas.geometry import dot_vectors
from compas.geometry import length_vector
from compas.geometry import length_vector_sqrd
from compas.geometry import normalize_vector
from compas.geometry import scale_vector
from compas.geometry import subtract_vectors
from compas.geometry import sum_vectors
from compas.itertools import pairwise
from compas.tolerance import TOL

__all__ = [
    'mesh_move_vertex_by',
    'mesh_move_by',
    'mesh_move_vertices_by',
    'mesh_move_vertex_to',
    'mesh_move_vertices_to',
    'trimesh_face_circle',
    'mesh_weld',
    'meshes_join',
    'meshes_join_and_weld',
]


def mesh_move_vertex_by(mesh, vector, vkey):
    """Move a mesh vertex by a vector.

    Parameters
    ----------
    mesh : Mesh
        A mesh.
    vector : list
        An XYZ vector.
    vkey : hashable
        A vertex key.
    """

    mesh.vertex[vkey]['x'] += vector[0]
    mesh.vertex[vkey]['y'] += vector[1]
    mesh.vertex[vkey]['z'] += vector[2]


def mesh_move_by(mesh, vector):
    """Move a mesh by a vector.

    Parameters
    ----------
    mesh : Mesh
        A mesh.
    vector : list
        An XYZ vector.
    """

    for vkey in mesh.vertices():
        mesh_move_vertex_by(mesh, vector, vkey)


def mesh_move_vertices_by(mesh, key_to_vector):
    """Move mesh vertices by different vectors.

    Parameters
    ----------
    mesh : Mesh
        A mesh.
    key_to_vector : dict
        A dictionary of vertex keys pointing to vectors.
    """

    for vkey, vector in key_to_vector.items():
        mesh_move_vertex_by(mesh, vector, vkey)


def mesh_move_vertex_to(mesh, point, vkey):
    """Move a mesh vertex to a point.

    Parameters
    ----------
    mesh : Mesh
        A mesh.
    point : list
        An XYZ point.
    vkey : hashable
        A vertex key.
    """

    mesh.vertex[vkey]['x'] = point[0]
    mesh.vertex[vkey]['y'] = point[1]
    mesh.vertex[vkey]['z'] = point[2]


def mesh_move_vertices_to(mesh, key_to_point):
    """Move mesh vertices to different points.

    Parameters
    ----------
    mesh : Mesh
        A mesh.
    key_to_point : dict
        A dictionary of vertex keys pointing to points.
    """

    for vkey, point in key_to_point.items():
        mesh_move_vertex_to(mesh, point, vkey)


def trimesh_face_circle(mesh, fkey):
    """Circumcircle of a triangular face.

    Parameters
    ----------
    mesh : Mesh
        A mesh.
    fkey : hashable
        A face key.

    Returns
    -------
    tuple or None
        ``(center, radius, normal)`` -- centre first, which is the ordering
        compas_singular relies on -- or None if the face is not a triangle.

    Notes
    -----
    This is deliberately not :meth:`compas.datastructures.Mesh.face_circle`.
    That one is a numpy least-squares *bestfit* through the face coordinates and
    returns a ``Circle`` for any face, never None. Callers here compare
    circumcentres by ``TOL.geometric_key`` equality and use one as a dict key
    that is looked up again later (see ``Skeleton.lines`` and
    ``SkeletonDecomposition.decomposition_polylines``), so the centre has to
    come out of the same closed-form arithmetic every time. A bestfit through
    three points is analytically the circumcircle but not bit-identical to the
    closed form, and a three-decimal geometric key is exactly where that
    difference surfaces -- as a KeyError, not as a slightly wrong number.
    """
    vertices = mesh.face_vertices(fkey)
    if len(vertices) != 3:
        return None

    u, v, w = vertices
    a = mesh.vertex_coordinates(u)
    b = mesh.vertex_coordinates(v)
    c = mesh.vertex_coordinates(w)

    ab = subtract_vectors(b, a)
    cb = subtract_vectors(b, c)
    ba = subtract_vectors(a, b)
    ca = subtract_vectors(a, c)
    ac = subtract_vectors(c, a)
    bc = subtract_vectors(c, b)

    normal = normalize_vector(cross_vectors(ab, ac))
    d = 2 * length_vector_sqrd(cross_vectors(ba, cb))
    A = length_vector_sqrd(cb) * dot_vectors(ba, ca) / d
    B = length_vector_sqrd(ca) * dot_vectors(ab, cb) / d
    C = length_vector_sqrd(ba) * dot_vectors(ac, bc) / d
    center = sum_vectors([scale_vector(a, A), scale_vector(b, B), scale_vector(c, C)])
    radius = length_vector(subtract_vectors(a, center))

    return center, radius, normal


def mesh_weld(mesh, precision=None, cls=None):
    """Weld vertices of a mesh within a precision distance, returning a new mesh.

    Parameters
    ----------
    mesh : Mesh
        A mesh.
    precision : int, optional
        Decimals of the geometric key used to match vertices.
        Defaults to the COMPAS tolerance precision.
    cls : type, optional
        The class of the resulting mesh. Defaults to the class of ``mesh``.

    Returns
    -------
    Mesh
        A new welded mesh.

    Notes
    -----
    This is deliberately not :meth:`compas.datastructures.Mesh.weld`. That one
    welds in place, returns None, keeps the original vertex keys, and leaves
    faces that collapse below three vertices in the mesh. Callers here need a
    *new* mesh (``SkeletonDecomposition.quadrangulate_polygonal_faces``) and
    need the collapse (the pseudo-quad pole bookkeeping in
    ``CoarsePseudoQuadMesh.densification`` relies on a pole quad ``[a, b, c, c]``
    becoming a triangle).
    """
    if cls is None:
        cls = type(mesh)

    key_xyz = {key: mesh.vertex_coordinates(key) for key in mesh.vertices()}
    gkey_key = {TOL.geometric_key(xyz, precision): key for key, xyz in key_xyz.items()}
    gkey_index = {gkey: index for index, gkey in enumerate(gkey_key)}

    vertices = [key_xyz[key] for gkey, key in gkey_key.items()]
    faces = [[gkey_index[TOL.geometric_key(key_xyz[key], precision)] for key in mesh.face_vertices(fkey)] for fkey in mesh.faces()]

    faces[:] = [[u for u, v in pairwise(face + face[:1]) if u != v] for face in faces]
    faces[:] = [face for face in faces if len(face) > 2]  # discard degenerate faces

    return cls.from_vertices_and_faces(vertices, faces)


def meshes_join(meshes, cls=None):
    """Join meshes without welding, returning a new mesh.

    Parameters
    ----------
    meshes : list
        A list of meshes.
    cls : type, optional
        The class of the resulting mesh. Defaults to the class of the first mesh.

    Returns
    -------
    Mesh
        A new joined mesh.

    Notes
    -----
    This is deliberately not :meth:`compas.datastructures.Mesh.join`, which
    joins one other mesh in place and returns None.
    """
    if cls is None:
        cls = type(meshes[0])

    vertices = []
    faces = []
    for mesh in meshes:
        key_index = {key: len(vertices) + i for i, key in enumerate(mesh.vertices())}
        vertices += [mesh.vertex_coordinates(key) for key in mesh.vertices()]
        faces += [[key_index[key] for key in mesh.face_vertices(fkey)] for fkey in mesh.faces()]

    return cls.from_vertices_and_faces(vertices, faces)


def meshes_join_and_weld(meshes, precision=None, cls=None):
    """Join and weld meshes within a precision distance, returning a new mesh.

    Parameters
    ----------
    meshes : list
        A list of meshes.
    precision : int, optional
        Decimals of the geometric key used to match vertices.
        Defaults to the COMPAS tolerance precision.
    cls : type, optional
        The class of the resulting mesh. Defaults to the class of the first mesh.

    Returns
    -------
    Mesh
        A new joined and welded mesh.
    """
    return mesh_weld(meshes_join(meshes, cls=cls), precision=precision)


# ==============================================================================
# Main
# ==============================================================================

if __name__ == '__main__':
    pass
