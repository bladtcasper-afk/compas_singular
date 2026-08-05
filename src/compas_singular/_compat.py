"""Backports of COMPAS functions used by compas_singular.

compas_singular was originally written against the COMPAS 0.x/1.x API. A number
of the free functions it relied on were removed or relocated in COMPAS 2.x.
This module re-implements or re-exports them so that only *compas_singular* has
to be adapted -- COMPAS itself is left untouched.

The implementations are ported from COMPAS 1.17.9 (the last 1.x release) with
the minimal changes required for COMPAS 2.x:

* ``geometric_key`` -> :meth:`compas.tolerance.TOL.geometric_key`.
* ``pairwise`` / ``window`` / ``linspace`` -> :mod:`compas.itertools`.
* ``adjacency_from_edges`` -> :func:`compas.topology.vertex_adjacency_from_edges`
  (identical behaviour, renamed).
* Edge methods take a single ``(u, v)`` tuple in COMPAS 2.x. compas_singular's
  own :class:`~compas_singular.datastructures.mesh.mesh.Mesh` re-exposes the old
  two-argument signatures, so the ported code keeps the original
  ``mesh.edge_midpoint(u, v)`` calling style.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import re
from math import cos
from math import sin

from compas.tolerance import TOL
from compas.itertools import pairwise
from compas.topology import connected_components
from compas.topology import vertex_adjacency_from_edges as adjacency_from_edges  # noqa: F401

from compas.geometry import subtract_vectors
from compas.geometry import cross_vectors
from compas.geometry import normalize_vector
from compas.geometry import dot_vectors
from compas.geometry import length_vector
from compas.geometry import length_vector_sqrd
from compas.geometry import scale_vector
from compas.geometry import sum_vectors


__all__ = [
    "geometric_key",
    "average",
    "circle_evaluate",
    "archimedean_spiral_evaluate",
    "circle_from_points",
    "adjacency_from_edges",
    "trimesh_face_circle",
    "mesh_substitute_vertex_in_faces",
    "mesh_insert_vertex_on_edge",
    "mesh_weld",
    "meshes_join",
    "meshes_join_and_weld",
    "mesh_unweld_edges",
    "mesh_explode",
    "network_polylines",
    "network_disconnected_nodes",
    "mesh_smooth_centroid",
    "mesh_smooth_area",
]


# =============================================================================
# utilities (compas.utilities in 0.x/1.x)
# =============================================================================

def geometric_key(xyz, precision=None, sanitize=True):
    """Convert XYZ coordinates to a string key.

    In COMPAS 2.x this lives on the global tolerance object. The 1.x API also
    accepted a string precision specifier such as ``'3f'``; that form is
    translated to the integer number of decimals expected by COMPAS 2.x.
    """
    if isinstance(precision, str):
        match = re.search(r"\d+", precision)
        precision = int(match.group()) if match else 0
    return TOL.geometric_key(xyz, precision=precision, sanitize=sanitize)


def average(values):
    """Compute the arithmetic mean of a sequence of numbers."""
    values = list(values)
    return sum(values) / float(len(values))


# =============================================================================
# geometry (compas.geometry in 0.x/1.x)
# =============================================================================

def circle_evaluate(t, r, z=0):
    """Evaluate a circle of radius ``r`` at parameter ``t``."""
    return [r * cos(t), r * sin(t), z]


def archimedean_spiral_evaluate(t, a, b, z=0):
    """Evaluate an archimedean spiral ``r = a + b * theta`` at parameter ``t``."""
    return [b * t * cos(t + a), b * t * sin(t + a), z]


def circle_from_points(a, b, c):
    """Construct the circumscribed circle of three points.

    Returns ``((center, normal), radius)`` (the COMPAS 1.x return shape).
    """
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
    Aa = scale_vector(a, A)
    Bb = scale_vector(b, B)
    Cc = scale_vector(c, C)
    center = sum_vectors([Aa, Bb, Cc])
    radius = length_vector(subtract_vectors(a, center))
    return (center, normal), radius


# =============================================================================
# mesh geometry (compas.datastructures in 0.x/1.x)
# =============================================================================

def trimesh_face_circle(mesh, fkey):
    """Circumcircle of a triangular face.

    Returns ``(center, radius, normal)`` -- the ``center`` first ordering that
    compas_singular relies on (``trimesh_face_circle(mesh, fkey)[0]`` is the
    circumcentre) -- or ``None`` if the face is not a triangle.
    """
    vertices = mesh.face_vertices(fkey)
    if len(vertices) != 3:
        return None
    u, v, w = vertices
    a = mesh.vertex_coordinates(u)
    b = mesh.vertex_coordinates(v)
    c = mesh.vertex_coordinates(w)
    (center, normal), radius = circle_from_points(a, b, c)
    return center, radius, normal


# =============================================================================
# mesh operations (compas.datastructures in 0.x/1.x)
# =============================================================================

def mesh_substitute_vertex_in_faces(mesh, old_vkey, new_vkey, fkeys=None):
    """Substitute a vertex by another one in all faces (or a given set)."""
    if fkeys is None:
        fkeys = list(mesh.faces())
    for fkey in fkeys:
        face_vertices = [new_vkey if key == old_vkey else key for key in mesh.face_vertices(fkey)]
        mesh.delete_face(fkey)
        mesh.add_face(face_vertices, fkey)
    return fkeys


def mesh_insert_vertex_on_edge(mesh, u, v, vkey=None):
    """Insert a vertex in the faces adjacent to an edge, between its vertices."""
    # add new vertex if there is none or if vkey not in vertices
    if vkey is None:
        vkey = mesh.add_vertex(attr_dict={attr: xyz for attr, xyz in zip(["x", "y", "z"], mesh.edge_midpoint(u, v))})
    elif vkey not in list(mesh.vertices()):
        vkey = mesh.add_vertex(
            key=vkey,
            attr_dict={attr: xyz for attr, xyz in zip(["x", "y", "z"], mesh.edge_midpoint(u, v))},
        )

    # insert vertex
    for fkey, halfedge in zip(mesh.edge_faces(u, v), [(u, v), (v, u)]):
        if fkey is not None:
            face_vertices = mesh.face_vertices(fkey)[:]
            face_vertices.insert(face_vertices.index(halfedge[-1]), vkey)
            mesh.delete_face(fkey)
            mesh.add_face(face_vertices, fkey)

    return vkey


def mesh_weld(mesh, precision=None, cls=None):
    """Weld vertices of a mesh within a precision distance, returning a new mesh."""
    if cls is None:
        cls = type(mesh)

    key_xyz = {key: mesh.vertex_coordinates(key) for key in mesh.vertices()}
    gkey_key = {geometric_key(xyz, precision): key for key, xyz in key_xyz.items()}
    gkey_index = {gkey: index for index, gkey in enumerate(gkey_key)}

    vertices = [key_xyz[key] for gkey, key in gkey_key.items()]
    faces = [[gkey_index[geometric_key(key_xyz[key], precision)] for key in mesh.face_vertices(fkey)] for fkey in mesh.faces()]

    faces[:] = [[u for u, v in pairwise(face + face[:1]) if u != v] for face in faces]
    faces[:] = [face for face in faces if len(face) > 2]  # discard degenerate faces

    return cls.from_vertices_and_faces(vertices, faces)


def meshes_join(meshes, cls=None):
    """Join meshes without welding."""
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
    """Join and weld meshes within a precision distance."""
    return mesh_weld(meshes_join(meshes, cls=cls), precision=precision)


def mesh_unweld_edges(mesh, edges):
    """Unweld a mesh along a set of edges (in place)."""
    # set of vertices touched by the edges to unweld
    vertices = set([i for edge in edges for i in edge])

    # store changes to apply them all at once
    vertex_changes = {}

    for vkey in vertices:

        # maps between mesh face index and local network node index
        old_to_new = {nbr: i for i, nbr in enumerate(mesh.vertex_faces(vkey))}
        new_to_old = {i: nbr for i, nbr in enumerate(mesh.vertex_faces(vkey))}

        # adjacency network of the faces around the vertex, excluding adjacency
        # through the edges to unweld
        network_edges = []
        for nbr in mesh.vertex_neighbors(vkey):
            if not mesh.is_edge_on_boundary(vkey, nbr) and (vkey, nbr) not in edges and (nbr, vkey) not in edges:
                network_edges.append((old_to_new[mesh.halfedge[vkey][nbr]], old_to_new[mesh.halfedge[nbr][vkey]]))

        adjacency = adjacency_from_edges(network_edges)
        for key, values in adjacency.items():
            adjacency[key] = {value: None for value in values}
        # include non-connected face nodes
        edge_vertices = list(set([i for edge in network_edges for i in edge]))
        for i in range(len(mesh.vertex_faces(vkey))):
            if i not in edge_vertices:
                adjacency[i] = {}

        # disconnected parts around the vertex due to unwelding
        vertex_changes[vkey] = [[new_to_old[key] for key in part] for part in connected_components(adjacency)]

    for vkey, changes in vertex_changes.items():
        # replace the vertex by a new one in the faces of each disconnected part
        for change in changes:
            mesh_substitute_vertex_in_faces(mesh, vkey, mesh.add_vertex(attr_dict=mesh.vertex[vkey]), change)
        # delete the old vertex
        mesh.delete_vertex(vkey)


def _mesh_disconnected_vertices(mesh):
    return connected_components(mesh.adjacency)


def _mesh_disconnected_faces(mesh):
    parts = _mesh_disconnected_vertices(mesh)
    return [set([fkey for vkey in part for fkey in mesh.vertex_faces(vkey)]) for part in parts]


def mesh_explode(mesh, cls=None):
    """Explode a mesh into its disconnected parts."""
    if cls is None:
        cls = type(mesh)

    parts = _mesh_disconnected_faces(mesh)
    exploded_meshes = []
    for part in parts:
        vertex_keys = list(set([vkey for fkey in part for vkey in mesh.face_vertices(fkey)]))
        vertices = [mesh.vertex_coordinates(vkey) for vkey in vertex_keys]
        key_to_index = {vkey: i for i, vkey in enumerate(vertex_keys)}
        faces = [[key_to_index[vkey] for vkey in mesh.face_vertices(fkey)] for fkey in part]
        exploded_meshes.append(cls.from_vertices_and_faces(vertices, faces))
    return exploded_meshes


# =============================================================================
# network operations (compas.datastructures in 0.x/1.x)
# =============================================================================

def network_disconnected_nodes(network):
    """Get the disconnected node groups in a network."""
    return connected_components(network.adjacency)


def network_polylines(network, splits=None):
    """Join network edges into polylines.

    The polylines stop at nodes with a valency different from 2, and at any
    optional split points.
    """
    if splits is None:
        splits = []
    stop_geom_keys = set([geometric_key(xyz) for xyz in splits])

    polylines = []
    edges_to_visit = set(network.edges())

    while len(edges_to_visit) > 0:
        polyline = list(edges_to_visit.pop())

        # extend the polyline until it is closed ...
        while polyline[0] != polyline[-1]:

            # ... or until both ends are not two-valent (or are split points)
            if len(network.neighbors(polyline[-1])) != 2 or geometric_key(network.node_coordinates(polyline[-1])) in stop_geom_keys:
                polyline = list(reversed(polyline))
                if len(network.neighbors(polyline[-1])) != 2 or geometric_key(network.node_coordinates(polyline[-1])) in stop_geom_keys:
                    break

            polyline.append([nbr for nbr in network.neighbors(polyline[-1]) if nbr != polyline[-2]][0])

        for u, v in pairwise(polyline):
            if (u, v) in edges_to_visit:
                edges_to_visit.remove((u, v))
            elif (v, u) in edges_to_visit:
                edges_to_visit.remove((v, u))

        polylines.append(polyline)

    return [[network.node_coordinates(vkey) for vkey in polyline] for polyline in polylines]


# =============================================================================
# mesh smoothing (compas.datastructures free functions -> Mesh methods)
# =============================================================================

def mesh_smooth_centroid(mesh, fixed=None, kmax=100, damping=0.5, callback=None, callback_args=None):
    """Smooth a mesh by moving each free vertex to the centroid of its neighbors."""
    return mesh.smooth_centroid(fixed=fixed, kmax=kmax, damping=damping, callback=callback, callback_args=callback_args)


def mesh_smooth_area(mesh, fixed=None, kmax=100, damping=0.5, callback=None, callback_args=None):
    """Smooth a mesh by moving each free vertex to the area-weighted centroid of its faces."""
    return mesh.smooth_area(fixed=fixed, kmax=kmax, damping=damping, callback=callback, callback_args=callback_args)
