"""**Delete a strip, welding the two sides together. Plus the queries that price one.**

Deleting a strip is not the mirror image of adding one. Adding can leave the
geometry untouched because the new vertices start on top of the old; deleting
MERGES vertices, and the survivor has to be placed somewhere. ``delete_strip``
therefore does move vertices -- but by a fixed rule, not by relaxation:

1. vertices the deletion disconnected, if any;
2. otherwise vertices that were on a boundary before it;
3. otherwise all of them.

Preferring the old boundary is what keeps a welded outline where it was instead
of pulling it inward.

**Pre-splitting to save a boundary is NOT done here.** A boundary represented by
fewer than three edges after a deletion collapses, and the remedy is to refine
the strips that survive on it (thesis 5.3.2, Fig 5.18). That is a policy, and it
belongs to whoever is driving:
:func:`strips_to_split_to_prevent_boundary_collapse` works out what to split and
:func:`~compas_singular.datastructures.mesh_quad.grammar.add_strip.split_strips`
performs it, so a caller that wants it does the two steps itself. It used to be a
``preserve_boundaries`` flag on this function, which meant the policy was decided
in two places at once -- here and in
:meth:`~compas_singular.editing.editor.MeshEditor._split_strips` -- and the
splitting it did silently smoothed the whole mesh.

Until 2026-09-18 a second, simpler implementation of this lived here while the
one below lived in ``grammar_pattern.py``. The simpler one merged to a plain
centroid and knew nothing about collateral deletions or poles.
"""
from compas.datastructures.mesh.operations.substitute import mesh_substitute_vertex_in_faces
from compas.geometry import centroid_points
from compas.itertools import pairwise
from compas.topology import connected_components

from compas_singular.datastructures.network import Network


__all__ = [
    'delete_strip',
    'delete_strips',
    'collateral_strip_deletions',
    'total_boundary_deletions',
    'strips_to_split_to_prevent_boundary_collapse',
]


def delete_strips(mesh, skeys):
    """Delete several strips.

    Strip keys are re-checked as the deletions go, because deleting one strip can
    delete others with it -- see :func:`collateral_strip_deletions`.
    """
    for skey in skeys:
        if skey in list(mesh.strips()):
            delete_strip(mesh, skey)


def delete_strip(mesh, skey):
    """**Delete the strip** ``skey``, welding the faces either side together.

    Parameters
    ----------
    mesh : QuadMesh
        A quad mesh.
    skey : hashable
        A strip key.

    Returns
    -------
    dict
        ``{old vertex: the vertex it was merged into}``. Empty if there was no
        such strip.

    Notes
    -----
    Strips entirely consumed by this one go with it, and ``attributes['face_pole']``
    is repointed at the surviving vertices. Vertices MOVE -- see the module note
    for the rule.
    """
    if skey not in list(mesh.strips()):
        return {}

    old_boundary_vertices = list(mesh.vertices_on_boundaries())

    # get strip data
    strip_edges = mesh.strip_edges(skey)
    strip_faces = mesh.strip_faces(skey)

    # collateral strip deletions
    collateral_deleted_strips = []
    for skey_2 in mesh.strips():
        if skey_2 == skey:
            continue
        faces_2 = mesh.strip_faces(skey_2)
        # ``faces_2`` non-empty: ``all([])`` is True, so a strip with no faces -- one
        # edge between two polygons -- would be "deleted" by every deletion.
        if faces_2 and all([fkey in strip_faces for fkey in faces_2]):
            collateral_deleted_strips.append(skey_2)

    # build network between vertices of the edges of the strip to delete to
    # get the disconnect parts of vertices to merge
    vertices = set([i for edge in strip_edges for i in edge])
    # maps between old and new indices
    old_to_new = {vkey: i for i, vkey in enumerate(vertices)}
    new_to_old = {i: vkey for i, vkey in enumerate(vertices)}
    # network
    vertex_coordinates = {i: mesh.vertex_coordinates(vkey) for i, vkey in enumerate(vertices)}
    edges = [(old_to_new[u], old_to_new[v]) for u, v in strip_edges]
    network = Network.from_nodes_and_edges(vertex_coordinates, edges)
    # disconnected parts
    parts = connected_components(network.adjacency)

    # delete strip faces
    for fkey in strip_faces:
        mesh.delete_face_in_strips(fkey)
    for fkey in strip_faces:
        mesh.delete_face(fkey)

    old_vkeys_to_new_vkeys = {}

    # merge strip edge vertices that are connected
    for part in parts:

        # move back from network vertices to mesh vertices
        vertices = [new_to_old[vkey] for vkey in part]

        # skip adding a vertex if all vertices of the part are disconnected
        if any(mesh.is_vertex_connected(vkey) for vkey in vertices):

            # get position based on disconnected vertices that used to be on
            # the boundary if any
            if any(not mesh.is_vertex_connected(vkey) for vkey in vertices):
                points = [mesh.vertex_coordinates(
                    vkey) for vkey in vertices if not mesh.is_vertex_connected(vkey)]
            # or based on old boundary vertices if any
            elif any(vkey in old_boundary_vertices for vkey in vertices):
                points = [mesh.vertex_coordinates(
                    vkey) for vkey in vertices if vkey in old_boundary_vertices]
            else:
                points = [mesh.vertex_coordinates(vkey) for vkey in vertices]

            # new vertex
            x, y, z = centroid_points(points)
            new_vkey = mesh.add_vertex(attr_dict={'x': x, 'y': y, 'z': z})
            old_vkeys_to_new_vkeys.update({old_vkey: new_vkey for old_vkey in vertices})

            # replace the old vertices
            for old_vkey in vertices:
                mesh.substitute_vertex_in_strips(old_vkey, new_vkey)
                mesh_substitute_vertex_in_faces(
                    mesh, old_vkey, new_vkey, mesh.vertex_faces(old_vkey))

        # delete the old vertices
        for old_vkey in vertices:
            mesh.delete_vertex(old_vkey)

    # delete data of deleted strip and collateral deleted strips
    del mesh.attributes['strips'][skey]
    for skey_2 in collateral_deleted_strips:
        del mesh.attributes['strips'][skey_2]

    if 'face_pole' in mesh.attributes:
        for fkey, pole in mesh.attributes['face_pole'].items():
            if fkey in mesh.attributes['face_pole']:
                if pole == mesh.attributes['face_pole'][fkey]:
                    if pole in old_vkeys_to_new_vkeys:
                        mesh.attributes['face_pole'][fkey] = old_vkeys_to_new_vkeys[pole]

    return old_vkeys_to_new_vkeys


def collateral_strip_deletions(mesh, skeys):
    """The strips that deleting ``skeys`` would delete as well. Mutates nothing."""
    deleted_fkeys = [fkey for skey in skeys for fkey in mesh.strip_faces(skey)]
    # A strip with no faces at all is not collateral -- see ``delete_strip``.
    return [skey for skey in mesh.strips() if skey not in skeys and mesh.strip_faces(skey)
            and all([fkey in deleted_fkeys for fkey in mesh.strip_faces(skey)])]


def total_boundary_deletions(mesh, skeys):
    """The boundaries that deleting ``skeys`` would collapse. Mutates nothing."""
    deleted_strips = list(skeys) + list(collateral_strip_deletions(mesh, skeys))
    deleted_boundaries = []
    deleted_edges = set([edge for skey in deleted_strips for edge in mesh.strip_edges(skey)])
    for boundary in mesh.boundaries():
        edges = [(u, v) for u, v in pairwise(boundary + boundary[:1])]
        if all([(u, v) in deleted_edges or (v, u) in deleted_edges for u, v in edges]):
            deleted_boundaries.append(boundary)
    return deleted_boundaries


def strips_to_split_to_prevent_boundary_collapse(mesh, skeys):
    """**What to split before deleting** ``skeys``, to keep every boundary alive.

    Thesis 5.3.2, Fig 5.18: a boundary collapses when fewer than three edges
    represent it after a deletion. One surviving strip on it, split it in three;
    two, split each in two, "to avoid any bias".

    Feed the result to
    :func:`~compas_singular.datastructures.mesh_quad.grammar.add_strip.split_strips`
    BEFORE calling :func:`delete_strip`. Mutates nothing.

    Returns
    -------
    dict or None
        **Tri-state, and all three states matter:**

        * ``None`` -- the collapse cannot be prevented: no strip survives on that
          boundary to refine;
        * ``{}`` -- nothing is at risk, so there is nothing to do;
        * ``{skey: n}`` -- refine these first.

        ``None`` and ``{}`` are NOT interchangeable. Flattening them together is
        what made asking to preserve boundaries refuse a deletion that risked
        none; :meth:`~compas_singular.editing.editor.MeshEditor._trial_delete`
        reads the difference as its ``can_preserve``.
    """
    to_split = {}
    for boundary in mesh.boundaries():
        non_deleted_strips = [mesh.edge_strip((u, v)) for u, v in pairwise(boundary + boundary[:1]) if mesh.edge_strip((u, v)) not in skeys]

        if len(non_deleted_strips) == 0:
            # Nothing survives on this boundary, so nothing can be refined to
            # save it. ``None``, NOT ``{}`` -- see the tri-state above.
            return None

        elif len(non_deleted_strips) == 1:
            skey = non_deleted_strips[0]
            if skey in to_split:
                if to_split[skey] == 3:
                    break
                if to_split[skey] == 2:
                    to_split[skey] = 3
            else:
                to_split[skey] = 3

        elif len(non_deleted_strips) == 2:
            for skey in non_deleted_strips:
                if skey in to_split:
                    break
            to_split.update({skey: 2 for skey in non_deleted_strips})

    return to_split
