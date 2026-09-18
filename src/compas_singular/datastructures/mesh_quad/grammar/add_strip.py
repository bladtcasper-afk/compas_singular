"""**Add a strip along a polyedge. Topology only -- nothing here moves a vertex.**

``add_strip`` splits every vertex of the polyedge in two and fills the band
between the copies, which is thesis 5.3.1. Both copies are created AT THE
POSITION of the vertex they replace, so the new strip has **zero width** until
something separates them. That separation is deliberately not done here: it is a
geometric decision, it differs per caller, and hiding it inside the grammar made
it impossible to add a strip without also moving the rest of the mesh.

Openers in the codebase, for reference:

* :meth:`~compas_singular.editing.CoarseEditor._open_strip` -- exact rule, the
  new pair goes at one third and two thirds of the span the polyedge crossed.
  Corners on the layout boundary are projected back onto the wall.
* :meth:`~compas_singular.editing.DenseMeshEditor.relax` -- centroid smoothing
  with every boundary ring held except the new pair.

**A caller that opens neither leaves coincident vertices**, which survive
``is_manifold`` but collapse a face to zero area -- and at float32 that is enough
to make a whole Rhino mesh fail to bake.

This module replaced a second implementation that lived in ``grammar_pattern.py``
until 2026-09-18. That one duplicated the whole polyedge in a single pass, keyed
by vertex, so a polyedge visiting the same vertex twice deleted it twice; it also
ended with a 20-iteration constrained smooth of the WHOLE mesh. The walk below
consumes the polyedge one vertex at a time and re-derives the remainder, which is
what makes U-turns and self-crossings tractable -- see ``add_strip``.
"""
from compas.topology import breadth_first_paths
from compas.datastructures.mesh.operations.substitute import mesh_substitute_vertex_in_faces
from compas.itertools import pairwise


__all__ = [
    'add_strip',
    'add_strips',
    'split_strip',
    'split_strips',
    'strip_polyedge_update',
    'is_polyedge_valid_for_strip_addition',
]


def add_strips(mesh, polyedges):
    """Add a strip along each polyedge, in order.

    The polyedges still to come are re-derived after every insertion: an
    insertion renumbers and replaces the very vertices they are written in terms
    of, so a polyedge collected before the first strip is meaningless after it.

    Parameters
    ----------
    mesh : QuadMesh
        A quad mesh, with ``attributes['strips']`` already collected.
    polyedges : list[list[int]]
        Polyedges, each a list of vertex keys.

    Returns
    -------
    list
        The new strip keys, in the order the strips were added.
    """
    pending = list(polyedges)
    new_skeys = []

    while pending:
        polyedge = pending.pop()
        # ``list(...)``: ``add_strip`` consumes the polyedge it is given.
        new_skey, old_to_new = add_strip(mesh, list(polyedge))
        new_skeys.append(new_skey)
        pending = [strip_polyedge_update(mesh, pending_polyedge, old_to_new)
                   for pending_polyedge in pending]

    return new_skeys


def add_strip(mesh, polyedge):
    """**Add a strip along** ``polyedge``. Topology only -- the strip has zero width.

    Each vertex ``Vi`` of the polyedge becomes two, ``Vi`` is substituted by the
    left copy in the faces on one side and by the right copy in those on the
    other, and the band between them is the new strip (thesis 5.3.1). Strip
    LABELS are preserved, so densities set per strip survive this.

    The polyedge is consumed. Pass a copy to keep yours -- and note that
    ``mesh.attributes['polyedges']`` hands out its own lists, so passing one
    straight in corrupts it.

    Parameters
    ----------
    mesh : QuadMesh
        A quad mesh. ``attributes['strips']`` must be non-empty --
        ``update_strip_data`` does ``max(...) + 1`` and raises on an empty dict.
    polyedge : list[int]
        Vertex keys in order, each joined to the next by an edge. Either closed,
        or with both ends on the boundary.

    Returns
    -------
    tuple[int, dict]
        The new strip key, and ``{old vertex: (left copy, right copy)}``. The
        left and right polyedges are the values of that map read in polyedge
        order; the pairs in it are exactly the vertices an opener has to
        separate.

    Notes
    -----
    The walk takes one vertex per iteration and rebuilds what is left of the
    polyedge through the two vertices just created
    (``polyedge_from_to_via_vertices``). That is what lets a polyedge revisit a
    vertex -- a U-turn or a self-crossing -- instead of deleting it twice.
    Support for those is INCOMPLETE: the walk handles them, but
    ``update_strip_data`` below does not, and a self-crossing polyedge such as
    ``[1, 2, 8, 2, 3]`` still raises there. Guard with
    ``is_polyedge_valid_for_strip_addition`` and keep polyedges simple.
    """
    full_updated_polyedge = []
    # store data
    left_polyedge = []
    right_polyedge = []
    new_faces = []

    # exception if closed
    is_closed = polyedge[0] == polyedge[-1]
    if is_closed:
        polyedge.pop()

    k = -1
    count = len(polyedge) * 2
    while count and len(polyedge) > 0:
        k += 1
        count -= 1

        # select u, v, w if not closed
        if not is_closed:
            # u
            if len(new_faces) != 0:
                u1, u2 = left_polyedge[-1], right_polyedge[-1]
            else:
                u1, u2 = None, None
            # v
            v = polyedge.pop(0)
            full_updated_polyedge.append(v)
            # w
            if len(polyedge) != 0:
                w = polyedge[0]
            else:
                w = None

        # select u, v, w if closed
        else:
            # u
            if len(new_faces) != 0:
                u1, u2 = left_polyedge[-1], right_polyedge[-1]
            else:
                u1 = polyedge[-1]  # artificial u1
            # v
            v = polyedge.pop(0)
            full_updated_polyedge.append(v)
            # w
            if len(polyedge) != 0:
                w = polyedge[0]
            else:
                w = left_polyedge[0]

        # add new vertices
        faces = sort_faces(mesh, u1, v, w)
        v1, v2 = mesh.add_vertex(attr_dict=mesh.vertex[v]), mesh.add_vertex(attr_dict=mesh.vertex[v])

        if type(faces[0]) == list:
            faces_1, faces_2 = faces
        else:
            # exception necessary for U-turns
            if faces[0] in mesh.vertex_faces(left_polyedge[-2]):
                faces_1 = faces
                faces_2 = []
            else:
                faces_1 = []
                faces_2 = faces
        mesh_substitute_vertex_in_faces(mesh, v, v1, faces_1)
        mesh_substitute_vertex_in_faces(mesh, v, v2, faces_2)
        mesh.delete_vertex(v)
        left_polyedge.append(v1)
        right_polyedge.append(v2)

        # add new faces, different if at the start, end or main part of the polyedge
        if len(new_faces) == 0:
            if not is_closed:
                new_faces.append(mesh.add_face([v1, w, v2]))
            else:
                new_faces.append(mesh.add_face([v1, v2, u1]))
                new_faces.append(mesh.add_face([v1, w, v2]))
        elif len(polyedge) == 0:
            if not is_closed:
                u1, u2 = left_polyedge[-2], right_polyedge[-2]
                face = new_faces.pop()
                mesh.delete_face(face)
                new_faces.append(mesh.add_face([u1, v1, v2, u2]))
            else:
                u1, u2 = left_polyedge[-2], right_polyedge[-2]
                face = new_faces.pop()
                mesh.delete_face(face)
                # ``[u1, v1, v2, u2]``, as every other rung of the strip. This was
                # ``[v1, u1, u2, v2]`` -- the same quad REVERSED -- so it shared a
                # same-direction halfedge with both its neighbours: the mesh still
                # passed ``is_manifold`` but every strip walk through that face went
                # wrong. The closing face below is right as it is: from the last
                # vertex to the first, ``v`` plays the role ``u`` plays here.
                new_faces.append(mesh.add_face([u1, v1, v2, u2]))
                face = new_faces.pop(0)
                mesh.delete_face(face)
                u1, u2 = left_polyedge[0], right_polyedge[0]
                new_faces.append(mesh.add_face([v1, u1, u2, v2]))
                # NOT followed by ``mesh_substitute_vertex_in_faces(mesh, v, v1/v2)``.
                # ``v`` was deleted above, and those calls pass no ``fkeys``, so they
                # default to EVERY face in the mesh and delete-and-re-add each one --
                # twice -- to substitute a vertex that no face references any more.
                # The ring is already closed by the two faces added above.
        else:
            face = new_faces.pop()
            mesh.delete_face(face)
            new_faces.append(mesh.add_face([u1, v1, v2, u2]))
            new_faces.append(mesh.add_face([v1, w, v2]))

        # update
        updated_polyedge = []
        via_vkeys = [v1, v2]
        for i, vkey in enumerate(polyedge):
            if vkey != v:
                updated_polyedge.append(vkey)
            else:
                from_vkey = polyedge[i - 1]
                to_vkey = polyedge[i + 1]
                updated_polyedge += polyedge_from_to_via_vertices(mesh, from_vkey, to_vkey, via_vkeys)[1:-1]
        polyedge = updated_polyedge

    old_vkeys_to_new_vkeys = {u0: (u1, u2) for u0, u1, u2 in zip(full_updated_polyedge, left_polyedge, right_polyedge)}

    n = update_strip_data(mesh, full_updated_polyedge, old_vkeys_to_new_vkeys, closed=is_closed)
    return n, old_vkeys_to_new_vkeys


def split_strip(mesh, skey, n=2):
    """Refine a strip into ``n`` strips. Topology only -- see the module note.

    Returns
    -------
    list
        The existing strip key, followed by the ``n - 1`` new ones.
    """
    return [skey] + [add_strip(mesh, list(mesh.strip_side_polyedges(skey)[0]))[0]
                     for _ in range(n - 1)]


def split_strips(mesh, skey_to_n):
    """Refine several strips. Topology only -- see the module note.

    Parameters
    ----------
    mesh : QuadMesh
        A quad mesh.
    skey_to_n : dict
        Strip keys pointing to the number of strips to refine each one into.

    Returns
    -------
    dict
        Each split strip key, pointing to the keys of the strips refining it.
    """
    return {skey: split_strip(mesh, skey, n) for skey, n in skey_to_n.items()}


def update_strip_data(mesh, full_updated_polyedge, old_vkeys_to_new_vkeys, closed=False):
    """Bring ``attributes['strips']`` up to date after a strip was added.

    ``closed`` matters: a closed polyedge arrives here without its repeated end
    vertex, so ``pairwise`` used to skip the closing edge. The strip crossing it
    was then treated as a PARALLEL strip and reached for a vertex the insertion
    had deleted -- ``KeyError`` on every closed polyedge, measured 3 of 3 on the
    rings around a hole in a densified plate.
    """
    sequence = full_updated_polyedge + full_updated_polyedge[:1] if closed else full_updated_polyedge

    # orthogonal strips
    orth_to_update = {}
    for old_u, old_v in pairwise(sequence):
        new_u = old_vkeys_to_new_vkeys[old_u][0]
        new_v = old_vkeys_to_new_vkeys[old_v][0]
        skey = mesh.edge_strip((old_u, old_v))
        edges = mesh.collect_strip(new_u, new_v)
        orth_to_update[skey] = edges
    for skey, edges in orth_to_update.items():
        mesh.attributes['strips'][skey] = edges

    # parallel strips
    paral_to_update = {}
    for skey, edges in mesh.attributes['strips'].items():
        if skey not in orth_to_update:
            for u, v in edges:
                if u in old_vkeys_to_new_vkeys:
                    u, v = v, u
                elif v in old_vkeys_to_new_vkeys:
                    u, v = u, v
                else:
                    continue
                new_v = [vkey for vkey in old_vkeys_to_new_vkeys[v] if vkey in mesh.halfedge[u]][0]
                new_edges = mesh.collect_strip(u, new_v)
                paral_to_update[skey] = new_edges
                break
    for skey, edges in paral_to_update.items():
        mesh.attributes['strips'][skey] = edges

    # self strip
    n = max(mesh.attributes['strips']) + 1
    strip_edges = [tuple(old_vkeys_to_new_vkeys[vkey]) for vkey in full_updated_polyedge]
    mesh.attributes['strips'][n] = strip_edges

    return n


def strip_polyedge_update(mesh, polyedge, vertex_modifications):
    """Rewrite ``polyedge`` in terms of the vertices an insertion left behind.

    Parameters
    ----------
    mesh : QuadMesh
        A quad mesh.
    polyedge : list[int]
        A polyedge, as a list of the OLD vertex keys.
    vertex_modifications : dict
        Old vertex keys pointing to the new ones that replaced them -- the
        second return value of :func:`add_strip`.

    Returns
    -------
    list
        The shortest path that visits the replacements in the original order.
    """
    closed = polyedge[0] == polyedge[-1]

    if closed:
        polyedge = polyedge[:-1]

    # update polyedge with candidate vertices
    polyedge_modifications = {vkey: (vertex_modifications[vkey] if vkey in vertex_modifications else [vkey]) for vkey in polyedge}
    # list all candidate vertices to form new polyedge
    candidate_vertices = tuple(set([vkey for vkeys in polyedge_modifications.values() for vkey in vkeys]))
    # adjacency restricted to candidate vertices
    adjacency = {vkey: [nbr for nbr in mesh.vertex_neighbors(vkey) if nbr in candidate_vertices] for vkey in mesh.vertices() if vkey in candidate_vertices}

    # for each combination of boundary vertex extremities, get the shortest
    # valid path through the modified vertices of the polyedges
    shortest_polyedge = None
    # start vertices on boundary
    for vkey_start in polyedge_modifications[polyedge[0]]:
        if mesh.is_vertex_on_boundary(vkey_start):
            # end vertices on boundary
            for vkey_end in polyedge_modifications[polyedge[-1]]:
                if mesh.is_vertex_on_boundary(vkey_end):
                    # iterate through all paths between start and end vertices
                    # starting by the shortest
                    for candidate_polyedge in breadth_first_paths(adjacency, vkey_start, vkey_end):
                        is_valid = True
                        # if was initailly closed, make sure that temporary end
                        # is adjacent to start
                        if closed:
                            if candidate_polyedge[0] not in mesh.vertex_neighbors(candidate_polyedge[-1]):
                                continue
                        # check that vertices in path come from the modified
                        # vertices of the polyedge in the same order
                        i = 0
                        for vkey in candidate_polyedge:
                            if vkey in polyedge_modifications[polyedge[i]]:
                                continue
                            elif vkey in polyedge_modifications[polyedge[i + 1]]:
                                i += 1
                            else:
                                is_valid = False
                                break
                        # update if shorter
                        if is_valid:
                            if shortest_polyedge is None or len(shortest_polyedge) > len(candidate_polyedge):
                                shortest_polyedge = candidate_polyedge
                            break

    if closed:
        shortest_polyedge.append(shortest_polyedge[0])

    return shortest_polyedge


def sort_faces(mesh, u, v, w):

    sorted_faces = [[], []]
    k = 0
    vertex_faces = mesh.vertex_faces(v, ordered=True, include_none=True)
    f0 = mesh.halfedge[w][v] if w is not None else None
    i0 = vertex_faces.index(f0)
    vertex_faces = vertex_faces[i0:] + vertex_faces[:i0]

    for face in vertex_faces:

        if face is not None:
            sorted_faces[k].append(face)

        if u is None:
            if face is None:
                k = 1 - k
        elif face == mesh.halfedge[v][u]:
            k = 1 - k

    # indeterminate exception if u == w
    if u == w:
        return [face for faces in sorted_faces for face in faces]
    else:
        return sorted_faces


def adjacency_from_to_via_vertices(mesh, from_vkey, to_vkey, via_vkeys):
    # get mesh adjacency constraiend to from_vkey and via_keys, via_keys and via_keys, and via_vkeys and to_vkey

    all_vkeys = set([from_vkey, to_vkey] + via_vkeys)
    adjacency = {}
    for vkey, nbrs in mesh.adjacency.items():
        if vkey not in all_vkeys:
            continue
        else:
            sub_adj = {}
            for nbr, face in nbrs.items():
                if nbr not in all_vkeys or (vkey == from_vkey and nbr == to_vkey) or (vkey == to_vkey and nbr == from_vkey):
                    continue
                else:
                    sub_adj.update({nbr: face})
            adjacency.update({vkey: sub_adj})
    return adjacency


def polyedge_from_to_via_vertices(mesh, from_vkey, to_vkey, via_vkeys):
    # return shortest polyedge from_vkey to_vkey via_vkeys

    adjacency = adjacency_from_to_via_vertices(mesh, from_vkey, to_vkey, via_vkeys)
    return next(breadth_first_paths(adjacency, from_vkey, to_vkey))


def is_polyedge_valid_for_strip_addition(mesh, polyedge):
    if len(polyedge) > 2:
        if polyedge[0] == polyedge[-1] or (mesh.is_vertex_on_boundary(polyedge[0]) and mesh.is_vertex_on_boundary(polyedge[-1])):
            return True
    return False
