"""Add a strip along a polyedge (thesis 5.3.1) and, by default, open it by an exact thirds rule.

With ``open_strip=False`` the new vertex pairs stay coincident until the caller separates them.
"""
from __future__ import annotations

from typing import Any
from typing import Callable
from typing import Iterable
from typing import TYPE_CHECKING

from compas.topology import breadth_first_paths
from compas.datastructures.mesh.operations.substitute import mesh_substitute_vertex_in_faces
from compas.geometry import distance_point_point
from compas.itertools import pairwise

if TYPE_CHECKING:
    from compas_singular.datastructures import QuadMesh


__all__ = [
    'add_strip',
    'add_strips',
    'split_strip',
    'split_strips',
    'strip_polyedge_update',
    'is_polyedge_valid_for_strip_addition',
    'polyedge_sides',
    'open_added_strip',
]


def add_strips(mesh: QuadMesh, polyedges: list[list[int]], open_strip: bool = True, project: Callable[[list[float]], list[float]] | None = None) -> list[int]:
    """Add a strip along each polyedge in order, re-deriving the remaining polyedges after each insertion.

    Parameters
    ----------
    mesh : QuadMesh
        A quad mesh, with ``attributes['strips']`` already collected.
    polyedges : list[list[int]]
        Polyedges, each a list of vertex keys.
    open_strip : bool, optional
        Open each new strip by ``open_added_strip``. Default ``True``.
    project : callable, optional
        Passed on to ``open_added_strip``.

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
        new_skey, old_to_new = add_strip(mesh, list(polyedge),
                                         open_strip=open_strip, project=project)
        new_skeys.append(new_skey)
        pending = [strip_polyedge_update(mesh, pending_polyedge, old_to_new)
                   for pending_polyedge in pending]

    return new_skeys


def add_strip(mesh: QuadMesh, polyedge: list[int], open_strip: bool = True, project: Callable[[list[float]], list[float]] | None = None) -> tuple[int, dict[int, tuple[int, int]]]:
    """Add a strip along ``polyedge``, opened by ``open_added_strip`` by default.

    The polyedge is consumed; self-crossing polyedges still raise in ``update_strip_data``.

    Parameters
    ----------
    mesh : QuadMesh
        A quad mesh. ``attributes['strips']`` must be non-empty --
        ``update_strip_data`` does ``max(...) + 1`` and raises on an empty dict.
    polyedge : list[int]
        Vertex keys in order, each joined to the next by an edge. Either closed,
        or with both ends on the boundary.
    open_strip : bool, optional
        Give the new strip its width with ``open_added_strip``. Default
        ``True``. ``False`` is topology only: both copies of each vertex stay on
        top of the vertex they replace.
    project : callable, optional
        ``project(xyz) -> xyz`` for the new pair at a BOUNDARY vertex, e.g. back
        onto a curved wall. Without it they stay on the chord across the vertex.
        Ignored when ``open_strip`` is ``False``.

    Returns
    -------
    tuple[int, dict]
        The new strip key, and ``{old vertex: (left copy, right copy)}``. The
        left and right polyedges are the values of that map read in polyedge
        order; the pairs in it are exactly the vertices an opener has to
        separate.
    """
    if open_strip:
        # Measured BEFORE the walk: it deletes every vertex of the polyedge, so
        # the neighbours that define each side have to be read off the mesh as
        # it still is.
        sides = polyedge_sides(mesh, polyedge)
        boundary = _boundary_vertex_set(mesh)
        on_boundary = set(vkey for vkey in polyedge if vkey in boundary)

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
    if open_strip:
        open_added_strip(mesh, old_vkeys_to_new_vkeys, sides, on_boundary, project=project)
    return n, old_vkeys_to_new_vkeys


def split_strip(mesh: QuadMesh, skey: int, n: int = 2, open_strip: bool = True, project: Callable[[list[float]], list[float]] | None = None) -> list[int]:
    """Refine a strip into ``n`` strips. ``open_strip``/``project`` as in ``add_strip``.

    Returns
    -------
    list
        The existing strip key, followed by the ``n - 1`` new ones.
    """
    return [skey] + [add_strip(mesh, list(mesh.strip_side_polyedges(skey)[0]),
                               open_strip=open_strip, project=project)[0]
                     for _ in range(n - 1)]


def split_strips(mesh: QuadMesh, skey_to_n: dict[int, int], open_strip: bool = True, project: Callable[[list[float]], list[float]] | None = None) -> dict[int, list[int]]:
    """Refine several strips. ``open_strip``/``project`` as in ``add_strip``.

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
    return {skey: split_strip(mesh, skey, n, open_strip=open_strip, project=project)
            for skey, n in skey_to_n.items()}


def polyedge_sides(mesh: QuadMesh, polyedge: list[int]) -> dict[int, tuple[list[int], list[int]] | None]:
    """``{vertex: (left neighbours, right neighbours)}`` across the polyedge, in XY."""
    closed = polyedge[0] == polyedge[-1]
    seq = polyedge[:-1] if closed else polyedge
    count = len(seq)
    out = {}
    for i, vkey in enumerate(seq):
        if closed:
            prev, nxt = seq[(i - 1) % count], seq[(i + 1) % count]
        else:
            prev = seq[i - 1] if i > 0 else None
            nxt = seq[i + 1] if i < count - 1 else None
        point = mesh.vertex_coordinates(vkey)
        a = mesh.vertex_coordinates(prev) if prev is not None else point
        b = mesh.vertex_coordinates(nxt) if nxt is not None else point
        dx, dy = b[0] - a[0], b[1] - a[1]
        if dx * dx + dy * dy < 1e-18:
            out[vkey] = None
            continue
        left, right = [], []
        for nbr in mesh.vertex_neighbors(vkey):
            if nbr == prev or nbr == nxt:
                continue
            q = mesh.vertex_coordinates(nbr)
            side = dx * (q[1] - point[1]) - dy * (q[0] - point[0])
            (left if side > 0.0 else right).append(nbr)
        out[vkey] = (left, right)
    return out


def open_added_strip(
    mesh: QuadMesh,
    old_to_new: dict[int, tuple[int, int]],
    sides: dict[int, tuple[list[int], list[int]] | None],
    on_boundary: Iterable[int] = (),
    project: Callable[[list[float]], list[float]] | None = None,
) -> int:
    """Give a new strip its width: place each new pair at thirds of the span across the old vertex.

    Only the new pairs move; boundary pairs go through ``project``. Returns the number repositioned.

    Parameters
    ----------
    mesh : QuadMesh
        The mesh after a topology-only ``add_strip``.
    old_to_new : dict
        ``{old vertex: (copy, copy)}``, the second return value of ``add_strip``.
    sides : dict
        ``polyedge_sides`` of the polyedge, measured BEFORE ``add_strip``.
    on_boundary : iterable, optional
        Old vertices that were on the boundary.
    project : callable, optional
        ``project(xyz) -> xyz`` applied to the new pair at a boundary vertex.

    Returns
    -------
    int
        The number of pairs actually repositioned. Planar (xy): z is set to 0.
    """
    on_boundary = set(on_boundary)
    opened = 0
    for old, pair in old_to_new.items():
        info = sides.get(old)
        if info is None or len(pair) != 2:
            continue
        left_keys, right_keys = info
        first, second = pair

        # Which copy took which side is MEASURED on the result rather than
        # assumed from the walk's left/right convention: the copy adjacent to a
        # left-hand neighbour is the left one.
        if any(k in mesh.halfedge.get(first, {}) for k in left_keys):
            low, high = first, second
        elif any(k in mesh.halfedge.get(second, {}) for k in left_keys):
            low, high = second, first
        elif any(k in mesh.halfedge.get(first, {}) for k in right_keys):
            low, high = second, first
        elif any(k in mesh.halfedge.get(second, {}) for k in right_keys):
            low, high = first, second
        else:
            continue

        here = mesh.vertex_coordinates(low)
        a = _centroid([mesh.vertex_coordinates(k) for k in left_keys], here)
        b = _centroid([mesh.vertex_coordinates(k) for k in right_keys], here)
        if distance_point_point(a, b) <= 1e-12:
            continue

        one = [a[0] + (b[0] - a[0]) / 3.0, a[1] + (b[1] - a[1]) / 3.0, 0.0]
        two = [a[0] + 2.0 * (b[0] - a[0]) / 3.0,
               a[1] + 2.0 * (b[1] - a[1]) / 3.0, 0.0]
        if project is not None and old in on_boundary:
            one, two = project(one), project(two)

        mesh.vertex_attributes(low, 'xyz', [one[0], one[1], 0.0])
        mesh.vertex_attributes(high, 'xyz', [two[0], two[1], 0.0])
        opened += 1
    return opened


def _centroid(points: list[list[float]], fallback: list[float]) -> list[float]:
    """The average of ``points``, or ``fallback`` when there are none."""
    if not points:
        return list(fallback)
    n = float(len(points))
    return [sum(p[0] for p in points) / n, sum(p[1] for p in points) / n, 0.0]


def _boundary_vertex_set(mesh: QuadMesh) -> set[int]:
    """Every vertex with a faceless halfedge -- holes included, no loop walk.

    The same rule as ``editing.editor.boundary_vertex_set``, repeated here
    because the grammar must not import from ``editing``.
    """
    out = set()
    for u, nbrs in mesh.halfedge.items():
        for v, fkey in nbrs.items():
            if fkey is None:
                out.add(u)
                out.add(v)
    return out


def update_strip_data(mesh: QuadMesh, full_updated_polyedge: list[int], old_vkeys_to_new_vkeys: dict[int, tuple[int, int]], closed: bool = False) -> int:
    """Bring ``attributes['strips']`` up to date after a strip was added, closing edge included."""
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


def strip_polyedge_update(mesh: QuadMesh, polyedge: list[int], vertex_modifications: dict[int, tuple[int, int]]) -> list[int]:
    """Rewrite ``polyedge`` in terms of the vertices an insertion left behind.

    Parameters
    ----------
    mesh : QuadMesh
        A quad mesh.
    polyedge : list[int]
        A polyedge, as a list of the OLD vertex keys.
    vertex_modifications : dict
        Old vertex keys pointing to the new ones that replaced them -- the
        second return value of ``add_strip``.

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


def sort_faces(mesh: QuadMesh, u: int | None, v: int, w: int | None) -> list[list[int]] | list[int]:

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


def adjacency_from_to_via_vertices(mesh: QuadMesh, from_vkey: int, to_vkey: int, via_vkeys: list[int]) -> dict[int, dict[int, Any]]:
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


def polyedge_from_to_via_vertices(mesh: QuadMesh, from_vkey: int, to_vkey: int, via_vkeys: list[int]) -> list[int]:
    # return shortest polyedge from_vkey to_vkey via_vkeys

    adjacency = adjacency_from_to_via_vertices(mesh, from_vkey, to_vkey, via_vkeys)
    return next(breadth_first_paths(adjacency, from_vkey, to_vkey))


def is_polyedge_valid_for_strip_addition(mesh: QuadMesh, polyedge: list[int]) -> bool:
    if len(polyedge) > 2:
        if polyedge[0] == polyedge[-1] or (mesh.is_vertex_on_boundary(polyedge[0]) and mesh.is_vertex_on_boundary(polyedge[-1])):
            return True
    return False
