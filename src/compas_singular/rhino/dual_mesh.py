"""The dual of a quad mesh, closed on the boundary with an extra layer of quads.

:func:`compas.datastructures.mesh_conway_dual` emits one dual face per primal
vertex, but only for the INTERIOR ones -- its comprehension reads
``if not mesh.is_vertex_on_boundary(vkey)``. That is right for a closed
polyhedron (its own doctest dualises a cube) and wrong for a plate or a folded
slab: the result shrinks inward by one full ring and the outline is gone.

This module keeps that function for the interior and closes the boundary
afterwards. A boundary vertex ``v`` gets the dual cell::

    [ midpoint(v, nbrs[0]) ] + [ centroid(f) for f in faces ] + [ midpoint(v, nbrs[-1]) ]

which works because COMPAS guarantees ``vertex_neighbors(v, ordered=True)``
starts and ends on the boundary, with ``vertex_faces(v, ordered=True)`` giving
the faces between them in the same angular order. The primal vertex itself is
inserted ONLY at a domain corner (valence 2), which is the one place the cell
would otherwise be a triangle:

===================  =======  =====  ==============================
primal bnd. vertex   valence  faces  dual cell
===================  =======  =====  ==============================
corner               2        1      mid + centroid + mid + v = quad
regular              3        2      mid + 2 centroids + mid = quad
boundary singularity 4        3      pentagon (unavoidable)
===================  =======  =====  ==============================

``Mesh.dual(include_boundary=True)`` also closes the boundary, and closes it on
the outline exactly -- but it inserts the primal vertex into EVERY closure face
(``duality.py:96``), so every regular boundary vertex becomes a pentagon. That
is the only thing this function does differently, and it is the whole point of
it. Measured:

===============================  ==============  ==============  ==========
input                            mesh.dual(bnd)  dual_mesh       primal area
===============================  ==============  ==============  ==========
5x4 planar grid                  16 quad, 14 pt  **30 quad**     20.0000
folded grid, floor + one wall    19 quad, 16 pt  **35 quad**     24.0000
annulus 4x16, curved boundary    48 quad, 32 pt  **80 quad**     64.2908
===============================  ==============  ==============  ==========

The cost of the all-quad boundary is paid on a CURVED wall, and the trade is not
removable: a mesh edge is a straight chord, so the boundary either runs midpoint
-> midpoint and cuts the corner at each primal vertex, or holds that vertex and
is a pentagon. On the annulus above, ``mesh.dual(include_boundary=True)`` keeps
the area exactly (64.2908) while this function loses 3.8% (61.8439). Where the
primal boundary is straight between its vertices -- a plate, a wall, a folded
slab -- the midpoints are collinear with the vertex they skip and nothing is
lost (grid: 20.0000 both ways, deviation 0.000000000000).

So: use this when every face must be a quad. Use
``Mesh.dual(include_boundary=True)`` when a curved outline matters more than the
face count.

Either way the raw dual leaves its boundary blocks HALF a block deep, because
its boundary sits on the primal's edge midpoints -- and a corner block, half in
both directions, a quarter. :func:`redistribute_blocks` slides every grid line
to even that out, and :func:`dual_mesh` runs it by default; see its docstring
for why the obvious alternative, smoothing, is wrong. Measured block area
max/min on a 5x4 grid: 4.00 raw, 1.00 redistributed.

One more caveat: **do not dualise a mesh with poles.** A pole is already a face
rather than a gathering joint, so duality trades it the wrong way and turns it
into a valence-3 joint. All-quad meshes only.

This module deliberately imports nothing from Rhino, so it can be exercised
outside the application.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from compas.datastructures import mesh_conway_dual
from compas.geometry import bestfit_plane
from compas.geometry import closest_point_on_segment
from compas.geometry import cross_vectors
from compas.geometry import distance_point_plane
from compas.geometry import distance_point_point
from compas.geometry import dot_vectors
from compas.geometry import length_vector
from compas.geometry import subtract_vectors

from compas.tolerance import TOL


__all__ = ["dual_mesh", "redistribute_blocks"]


def dual_mesh(mesh, redistribute=True):
    """The dual of a mesh, with its boundary closed by a layer of quads.

    Parameters
    ----------
    mesh : :class:`compas.datastructures.Mesh`

    Returns
    -------
    :class:`compas.datastructures.Mesh`
        The dual, of the same class as ``mesh``.

    """
    dual = mesh_conway_dual(mesh)

    # renumbers the keys 0..n-1, which is what from_vertices_and_faces wants
    vertices, faces = dual.to_vertices_and_faces()

    index = {TOL.geometric_key(xyz): i for i, xyz in enumerate(vertices)}

    #Extent to the boundary
    def add(xyz):
        gkey = TOL.geometric_key(xyz)
        if gkey not in index:
            index[gkey] = len(vertices)
            vertices.append([float(c) for c in xyz])
        return index[gkey]

    for vkey in mesh.vertices():
        if not mesh.is_vertex_on_boundary(vkey):
            continue
        nbrs = mesh.vertex_neighbors(vkey, ordered=True)
        if not nbrs:
            continue
        ring = [add(mesh.edge_midpoint((vkey, nbrs[0])))]
        ring += [add(mesh.face_centroid(fkey))
                 for fkey in mesh.vertex_faces(vkey, ordered=True)]
        ring.append(add(mesh.edge_midpoint((vkey, nbrs[-1]))))
        if len(nbrs) == 2:
            ring.append(add(mesh.vertex_coordinates(vkey)))
        # mesh_conway_dual winds its faces with reversed(vertex_faces(ordered)),
        # and the ring above follows vertex_faces. A face added the other way
        # round does not raise, it overwrites a halfedge -- so reverse it here
        # rather than leaving unify_cycles to repair a broken map.
        faces.append(ring[::-1])

    out = type(dual).from_vertices_and_faces(vertices, faces)
    
    # There is already a way to extent to the boundary in compas, but I have noticed that it gives weird corners
    # out = mesh.dual(include_boundary=True)
    # _vertices, faces = out.to_vertices_and_faces()

    # Both halves wind the same way, so the result is already oriented -- and it
    # comes out with the primal's own orientation, normals dotting +1.00 with it.
    # unify_cycles is only the repair for a primal that arrived mis-wound, and it
    # is QUADRATIC: on 2880 primal faces it took 39.8 s, against 0.06 s for the
    # conway dual it repairs. So pay for it only when the check says to.
    if not _is_oriented(faces):
        out.unify_cycles()
    # Redistribute
    if redistribute:
        redistribute_blocks(out, mesh)

    return out


def _is_oriented(faces):
    """True if no directed halfedge is claimed twice -- the orientability test.

    Two faces claiming ``(u, v)`` the same way round are wound against each
    other. ``Mesh.add_face`` does not raise on that, it overwrites the halfedge,
    so the damage is silent and this is what catches it.
    """
    seen = set()
    for face in faces:
        for u, v in zip(face, face[1:] + face[:1]):
            if (u, v) in seen:
                return False
            seen.add((u, v))
    return True


# ----------------------------------------------------------------------------
# evening out the block sizes
# ----------------------------------------------------------------------------

def redistribute_blocks(dual, primal, project=True):
    """Even out the block sizes of a dual, in place.

    The raw dual puts its boundary on the primal's edge MIDPOINTS, so the outer
    ring of blocks is half a block deep and a corner block -- half-size in both
    directions -- is a quarter. Nothing else is wrong: the interior blocks
    already span centroid to centroid, a full primal cell. So the fix is to
    slide every grid line, not to inflate the edge ring.

    Every line to slide lies on a CHAIN with two fixed ends, and on a chain the
    problem is one-dimensional and has a closed form. Two families of chain
    cover the whole mesh:

    * a primal **strip** -- the run of quads reached by stepping across opposite
      edges -- whose chain is ``[boundary midpoint, centroid, ..., boundary
      midpoint]``. Strips are well defined through singularities, so this needs
      no special case for them;
    * a primal **boundary arc** between two domain corners, whose chain is
      ``[corner, edge midpoint, ..., corner]``.

    On both, the chain's segments ARE the block widths along it, and the two end
    segments are the half blocks. The rule is therefore the whole algorithm:

        double the two end segments, then rescale all of them to the chain length

    On a uniform strip of n cells that turns widths ``0.5, 1, ..., 1, 0.5`` into
    ``n+1`` equal ones, which is the intended result; on a graded strip it keeps
    the grading, because every interior width is carried through untouched. A
    face lies on two strips, so it takes the sum of the two displacements --
    exact on a grid, where they are orthogonal.

    Chains with no end are skipped, and correctly: a strip that closes into a
    loop and a boundary ring with no corner -- a smooth hole -- have no half
    block to fix. Their DEPTH is fixed by the strips that cross them.

    Measured, block area max/min:

    =============================  =====  =============  =============
    input                          raw    redistributed  primal itself
    =============================  =====  =============  =============
    5x4 planar grid                4.00   **1.00**       1.00
    folded grid, floor + one wall  4.00   **1.35**       1.00
    annulus 4x16, curved           3.89   **2.04**       1.95
    =============================  =====  =============  =============

    The annulus does not reach 1.00 and should not: its own primal blocks
    measure 1.95, because a ring at r=5 is genuinely wider than one at r=2. What
    is left after redistribution is the shape, not the artefact.

    Note what this is NOT. Smoothing the dual towards uniform blocks looks like
    the obvious alternative and is wrong on any curved domain: the harmonic
    equilibrium of an annulus piles its rings against the inner boundary, giving
    max/min **34.50** -- far worse than the 3.89 it started from. That is
    Laplacian shrinkage and no amount of tuning removes it.

    Parameters
    ----------
    dual : :class:`compas.datastructures.Mesh`
        The dual to modify, as returned by :func:`dual_mesh`. Modified in place.
    primal : :class:`compas.datastructures.Mesh`
        The mesh it was built from. Supplies the outline, the strips and the
        surface.
    project : bool, optional
        Pull the moved vertices back onto the primal surface afterwards. A chain
        segment that spans a crease is a chord, so a vertex placed along it
        would otherwise float off a folded slab. Skipped on a planar primal,
        where it can do nothing.

    Returns
    -------
    :class:`compas.datastructures.Mesh`
        ``dual``, modified in place.

    """
    gkey_vertex = {TOL.geometric_key(dual.vertex_coordinates(v)): v
                   for v in dual.vertices()}

    def at(xyz):
        return gkey_vertex.get(TOL.geometric_key(xyz))

    # A boundary chain rides the primal OUTLINE, not its own chords: on a wall
    # that bends, midpoint -> midpoint cuts the corner and the dual boundary
    # would leave the outline it is supposed to trace (measured 3.7e-2 off a
    # folded slab). A strip chain has no such guide and uses its own points --
    # what it chords is the crease, which the projection below pulls back.
    chains = []
    for chain, polyline in _boundary_chains(primal, gkey_vertex):
        lengths = [distance_point_point(a, b) for a, b in zip(polyline, polyline[1:])]
        stops, running = [0.0], 0.0
        for length in lengths:
            stops.append(running + length / 2.0)  # where the edge midpoint sits
            running += length
        stops.append(running)
        chains.append((chain, polyline, lengths, stops))

    for faces, ends in _strips(primal):
        chain = ([at(ends[0])]
                 + [at(primal.face_centroid(f)) for f in faces]
                 + [at(ends[1])])
        if any(vkey is None for vkey in chain):
            continue
        points = [dual.vertex_coordinates(v) for v in chain]
        lengths = [distance_point_point(a, b) for a, b in zip(points, points[1:])]
        stops, running = [0.0], 0.0
        for length in lengths:
            running += length
            stops.append(running)
        chains.append((chain, points, lengths, stops))

    # every chain is measured against the ORIGINAL dual and the displacements
    # are summed, never applied in sequence. A boundary midpoint slides ALONG
    # the wall and the strip it caps runs INTO the wall, so the two moves are
    # orthogonal and independent -- but measure the strip after the boundary has
    # already moved and its end segment reads as the diagonal of the two, which
    # skews the whole chain (measured: block spread 1.00 -> 1.07 on a grid).
    moves = {}
    for chain, polyline, lengths, stops in chains:
        if stops[-1] <= 0.0:
            continue
        for vkey, stop in zip(chain[1:-1], _respaced(stops)):
            was = dual.vertex_coordinates(vkey)
            xyz = _point_at(polyline, lengths, stop)
            delta = moves.setdefault(vkey, [0.0, 0.0, 0.0])
            for k in range(3):
                delta[k] += xyz[k] - was[k]

    home = {}
    for fkey in primal.faces():
        vkey = at(primal.face_centroid(fkey))
        if vkey is not None:
            home[vkey] = fkey

    for vkey, delta in moves.items():
        xyz = dual.vertex_coordinates(vkey)
        dual.vertex_attributes(vkey, "xyz", [xyz[k] + delta[k] for k in range(3)])

    # Only a vertex near a CREASE can be off the surface. Everywhere else the
    # mesh is locally flat over the whole half-cell a vertex can travel, the
    # chain segments lie in that plane, and projecting is a no-op that costs
    # nine face tests a vertex. On a folded 8 000-face slab only 160 faces are
    # creased, and skipping the rest took the run from 11.9 s to 1.3 s.
    if project and not _is_planar(primal):
        creased = _creased_faces(primal)
        around = {}
        for vkey in moves:
            fkey = home.get(vkey)
            if fkey is None or fkey not in creased:
                continue
            if fkey not in around:
                around[fkey] = list(set(
                    [fkey] + [g for v in primal.face_vertices(fkey)
                              for g in primal.vertex_faces(v)]))
            xyz, _ = _project(dual.vertex_coordinates(vkey), primal, around[fkey])
            dual.vertex_attributes(vkey, "xyz", xyz)

    return dual


def _respaced(stops):
    """New positions along a chain for its interior points, evening the gaps.

    ``stops`` are the arc positions of the chain's points, so the gaps between
    them are the widths of the blocks along the chain -- and the first and last
    gap are the half blocks the dual is born with. Doubling those two and
    rescaling all of them back to the chain length is the whole redistribution:
    every interior width is carried through untouched, which is what keeps a
    graded strip graded.
    """
    widths = [b - a for a, b in zip(stops, stops[1:])]
    span = stops[-1] - stops[0]
    widths[0] *= 2.0
    widths[-1] *= 2.0
    scale = span / sum(widths)

    out, running = [], stops[0]
    for width in widths[:-1]:
        running += width * scale
        out.append(running)
    return out


def _point_at(polyline, lengths, distance):
    """The point at ``distance`` along a polyline of known segment lengths."""
    for i, length in enumerate(lengths):
        if distance <= length or i == len(lengths) - 1:
            t = distance / length if length else 0.0
            a, b = polyline[i], polyline[i + 1]
            return [a[k] + t * (b[k] - a[k]) for k in range(3)]
        distance -= length
    return list(polyline[-1])


def _boundary_chains(primal, gkey_vertex):
    """Each primal boundary arc, as ``(dual vertices, primal outline points)``.

    The dual vertices run corner, edge midpoint, ..., edge midpoint, corner. The
    outline is the primal boundary polyline over the same arc, which is what the
    dual vertices get placed along -- their own chords would cut every bend.
    """
    chains = []
    for ring in primal.vertices_on_boundaries():
        if ring[0] == ring[-1]:
            ring = ring[:-1]
        n = len(ring)
        if n < 3:
            continue

        # dual_mesh inserts the primal vertex itself only at a valence-2 domain
        # corner, so a primal boundary vertex that IS a dual vertex is a corner
        corner = [gkey_vertex.get(TOL.geometric_key(primal.vertex_coordinates(v)))
                  for v in ring]
        marks = [i for i, vkey in enumerate(corner) if vkey is not None]
        if not marks:
            continue  # a smooth ring has no half block along itself

        mids = [gkey_vertex.get(
            TOL.geometric_key(primal.edge_midpoint((ring[i], ring[(i + 1) % n]))))
            for i in range(n)]

        for a, b in zip(marks, marks[1:] + marks[:1]):
            k = (b - a) % n
            chain = [corner[a]] + [mids[(a + j) % n] for j in range(k)] + [corner[b]]
            if len(chain) > 2 and all(vkey is not None for vkey in chain):
                outline = [primal.vertex_coordinates(ring[(a + j) % n])
                           for j in range(k + 1)]
                chains.append((chain, outline))
    return chains


def _strips(primal):
    """Open primal strips, as ``(faces, (start midpoint, end midpoint))``.

    A strip is the run of quads reached by stepping across opposite edges. Only
    the ones that END on the boundary are returned -- a closed strip has no half
    block to fix, and no fixed point to redistribute against.
    """
    strips, seen = [], set()
    # hoisted: Mesh.number_of_faces() is len(list(self.faces())), so leaving it
    # in the walk below makes the whole pass quadratic -- 14 s of a 15 s run on
    # an 8 000-face mesh, in 128 million generator steps
    cap = primal.number_of_faces()
    for ring in primal.edges_on_boundaries():
        for edge in ring:
            if frozenset(edge) in seen:
                continue
            faces, walk, entry, exit_edge = [], _face_of(primal, edge), edge, None
            while walk is not None:
                faces.append(walk)
                exit_edge = _opposite_edge(primal, walk, entry)
                if exit_edge is None:
                    break
                nxt = None
                for candidate in primal.edge_faces(exit_edge):
                    if candidate is not None and candidate != walk:
                        nxt = candidate
                        break
                entry, walk = exit_edge, nxt
                if len(faces) > cap:
                    exit_edge = None  # a strip that will not close: give up
                    break
            if not faces or exit_edge is None:
                continue
            seen.add(frozenset(edge))
            seen.add(frozenset(exit_edge))
            strips.append((faces, (primal.edge_midpoint(edge),
                                   primal.edge_midpoint(exit_edge))))
    return strips


def _face_of(mesh, edge):
    for fkey in mesh.edge_faces(edge):
        if fkey is not None:
            return fkey
    return None


def _opposite_edge(mesh, fkey, edge):
    """The edge of a quad facing ``edge``; None if the face is not a quad."""
    vertices = mesh.face_vertices(fkey)
    if len(vertices) != 4:
        return None
    pair = set(edge)
    for i in range(4):
        if {vertices[i], vertices[(i + 1) % 4]} == pair:
            return (vertices[(i + 2) % 4], vertices[(i + 3) % 4])
    return None


def _creased_faces(mesh, tol=1e-9):
    """Faces with a non-coplanar face in their 1-ring -- where a fold shows."""
    normals = {fkey: mesh.face_normal(fkey) for fkey in mesh.faces()}
    creased = set()
    for fkey in mesh.faces():
        ring = set(g for vkey in mesh.face_vertices(fkey)
                   for g in mesh.vertex_faces(vkey))
        for g in ring:
            if dot_vectors(normals[fkey], normals[g]) < 1.0 - tol:
                creased.add(fkey)
                break
    return creased


def _is_planar(mesh, tol=None):
    points = [mesh.vertex_coordinates(v) for v in mesh.vertices()]
    if len(points) < 4:
        return True
    if tol is None:
        lengths = [mesh.edge_length(edge) for edge in mesh.edges()]
        tol = 1e-6 * sum(lengths) / float(len(lengths))
    plane = bestfit_plane(points)
    return all(distance_point_plane(point, plane) <= tol for point in points)


def _project(point, mesh, candidates):
    """Closest point on a local patch of faces; also returns the winning face."""
    best, best_distance, best_face = point, None, candidates[0]
    for fkey in candidates:
        polygon = mesh.face_coordinates(fkey)
        centroid = mesh.face_centroid(fkey)
        for i in range(len(polygon)):
            near = _closest_on_triangle(
                point, polygon[i], polygon[(i + 1) % len(polygon)], centroid)
            distance = distance_point_point(point, near)
            if best_distance is None or distance < best_distance:
                best, best_distance, best_face = near, distance, fkey
    return best, best_face


def _closest_on_triangle(point, a, b, c):
    """Closest point of triangle abc to point, in its interior or on an edge."""
    ab = subtract_vectors(b, a)
    ac = subtract_vectors(c, a)
    normal = cross_vectors(ab, ac)
    norm = length_vector(normal)
    if norm > 1e-12:
        normal = [n / norm for n in normal]
        ap = subtract_vectors(point, a)
        inplane = [point[k] - dot_vectors(ap, normal) * normal[k] for k in range(3)]
        u, v = _barycentric(inplane, a, ab, ac)
        if u >= 0.0 and v >= 0.0 and u + v <= 1.0:
            return inplane
    best, best_distance = None, None
    for start, end in ((a, b), (b, c), (c, a)):
        near = closest_point_on_segment(point, (start, end))
        distance = distance_point_point(point, near)
        if best_distance is None or distance < best_distance:
            best, best_distance = near, distance
    return best


def _barycentric(point, a, ab, ac):
    d00 = dot_vectors(ab, ab)
    d01 = dot_vectors(ab, ac)
    d11 = dot_vectors(ac, ac)
    ap = subtract_vectors(point, a)
    d20 = dot_vectors(ap, ab)
    d21 = dot_vectors(ap, ac)
    denom = d00 * d11 - d01 * d01
    if abs(denom) < 1e-18:
        return -1.0, -1.0
    return (d11 * d20 - d01 * d21) / denom, (d00 * d21 - d01 * d20) / denom


