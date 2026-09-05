from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

from math import floor
from operator import itemgetter

from compas.geometry import centroid_points
from compas.geometry import Polyline
from compas.itertools import pairwise

from compas_singular.utilities import list_split

from ..mesh import Mesh


__all__ = ['QuadMesh']


class QuadMesh(Mesh):

    def __init__(self, *args, **kwargs):
        super(QuadMesh, self).__init__(*args, **kwargs)
        self.attributes['strips'] = {}
        self.attributes['polyedges'] = {}

    def strips(self, data=False):
        if not bool(self.attributes['strips']):
            self.collect_strips()
        else:
            pass
            # print("Using earlier collected strips. Pay attention that the mesh has not changed since.")
        for skey in self.attributes['strips']:
            if data:
                yield skey, self.attributes['strips'][skey]
            else:
                yield skey

    def polyedges(self, data=False):
        if not bool(self.attributes['polyedges']):
            self.collect_polyedges()
        else:
            pass
            # print("Using earlier collected polyedges. Pay attention that the mesh has not changed since.")
        for key in self.attributes['polyedges']:
            if data:
                yield key, self.attributes['polyedges'][key]
            else:
                yield key

    # --------------------------------------------------------------------------
    # opposite elements
    # --------------------------------------------------------------------------

    def face_opposite_edge(self, u, v):
        """Returns the opposite edge in the quad face.

        Parameters
        ----------
        u : int
            The identifier of the edge start.
        v : int
            The identifier of the edge end.

        Returns
        -------
        (w, x) : tuple, None
            The opposite edge.
            None if (u, v) is a boundary halfedge, i.e. has no face.

        """

        fkey = self.halfedge[u][v]
        if fkey is None:
            return None
        w = self.face_vertex_descendant(fkey, v)
        x = self.face_vertex_descendant(fkey, w)
        return (w, x)

    def vertex_opposite_vertex(self, u, v, strict=False):
        """Returns the opposite vertex to u accross vertex v.

        Parameters
        ----------
        u : hashable
            A vertex key.
        v : hashable
            A vertex key.
        strict : bool, optional
            Use the stricter crossing rules. Default is False, i.e. the historical
            behaviour.

            The default rules decide on vertices alone, which mis-steers a walk in
            two cases. They turn an interior polyedge onto a boundary when ``u``
            merely *is* a boundary vertex without the edge ``(u, v)`` being a
            boundary edge, and they cannot tell a boundary loop from a boundary.
            With ``strict=True`` the crossing is decided on edges instead: step
            across ``v`` only if ``v`` is a regular interior vertex, or if both
            ``(u, v)`` and the edge walked to are boundary edges.

            This is the engine under :meth:`collect_polyedge` and therefore under
            the whole coarse-layout pipeline, so switching it changes decompositions.

        Returns
        -------
        hashable, None
            The opposite vertex.
            None if v is a singularity or if (u, v) leads outwards.

        """

        if strict:
            nbrs = self.vertex_neighbors(v, ordered=True)
            n = len(nbrs)

            # regular interior vertex: cross straight over
            if n == 4 and not self.is_vertex_on_boundary(v):
                return nbrs[nbrs.index(u) - 2]

            # regular boundary vertex reached along the boundary: follow the boundary
            if n == 3 and self.is_edge_on_boundary(u, v):
                for nbr in nbrs:
                    if nbr != u and self.is_edge_on_boundary(v, nbr):
                        return nbr

            return None

        if self.is_vertex_singular(v):
            return None

        elif self.is_vertex_on_boundary(v):

            if not self.is_vertex_on_boundary(u):
                return None

            else:
                return [nbr for nbr in self.vertex_neighbors(v) if nbr != u and self.is_vertex_on_boundary(nbr)][0]

        else:
            nbrs = self.vertex_neighbors(v, ordered=True)
            return nbrs[nbrs.index(u) - 2]

    # --------------------------------------------------------------------------
    # singularities
    # --------------------------------------------------------------------------

    def is_vertex_singular(self, vkey):
        """Output whether a vertex is quad mesh singularity.

        Parameters
        ----------
        vkey : int
            The vertex key.

        Returns
        -------
        bool
            True if the vertex is a quad mesh singularity. False otherwise.

        """

        if (self.is_vertex_on_boundary(vkey) and self.vertex_degree(vkey) != 3) or (not self.is_vertex_on_boundary(vkey) and self.vertex_degree(vkey) != 4):
            return True

        else:
            return False

    def singularities(self):
        """Returns all the singularity indices in the quad mesh.

        Returns
        -------
        list
            The list of vertex indices that are quad mesh singularities.

        """
        return [vkey for vkey in self.vertices() if self.is_vertex_singular(vkey)]

    def vertex_topo_index(self, vkey):
        """Compute vertex index.

        Parameters
        ----------
        vkey : int
            The vertex key.

        Returns
        -------
        int
            Vertex index.

        """

        if self.vertex_degree(vkey) == 0:
            return 0

        regular_valency = 4 if not self.is_vertex_on_boundary(vkey) else 3

        return (regular_valency - self.vertex_degree(vkey)) / 4

    # --------------------------------------------------------------------------
    # polyedges
    # --------------------------------------------------------------------------

    def collect_polyedge(self, u0, v0, both_sides=True, oriented=False, strict=False):
        """Collect all the edges in the polyedge of the input edge.

        Parameters
        ----------
        u : int
            The identifier of the edge start.
        v : int
            The identifier of the edge end.
        both_sides : bool, optional
            Whether to walk in both directions from the seed halfedge. Default is
            True. With False the walk stops at the first extremity, which is what
            a directional tracer wants.
        oriented : bool, optional
            Whether to return the polyedge in the direction of the seed halfedge.
            Default is False, which returns the reversed polyedge whenever the walk
            reaches an extremity on the first side, i.e. the seed direction is lost.
        strict : bool, optional
            Passed to :meth:`vertex_opposite_vertex`. Default is False.

        Returns
        -------
        polyedge : list
            The list of the vertices in polyedge.
        """

        flipped = False
        polyedge = [u0, v0]

        while len(polyedge) <= self.number_of_vertices():

            # end if closed loop
            if polyedge[0] == polyedge[-1]:
                break

            # get next vertex accros four-valent vertex
            w = self.vertex_opposite_vertex(*polyedge[-2:], strict=strict)

            # flip if end of first extremity
            if w is None:
                if not both_sides:
                    break
                polyedge = list(reversed(polyedge))
                flipped = True
                # stop if end of second extremity
                w = self.vertex_opposite_vertex(*polyedge[-2:], strict=strict)
                if w is None:
                    break

            # add next vertex
            polyedge.append(w)

        if oriented and flipped:
            polyedge = list(reversed(polyedge))

        return polyedge

    def collect_polyedges(self, strict=False):
        """Collect the polyedges accross four-valent vertices between boundaries and/or singularities and store it in the mesh data attributes.

        Parameters
        ----------
        strict : bool, optional
            Passed to :meth:`vertex_opposite_vertex`. Default is False.

        Returns
        -------
        polyedges : list
            List of quad polyedges as list of vertices.

        """

        # the list fixes the seed order (deterministic, so polyedge keys are stable
        # across runs); the set is only for O(1) membership and removal.
        edges = list(self.edges())
        remaining = set(edges)

        nb_polyedges = -1
        for u0, v0 in reversed(edges):

            if (u0, v0) not in remaining:
                # already consumed by an earlier polyedge
                continue

            nb_polyedges += 1

            # collect new polyedge
            polyedge = self.collect_polyedge(u0, v0, strict=strict)
            self.attributes['polyedges'].update({nb_polyedges: polyedge})

            # remove collected edges
            for u, v in pairwise(polyedge):
                remaining.discard((u, v))
                remaining.discard((v, u))

        return self.polyedges(data=True)

    def is_polyedge_closed(self, pkey):
        """Output whether a polyedge is closed.

        Parameters
        ----------
        pkey : hashable
            A strip key.

        Returns
        -------
        bool
            True if the polyedge is closed. False otherwise.
        """

        return self.attributes['polyedges'][pkey][0] == self.attributes['polyedges'][pkey][-1]

    def number_of_polyedges(self):
        """Count the number of polyedges in the mesh."""
        return len(list(self.polyedges()))

    def polyedge_vertices(self, pkey):
        """Return the vertices of a polyedge.

        Parameters
        ----------
        pkey : hashable
            A polyedge key.

        Returns
        -------
        list
            The vertices of the polyedge.
        """

        return self.attributes['polyedges'][pkey]

    def polyedge_edges(self, pkey):
        """Return the edges of a polyedge.

        Parameters
        ----------
        pkey : hashable
            A polyedge key.

        Returns
        -------
        list
            The edges of the polyedge, as pairs of vertex keys.
        """

        return list(pairwise(self.polyedge_vertices(pkey)))

    def polyedge_midpoint(self, pkey):
        """Return the point at mid-length of a polyedge.

        Parameters
        ----------
        pkey : hashable
            A polyedge key.

        Returns
        -------
        Point
            The midpoint.
        """

        return Polyline(self.polyline(pkey)).point_at(0.5)

    def polyedge_length(self, pkey):
        """Return the length of a polyedge.

        Parameters
        ----------
        pkey : hashable
            A polyedge key.

        Returns
        -------
        float
            The sum of the lengths of the edges of the polyedge.
        """

        return sum([self.edge_length(u, v) for u, v in self.polyedge_edges(pkey)])

    def singularity_polyedges(self):
        """Collect the polyedges connected to singularities.

        Returns
        -------
        list
            The polyedges connected to singularities.

        """

        # keep only polyedges connected to singularities or along the boundary
        polyedges = [polyedge for key, polyedge in self.polyedges(data=True) if self.is_vertex_singular(
            polyedge[0]) or self.is_vertex_singular(polyedge[-1]) or self.is_edge_on_boundary(polyedge[0], polyedge[1])]

        # get intersections between polyedges for split
        vertices = [vkey for polyedge in polyedges for vkey in set(polyedge)]
        split_vertices = [vkey for vkey in self.vertices() if vertices.count(vkey) > 1]

        # split singularity polyedges
        return [split_polyedge for polyedge in polyedges for split_polyedge in list_split(polyedge, [polyedge.index(vkey) for vkey in split_vertices if vkey in polyedge])]

    def singularity_polyedge_decomposition(self, strict=False):
        """Returns a quad patch decomposition of the mesh based on the singularity polyedges, including boundaries and additionnal splits on the boundaries.

        Parameters
        ----------
        strict : bool, optional
            Passed to :meth:`vertex_opposite_vertex`. Default is False.

        Returns
        -------
        list
            The polyedges forming the decomposition.

        """
        if self.attributes['polyedges'] == {}:
            self.collect_polyedges(strict=strict)

        polyedges = [polyedge for key, polyedge in self.polyedges(data=True) if (self.is_vertex_singular(
            polyedge[0]) or self.is_vertex_singular(polyedge[-1])) and not self.is_edge_on_boundary(polyedge[0], polyedge[1])]

        # split boundaries
        all_splits = list(set([vkey for polyedge in polyedges for vkey in polyedge] + self.singularities()))

        for boundary in self.boundaries():
            splits = [vkey for vkey in boundary if vkey in all_splits]
            new_splits = []

            if len(splits) == 0:
                new_splits += [vkey for vkey in list(itemgetter(0, int(floor(len(boundary) / 3)), int(floor(len(boundary) * 2 / 3)))(boundary))]

            elif len(splits) == 1:
                i = boundary.index(splits[0])
                new_splits += list(itemgetter(i - int(floor(len(boundary) * 2 / 3)), i - int(floor(len(boundary) / 3)))(boundary))

            elif len(splits) == 2:
                one, two = list_split(boundary + boundary[:1], [boundary.index(vkey) for vkey in splits])
                half = one if len(one) > len(two) else two
                new_splits.append(half[int(floor(len(half) / 2))])

            for vkey in new_splits:
                for nbr in self.vertex_neighbors(vkey):
                    if not self.is_edge_on_boundary(vkey, nbr):
                        new_polyedge = self.collect_polyedge(vkey, nbr, strict=strict)
                        polyedges.append(new_polyedge)
                        all_splits = list(set(all_splits + new_polyedge))
                        break

        # add boundaries
        polyedges += [polyedge for key, polyedge in self.polyedges(data=True) if self.is_edge_on_boundary(polyedge[0], polyedge[1])]

        # get intersections between polyedges for split
        vertices = [vkey for polyedge in polyedges for vkey in set(polyedge)]
        split_vertices = [vkey for vkey in self.vertices() if vertices.count(vkey) > 1]

        # split singularity polyedges
        return [
            split_polyedge for polyedge in polyedges
            for split_polyedge in list_split(polyedge, [polyedge.index(vkey) for vkey in split_vertices if vkey in polyedge])]

    # --------------------------------------------------------------------------
    # polylines
    # --------------------------------------------------------------------------

    def polyedge_graph(self, legacy=True):
        """Compute the vertices and edges of the graph representing the polyedge connectivity,
        where each graph vertex is a mesh polyedge and each graph edge a non-compas_singular mesh vertex representing the crossing of two polyedges.
        Polyedges connected by their extremities, which are singularities, do not count as overlapping.

        Parameters
        ----------
        legacy : bool, optional
            Use the historical implementation. Default is True.

            The legacy graph is quadratic in the number of polyedges and emits, for
            every non-singular vertex, an edge to the *first* polyedge containing
            that vertex -- which is usually the polyedge itself. Roughly half of the
            edges it returns are therefore self-loops (u, u). With ``legacy=False``
            each non-singular vertex contributes one edge between the two polyedges
            that actually cross there, and vertices that are not a crossing are
            skipped.

        Returns
        -------
        tuple
            A tuple of two objects, the dictionary of mesh polyedge indices pointing to their centroid coordinates, and the list of edges between graph vertices.
        """

        if legacy:
            vertices = {key: centroid_points([self.vertex_coordinates(vkey) for vkey in polyedge]) for key, polyedge in self.polyedges(data=True)}
            edges = []
            for key, polyedge in self.polyedges(data=True):
                for vkey in polyedge:
                    if not self.is_vertex_singular(vkey):
                        for key_2, polyedge_2 in self.polyedges(data=True):
                            if vkey in polyedge_2:
                                edges.append((key, key_2))
                                break
            return vertices, edges

        vertices = {pkey: centroid_points(self.polyline(pkey)) for pkey in self.polyedges()}

        vkey_to_pkeys = {vkey: set() for vkey in self.vertices()}
        for pkey, polyedge in self.polyedges(data=True):
            for vkey in polyedge:
                if not self.is_vertex_singular(vkey):
                    vkey_to_pkeys[vkey].add(pkey)

        # a vertex that two polyedges pass through is a crossing; anything else
        # (a singularity, an extremity) is not an adjacency and is skipped
        edges = [tuple(pkeys) for pkeys in vkey_to_pkeys.values() if len(pkeys) == 2]

        return vertices, edges

    # --------------------------------------------------------------------------
    # polylines
    # --------------------------------------------------------------------------

    def polyline(self, pkey):
        """Return the coordinates of the vertices of a polyedge.

        Parameters
        ----------
        pkey : hashable
            A polyedge key.

        Returns
        -------
        list
            The polyline as a list of XYZ points.
        """

        return [self.vertex_coordinates(vkey) for vkey in self.polyedge_vertices(pkey)]

    def polylines(self):
        """Return the polylines of the quad mesh.

        Returns
        -------
        list
            The polylines.
        """

        return [[self.vertex_coordinates(vkey) for vkey in polyedge] for key, polyedge in self.polyedges(data=True)]

    def singularity_polylines(self):
        """Return the polylines connected to singularities.

        Returns
        -------
        list
            The polylines connected to singularities.

        """
        return [[self.vertex_coordinates(vkey) for vkey in polyedge] for polyedge in self.singularity_polyedges()]

    def singularity_polyline_decomposition(self):
        """Return the polylines forming a quad patch decomposition of the mesh.

        Returns
        -------
        list
            The polylines connected to singularities.

        """
        return [[self.vertex_coordinates(vkey) for vkey in polyedge] for polyedge in self.singularity_polyedge_decomposition()]

    # --------------------------------------------------------------------------
    # strips
    # --------------------------------------------------------------------------

    def number_of_strips(self):
        """Count the number of strips in the mesh."""
        return len(list(self.strips()))

    def collect_strip(self, u0, v0, both_sides=True):
        """Returns all the edges in the strip of the input edge.

        Parameters
        ----------
        u : int
            The identifier of the edge start.
        v : int
            The identifier of the edge end.
        both_sides : bool, optional
            Whether to walk in both directions from the seed halfedge. Default is
            True. With False the walk stops at the first extremity.

        Returns
        -------
        strip : list
            The list of the edges in strip.
        """

        if self.halfedge[u0][v0] is None:
            if not both_sides:
                return [(u0, v0)]
            u0, v0 = v0, u0

        edges = [(u0, v0)]

        count = self.number_of_edges()
        while count > 0:
            count -= 1

            u, v = edges[-1]
            opposite = self.face_opposite_edge(u, v)
            if opposite is None:
                break
            w, x = opposite

            if (x, w) == edges[0]:
                break

            edges.append((x, w))

            if w not in self.halfedge[x] or self.halfedge[x][w] is None:
                if not both_sides:
                    break
                edges = [(v, u) for u, v in reversed(edges)]
                u, v = edges[-1]
                if v not in self.halfedge[u] or self.halfedge[u][v] is None:
                    break

        return edges

    def collect_strips(self):
        """Collect the strip data and store it in the mesh data attributes.

        Returns
        -------
        strips : dict
            The strip data.
        """

        # CLEAR FIRST. ``update`` only adds and overwrites, so re-collecting a mesh
        # that has LOST a strip left the old keys behind -- and worse, a key that was
        # deleted and is now re-assigned lands at the END of the dict's insertion
        # order. Measured on a 4x4 grid: delete strip 3, re-collect, and the order is
        # [0, 1, 2, 4, 5, 6, 7, 3] -- so ``list(strips())[-1]`` is 3 while the maximum
        # is 7, and ``grammar_pattern.add_strip``'s ``last + 1`` then names strip 4,
        # which already exists and is silently overwritten. Stale entries are also
        # read as real by ``is_strip_closed``, which only ever looks at
        # ``strips[skey][0]``.
        self.attributes['strips'].clear()

        # see collect_polyedges: list for a stable seed order, set for O(1) removal
        edges = [(u, v) if self.halfedge[u][v] is not None else (v, u) for u, v in self.edges()]
        remaining = set(edges)

        nb_strip = -1
        for u0, v0 in reversed(edges):

            if (u0, v0) not in remaining:
                continue

            nb_strip += 1

            strip_edges = self.collect_strip(u0, v0)
            self.attributes['strips'].update({nb_strip: strip_edges})
            for u, v in strip_edges:
                remaining.discard((u, v))
                remaining.discard((v, u))

        return self.strips(data=True)

    def is_strip_closed(self, skey):
        """Output whether a strip is closed.

        Parameters
        ----------
        skey : hashable
            A strip key.

        Returns
        -------
        bool
            True if the strip is closed. False otherwise.
        """

        return not self.is_edge_on_boundary(*self.strip_edges(skey)[0])

    def strip_edges(self, skey):
        """Return the edges of a strip.

        Parameters
        ----------
        skey : hashable
            A strip key.
        Returns
        -------
        list
            The edges of the strip.

        """

        return self.attributes['strips'][skey]

    def edge_strip(self, edge):
        """Return the strip of an edge.

        Parameters
        ----------
        edge : tuple
            An edge as two vertex keys.

        Returns
        -------
        strip
            The strip of the edge.
        """

        for skey, edges in self.strips(data=True):
            if edge in edges or tuple(reversed(edge)) in edges:
                return skey

    def strip_faces(self, skey):
        """Return the faces of a strip.

        Parameters
        ----------
        skey : hashable
            A strip key.

        Returns
        -------
        list
            The faces of the strip.

        """

        return [self.halfedge[u][v] for u, v in self.strip_edges(skey) if self.halfedge[u][v] is not None]

    def face_strips(self, fkey):
        """Return the two strips of a face.

        Parameters
        ----------
        fkey : hashable

        Returns
        -------
        list
            The two strips of the face.
        """

        return [self.edge_strip((u, v)) for u, v in list(self.face_halfedges(fkey))[:2]]

    # --------------------------------------------------------------------------
    # strip data operations
    # --------------------------------------------------------------------------

    def substitute_vertex_in_strips(self, old_vkey, new_vkey, strips=None):
        """Substitute a vertex by another one.

        Parameters
        ----------
        old_vkey : hashable
            The old vertex key.
        new_vkey : hashable
            The new vertex key.
        strips : list
            List of specific strip keys. Per default None, i.e. all.

        """

        if strips is None:
            strips = list(self.strips())
        self.attributes['strips'].update({skey: [tuple([new_vkey if vkey == old_vkey else vkey for vkey in list(edge)])
                                                         for edge in self.strip_edges(skey)] for skey in strips})

    def delete_face_in_strips(self, fkey):
        """Delete face in strips.

        Parameters
        ----------
        old_vkey : hashable
            The old vertex key.
        new_vkey : hashable
            The new vertex key.

        """

        self.attributes['strips'] = {skey: [(u, v) for u, v in self.strip_edges(skey) if self.halfedge[u][v] != fkey] for skey in self.strips()}

    # --------------------------------------------------------------------------
    # strip graph
    # --------------------------------------------------------------------------

    def strip_graph(self):
        """Compute the vertices and edges of the graph representing the strip connectivity,
        where each graph vertex is a mesh strip and each graph edge a mesh face representing the crossing of two strips.
        Potentially includes loop edges (u, u) or multiple arallel edges (u, v) and/or (v, u).

        Returns
        -------
        tuple
            A tuple of two objects, the dictionary of mesh strip keys pointing to graph vertex coordinates,
            and the list of edges between graph vertices.
        """

        vertices = {skey: centroid_points(self.strip_edge_midpoint_polyline(skey) if not self.is_strip_closed(skey)
                                          else self.strip_edge_midpoint_polyline(skey)[:-1]) for skey in self.strips()}
        edges = [tuple(self.face_strips(fkey)) for fkey in self.faces()]
        return vertices, edges

    # --------------------------------------------------------------------------
    # strip polyedges
    # --------------------------------------------------------------------------

    def strip_side_polyedges(self, skey):
        """Return the two side polyedges of a strip.

        Parameters
        ----------
        skey : hashable
            A strip key.

        Returns
        -------
        tuple
            The pair of polyedges on the side of the strip.
        """

        strip_edges = self.strip_edges(skey)

        starts = [edge[0] for edge in strip_edges]
        ends = [edge[1] for edge in strip_edges]

        if self.is_strip_closed(skey):
            starts += starts[:1]
            ends += ends[:1]

        return (starts, ends)

    # --------------------------------------------------------------------------
    # strip polylines
    # --------------------------------------------------------------------------

    def strip_edge_midpoint_polyline(self, skey):
        """Return the strip polyline connecting edge midpoints.

        Parameters
        ----------
        skey : hashable
            A strip key.

        Returns
        -------
        list
            The edge midpoint polyline.
        """

        polyline = [self.edge_midpoint(u, v) for u, v in self.strip_edges(skey)]

        if self.is_strip_closed(skey):
            return polyline + polyline[: 1]

        else:
            return polyline

    def strip_face_centroid_polyline(self, skey):
        """Return the strip polyline connecting face centroids.

        Parameters
        ----------
        skey : hashable
            A strip key.

        Returns
        -------
        list
            The face centroid polyline.
        """

        polyline = [self.face_centroid(fkey) for fkey in self.strip_faces(skey)]

        if self.is_strip_closed(skey):
            return polyline + polyline[: 1]

        else:
            return polyline

    def strip_side_polylines(self, skey):
        """Return the two side polylines of a strip.

        Parameters
        ----------
        skey : hashable
            A strip key.

        Returns
        -------
        tuple
            The pair of polylines on the side of the strip.
        """

        starts, ends = self.strip_side_polyedges(skey)
        return ([self.vertex_coordinates(vkey) for vkey in starts], [self.vertex_coordinates(vkey) for vkey in ends])


# ==============================================================================
# Main
# ==============================================================================

if __name__ == '__main__':
    pass

    # import compas
    # from compas_plotters.meshplotter import MeshPlotter

    # # mesh = QuadMesh.from_obj(compas.get('faces.obj'))
    # # mesh = QuadMesh.from_json('/Users/Robin/Desktop/json/debug.json')

    # # mesh.collect_strips()
    # # mesh.collect_polyedges()

    # # print(mesh.singularities())
    # # print(len(list(mesh.strips())))
    # # print(len(list(mesh.polyedges())))

    # # print(len(mesh.singularity_polyedge_decomposition()))

    # # print(mesh.strip_graph())
    # # print(mesh.polyedge_graph())

    # #plotter = MeshPlotter(mesh, figsize=(20, 20))
    # #plotter.draw_vertices(radius=0.4, text='key')
    # # plotter.draw_edges()
    # # plotter.draw_faces()
    # # plotter.show()
