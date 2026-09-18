from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

import os

from compas.datastructures import Mesh
from compas.geometry import centroid_points
from compas.geometry import angle_points
from compas.geometry import Point


__all__ = ['Mesh']


class Mesh(Mesh):

    def __init__(self, *args, **kwargs):
        super(Mesh, self).__init__(*args, **kwargs)

    # ------------------------------------------------------------------------
    # COMPAS 2.x compatibility
    #
    # In COMPAS 2.x the edge methods take a single ``(u, v)`` tuple, whereas
    # compas_singular was written against the COMPAS 0.x/1.x API where they
    # took ``u`` and ``v`` as two separate arguments. These thin overrides
    # accept both calling styles (the tuple form is what COMPAS 2.x uses
    # internally) so the rest of compas_singular can keep the old signatures.
    # ------------------------------------------------------------------------

    @staticmethod
    def _as_edge(u, v):
        return u if v is None else (u, v)

    def edge_midpoint(self, u, v=None):
        return super(Mesh, self).edge_midpoint(self._as_edge(u, v))

    def edge_length(self, u, v=None):
        return super(Mesh, self).edge_length(self._as_edge(u, v))

    def edge_faces(self, u, v=None):
        return super(Mesh, self).edge_faces(self._as_edge(u, v))

    def is_edge_on_boundary(self, u, v=None):
        return super(Mesh, self).is_edge_on_boundary(self._as_edge(u, v))

    def edge_point(self, u, v=None, t=0.5):
        # old style: edge_point(u, v, t); new style: edge_point((u, v), t)
        if isinstance(u, (list, tuple)):
            return super(Mesh, self).edge_point(u, 0.5 if v is None else v)
        return super(Mesh, self).edge_point((u, v), t)

    # ------------------------------------------------------------------------
    # JSON -- the state a CAD round trip cannot carry
    # ------------------------------------------------------------------------
    #
    # A mesh baked into a document is vertices and faces. Everything a
    # ``CoarseQuadMesh`` KNOWS lives in ``attributes`` -- ``strips``,
    # ``strips_density``, ``polyedges``, ``face_pole``, the coarse-to-dense maps
    # -- and a bake drops all of it. compas serialises that dict, so these two
    # are the one mechanism that moves a mesh's knowledge and not just its
    # shape. Measured on a 9-strip plate with a point feature: 9 strips,
    # per-strip densities including a hand-set 9, and 3 ``face_pole`` entries,
    # in 2.3 kB.

    def save_to_json(self, filepath, pretty=False):
        """Write the mesh and everything in :attr:`attributes` to ``filepath``.

        ``to_json`` with two differences worth having: the parent directory is
        created if it does not exist (``to_json`` raises ``FileNotFoundError``),
        and the path is returned so a caller can report where it went.

        Extra state goes in :attr:`attributes` rather than a second file -- it
        is a plain dict and it round-trips::

            coarse.attributes['route'] = 'skeleton'
            coarse.attributes['solve'] = {'target_length': 0.5, 'mode': 'tangent'}
            coarse.save_to_json(path)

        Parameters
        ----------
        filepath : str
        pretty : bool, optional
            Indent the output. Roughly doubles the size; useful while debugging.

        Returns
        -------
        str
            ``filepath``.
        """
        folder = os.path.dirname(os.path.abspath(filepath))
        if folder and not os.path.isdir(folder):
            os.makedirs(folder)
        self.to_json(filepath, pretty=pretty)
        return filepath

    @classmethod
    def load_from_json(cls, filepath, default=None):
        """**Construct a mesh from what** :meth:`save_to_json` **wrote.**

        The class comes from the FILE, not from ``cls``: compas stores a
        ``dtype`` and dispatches on it, so a ``CoarsePseudoQuadMesh`` comes back
        as one however it is asked for -- ``Mesh.load_from_json`` included.
        ``cls`` still acts as an assertion, because ``Data.from_json`` refuses
        anything that is not an instance of it ("The data in the file is not a
        ..."), so asking through the class you expect is worth doing.

        **Strip keys are repaired here.** JSON object keys are strings, so
        ``attributes['strips']``, ``['strips_density']`` and ``['polyedges']``
        come back keyed ``'0'``, ``'1'``, ... The mesh still densifies --
        everything internal reads whatever keys it finds, verified identical
        face counts -- but ``get_strip_density(0)`` raises ``KeyError`` and
        ``sorted(strips())`` puts ``'10'`` before ``'2'``. Both are silent in
        the ways that matter, so digit-string keys become ints again on the way
        in and the round trip is invisible.

        Parameters
        ----------
        filepath : str
        default : optional
            Returned when the file does not exist, instead of raising. For a
            caller resuming from a cache that may not have been written yet.

        Returns
        -------
        Mesh
            Of whatever class the file names.

        Raises
        ------
        TypeError
            From ``Data.from_json``, if the file holds a mesh that is not an
            instance of ``cls``.
        """
        if not os.path.isfile(filepath):
            return default

        mesh = cls.from_json(filepath)

        for key in ('strips', 'strips_density', 'polyedges', 'edges_to_curves', 'dense_pattern'):
            table = mesh.attributes.get(key)
            if not isinstance(table, dict):
                continue
            mesh.attributes[key] = {
                int(k) if isinstance(k, str) and k.lstrip('-').isdigit() else k: v
                for k, v in table.items()}
        return mesh

    def to_vertices_and_faces(self, keep_keys=True):

        if keep_keys:
            vertices = {vkey: self.vertex_coordinates(vkey) for vkey in self.vertices()}
            faces = {fkey: self.face_vertices(fkey) for fkey in self.faces()}
        else:
            vertex_index = self.vertex_index()
            vertices = [self.vertex_coordinates(key) for key in self.vertices()]
            faces = [[vertex_index[key] for key in self.face_vertices(fkey)] for fkey in self.faces()]
        return vertices, faces

    def boundaries(self):
        """Collect the mesh boundaries as lists of vertices.

        Parameters
        ----------
        mesh : Mesh
            Mesh.

        Returns
        -------
        boundaries : list
            List of boundaries as lists of vertex keys.

        """

        boundary_edges = {}
        for u, v in self.edges():
            if self.halfedge[u][v] is None:
                boundary_edges[u] = v
            elif self.halfedge[v][u] is None:
                boundary_edges[v] = u

        boundaries = []
        boundary = list(boundary_edges.popitem())
        while len(boundary_edges) > 0:
            w = boundary_edges.pop(boundary[-1])
            if w == boundary[0]:
                boundaries.append(boundary)
                if len(boundary_edges) > 0:
                    boundary = list(boundary_edges.popitem())
            else:
                boundary.append(w)

        return boundaries

    def is_boundary_vertex_kink(self, vkey, threshold_angle):
        """Return whether there is a kink at a boundary vertex according to a threshold angle.

        Parameters
        ----------
        vkey : Key
            The boundary vertex key.
        threshold_angle : float
            Threshold angle in rad.

        Returns
        -------
        bool
            True if vertex is on the boundary and has an angle larger than the threshold angle. False otherwise.
        """

        # check if vertex is on boundary
        if not self.is_vertex_on_boundary(vkey):
            return False

        # get the two adjacent boundary vertices (exactly two for manifold meshes)
        ukey, wkey = [nbr for nbr in self.vertex_neighbors(vkey) if self.is_edge_on_boundary(vkey, nbr)]

        # compare boundary angle with threshold angle
        return angle_points(self.vertex_coordinates(ukey), self.vertex_coordinates(vkey), self.vertex_coordinates(wkey)) > threshold_angle

    def boundary_kinks(self, threshold_angle):
        """Return the boundary vertices with kinks.

        Parameters
        ----------
        threshold_angle : float
            Threshold angle in rad.

        Returns
        -------
        list
            The list of the boundary vertices at kink angles higher than the threshold value.

        """

        return [vkey for bdry in self.vertices_on_boundaries() for vkey in bdry if self.is_boundary_vertex_kink(vkey, threshold_angle)]

    def vertex_centroid(self):
        """Calculate the centroid of the mesh vertices.

        Parameters
        ----------

        Returns
        -------
        list
            The coordinates of the centroid of the mesh vertices.
        """

        return centroid_points([self.vertex_coordinates(vkey) for vkey in self.vertices()])

    def vertex_map(self, view=False):
        vkeys, _ = self.to_vertices_and_faces()

        vertices = []
        for vkey in vkeys:
            vertex = Point(*self.vertex_coordinates(vkey))
            vertices.append((vkey, vertex))

        if view:
            from compas_viewer.viewer import Viewer
            from compas_viewer.scene.tagobject import Tag

            viewer = Viewer()
            group = viewer.scene.add_group("Vertex map")

            for vertex in vertices:
                vkey = vertex[0]
                vertex = vertex[1]

                group.add(vertex, name="Vertex: " + str(vkey))
                group.add(Tag(text=str(vkey), position=vertex))

            viewer.show()
            viewer.scene.clear()

        return vertices

# ==============================================================================
# Main
# ==============================================================================

if __name__ == '__main__':
    pass
