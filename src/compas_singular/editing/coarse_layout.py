"""**Hand-editing a coarse quad layout: the mesh, with no Rhino in it.**

:class:`CoarseLayoutEditor` holds a coarse layout being edited and offers the
four things that can be done to one -- move a corner, cut it with a drawn
curve, delete a strip, and commit the result back to the decomposition it came
from. Nothing here picks, prompts, draws or prints. A refusal is recorded in
:attr:`~CoarseLayoutEditor.last_reason` and returned as ``False``; whoever is
driving decides how to say so.

**The layout is edited as a MESH, not as the polylines baked for it.** Both are
legal inputs to :meth:`FieldDecomposition.edit_coarse` -- ``edit.faces_from_geometry``
accepts either -- but only one of them is safe to move a corner in:

* **polylines duplicate every shared corner.** ``edit.face_polylines`` writes
  one closed polyline per patch, so an interior corner where four patches meet
  is four coincident points. Moving one of them and leaving the other three is
  a tear, and the weld in ``edit.mesh_from_faces`` rounds at 3 decimals, so the
  moved copy silently becomes a NEW vertex: the layout comes back with a slit
  in it and nothing reports an error.
* **the mesh shares them.** ``vertex_attributes(vkey, 'xyz', ...)`` moves the
  one corner that all four patches reference, which is what the user meant.

It also keeps the strip-index trap in ``RHINO_PLUGIN.md`` §4.2 out of the edit
ITSELF -- indices permute because a bake/read-back rebuilds the topology from
geometry, and the mesh here never leaves Python, so keys are stable for as long
as corners are being moved. **They are not stable across** :meth:`CoarseLayoutEditor.commit`:
``edit_coarse`` welds and repairs, which is the same rebuild, and it was
measured to leave 2 of 8 vertex keys pointing at the same corner on a disc. So
a density map held across a commit must be keyed geometrically -- rounded edge
midpoints, or one remembered ``(gkey(u), gkey(v))`` -- exactly as §4.2 says.

**The edit is not finished until** :meth:`~CoarseLayoutEditor.commit` **runs.**
Moving vertices in place produces a layout whose edges no longer match any
traced separatrix by geometric key. Densifying THAT directly gives straight
chords and throws away the field alignment the whole front end exists to
produce. ``edit_coarse`` re-matches each edge to the separatrix it came from
and warps that curve onto the moved corners (``edit.warp_polyline``), so the
nudge costs the nudge and not the curvature. It also snaps boundary corners
onto the domain walls, which is load-bearing rather than cosmetic -- see
``framefield/edit.py``.

**A DRAWN CURVE IS A CUT, AND A CUT HAS TO LEAVE QUADS**
--------------------------------------------------------

:meth:`~CoarseLayoutEditor.insert_curve` takes a line, an arc or a polyline as
a point list and makes it part of the layout: every coarse edge it crosses is
split, and every patch it passes through is split in two. ``Mesh.split_edge``
puts the new corner into BOTH faces of the edge, so nothing is ever left
hanging -- but that is also why the cut cannot stop wherever it likes. A curve
that ends halfway along an interior edge leaves the patch on the FAR side of
that edge with five corners, and ``solve_non_quad_faces`` then fans the
pentagon into a quad plus a triangle and keeps the triangle as a pole. Measured
on a two-quad test layout, one such line: ``{5: 2, 3: 1}`` in, three poles out
-- singularities nobody drew.

So the cut is planned first and checked before anything is mutated, and the
rule it is checked against is simply **every patch of the result has four
sides** (pseudo-quad poles that were already there excepted). Two things follow
from that rule rather than being imposed on top of it, and
:attr:`~CoarseLayoutEditor.last_reason` says which one was hit:

* a patch must be entered and left through OPPOSITE sides -- through adjacent
  ones, or through a corner, the split makes a triangle. Nothing can rescue
  this one, and it is refused;
* every edge on the way must be crossed THROUGH, which puts the two ends of the
  curve on the layout's own boundary. A cut from one wall to another is
  therefore what gets accepted, and it is the same condition
  ``is_polyedge_valid_for_strip_addition`` puts on adding a strip.

**The second one does not have to be drawn, only finished.** A curve stopping
in the middle of the layout is not wrong, it is short: the patch behind that
edge has to be split too, and the one behind that, until the run reaches a wall
-- which is to say the cut follows a STRIP, and the continuation is
``face_opposite_edge`` at the same parameter. :meth:`~CoarseLayoutEditor._extend_plan`
does exactly that, and ``insert_curve``'s ``extend`` argument offers it whenever
a cut comes up short, so cutting one patch across the middle works and costs
the patches downstream of it rather than the layout's quads. The part the user
drew keeps its shape; the continuation is straight, because it should look like
what it is.

The mesh is only replaced once the whole cut is known to be legal: the plan is
computed against the untouched layout, applied to a COPY, and the copy is
adopted only if it passes. A refusal costs nothing and leaves the layout
exactly as it was.

**The drawn SHAPE lives in the curve map, not in the mesh.** A coarse edge is
always a straight chord -- curvature reaches the mesh through
``FieldDecomposition.edges_to_curves``, which hands one polyline per coarse
edge to the densifier. So an arc inserted here is stored in
:attr:`CoarseLayoutEditor.curves` and published to ``decomposition.user_curves``
on commit, where ``edges_to_curves``' exact-match branch finds it. The same map
holds the two halves of any curved edge the cut crossed: without them the
halves reach the warp branch, which would drag the WHOLE separatrix onto a
half-length edge and fold it.

**DELETING IS A STRIP, NOT AN EDGE**
------------------------------------

A quad layout cannot lose a single edge: removing one merges two patches into a
hexagon. The only removal that keeps every patch four-sided runs the full width
of the layout, which is what a strip is -- so :meth:`~CoarseLayoutEditor.delete_strip`
takes an edge and deletes the strip through it, exactly as
``CMD_edit_quad_mesh``'s ``remove_line`` does on the dense mesh. It is
``grammar_pattern.delete_strip``, the library's only deletion grammar.

Note what that is and is not. It **collapses a band** of patches and welds the
band's two sides together; it does not *dissolve* the picked line and merge the
patches either side of it. Dissolving is the exact inverse of a cut, and no
such operation exists in ``compas_singular``. Before a commit,
:meth:`~CoarseLayoutEditor.reset` undoes a cut; after one, nothing does.

Deletion is not local, and :meth:`~CoarseLayoutEditor.plan_strip_deletion`
reports the three ways it reaches further than the picked strip -- collateral
strips, welded neighbours and collapsed boundaries -- so a caller can put the
consequences in front of the user before anything happens.

**And some strips simply cannot go, which nothing about the strip predicts.**
Welding a strip's two sides together can put three patches on one edge. Swept
over every strip of four layouts: square+square-hole 8 of 8 fine, L-plate 4 of
4, square+two-circular-holes 23 of 24 (the 24th raised ``KeyError`` from inside
the grammar) -- but **square+one-circular-hole, 10 of 20 non-manifold**. Poles
are not the discriminator, which was the obvious guess and is wrong: in that
layout all 7 strips running into a pole deleted cleanly, and every one of the
10 failures had no pole edge at all. So there is no cheap test, and
``plan_strip_deletion`` does not attempt one -- **it performs the deletion on a
copy and reports what came out**. ``edit_coarse`` refuses a non-manifold layout
anyway, so checking here turns a failure several commands later into a refusal
that costs nothing.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from compas.geometry import closest_point_on_segment_xy
from compas.geometry import distance_point_point
from compas.geometry import intersection_segment_segment_xy
from compas.geometry import is_point_in_polygon_xy

# From grammar_pattern, NOT grammar.delete_strip: the two modules both define
# ``delete_strip`` and ``datastructures/__init__`` star-imports ``grammar``
# last, so the package-level name resolves to the MODULE and calling it raises
# ``TypeError: 'module' object is not callable``. This is also the version that
# merges the two sides, handles collateral deletions and repairs 'face_pole'.
from ..datastructures.mesh_quad.grammar_pattern import collateral_strip_deletions
from ..datastructures.mesh_quad.grammar_pattern import delete_strip as _grammar_delete_strip
from ..datastructures.mesh_quad.grammar_pattern import total_boundary_deletions

# ``coarse_curves`` lives under ``rhino`` for historical reasons but touches
# neither ``rhinoscriptsyntax`` nor ``Rhino`` -- its own docstring says so. This
# is the same wall projection ``CMD_edit_quad_mesh.project_to_wall`` uses, so
# both editors hold a moved boundary vertex on the outline the same way.
from ..rhino.coarse_curves import BoundaryLoop


__all__ = ['CoarseLayoutEditor', 'PRECISION']


#: Decimals the curve map is keyed at. The resolution ``geometric_key`` and
#: ``from_polylines``' endpoint matching already use, so a key made here still
#: matches after a round trip through either.
PRECISION = 3


class CoarseLayoutEditor(object):
    """A coarse quad layout being hand-edited. No Rhino, no prompts, no prints.

    Parameters
    ----------
    decomposition : :class:`FieldDecomposition`
        The one that produced the layout. Kept because :meth:`commit` hands the
        edited mesh back to its :meth:`~FieldDecomposition.edit_coarse`, which
        needs the separatrices and the domain loops that only it has.
    coarse : :class:`CoarsePseudoQuadMesh`, optional
        The layout to edit. Defaults to ``decomposition.mesh``.
    loops : list, optional
        The domain walls, outer first, as point lists or closed polylines. A
        boundary corner that moves is projected onto the nearest one, so that a
        caller's preview shows where the corner will actually end up. Defaults
        to the decomposition's own loops.
    snap_tol : float, optional
        Passed straight through to ``edit_coarse`` on :meth:`commit`.

    Attributes
    ----------
    mesh : :class:`CoarsePseudoQuadMesh`
        The layout as it currently stands. Replaced wholesale by every
        operation that succeeds, so hold the editor, not the mesh.
    curves : dict
        ``(rounded pa, rounded pb) -> [[x, y, z], ...]`` -- the SHAPE of a
        coarse edge that has one: a curve the user drew, or a half of a curved
        edge that a cut crossed. Keyed geometrically and never by vertex key,
        because :meth:`commit` renumbers them.
    last_reason : str
        Why the last operation was refused. Empty when nothing was.

    Examples
    --------
    >>> decomposition = FieldDecomposition.from_boundary(outer, inners)
    >>> decomposition.decomposition_mesh()
    >>> editor = CoarseLayoutEditor(decomposition)
    >>> editor.insert_curve([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]], extend=True)
    True
    >>> ok, notes = editor.commit()
    >>> mesh = decomposition.quad_mesh()
    """

    def __init__(self, decomposition, coarse=None, loops=None, snap_tol=None):
        self.decomposition = decomposition
        self.mesh = coarse if coarse is not None else decomposition.mesh
        if self.mesh is None:
            raise ValueError('no coarse layout -- call decomposition_mesh() first')
        # ``_loops`` is private but it is the only accessor for outer + inners,
        # and ``edit_coarse`` uses the same one internally. Taking it from the
        # decomposition by default keeps the preview projection and the commit
        # snap measuring against the same walls.
        if loops is None:
            loops = decomposition._loops()
        self.loops = list(loops)
        self._walls = [BoundaryLoop(loop) for loop in self.loops
                       if len(list(loop)) >= 3]
        self.snap_tol = snap_tol
        self.curves = {}
        #: The layout as it was before the current round of edits. A whole mesh
        #: rather than a coordinate map: a cut changes the TOPOLOGY, and putting
        #: coordinates back would leave the new patches in place.
        self._snapshot = self.mesh.copy()
        self._snapshot_curves = {}
        self.last_reason = ''
        #: What the last successful cut / deletion did, for a caller to report.
        self.last_cut = {}
        self.last_deletion = {}

    # ------------------------------------------------------------------
    # tolerances -- everything here is relative to the layout
    # ------------------------------------------------------------------

    def mean_edge(self):
        """Average coarse edge length. The scale every tolerance here is in."""
        lengths = [self.mesh.edge_length(edge) for edge in self.mesh.edges()]
        lengths = [length for length in lengths if length > 0.0]
        return sum(lengths) / len(lengths) if lengths else 1.0

    def corner_tol(self):
        """How close to a corner a point is TAKEN AS that corner.

        Not politeness: splitting an edge a hair from its end makes a coarse
        edge a hundredth of its neighbour's length, and that strip densifies
        into slivers however sensible the density is. It also stays clear of
        the 3-decimal weld in ``edit.mesh_from_faces``, which would otherwise
        merge the two corners at commit and silently drop the patch as
        degenerate.
        """
        return 0.05 * self.mean_edge()

    def on_tol(self):
        """How far off an edge a point may be and still count as on it.

        Tight, because a point that reaches here was picked with the cursor
        constrained onto a layout edge -- this absorbs float noise, not aim.
        """
        return 1e-3 * self.mean_edge()

    # ------------------------------------------------------------------
    # queries a front end needs
    # ------------------------------------------------------------------

    def locate(self, point):
        """Where a point sits on the layout: on a corner, on an edge, or not.

        Returns ``('vertex', vkey)``, ``('edge', (u, v), xyz)``, or ``None``.

        Corners win over edges, and by a whole :meth:`corner_tol` -- see there.
        ``None`` means the point is on the layout nowhere, which is a refusal
        rather than something to guess at.
        """
        mesh = self.mesh
        tol = self.corner_tol()
        best = None
        for vkey in mesh.vertices():
            d = distance_point_point(point, mesh.vertex_coordinates(vkey))
            if d <= tol and (best is None or d < best[0]):
                best = (d, vkey)
        if best is not None:
            return ('vertex', best[1])

        tol = self.on_tol()
        best = None
        for u, v in mesh.edges():
            segment = (mesh.vertex_coordinates(u), mesh.vertex_coordinates(v))
            q = closest_point_on_segment_xy(point, segment)
            d = distance_point_point(point, q)
            if d <= tol and (best is None or d < best[0]):
                best = (d, (u, v), [q[0], q[1], 0.0])
        if best is not None:
            return ('edge', best[1], best[2])
        return None

    def project(self, xyz, boundary):
        """A boundary corner goes onto the nearest wall; an interior one does not.

        Eligibility is TOPOLOGICAL, exactly as in ``edit.snap_to_loops``, and
        for the same reason: projecting anything merely near a wall would drag
        an interior corner that legitimately sits in a narrow slot onto it and
        collapse a patch nobody touched.

        Doing it here as well as in ``edit_coarse`` is not redundant -- it is
        what lets a front end preview the position the commit will actually
        produce.
        """
        if not boundary or not self._walls:
            return list(xyz)
        best_d, best_p = float('inf'), None
        for wall in self._walls:
            d, _s, q = wall.project(xyz)
            if d < best_d:
                best_d, best_p = d, q
        return best_p if best_p is not None else list(xyz)

    def is_vertex_on_boundary(self, vkey):
        """Whether a corner is on the layout's boundary, and so wall-held."""
        return vkey in set(self.mesh.vertices_on_boundary())

    def edge_shape(self, u, v):
        """The drawn SHAPE this edge carries, or ``None`` if it is a chord.

        For display: an edge that carries an inserted arc should be drawn as
        that arc, or a cut that worked looks like one that snapped straight.
        """
        return self._own_curve(self.mesh.vertex_coordinates(u),
                               self.mesh.vertex_coordinates(v))

    # ------------------------------------------------------------------
    # move a corner
    # ------------------------------------------------------------------

    def move_vertex(self, vkey, xyz, project=True):
        """Put one corner at ``xyz``. ``True`` if it moved.

        A corner on the layout's boundary is projected onto its domain wall
        unless ``project`` is off -- so a caller that already previewed the
        projected point may pass it straight back in, the projection being
        idempotent for a point already on the wall.
        """
        if vkey not in list(self.mesh.vertices()):
            return self._refuse('that corner is not part of the layout')
        if project:
            xyz = self.project(xyz, self.is_vertex_on_boundary(vkey))
        self.mesh.vertex_attributes(vkey, 'xyz', [xyz[0], xyz[1], 0.0])
        self.last_reason = ''
        return True

    # ------------------------------------------------------------------
    # the curve map
    # ------------------------------------------------------------------

    @staticmethod
    def _curve_key(pa, pb):
        """Geometric key for the curve map. 3 decimals, as everything else."""
        return (round(pa[0], PRECISION), round(pa[1], PRECISION),
                round(pb[0], PRECISION), round(pb[1], PRECISION))

    def _own_curve(self, pa, pb):
        """The shape registered for this edge, either way round, or ``None``."""
        curve = self.curves.get(self._curve_key(pa, pb))
        if curve is not None:
            return [list(p) for p in curve]
        curve = self.curves.get(self._curve_key(pb, pa))
        if curve is not None:
            return [list(p) for p in reversed(curve)]
        return None

    def _edge_curve(self, pa, pb):
        """The shape this coarse edge follows: drawn, or traced.

        The traced half matters as much as the drawn one. Cutting across an
        interior edge splits the separatrix it followed, and the two halves
        have to be told what they are -- left to itself, ``edges_to_curves``
        reaches its warp branch, finds the WHOLE separatrix still anchored at
        the far corner, and warps all of it onto a half-length edge, which
        folds.
        """
        curve = self._own_curve(pa, pb)
        if curve is not None:
            return curve

        tol = 1e-3
        for polyline in (getattr(self.decomposition, 'polylines', None) or []):
            points = [list(p)[:3] for p in polyline]
            # Two points is a chord and carries no shape worth splitting.
            if len(points) < 3:
                continue
            if (distance_point_point(points[0], pa) < tol
                    and distance_point_point(points[-1], pb) < tol):
                return points
            if (distance_point_point(points[0], pb) < tol
                    and distance_point_point(points[-1], pa) < tol):
                return list(reversed(points))
        return None

    @staticmethod
    def _clean(points):
        """Drop consecutive duplicates -- ``Polyline.point_at`` divides by them."""
        out = [list(points[0])]
        for point in points[1:]:
            if distance_point_point(out[-1], point) > 1e-9:
                out.append(list(point))
        return out

    def _split_curve(self, curve, point):
        """``(point ON the curve, head, tail)`` for the curve cut at ``point``.

        The returned point is the one to put the new corner at, and it is not
        the point handed in: a corner that a cut puts halfway along a curved
        edge belongs ON that curve, not on the chord the layout draws for it.
        ``(None, None, None)`` when the split leaves a stub.
        """
        best = None
        for i in range(len(curve) - 1):
            q = closest_point_on_segment_xy(point, (curve[i], curve[i + 1]))
            d = distance_point_point(point, q)
            if best is None or d < best[0]:
                best = (d, i, [q[0], q[1], 0.0])
        if best is None:
            return None, None, None

        _d, i, q = best
        head = self._clean([list(p) for p in curve[:i + 1]] + [q])
        tail = self._clean([q] + [list(p) for p in curve[i + 1:]])
        if len(head) < 2 or len(tail) < 2:
            return None, None, None
        return q, head, tail

    # ------------------------------------------------------------------
    # planning a cut
    # ------------------------------------------------------------------

    def _node_point(self, node):
        """The xyz a plan node sits at."""
        if node[0] == 'vertex':
            return self.mesh.vertex_coordinates(node[1])
        return node[2]

    @staticmethod
    def _direction(points, index, start):
        """Curve direction leaving ``start`` on segment ``index``."""
        for i in range(index, len(points) - 1):
            dx = points[i + 1][0] - start[0]
            dy = points[i + 1][1] - start[1]
            if (dx * dx + dy * dy) > 1e-18:
                return [dx, dy, 0.0]
            start = points[i + 1]
        return None

    def _face_across(self, edge, direction):
        """The face on the side of ``edge`` the curve is heading into.

        Decided by which side of the edge the candidate's own CENTROID is on,
        rather than by the halfedge convention, so it does not depend on the
        cycles having been unified one way round. ``None`` when the curve runs
        along the edge instead of across it, or when that side is open.
        """
        mesh = self.mesh
        u, v = edge
        pu = mesh.vertex_coordinates(u)
        pv = mesh.vertex_coordinates(v)
        ex, ey = pv[0] - pu[0], pv[1] - pu[1]
        side = ex * direction[1] - ey * direction[0]
        if abs(side) <= 1e-12:
            return None
        for fkey in (mesh.halfedge[u][v], mesh.halfedge[v][u]):
            if fkey is None:
                continue
            centroid = mesh.face_centroid(fkey)
            cs = ex * (centroid[1] - pu[1]) - ey * (centroid[0] - pu[0])
            if cs * side > 0.0:
                return fkey
        return None

    def _face_around(self, vkey, points, index, start):
        """The face at a corner the curve heads into -- containment of a probe."""
        direction = self._direction(points, index, start)
        if direction is None:
            return None
        length = (direction[0] ** 2 + direction[1] ** 2) ** 0.5
        step = 0.25 * self.mean_edge() / length
        probe = [start[0] + direction[0] * step, start[1] + direction[1] * step, 0.0]
        for fkey in self.mesh.vertex_faces(vkey):
            polygon = [self.mesh.vertex_coordinates(w)
                       for w in self.mesh.face_vertices(fkey)]
            if is_point_in_polygon_xy(probe, polygon):
                return fkey
        return None

    def _exit_of(self, fkey, points, index, start, end):
        """Where the curve leaves patch ``fkey``, entering it at ``start``.

        Returns ``(node, index, subcurve, reason)``. ``node`` is ``None`` on a
        refusal and ``reason`` says why. ``subcurve`` is the piece of the drawn
        curve inside this patch, which is what the new coarse edge will follow.
        """
        mesh = self.mesh
        tol = self.on_tol()
        skip = max(self.corner_tol(), tol)
        halfedges = mesh.face_halfedges(fkey)

        cursor = list(start)
        for i in range(index, len(points) - 1):
            a = cursor if i == index else points[i]
            b = points[i + 1]
            best = None
            for u, v in halfedges:
                segment = (mesh.vertex_coordinates(u), mesh.vertex_coordinates(v))
                x = intersection_segment_segment_xy((a, b), segment, tol)
                if x is None:
                    continue
                # The way IN is an intersection too, and it is not the way out.
                if distance_point_point(x, start) <= skip:
                    continue
                t = distance_point_point(a, x)
                if best is None or t < best[0]:
                    best = (t, [x[0], x[1], 0.0], (u, v))
            if best is None:
                continue
            _t, x, edge = best
            sub = self._clean([list(start)] + [list(p) for p in points[index + 1:i + 1]] + [x])
            return self._node_at(fkey, edge, x), i, sub, ''

        # No crossing left: the curve stops inside this patch. That is only
        # legal if it stopped ON this patch's own rim -- which the exact
        # arithmetic above can miss by a float when the last point IS the exit.
        point = self._node_point(end)
        corners = set(mesh.face_vertices(fkey))
        on_rim = ((end[0] == 'vertex' and end[1] in corners)
                  or (end[0] == 'edge' and (tuple(end[1]) in halfedges
                                            or tuple(reversed(end[1])) in halfedges)))
        if on_rim and distance_point_point(point, start) > skip:
            sub = self._clean([list(start)]
                              + [list(p) for p in points[index + 1:-1]] + [list(point)])
            node = end if end[0] == 'vertex' else self._node_at(fkey, end[1], point)
            return node, len(points) - 2, sub, ''

        return None, index, None, (
            'the line stops inside a patch. It has to reach a patch edge -- and '
            'to leave the layout all-quad, a boundary')

    def _node_at(self, fkey, edge, point):
        """An exit node, collapsed onto a corner when it lands on one."""
        for vkey in self.mesh.face_vertices(fkey):
            if distance_point_point(point, self.mesh.vertex_coordinates(vkey)) <= self.corner_tol():
                return ('vertex', vkey)
        return ('edge', tuple(edge), point)

    def _plan_cut(self, points):
        """**Walk the drawn curve across the layout, without touching it.**

        Produces ``{'nodes': [...], 'faces': [(fkey, i, j, subcurve), ...]}``:
        every corner the cut needs and every patch it splits, entry node ``i``
        to exit node ``j``. Nothing is mutated here, so every refusal below
        costs the user nothing but the message.
        """
        points = [[float(p[0]), float(p[1]), 0.0] for p in points]
        points = self._clean(points)
        if len(points) < 2:
            return None, 'the drawn curve has no length'

        start = self.locate(points[0])
        if start is None:
            return None, ('the line does not START on the layout -- pick its '
                          'first point on a coarse edge')
        end = self.locate(points[-1])
        if end is None:
            return None, ('the line does not END on the layout -- pick its last '
                          'point on a coarse edge')
        points[0] = list(self._node_point(start))
        points[-1] = list(self._node_point(end))
        if distance_point_point(points[0], points[-1]) <= self.corner_tol():
            return None, 'the line starts and ends at the same corner'

        nodes = [start]
        faces = []
        index, entry = 0, 0
        while True:
            here = self._node_point(nodes[entry])
            if nodes[entry][0] == 'vertex':
                fkey = self._face_around(nodes[entry][1], points, index, here)
            else:
                direction = self._direction(points, index, here)
                fkey = (self._face_across(nodes[entry][1], direction)
                        if direction is not None else None)
            if fkey is None:
                return None, ('the line leaves the layout, or runs along a patch '
                              'edge instead of across it')
            if any(item[0] == fkey for item in faces):
                return None, ('the line crosses the same patch twice. Draw it as '
                              'two cuts')
            if self.mesh.is_face_pseudo_quad(fkey):
                return None, ('the line crosses a pseudo-quad (a patch with a '
                              'collapsed corner, drawn as a triangle). Its fourth '
                              'side is a pole, so there is nothing there to cut')

            node, index, sub, reason = self._exit_of(fkey, points, index, here, end)
            if node is None:
                return None, reason
            nodes.append(node)
            faces.append((fkey, entry, len(nodes) - 1, sub))
            entry = len(nodes) - 1

            if distance_point_point(self._node_point(node), points[-1]) <= self.on_tol():
                break
            if len(faces) > self.mesh.number_of_faces():
                return None, 'the line does not terminate on the layout'

        return {'nodes': nodes, 'faces': faces}, ''

    def _open_ends(self, plan):
        """The ends of a planned cut that stop on an INTERIOR edge.

        Each one is a patch left with five sides, and they are the only reason
        a cut that is otherwise perfectly good gets refused -- which is what
        :meth:`_extend_plan` is for. Returned as ``(node index, the patch the
        cut came from)`` so the extension knows which way to go.
        """
        out = []
        if not plan['faces']:
            return out
        for index, fkey in ((plan['faces'][0][1], plan['faces'][0][0]),
                            (plan['faces'][-1][2], plan['faces'][-1][0])):
            node = plan['nodes'][index]
            if node[0] != 'edge':
                continue
            if not self.mesh.is_edge_on_boundary(node[1]):
                out.append((index, fkey))
        return out

    def _extend_plan(self, plan, index, fkey):
        """**Carry a cut on to the boundary, patch by patch.**

        A cut that stops on an interior edge is not wrong, it is unfinished:
        the patch behind that edge has to be split too, and then the one behind
        that, until the run reaches a wall. Which is the same thing as saying a
        cut follows a STRIP -- each patch is entered and left through opposite
        sides, so the continuation is ``face_opposite_edge`` at the same
        parameter, and the layout stays all-quad the whole way.

        Straight lines, deliberately: this is the part the user did not draw,
        so it should look like what it is. Their own curve keeps its shape.
        """
        mesh = self.mesh
        while True:
            node = plan['nodes'][index]
            edge, point = node[1], node[2]
            nxt = None
            for candidate in (mesh.halfedge[edge[0]][edge[1]],
                              mesh.halfedge[edge[1]][edge[0]]):
                if candidate is not None and candidate != fkey:
                    nxt = candidate
            if nxt is None:
                return True, ''
            if any(item[0] == nxt for item in plan['faces']):
                return False, ('extending the cut runs back into a patch it has '
                               'already crossed -- draw it as two cuts')
            if len(mesh.face_vertices(nxt)) != 4 or mesh.is_face_pseudo_quad(nxt):
                return False, ('the cut would have to continue through a patch '
                               'that is not a quad, and there is no opposite '
                               'side there to continue to')

            # The halfedge as THIS patch sees it, or ``face_opposite_edge``
            # answers for the patch on the other side.
            uv = None
            for u, v in mesh.face_halfedges(nxt):
                if set((u, v)) == set(edge):
                    uv = (u, v)
            if uv is None:
                return False, 'the cut lost track of the patch it was crossing'
            u, v = uv
            w, x = mesh.face_opposite_edge(u, v)

            pu, pv = mesh.vertex_coordinates(u), mesh.vertex_coordinates(v)
            span = distance_point_point(pu, pv)
            if span <= 0.0:
                return False, 'the cut would continue across a zero-length edge'
            t = distance_point_point(pu, closest_point_on_segment_xy(point, (pu, pv))) / span
            # (u, v, w, x) go round the quad, so u pairs with x and v with w:
            # the point t of the way from u to v is t of the way from x to w.
            pw, px = mesh.vertex_coordinates(w), mesh.vertex_coordinates(x)
            q = [px[0] + t * (pw[0] - px[0]), px[1] + t * (pw[1] - px[1]), 0.0]

            plan['nodes'].append(self._node_at(nxt, (w, x), q))
            plan['faces'].append((nxt, index, len(plan['nodes']) - 1,
                                  [list(point), list(q)]))
            index = len(plan['nodes']) - 1
            fkey = nxt
            if plan['nodes'][index][0] != 'edge':
                return True, ''
            if mesh.is_edge_on_boundary(plan['nodes'][index][1]):
                return True, ''

    # ------------------------------------------------------------------
    # applying a cut
    # ------------------------------------------------------------------

    def _current_edge(self, work, u, v, point):
        """The piece of edge ``(u, v)`` that ``point`` is on, after any splits."""
        if work.has_edge((u, v)):
            return (u, v)
        pu = work.vertex_coordinates(u)
        pv = work.vertex_coordinates(v)
        tol = self.on_tol()
        for a, b in work.edges():
            pa = work.vertex_coordinates(a)
            pb = work.vertex_coordinates(b)
            if (distance_point_point(pa, closest_point_on_segment_xy(pa, (pu, pv))) > tol
                    or distance_point_point(pb, closest_point_on_segment_xy(pb, (pu, pv))) > tol):
                continue
            if distance_point_point(point, closest_point_on_segment_xy(point, (pa, pb))) <= tol:
                return (a, b)
        return None

    def _split_edge_at(self, work, edge, point, pending):
        """Split one coarse edge, putting the new corner where it belongs.

        Three things happen that a bare ``split_edge`` does not do:

        * the corner is placed on the edge's own CURVE where it has one, and
          projected onto the domain wall where the edge is a naked one, so the
          layout does not lose the area between chord and arc and
          ``edit_coarse``'s own snap has nothing left to do -- which is what
          keeps the registered curve's endpoint an exact match at commit;
        * both halves of a curved edge inherit their piece of it;
        * an edge already split by an earlier crossing is resolved to the piece
          the point is actually on.
        """
        found = self._current_edge(work, edge[0], edge[1], point)
        if found is None:
            return None
        u, v = found
        pu = work.vertex_coordinates(u)
        pv = work.vertex_coordinates(v)

        curve = self._edge_curve(pu, pv)
        exact, head, tail = (self._split_curve(curve, point)
                             if curve is not None else (None, None, None))
        if exact is None:
            exact = self.project(list(point), work.is_edge_on_boundary((u, v)))
            head = tail = None

        span = distance_point_point(pu, pv)
        if span <= 0.0:
            return None
        q = closest_point_on_segment_xy(point, (pu, pv))
        t = distance_point_point(pu, q) / span
        t = min(max(t, 1e-6), 1.0 - 1e-6)

        w = work.split_edge((u, v), t, allow_boundary=True)
        if w is None or w in (u, v):
            return None
        work.vertex_attributes(w, 'xyz', [exact[0], exact[1], 0.0])

        if head is not None and tail is not None:
            pending[self._curve_key(pu, exact)] = head
            pending[self._curve_key(exact, pv)] = tail
        return w

    def _apply_cut(self, work, plan):
        """Perform a planned cut on ``work``. ``(ok, reason, curves)``."""
        pending = {}
        vertex_of = {}
        for index, node in enumerate(plan['nodes']):
            if node[0] == 'vertex':
                vertex_of[index] = node[1]
                continue
            w = self._split_edge_at(work, node[1], node[2], pending)
            if w is None:
                return False, 'a coarse edge could not be split there', {}
            vertex_of[index] = w

        for fkey, i, j, _sub in plan['faces']:
            a, b = vertex_of[i], vertex_of[j]
            if a == b:
                return False, 'the line enters and leaves a patch at one corner', {}
            try:
                work.split_face(fkey, a, b)
            except ValueError:
                # ``split_face`` refuses neighbouring corners, which is exactly
                # the cut that would slice a triangle off the patch.
                return False, ('the line cuts a corner off a patch instead of '
                               'crossing it: it enters and leaves through sides '
                               'that meet'), {}

        ok, reason = self._check_quads(work)
        if not ok:
            return False, reason, {}

        for fkey, i, j, sub in plan['faces']:
            pa = work.vertex_coordinates(vertex_of[i])
            pb = work.vertex_coordinates(vertex_of[j])
            curve = self._clean([list(pa)] + [list(p) for p in sub[1:-1]] + [list(pb)])
            # Two points IS the chord, which is what an unregistered edge
            # densifies as. Registering it would only add noise to the map.
            if len(curve) > 2:
                pending[self._curve_key(pa, pb)] = curve
        return True, '', pending

    def _check_quads(self, work):
        """**The gate: the edit has to leave every patch four-sided.**

        Run on the copy, before it is adopted. A pseudo-quad that was already
        there stays legal -- it is a quad with one side collapsed and
        ``densification`` handles it -- but one that the edit turned into
        something else does not, because its ``face_pole`` entry now describes
        a face that no longer exists.
        """
        face_pole = work.attributes.get('face_pole') or {}
        for fkey in work.faces():
            n = len(work.face_vertices(fkey))
            if n == 4 and fkey not in face_pole:
                continue
            if n == 3 and fkey in face_pole:
                continue
            x, y, _z = work.face_centroid(fkey)
            if fkey in face_pole:
                return False, ('the edit runs into the pseudo-quad near '
                               '({:.2f}, {:.2f}) and would break its pole'.format(x, y))
            return False, (
                'it would leave the patch near ({:.2f}, {:.2f}) with {} sides. A '
                'cut has to enter and leave each patch through OPPOSITE sides, '
                'and every edge it crosses has to be crossed through -- so both '
                'its ends belong on the layout boundary. Continue the line to a '
                'wall, or start it on one'.format(x, y, n))
        return True, ''

    def _refuse(self, reason):
        """Record why an operation was not performed. Always ``False``.

        Recorded rather than printed: this module has no idea whether it is
        driving a Rhino dialog, a log or a test.
        """
        self.last_reason = reason
        return False

    def insert_curve(self, points, extend=None):
        """**Cut the layout with a drawn curve.** ``True`` if it was.

        Plan against the untouched layout, apply to a copy, adopt the copy only
        if it passes :meth:`_check_quads`. A refusal records why in
        :attr:`last_reason` and leaves the mesh exactly as it was, which is the
        whole reason for the copy.

        Parameters
        ----------
        points : list[[x, y, z]]
            The curve, already sampled. The samples ARE the shape from here on:
            an arc handed over as two points is a chord.
        extend : bool or callable, optional
            What to do when the curve stops on an interior edge -- the case
            that leaves a five-sided patch behind. ``None`` refuses, ``True``
            carries the cut on to the boundary (:meth:`_extend_plan`), and a
            callable is asked, with the number of ends that need it, and
            answers.
        """
        plan, reason = self._plan_cut(points)
        if plan is None:
            return self._refuse(reason)

        ends = self._open_ends(plan)
        if ends:
            wanted = extend(len(ends)) if callable(extend) else bool(extend)
            if not wanted:
                return self._refuse(
                    'the line stops on an interior patch edge, which would '
                    'leave the patch behind it with five sides. Either draw it '
                    'out to a boundary, or let the command carry it there')
            for index, fkey in ends:
                ok, reason = self._extend_plan(plan, index, fkey)
                if not ok:
                    return self._refuse(reason)

        work = self.mesh.copy()
        ok, reason, curves = self._apply_cut(work, plan)
        if not ok:
            return self._refuse(reason)

        self.mesh = work
        self.curves.update(curves)
        self.last_reason = ''
        self.last_cut = {
            'faces': len(plan['faces']),
            'corners': sum(1 for node in plan['nodes'] if node[0] == 'edge'),
            'faces_out': work.number_of_faces(),
        }
        return True

    # ------------------------------------------------------------------
    # deleting a strip
    # ------------------------------------------------------------------

    def _trial_delete(self, edge):
        """Delete the strip through ``edge`` on a COPY. ``(mesh or None, info)``.

        One code path for both :meth:`plan_strip_deletion` and
        :meth:`delete_strip`, so what the user is told and what then happens
        cannot drift apart. Nothing here touches ``self``.
        """
        info = {'skey': None, 'faces': 0, 'collateral': 0, 'boundaries_lost': 0,
                'ok': False, 'reason': ''}

        work = self.mesh.copy()
        work.collect_strips()
        skey = work.edge_strip(tuple(edge))
        if skey is None:
            info['reason'] = 'that edge belongs to no strip -- nothing to remove'
            return None, info

        info['skey'] = skey
        info['faces'] = len(work.strip_faces(skey))
        info['collateral'] = len(collateral_strip_deletions(work, [skey]))
        info['boundaries_lost'] = len(total_boundary_deletions(work, [skey]))

        try:
            _grammar_delete_strip(work, skey)
        except Exception as e:
            # Measured: 1 strip of 24 on a two-hole plate raises ``KeyError``
            # from inside the grammar. Catching it here costs the user a
            # refusal instead of a traceback.
            info['reason'] = 'the grammar could not delete that strip ({}: {})'.format(
                type(e).__name__, e)
            return None, info

        if work.number_of_faces() == 0:
            info['reason'] = 'that would delete every patch of the layout'
            return None, info

        ok, reason = self._check_quads(work)
        if not ok:
            info['reason'] = reason
            return None, info

        # **The gate that actually decides it, and it is not about poles.**
        # ``delete_strip`` welds the two sides of the strip together, and on a
        # layout where the two sides are not simply opposite that weld can put
        # three patches on one edge. Measured over four layouts, sweeping every
        # strip: square+square-hole 8/8 fine, L-plate 4/4 fine,
        # square+two-holes 23/24 fine, but square+one-circular-hole **10 of 20
        # non-manifold**. Poles do NOT predict it -- in that layout all 7
        # strips WITH pole edges deleted cleanly and every non-manifold one had
        # none. So there is nothing to refuse in advance; run it and look.
        # ``edit_coarse`` refuses a non-manifold layout anyway, so catching it
        # here turns a failure at commit into a refusal that costs nothing.
        if not work.is_manifold():
            info['reason'] = (
                'deleting that strip would leave the layout non-manifold -- '
                'welding its two sides together puts more than two patches on '
                'one edge. Nothing has been changed; try a different strip')
            return None, info

        info['ok'] = True
        return work, info

    def plan_strip_deletion(self, edge):
        """**What deleting the strip through** ``edge`` **would cost.** No mutation.

        Deleting a strip is not a local edit, and a front end has to be able to
        say so before it happens. Three things reach past the picked strip, and
        all three are predicted here, on a copy of the untouched layout:

        * it MERGES the two sides of the strip -- for each disconnected part of
          the strip's edge network one new corner is added at the centroid and
          every old corner is substituted for it, so surrounding corners move;
        * it can take other strips with it: any strip whose faces all lie
          inside the deleted one goes too;
        * it can COLLAPSE a boundary, which on a coarse layout means losing a
          hole the domain still has.

        **The prediction is the operation.** Rather than guess from the strip's
        shape, this performs the deletion on a copy and reports what came out,
        so ``ok`` is the truth and not an estimate -- see :meth:`_trial_delete`
        for why nothing cheaper works.

        Returns
        -------
        dict
            ``skey``, ``faces``, ``collateral``, ``boundaries_lost``, ``ok``
            and ``reason``. ``ok`` False means :meth:`delete_strip` will refuse,
            and ``reason`` says why.
        """
        _work, info = self._trial_delete(edge)
        self.last_reason = '' if info['ok'] else info['reason']
        return info

    def delete_strip(self, edge):
        """**Delete the strip through** ``edge``. ``True`` if it was.

        The library's ``grammar_pattern.delete_strip``, applied to a copy and
        adopted only if the copy survives every gate -- the same discipline
        :meth:`insert_curve` uses, so a refusal costs nothing and leaves the
        layout exactly as it was.

        This is what ``CMD_edit_quad_mesh``'s ``remove_line`` does to the dense
        mesh, and it means the same thing here: the band of patches collapses
        and its two sides weld together. It does NOT dissolve the picked line
        and merge the patches either side of it -- see the module docstring.

        Note what is deliberately NOT offered: ``preserve_boundaries``. It
        works by ``split_strips`` -> ``add_strip`` -> ``func_1``, a
        20-iteration constrained smooth, and on a COARSE layout that smoothing
        slides corners along the chorded boundary and clusters them -- measured
        gaps of 0.096 / 0.071 / 0.130 ... 2.782 on a 4.389 loop, min angle
        0.28 degrees. Predict the collapse, say so, and let the user decide.
        """
        work, info = self._trial_delete(edge)
        if work is None:
            return self._refuse(info['reason'])

        faces_in = self.mesh.number_of_faces()
        self.mesh = work
        self.last_reason = ''
        self.last_deletion = dict(info, faces_in=faces_in,
                                  faces_out=work.number_of_faces())
        return True

    # ------------------------------------------------------------------
    # undo, and handing it back
    # ------------------------------------------------------------------

    def reset(self):
        """Put the layout back as it was before this round of edits.

        Restores the whole MESH, not just the corner positions: a cut or a
        deletion changes the topology, and putting coordinates back would leave
        the patches it created in place with nothing to say where they came
        from. The snapshot is taken at construction and again after every
        successful :meth:`commit`, so reset undoes what has not been committed.
        """
        self.mesh = self._snapshot.copy()
        self.curves = dict(self._snapshot_curves)
        self.last_reason = ''

    def _rekey_curves(self, mesh):
        """Re-key the curve map onto a layout whose vertex keys were renumbered.

        ``edit_coarse`` welds and repairs, so it returns a mesh with different
        keys for the same corners -- 2 of 8 unchanged, measured. The map was
        keyed on rounded coordinates precisely so this is a re-match and not a
        loss; a curve whose edge did not survive the repair is dropped and
        counted, because the alternative is publishing a polyline that
        ``edges_to_curves`` will hand to an edge it does not belong to.
        """
        tol = 1e-3
        kept = {}
        lost = 0
        edges = [(mesh.vertex_coordinates(u), mesh.vertex_coordinates(v))
                 for u, v in mesh.edges()]
        for curve in self.curves.values():
            match = None
            for pa, pb in edges:
                if (distance_point_point(curve[0], pa) < tol
                        and distance_point_point(curve[-1], pb) < tol):
                    match = (pa, pb)
                    break
                if (distance_point_point(curve[0], pb) < tol
                        and distance_point_point(curve[-1], pa) < tol):
                    match = (pb, pa)
                    break
            if match is None:
                lost += 1
                continue
            curve = [list(p) for p in curve]
            curve[0], curve[-1] = list(match[0]), list(match[1])
            kept[self._curve_key(curve[0], curve[-1])] = curve
        self.curves = kept
        return lost

    def commit(self):
        """**Hand the edited layout to** :meth:`FieldDecomposition.edit_coarse`.

        This is the step that makes the edit worth anything: the moved edges
        are re-matched to the separatrices they were traced from and those
        curves are warped onto the new corners, so the layout keeps its field
        alignment. Skipping it and densifying the mutated mesh directly gives
        straight chords between the moved corners.

        ``edit_coarse`` deliberately RAISES on a layout that will not densify
        rather than falling back, so the failure is caught here and returned
        while the user can still fix it. On failure the mesh is left as edited
        so the offending corner can be moved again.

        Returns
        -------
        (bool, dict)
            Whether the layout was adopted, and ``edit_coarse``'s own notes --
            ``faces_in``, ``faces_out``, ``snapped``, ``off_wall``, ``sides``
            -- plus ``lost_curves``, the number of drawn shapes that no longer
            match a coarse edge. On failure, ``{'error': ...}``.
        """
        kwargs = {} if self.snap_tol is None else {'snap_tol': self.snap_tol}
        try:
            mesh = self.decomposition.edit_coarse(self.mesh, **kwargs)
        except ValueError as e:
            self.last_reason = str(e)
            return False, {'error': str(e)}

        # ``edit_coarse`` welds and repairs, so it returns a DIFFERENT mesh
        # object whose vertex keys are RENUMBERED -- measured 2 of 8 unchanged
        # on a disc. Rebind, or the next move edits an orphan the decomposition
        # no longer knows about, and re-key anything held by index.
        self.mesh = mesh
        lost = self._rekey_curves(mesh)
        # Publish the drawn shapes. This is the whole reason an inserted arc is
        # an arc in the final mesh rather than a chord: it puts them in the
        # polyline network, where ``edges_to_curves`` finds each added edge by
        # its endpoints. It happens HERE, after the re-key, because the
        # corrected endpoints are only known now -- and it goes through
        # ``set_user_curves`` rather than assigning the attribute, because
        # ``edit_coarse`` has already built ``polylines`` and would not see it.
        self.decomposition.set_user_curves(self.curves.values())
        self._snapshot = mesh.copy()
        self._snapshot_curves = dict(self.curves)
        self.last_reason = ''

        notes = dict(self.decomposition.edit_notes or {})
        notes['lost_curves'] = lost
        return True, notes
