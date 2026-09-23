"""**What every hand-editor of a quad mesh needs, with no CAD in it.**

:class:`MeshEditor` is the half that ``CoarseEditor`` and
``DenseMeshEditor`` had written out twice: holding a mesh being edited, holding
the domain walls a boundary vertex is kept on, moving one vertex, undoing, and
deleting the strip through an edge. Nothing here picks, prompts, draws or prints
-- a refusal is recorded in :attr:`MeshEditor.last_reason` and returned, and
whoever is driving decides how to say so.

**The mesh is edited on a COPY, and the caller's object is untouched until a
commit.** :attr:`MeshEditor.target` is what was handed in and
:attr:`MeshEditor.mesh` is the working copy. That is a tightening: previously
``move_vertex`` wrote straight into the caller's mesh while the topological
operations replaced ``self.mesh`` with a copy, so half an edit was visible
outside the editor before anything was committed -- which contradicts the promise
the Rhino command already makes, that nothing is written until commit. It also
makes :meth:`MeshEditor.reset` mean what it says for a move.

**Refusals return** ``(False, {'error': reason})`` **and successes**
``(True, notes)``, and both subclasses publish the same shape. ``DenseMeshEditor``
used to return bare booleans; it moved to this convention when its Rhino command
was rebuilt on it.

**Deleting a strip is a template, because the two subclasses gate it
differently.** The walk is identical -- resolve the edge to a strip, optionally
pre-split to save a boundary, run the grammar on a copy, check the copy -- but
what counts as an acceptable result is not, so :meth:`MeshEditor._gate` is left
to the subclass. On a coarse layout that gate is all-quad **plus manifold**, and
the manifold half is the one that actually decides it: swept over four layouts, a
square with one circular hole had **10 of 20** strips come out non-manifold,
while every one of its 7 pole-bearing strips deleted cleanly. Nothing about a
strip predicts it, which is why :meth:`MeshEditor.plan_strip_deletion` performs
the deletion on a copy and reports what came out instead of estimating.

**The plan and the deletion are the same code path.** They have to be: a front
end puts the numbers from the plan in front of a user, and if the deletion then
did something else the confirmation would be a lie.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from copy import deepcopy

from compas_singular.datastructures.mesh_quad.grammar.add_strip import split_strips
from compas_singular.datastructures.mesh_quad.grammar.delete_strip import collateral_strip_deletions
from compas_singular.datastructures.mesh_quad.grammar.delete_strip import delete_strip as _grammar_delete_strip
from compas_singular.datastructures.mesh_quad.grammar.delete_strip import (
    strips_to_split_to_prevent_boundary_collapse)
from compas_singular.datastructures.mesh_quad.grammar.delete_strip import total_boundary_deletions
from compas_singular.datastructures.mesh_quad_coarse.coarse_curves import BoundaryLoop
from compas_singular.datastructures.mesh_quad_coarse.coarse_curves import mean_edge_length


__all__ = ['MeshEditor', 'boundary_vertex_set']


def boundary_vertex_set(mesh):
    """Every vertex with a faceless halfedge. See :meth:`MeshEditor.boundary_vertices`."""
    out = set()
    for u, nbrs in mesh.halfedge.items():
        for v, fkey in nbrs.items():
            if fkey is None:
                out.add(u)
                out.add(v)
    return out


class MeshEditor(object):
    """A quad mesh being hand-edited. No CAD, no prompts, no prints.

    Parameters
    ----------
    mesh : :class:`compas_singular.datastructures.QuadMesh`
        The mesh to edit. Kept as :attr:`target`; the editing happens on a copy.
    walls : list, optional
        The domain walls, as point lists, closed polylines or
        :class:`BoundaryLoop` instances. A boundary vertex that moves is projected
        back onto the nearest one, so nudging the edge of the mesh does not eat
        the outline. Without walls a moved boundary vertex is left where it was
        put.
    work_on_copy : bool, optional
        Edit a copy and leave ``mesh`` untouched until a commit transplants into
        it. Default ``True``; ``False`` edits in place, for a caller that has
        already made its own copy.

    Attributes
    ----------
    target : mesh
        The object handed in. Only a commit writes to it.
    mesh : mesh
        The working copy. Replaced wholesale by every topological operation that
        succeeds, so hold the editor, not the mesh.
    walls : list
        The domain walls as :class:`BoundaryLoop` instances.
    last_reason : str
        Why the last operation was refused. Empty when nothing was.
    """

    #: Whether a moved vertex is flattened onto z = 0. A coarse layout is planar
    #: by construction and its subclass turns this on; a dense mesh need not be.
    PLANAR = False

    #: ``attributes`` keys that describe the layout as it was BEFORE this edit --
    #: a DENSE mesh built from it, or geometry keyed to its vertices. They go
    #: stale the moment the topology changes, and carrying them across a rebuild
    #: is how a stale dense mesh gets served for an edited layout.
    #:
    #: ``edges_to_curves`` is here for the second reason and is worth spelling
    #: out: it maps ``(u, v)`` to the polyline that edge densifies along, so a
    #: rebuild that renumbers points it at the wrong edges, and a MOVED corner
    #: leaves it anchored where the corner used to be -- the two patches sharing
    #: that edge would then be handed different endpoints and fail to weld.
    #: Empty on any layout that did not come from ``from_coarse_polylines``, so
    #: clearing it is a no-op everywhere else.
    DERIVED_ATTRIBUTES = ('quad_mesh', 'polygonal_mesh',
                          'vertex_coarse_to_dense', 'edge_coarse_to_dense',
                          'edges_to_curves')

    #: ``attributes`` keys that are keyed BY STRIP INDEX, and so are meaningless
    #: once a rebuild has renumbered the strips.
    STRIP_KEYED_ATTRIBUTES = ('strips', 'polyedges', 'strips_density')

    def __init__(self, mesh, walls=None, work_on_copy=True):
        if mesh is None:
            raise ValueError('no mesh to edit')
        self.target = mesh
        self.mesh = mesh.copy() if work_on_copy else mesh
        self.walls = self._as_walls(walls)
        #: The mesh as it was before the current round of edits. A whole mesh
        #: rather than a coordinate map: the topological operations change the
        #: connectivity, and putting coordinates back would leave the faces they
        #: created in place with nothing to say where they came from.
        self._snapshot = self.mesh.copy()
        self.last_reason = ''
        self.last_deletion = {}
        #: Whether anything has been changed since this editor was built. What
        #: the curvature branch keys off: on an UNEDITED layout every coarse edge
        #: still matches the separatrix it was traced from exactly, so there is
        #: nothing to warp and trying would only risk a wrong match.
        self.edited = False
        #: One entry per completed edit, oldest first -- see :meth:`push_undo`.
        self._undo_stack = []

    # ------------------------------------------------------------------
    # walls
    # ------------------------------------------------------------------

    @staticmethod
    def _as_walls(walls):
        """``BoundaryLoop`` per wall. Accepts loops already built as one."""
        out = []
        for wall in (walls or []):
            if isinstance(wall, BoundaryLoop):
                out.append(wall)
                continue
            points = list(wall)
            if len(points) >= 3:
                out.append(BoundaryLoop(points))
        return out

    def project_to_wall(self, xyz):
        """The nearest point of the nearest wall, or ``xyz`` if there is no wall.

        Only ever called for a vertex that is TOPOLOGICALLY on the mesh boundary.
        Projecting anything merely near a wall would drag an interior vertex that
        legitimately sits in a narrow slot onto it -- the same rule
        ``rebuild.snap_to_loops`` keeps, and for the same reason.
        """
        if not self.walls:
            return list(xyz)
        best_d, best_p = float('inf'), None
        for wall in self.walls:
            d, _s, q = wall.project(xyz)
            if d < best_d:
                best_d, best_p = d, q
        return list(best_p) if best_p is not None else list(xyz)

    # ------------------------------------------------------------------
    # measuring the mesh
    # ------------------------------------------------------------------

    def mean_edge(self, mesh=None):
        """Mean edge length -- the scale every tolerance here is expressed in."""
        return mean_edge_length(self.mesh if mesh is None else mesh)

    def boundary_vertices(self, mesh=None):
        """Every vertex on every boundary -- holes included -- as a set.

        Read straight off the halfedges: a vertex is on a boundary exactly when
        one of its halfedges, either way round, has no face. Two ways of getting
        this wrong came first:

        * ``vertices_on_boundary`` -- singular -- returns only the LONGEST
          boundary, so on a mesh with a hole every vertex of that hole was
          reported as interior. Measured on an annulus: 12 of 16 seen;
        * ``vertices_on_boundaries`` -- plural -- walks each boundary as a loop,
          and on a mesh where two faces only touch at a corner that walk can
          NEVER END. Measured on a dense mesh after random removals: it spun
          until killed, from inside ``move_vertex``. A set needs no walk.
        """
        mesh = self.mesh if mesh is None else mesh
        return boundary_vertex_set(mesh)

    def is_vertex_on_boundary(self, vkey, mesh=None):
        """Whether a vertex is on ANY boundary of the mesh -- holes included."""
        return vkey in self.boundary_vertices(mesh)

    def check_strip_count(self, mesh=None):
        """**Thesis Eq 5.6:** the number of OPEN strips is ``E - 2F``.

        A free post-condition on the strip data. Every boundary edge is one of the
        two extremity edges of an open strip, and summing four edges per face
        counts boundary edges once and interior edges twice, which gives
        ``S_open = E - 2F`` without collecting anything. If the stored strips
        disagree, an edit corrupted them -- so this catches the class of bug where
        a strip table quietly stops describing the mesh it belongs to.

        Returns
        -------
        (bool, int, int)
            Whether it holds, the counted number of open strips, and ``E - 2F``.

        Collected on a COPY, so calling it never leaves strip data on a mesh the
        caller did not ask to have it.
        """
        mesh = self.mesh if mesh is None else mesh
        work = mesh.copy()
        work.collect_strips()
        counted = sum(1 for skey in work.attributes['strips']
                      if not work.is_strip_closed(skey))
        expected = work.number_of_edges() - 2 * work.number_of_faces()
        return counted == expected, counted, expected

    # ------------------------------------------------------------------
    # returning
    # ------------------------------------------------------------------

    def _refuse(self, reason):
        """Record why an operation was not performed. Always ``(False, notes)``.

        Recorded rather than printed: this module has no idea whether it is
        driving a Rhino dialog, a log or a test. :attr:`last_reason` and
        ``notes['error']`` carry the SAME string, because a front end that shows
        one and a test that asserts the other must not be able to disagree.
        """
        self.last_reason = reason
        return False, {'error': reason}

    def _accept(self, **notes):
        """Record a success. Always ``(True, notes)``."""
        self.last_reason = ''
        return True, notes

    # ------------------------------------------------------------------
    # move one vertex
    # ------------------------------------------------------------------

    def move_vertex(self, vkey, xyz, project=True):
        """Move one vertex. ``(ok, notes)``.

        A vertex on the mesh boundary is projected back onto the nearest domain
        wall unless ``project`` is off -- so a caller that already previewed the
        projected point may pass it straight back in, the projection being
        idempotent for a point already on the wall. An interior vertex is never
        projected, whatever its distance.

        Topology is unchanged, so this is the one operation that stays available
        on a mesh with poles.
        """
        if vkey not in self.mesh.vertex:
            return self._refuse('vertex {} is not in the mesh'.format(vkey))

        point = self._as_point(xyz)
        projected = False
        if project and self.is_vertex_on_boundary(vkey):
            point = self._as_point(self.project_to_wall(point))
            projected = True

        self.mesh.vertex_attributes(vkey, 'xyz', point)
        self.edited = True
        return self._accept(moved=vkey, xyz=point, projected=projected)

    def _as_point(self, xyz):
        """``xyz`` as three floats, flattened onto z = 0 when :attr:`PLANAR`."""
        point = list(xyz)
        while len(point) < 3:
            point.append(0.0)
        return [float(point[0]), float(point[1]),
                0.0 if self.PLANAR else float(point[2])]

    # ------------------------------------------------------------------
    # undo
    # ------------------------------------------------------------------

    def snapshot(self):
        """Make the current mesh the state :meth:`reset` returns to."""
        self._snapshot = self.mesh.copy()

    def reset(self):
        """Put the mesh back as it was at construction, or at the last snapshot.

        Restores the whole MESH, not just the coordinates: the topological
        operations change the connectivity, and putting coordinates back would
        leave the faces they created in place with nothing to say where they came
        from.
        """
        self.mesh = self._snapshot.copy()
        self.last_deletion = {}
        self.edited = False
        return self._accept(faces=self.mesh.number_of_faces())

    # ------------------------------------------------------------------
    # undo -- one step back, not all the way to the start
    # ------------------------------------------------------------------
    #
    # :meth:`reset` answers "throw away this whole round of edits"; this
    # answers "that last one, not the others". A front end that acts on a pick
    # immediately -- deleting the strip the moment it is clicked, say, rather
    # than asking first -- needs this to make an accidental pick cheap to walk
    # back, without losing everything edited before it.

    def _state(self):
        """Everything :meth:`undo` needs to put back. A subclass with more
        state than the mesh (a curve map, a dirty flag) extends this and
        :meth:`_restore` together, never one without the other."""
        return {'mesh': self.mesh.copy(), 'edited': self.edited,
                'last_deletion': dict(self.last_deletion)}

    def _restore(self, state):
        """The inverse of :meth:`_state`."""
        self.mesh = state['mesh']
        self.edited = state['edited']
        self.last_deletion = state['last_deletion']

    def push_undo(self):
        """Remember the mesh as it is now, before the change about to happen.

        Call this immediately before an operation that may mutate ``self.mesh``
        -- one entry per attempt, whether or not it succeeds. A refused
        operation never touches the mesh, so its snapshot is a harmless no-op
        to undo back through; :meth:`discard_last_undo` is there for a caller
        that would rather not leave one.
        """
        self._undo_stack.append(self._state())

    def discard_last_undo(self):
        """Drop the most recent :meth:`push_undo` snapshot -- nothing changed."""
        if self._undo_stack:
            self._undo_stack.pop()

    def undo(self):
        """Put back the mesh as it was before the last :meth:`push_undo`. ``(ok, notes)``.

        Refuses with nothing to restore rather than silently doing nothing, so
        a front end can tell "undid something" from "there was nothing left".
        """
        if not self._undo_stack:
            return self._refuse('nothing to undo')
        self._restore(self._undo_stack.pop())
        self.last_reason = ''
        return self._accept(faces=self.mesh.number_of_faces(),
                            remaining=len(self._undo_stack))

    # ------------------------------------------------------------------
    # deleting a strip -- the template
    # ------------------------------------------------------------------

    def _split_strips(self, work, to_split):
        """Refine the strips that would otherwise let a boundary collapse.

        Thesis 5.3.2, Fig 5.18: a boundary collapses when fewer than three edges
        represent it after a deletion, and the remedy is to REFINE the strips that
        survive on it -- one left, split it in three; two left, split each in two,
        "to avoid any bias". ``strips_to_split_to_prevent_boundary_collapse``
        works that out; this performs it.

        **The strips this adds have ZERO WIDTH**, because the grammar's
        ``add_strip`` only ever does topology. A subclass MUST open them, and the
        two editors do it differently -- exact thirds on a coarse layout, centroid
        relaxation on a dense mesh -- which is why this is a hook and why the base
        deliberately does not pick one. Leaving them closed welds coincident
        vertices into the result, which survives ``is_manifold`` but bakes as a
        broken mesh.
        """
        return split_strips(work, to_split)

    def _gate(self, work):
        """**Is this result acceptable?** ``(ok, reason)``. Subclasses decide.

        Called on the COPY, before it is adopted. The base accepts anything, so a
        subclass that does not override gets the grammar's own guarantees and
        nothing more.
        """
        return True, ''

    def strip_through(self, edge):
        """The strip key through ``edge``, resolved on a copy. ``None`` if none.

        Collected on a COPY so the live mesh keeps no state the caller did not ask
        for, and re-collected every time because an earlier edit renumbered. The
        key is only valid on the copy it came from, which is why nothing here
        returns it for storage.
        """
        work = self.mesh.copy()
        work.collect_strips()
        return work.edge_strip(tuple(edge))

    def _trial_delete(self, edge=None, preserve_boundaries=False, skey=None):
        """Delete a strip on a COPY. ``(mesh or None, info)``. Touches no state.

        One code path for both :meth:`plan_strip_deletion` and
        :meth:`remove_strip`, so what the user is told and what then happens
        cannot drift apart.

        Addressed by ``edge`` or directly by ``skey`` -- a front end that lets the
        user pick the strip itself has a key already and should not have to invent
        an edge to pass it.
        """
        info = {'ok': False, 'reason': '', 'skey': None, 'faces': 0,
                'collateral': [], 'boundaries_lost': [], 'to_split': {},
                'can_preserve': True,
                'faces_before': self.mesh.number_of_faces(),
                'faces_after': None,
                'vertices_before': self.mesh.number_of_vertices(),
                'vertices_after': None}

        work = self.mesh.copy()
        work.collect_strips()

        if skey is None:
            if edge is None:
                info['reason'] = 'no edge or strip given'
                return None, info
            skey = work.edge_strip(tuple(edge))
        if skey is None or skey not in work.attributes['strips']:
            info['reason'] = 'that edge belongs to no strip -- nothing to remove'
            return None, info

        info['skey'] = skey
        info['faces'] = len(work.strip_faces(skey))
        info['collateral'] = list(collateral_strip_deletions(work, [skey]))
        # ``total_boundary_deletions`` returns the BOUNDARIES themselves -- a list
        # of vertex rings -- despite the name. Formatting it with ``{}`` prints a
        # list of rings where a number was meant, so callers take ``len``.
        info['boundaries_lost'] = list(total_boundary_deletions(work, [skey]))

        # **Thesis Eq 5.11, and this helper is the only thing that implements
        # it.** A boundary collapses when fewer than THREE edges represent it
        # after the deletion. ``total_boundary_deletions`` above asks a stricter
        # and different question -- whether EVERY one of its edges is consumed --
        # so it cannot stand in for this, and gating the pre-split on it meant
        # never pre-splitting in the cases the thesis is actually about.
        #
        # The helper is tri-state and all three states matter:
        #
        #   ``None``   the collapse cannot be prevented: no strip survives on that
        #              boundary to refine;
        #   ``{}``     nothing is at risk, so there is nothing to do;
        #   ``{...}``  refine these first -- one survivor, split it in three; two
        #              survivors, split each in two, "to avoid any bias" (5.3.2).
        #
        # Flattening ``None`` and ``{}`` together is what made asking to preserve
        # boundaries refuse a deletion that risked none.
        at_risk = strips_to_split_to_prevent_boundary_collapse(work, [skey])
        # Reported whether or not it was asked for: a front end wants to say "this
        # would cost you a hole, and here is what I would have to split to save
        # it" before the user decides.
        info['to_split'] = {} if at_risk is None else dict(at_risk)
        info['can_preserve'] = at_risk is not None

        if preserve_boundaries:
            if at_risk is None:
                info['reason'] = ('this strip cannot go without collapsing a '
                                  'boundary, and there is no strip left on that '
                                  'boundary to split to prevent it')
                return None, info
            if at_risk:
                self._split_strips(work, at_risk)
                # ``to_split`` never contains the strip being deleted, but the
                # split rebuilt the strip data, so re-resolve rather than trust
                # the key.
                if edge is not None:
                    skey = work.edge_strip(tuple(edge))
                    if skey is None:
                        info['reason'] = ('lost track of the strip after the '
                                          'pre-split')
                        return None, info
                    info['skey'] = skey

        try:
            _grammar_delete_strip(work, skey)
        except Exception as exc:
            # Measured: 1 strip of 24 on a two-hole plate raises ``KeyError`` from
            # inside the grammar. Catching it costs the user a refusal instead of
            # a traceback.
            info['reason'] = 'the grammar could not delete that strip ({}: {})'.format(
                type(exc).__name__, exc)
            return None, info

        if work.number_of_faces() == 0:
            info['reason'] = 'that would delete every face of the mesh'
            return None, info

        ok, reason = self._gate(work)
        if not ok:
            info['reason'] = reason
            return None, info

        info['ok'] = True
        info['faces_after'] = work.number_of_faces()
        info['vertices_after'] = work.number_of_vertices()
        return work, info

    def plan_strip_deletion(self, edge=None, preserve_boundaries=False, skey=None):
        """**What deleting the strip would cost.** No mutation.

        Deleting a strip is not a local edit, and a front end has to be able to
        say so before it happens. Three things reach past the picked strip, and
        all three are predicted here on a copy of the untouched mesh:

        * it MERGES the two sides of the strip -- for each disconnected part of
          the strip's edge network one new vertex is added at the centroid and
          every old vertex is substituted for it, so surrounding vertices move;
        * it can take other strips with it: any strip whose faces all lie inside
          the deleted one goes too;
        * it can COLLAPSE a boundary, which on a coarse layout means losing a hole
          the domain still has.

        **The prediction is the operation** -- see :meth:`_trial_delete`.

        Returns
        -------
        dict
            ``ok``, ``reason``, ``skey``, ``faces``, ``collateral`` (a LIST),
            ``boundaries_lost`` (a list of vertex rings), ``to_split``, and the
            before/after face and vertex counts. A subclass that publishes
            different names or shapes adapts this at its own boundary.
        """
        _work, info = self._trial_delete(
            edge, preserve_boundaries=preserve_boundaries, skey=skey)
        self.last_reason = '' if info['ok'] else info['reason']
        return info

    def remove_strip(self, edge=None, preserve_boundaries=False, skey=None):
        """**Delete the strip through** ``edge`` (or ``skey``). ``(ok, notes)``.

        Applied to a copy and adopted only if the copy survives :meth:`_gate`, so
        a refusal costs nothing and leaves the mesh exactly as it was.

        It COLLAPSES a band of faces and welds the band's two sides together; it
        does not *dissolve* the picked line and merge the faces either side of it.
        Dissolving is the exact inverse of a cut and no such operation exists in
        ``compas_singular``.
        """
        work, info = self._trial_delete(
            edge, preserve_boundaries=preserve_boundaries, skey=skey)
        if work is None:
            return self._refuse(info['reason'])

        self.mesh = work
        self.edited = True
        self.last_deletion = dict(info)
        return self._accept(**info)

    # ------------------------------------------------------------------
    # handing the result back to the caller's own object
    # ------------------------------------------------------------------

    def _transplant(self, source, carry_strip_data=False):
        """**Make** :attr:`target` **become** ``source``, keeping its identity.

        The caller holds a reference to the mesh it handed in, so a commit that
        rebound a new object would leave them looking at the pre-edit layout. No
        helper for this exists in compas: ``Mesh.clear()`` resets the topology and
        ``_max_vertex`` / ``_max_face`` but **does not touch** ``attributes``,
        ``__data__`` is a read-only property whose getter stringifies every key,
        and ``join`` adds rather than replaces. So it is done by hand.

        Vertex and face keys are carried across, so a caller holding keys keeps
        them **whenever the source was built by mutating this mesh**. They do not
        survive a weld-and-repair rebuild, which renumbers -- that is the caller's
        cue to re-key geometrically, and why the commit notes say which path ran.

        ``carry_strip_data`` decides the interesting half. The strip grammar
        maintains ``attributes['strips']`` incrementally and PRESERVES strip
        labels (thesis 5.3.3), so after an add or a delete the densities the user
        set are still about the right bands and travel across. A rebuild renumbers
        them, and carrying them there would silently apply the wrong numbers to
        the wrong strips.
        """
        target = self.target
        route = target.attributes.get('route')
        # A symmetric unit's seams, group and tolerances (``compas_singular.symmetry``)
        # describe the DOMAIN the layout lives in, like the route -- a rebuild that
        # renumbers the layout does not change them, and dropping them would leave
        # a unit that can no longer be expanded.
        symmetry = target.attributes.get('symmetry')

        target.clear()
        # ``clear()`` deliberately does not do this, and the omission is the trap:
        # a stale 'face_pole' or 'quad_mesh' would survive the transplant.
        target.attributes.clear()
        target.attributes.update(deepcopy(source.attributes))

        for vkey in source.vertices():
            target.add_vertex(key=vkey, attr_dict=dict(source.vertex[vkey]))
        for fkey in source.faces():
            target.add_face(source.face_vertices(fkey), fkey=fkey)
        for fkey, data in getattr(source, 'facedata', {}).items():
            target.facedata[fkey] = deepcopy(data)
        for edge, data in getattr(source, 'edgedata', {}).items():
            target.edgedata[edge] = deepcopy(data)

        # The dense mesh and the coarse-to-dense maps describe the layout as it
        # was BEFORE this edit. Keeping them is how a stale mesh gets served.
        for key in self.DERIVED_ATTRIBUTES:
            if key in target.attributes:
                target.attributes[key] = None if 'mesh' in key else {}
        if not carry_strip_data:
            for key in self.STRIP_KEYED_ATTRIBUTES:
                if key in target.attributes:
                    target.attributes[key] = {}
        # The route is a property of the layout, not of the rebuild that produced
        # it, so it comes from the object being written into rather than the copy.
        if route is not None:
            target.attributes['route'] = route
        if symmetry:
            target.attributes['symmetry'] = deepcopy(symmetry)
        return target
