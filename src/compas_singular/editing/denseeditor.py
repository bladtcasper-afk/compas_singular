"""**Hand-editing the FINAL dense quad mesh: the mesh, with no Rhino in it.**

:class:`DenseMeshEditor` is the counterpart of ``CoarseEditor`` one stage
further down the pipeline. It holds a dense quad mesh being edited and offers the
three things that can be done to one -- move a vertex, add a line, remove a line
-- plus the relaxation a newly added line needs to become visible. Nothing here
picks, prompts, draws or prints. A refusal is recorded in
:attr:`~DenseMeshEditor.last_reason` and returned as ``False``; whoever is
driving decides how to say so.

**What is shared with the coarse editor now lives in**
:class:`~compas_singular.editing.editor.MeshEditor`: holding the mesh and the
domain walls, moving one vertex, undo, and the strip-deletion walk. This class is
what only the DENSE stage needs -- ``add_line``, ``relax``, the polyedge lookup,
and the rule that poles disable the line operations.

**This class publishes bare booleans, and the base publishes** ``(ok, notes)``.
Deliberately: its callers -- ``CMD_edit_quad_mesh`` and the agent tools -- read
``True``/``False`` and take the reason from :attr:`last_reason`, and moving them
to tuples is a separate change from moving the code. The wrappers below are the
one place the two conventions meet.

**A LINE IS NEVER ONE EDGE**
----------------------------

A quad mesh cannot gain or lose a single edge: removing one merges two quads
into a hexagon, and adding one leaves a five-sided face. The only operations
that keep every face four-sided run the full width of the mesh -- which is what
a strip is. So the pick is an EDGE and the unit of change is the strip through
it, and the two operations are different things despite both being called
"adding and removing lines":

* :meth:`~DenseMeshEditor.add_line` takes the whole POLYEDGE through the edge
  and grows a strip beside it (``grammar.add_strip``). One new row of elements,
  wall to wall or all the way round;
* :meth:`~DenseMeshEditor.remove_line` deletes the whole STRIP through the edge
  (``grammar_pattern.delete_strip``) and welds the two sides together.

Both are planned on a COPY and adopted only if the copy is still all-quad, so a
refusal costs nothing and leaves the mesh exactly as it was.

**STRIP OPERATIONS NEED AN ALL-QUAD MESH, AND A MESH WITH POLES IS NOT ONE**
----------------------------------------------------------------------------

``collect_strip`` walks ``face_opposite_edge``, and the triangle fan around a
pole has none, so a strip through a pole is not defined. Rather than produce a
quietly wrong strip, :meth:`~DenseMeshEditor.strips_available` is ``False`` as
soon as the mesh has a non-quad face, and both line operations refuse.
:meth:`~DenseMeshEditor.move_vertex` still works, because it does not touch
topology.

**NOTHING IS KEYED BY INDEX ACROSS AN EDIT**
--------------------------------------------

Strip and polyedge keys come from popping ``self.edges()``, so they follow
vertex insertion order -- and every edit here renumbers. They are therefore
re-collected from the live mesh immediately before each use and never stored
between operations. That is the only safe way to use them at all, and it is why
every operation starts from an EDGE rather than from a remembered key. A caller
holding a selection across an edit must hold it geometrically.

**A NEW STRIP HAS ZERO WIDTH UNTIL THE MESH IS RELAXED**
--------------------------------------------------------

``add_strip`` creates its two vertices ON TOP of the vertex they replace, so
the added line is invisible and degenerate until something separates them.
:meth:`~DenseMeshEditor.add_line` therefore relaxes as part of the operation,
with ``old_to_new``'s vertices as the only boundary vertices allowed to move.
That is not cosmetic: the pair has to separate, and every other boundary vertex
has to stay where it is or the outline of the mesh moves.

**DELETING IS NOT LOCAL, AND THE CALLER HAS TO BE TOLD BEFORE IT HAPPENS**
--------------------------------------------------------------------------

Three things follow from what ``delete_strip`` does, and
:meth:`~DenseMeshEditor.plan_strip_deletion` predicts all three on the
untouched mesh so a driver can put them in front of a user, or an automated
caller can gate on them:

* it MERGES the two sides of the strip -- for each disconnected part of the
  strip's edge network one new vertex is added at the centroid and every old
  vertex is substituted for it, so surrounding vertices move;
* it can take other strips with it: any strip whose faces all lie inside the
  deleted one goes too;
* it can COLLAPSE a boundary that would be left with too few splits.

The remedy for the third is to pre-split the strips that would otherwise
collapse, which is what ``preserve_boundaries=True`` does. It goes through
``strips_to_split_to_prevent_boundary_collapse`` + ``split_strips`` rather than
through ``delete_strip``'s own ``preserve_boundaries`` flag, which calls a
helper -- ``boundary_strip_preserve`` -- that is not defined anywhere in the
package and raises ``NameError``.

**One correction came with the move to the base class.** ``preserve_boundaries``
used to be acted on whenever the flag was set, and the helper returns an EMPTY
dict when nothing needs splitting as well as ``None`` when a collapse cannot be
prevented. Treating those two the same way meant asking to preserve boundaries
REFUSED a deletion that would not have cost a boundary at all, with a message
saying one would collapse. It is now guarded on a boundary actually being at
risk. Nothing in the test suite exercised it; the agent's ``remove_line`` tool
is the one caller that can reach it.

And because no cheap test predicts every failure,
:meth:`~DenseMeshEditor.plan_strip_deletion` does not attempt one: it performs
the deletion on a copy and reports what came out, exactly as its counterpart on
the coarse layout does.

**THIS IS THE END OF THE PIPELINE, AND THAT IS DELIBERATE**
-----------------------------------------------------------

A hand-edited dense mesh has no coarse layout that produces it, so re-densifying
regenerates from the layout and throws the edit away. Nothing here writes back
to a layout and nothing tries to replay an edit onto a new one -- the
correspondence is not well defined once the layout changes. If an edit can be
expressed on the COARSE layout instead, do it there: ``CoarseEditor``
keeps the strip grammar, the field alignment and the boundary curvature, and a
mistake there costs one ``reset``.

That is also why this editor works **in place** rather than on a copy: it has no
``commit``, so a working copy would be a way to lose the edit rather than a way
to stage it. ``editor.mesh is the_mesh_you_passed`` holds here, and does not on
the coarse editor.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from compas.itertools import pairwise

# compas 2.x dropped ``mesh_smooth_centroid`` from the ``compas.datastructures``
# namespace, but not from the module it lives in. This is the import
# ``lizard.py`` already uses.
from compas.datastructures.mesh.smoothing import mesh_smooth_centroid

from ..datastructures.mesh_quad.grammar.add_strip import add_strip as _grammar_add_strip
from ..datastructures.mesh_quad.grammar.add_strip import is_polyedge_valid_for_strip_addition
from .editor import MeshEditor


__all__ = ['DenseMeshEditor', 'RELAX_ITERATIONS']


#: Smoothing passes used to open a newly added strip. ``add_strip`` creates its
#: two vertices on top of the one they replace, so without this the new strip
#: has zero width.
RELAX_ITERATIONS = 50


class DenseMeshEditor(MeshEditor):
    """A dense quad mesh being hand-edited. No Rhino, no prompts, no prints.

    Parameters
    ----------
    mesh : :class:`compas_singular.datastructures.QuadMesh`
        The mesh to edit. A ``QuadMesh``, not a plain
        :class:`compas.datastructures.Mesh`: the polyedges and strips are what
        the line operations select from. Edited **in place** by
        :meth:`move_vertex`, and REPLACED wholesale by :meth:`add_line` and
        :meth:`remove_line` -- so hold the editor, not the mesh.
    walls : list, optional
        The domain walls, as point lists, closed polylines or
        :class:`BoundaryLoop` instances. A boundary vertex that moves is
        projected onto the nearest one, so nudging the edge of the mesh does
        not eat the outline. Without walls a moved boundary vertex is left
        where it was put.
    relax_iterations : int, optional
        Smoothing passes used to open a newly added strip. Defaults to
        :data:`RELAX_ITERATIONS`.

    Attributes
    ----------
    mesh : :class:`compas_singular.datastructures.QuadMesh`
        The mesh as it currently stands. The object that was passed in.
    last_reason : str
        Why the last operation was refused. Empty when nothing was.
    last_addition : dict
        What the last successful :meth:`add_line` did.
    last_deletion : dict
        What the last successful :meth:`remove_line` did.

    Examples
    --------
    >>> editor = DenseMeshEditor(mesh, walls=loops)
    >>> plan = editor.plan_strip_deletion(edge)
    >>> plan['ok']
    True
    >>> editor.remove_line(edge)
    True
    """

    def __init__(self, mesh, walls=None, relax_iterations=None):
        if not hasattr(mesh, 'collect_polyedges'):
            raise TypeError(
                'DenseMeshEditor needs a QuadMesh -- the line operations select '
                'from its polyedges and strips, which a plain Mesh does not '
                'have. Got {}.'.format(type(mesh).__name__))
        # ``work_on_copy=False``: this stage has no commit to transplant a copy
        # back through, so a copy would silently strand the edit. See the module
        # docstring.
        super(DenseMeshEditor, self).__init__(mesh, walls=walls,
                                              work_on_copy=False)
        self.relax_iterations = (RELAX_ITERATIONS if relax_iterations is None
                                 else relax_iterations)
        self.last_addition = {}

    # ------------------------------------------------------------------
    # the two return conventions meet here
    # ------------------------------------------------------------------

    def _no(self, reason):
        """Record a refusal and return **this class's** value for one: ``False``.

        The base returns ``(False, notes)``. Every method below that publishes a
        bare boolean goes through here rather than through ``_refuse`` directly,
        so the two conventions cannot be mixed up by accident.
        """
        self._refuse(reason)
        return False

    def move_vertex(self, vkey, xyz, project=True):
        """Move one vertex. ``True`` if it moved.

        A vertex on the mesh boundary is projected back onto the nearest domain
        wall, so nudging the edge of the mesh does not eat the outline. An
        interior vertex is never projected, whatever its distance.

        Topology is unchanged, so this is the one operation that stays
        available on a mesh with poles.
        """
        ok, _notes = super(DenseMeshEditor, self).move_vertex(
            vkey, xyz, project=project)
        return ok

    # ------------------------------------------------------------------
    # queries -- none of these mutate anything
    # ------------------------------------------------------------------

    def _not_quads(self, mesh=None):
        mesh = self.mesh if mesh is None else mesh
        return [fkey for fkey in mesh.faces()
                if len(mesh.face_vertices(fkey)) != 4]

    def non_quad_faces(self):
        """Every face that is not four-sided. Empty on a clean quad mesh."""
        return self._not_quads()

    def strips_available(self):
        """Whether the line operations are defined on this mesh.

        ``False`` as soon as there is a non-quad face -- a pole fan, most
        likely. See the module docstring.
        """
        return not self._not_quads()

    def _pole_refusal(self):
        return self._no(
            '{} face(s) are not quads -- a strip through a pole is not defined, '
            'because collect_strip walks face_opposite_edge and a triangle has '
            'none. move_vertex still works.'.format(len(self._not_quads())))

    def polyedge_through(self, edge, mesh=None):
        """``(pkey, polyedge)`` for the polyedge containing ``edge``, as a new list.

        ``collect_polyedges`` is re-run every time rather than cached: an added
        or deleted strip renumbers the keys, and a remembered pkey would then
        name a different polyedge without erroring. The list is a COPY because
        ``add_strip`` pops from the list it is given, which would corrupt
        ``attributes['polyedges']`` if it were handed the stored one.
        """
        mesh = self.mesh if mesh is None else mesh
        u, v = edge
        for pkey, polyedge in dict(mesh.collect_polyedges()).items():
            for a, b in pairwise(polyedge):
                if (a, b) == (u, v) or (a, b) == (v, u):
                    return pkey, list(polyedge)
        return None, None

    # ------------------------------------------------------------------
    # the gate the base's strip-deletion template asks for
    # ------------------------------------------------------------------

    def _gate(self, work):
        """All-quad, and nothing more. The coarse editor also demands manifold."""
        bad = self._not_quads(work)
        if bad:
            return False, ('the result would have {} face(s) that are not '
                           'quads'.format(len(bad)))
        return True, ''

    # ------------------------------------------------------------------
    # add a line -- a strip along the polyedge through the picked edge
    # ------------------------------------------------------------------

    def add_line(self, edge, relax=True):
        """**Grow a strip along the polyedge through** ``edge``. ``True`` if grown.

        Planned on a copy and adopted only if the result is still all-quad, so
        a refusal costs nothing.

        Parameters
        ----------
        edge : tuple[int, int]
            Any edge of the polyedge to add a strip along.
        relax : bool, optional
            Separate the pair of vertices ``add_strip`` created on top of each
            other. On by default; the new strip has zero width without it. See
            the module docstring.
        """
        if not self.strips_available():
            return self._pole_refusal()

        pkey, polyedge = self.polyedge_through(edge)
        if polyedge is None:
            return self._no(
                'edge {} belongs to no polyedge -- nothing to add along'.format(
                    tuple(edge)))

        closed = polyedge[0] == polyedge[-1]
        # >2 vertices, and closed OR both ends on a boundary. A strip has to run
        # the full width of the mesh, for the same reason a cut on the coarse
        # layout has to reach a wall: anything less leaves a face with five
        # sides.
        if not is_polyedge_valid_for_strip_addition(self.mesh, polyedge):
            return self._no(
                'the line through edge {} has {} vertices, is not closed, and '
                'does not end on the boundary at both ends. A strip has to run '
                'the full width of the mesh -- wall to wall, or all the way '
                'round.'.format(tuple(edge), len(polyedge)))

        work = self.mesh.copy()
        # A PRECONDITION, not housekeeping: add_strip writes into
        # attributes['strips'] and fails without it.
        work.collect_strips()
        # Re-resolved on the copy: mesh.copy() preserves keys today, but relying
        # on that would make this break silently if it ever stopped being true.
        _pkey, work_polyedge = self.polyedge_through(edge, mesh=work)
        if work_polyedge is None:
            return self._no(
                'lost track of the line through edge {} on the working '
                'copy'.format(tuple(edge)))

        try:
            skey, old_to_new = _grammar_add_strip(work, list(work_polyedge))
        except Exception as exc:
            return self._no('add_strip failed: {}: {}'.format(
                type(exc).__name__, exc))

        ok, reason = self._gate(work)
        if not ok:
            return self._no(reason)

        # add_strip creates its two vertices ON TOP of the vertex they replace,
        # so the new strip has zero width and is invisible until the mesh is
        # relaxed. old_to_new maps old vkey -> the pair that has to separate,
        # and those are the only boundary vertices allowed to move.
        new_vertices = set()
        for pair in old_to_new.values():
            new_vertices.update(pair)
        if relax:
            self.relax(free=new_vertices, mesh=work)

        self.last_addition = {
            'strip': skey,
            'polyedge': pkey,
            'polyedge_vertices': len(work_polyedge),
            'closed': closed,
            'faces_before': self.mesh.number_of_faces(),
            'faces_after': work.number_of_faces(),
            'new_vertices': sorted(new_vertices),
            'relaxed': bool(relax),
        }
        self.mesh = work
        self.last_reason = ''
        return True

    # ------------------------------------------------------------------
    # remove a line -- the strip through the picked edge
    # ------------------------------------------------------------------

    @staticmethod
    def _as_dense_plan(info):
        """The base's canonical info dict in this class's published shape.

        Two renames and one reshape, and they are the whole difference: the base
        calls the key ``skey`` where this class has always published ``strip``,
        and it returns the collapsed boundaries THEMSELVES where this class
        publishes a count plus the rings under a separate name. Formatting a list
        of vertex rings with ``{}`` where a number was meant is exactly the bug
        the split name exists to prevent.
        """
        plan = dict(info)
        plan['strip'] = info['skey']
        plan['boundaries_lost'] = len(info['boundaries_lost'])
        plan['boundaries_lost_vertices'] = list(info['boundaries_lost'])
        return plan

    def plan_strip_deletion(self, edge, preserve_boundaries=False):
        """**What deleting the strip through** ``edge`` **would cost.** Mutates nothing.

        The three ways a deletion reaches further than the picked strip are
        predicted on the untouched mesh -- collateral strips, welded neighbours
        and collapsed boundaries -- and then the deletion is actually PERFORMED
        on a copy, because no cheap test predicts every failure. See the module
        docstring.

        Returns
        -------
        dict
            ``ok`` -- whether :meth:`remove_line` would succeed with these
            arguments -- and ``reason`` when it would not. Plus ``strip``,
            ``faces``, ``collateral``, ``boundaries_lost`` (a COUNT),
            ``boundaries_lost_vertices`` (which ones), ``to_split``,
            ``faces_before`` / ``faces_after`` and
            ``vertices_before`` / ``vertices_after``. The counts after are
            ``None`` when the deletion would not go through.
        """
        if not self.strips_available():
            info = {'ok': False, 'skey': None, 'faces': 0, 'collateral': [],
                    'boundaries_lost': [], 'to_split': [],
                    'faces_before': self.mesh.number_of_faces(),
                    'faces_after': None,
                    'vertices_before': self.mesh.number_of_vertices(),
                    'vertices_after': None,
                    'reason': '{} face(s) are not quads -- a strip through a '
                              'pole is not defined'.format(len(self._not_quads()))}
            return self._as_dense_plan(info)

        info = super(DenseMeshEditor, self).plan_strip_deletion(
            edge=edge, preserve_boundaries=preserve_boundaries)
        return self._as_dense_plan(info)

    def remove_line(self, edge, preserve_boundaries=False):
        """**Delete the strip through** ``edge``. ``True`` if deleted.

        Call :meth:`plan_strip_deletion` first if the consequences matter to the
        caller -- this method performs the same trial and simply adopts the
        result, so a refusal here costs nothing either.

        Parameters
        ----------
        edge : tuple[int, int]
            Any edge of the strip to delete.
        preserve_boundaries : bool, optional
            Pre-split the strips that would otherwise leave a boundary with too
            few splits. Only acts when a boundary is actually at risk, and
            refuses when one is and there is no strip left to split.
        """
        if not self.strips_available():
            return self._pole_refusal()

        faces_before = self.mesh.number_of_faces()
        vertices_before = self.mesh.number_of_vertices()
        ok, notes = super(DenseMeshEditor, self).remove_strip(
            edge=edge, preserve_boundaries=preserve_boundaries)
        if not ok:
            # ``remove_strip`` already recorded the reason; republish it as this
            # class's bare ``False``.
            return False

        plan = self._as_dense_plan(notes)
        self.last_deletion = {
            'strip': plan['strip'],
            'faces': plan['faces'],
            'collateral': plan['collateral'],
            'boundaries_lost': plan['boundaries_lost'],
            'pre_split': plan['to_split'] if preserve_boundaries else [],
            'faces_before': faces_before,
            'faces_after': self.mesh.number_of_faces(),
            'vertices_before': vertices_before,
            'vertices_after': self.mesh.number_of_vertices(),
        }
        return True

    # ------------------------------------------------------------------
    # relaxation
    # ------------------------------------------------------------------

    def relax(self, free=(), iterations=None, mesh=None):
        """Smooth the interior, holding the boundary. ``free`` vertices are released.

        ``vertices_on_boundaries()`` -- PLURAL, via
        :meth:`~compas_singular.editing.editor.MeshEditor.boundary_vertices`. The
        singular form returns only the longest boundary, so every hole would be
        left free and smoothing would round it off. ``free`` exists for the
        vertices a newly added strip has to separate: they sit on the boundary
        and must be allowed to move apart, and they are the only boundary
        vertices that may.

        Deliberately not ``framefield.relax.relax_mesh``: that one holds seams
        and the outline on their own curves and keeps a result only if it
        improves three metrics together. This one is the passes ``add_line``
        needs to open a strip.
        """
        mesh = self.mesh if mesh is None else mesh
        iterations = self.relax_iterations if iterations is None else iterations
        free = set(free)
        fixed = [vkey for vkey in self.boundary_vertices(mesh)
                 if vkey not in free]
        mesh_smooth_centroid(mesh, kmax=iterations, fixed=fixed)
        return mesh

    # ------------------------------------------------------------------
    # undo
    # ------------------------------------------------------------------

    def reset(self):
        """Put the mesh back as it was at construction, or at the last snapshot.

        Restores the whole MESH, not just the coordinates: the line operations
        change the topology, and putting coordinates back would leave the faces
        they created in place with nothing to say where they came from.
        """
        super(DenseMeshEditor, self).reset()
        self.last_addition = {}
