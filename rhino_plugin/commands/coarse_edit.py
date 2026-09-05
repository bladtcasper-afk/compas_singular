"""**Step 5 in Rhino: move the coarse layout's corners, then hand it to
``edit_coarse``.**

The layout is edited as a MESH held in Python, not as the polylines baked on
``ff::coarse``. Both are legal inputs to
:meth:`FieldDecomposition.edit_coarse` -- ``edit.faces_from_geometry`` accepts
either -- but only one of them is safe to drag a corner in:

* **polylines duplicate every shared corner.** ``edit.face_polylines`` writes one
  closed polyline per patch, so an interior corner where four patches meet is
  four coincident points. Moving one of them and leaving the other three is a
  tear, and the weld in ``edit.mesh_from_faces`` rounds at 3 decimals, so the
  moved copy silently becomes a NEW vertex: the layout comes back with a slit in
  it and nothing reports an error.
* **the mesh shares them.** ``vertex_attributes(vkey, 'xyz', ...)`` moves the one
  corner that all four patches reference, which is what the user meant.

It also keeps the strip-index trap in ``RHINO_PLUGIN.md`` §4.2 out of the edit
ITSELF -- indices permute because a bake/read-back rebuilds the topology from
geometry, and the mesh here never leaves Python, so keys are stable for as long
as corners are being dragged. **They are not stable across**
:meth:`CoarseLayoutObject.commit`: ``edit_coarse`` welds and repairs, which is
the same rebuild, and it was measured to leave 2 of 8 vertex keys pointing at
the same corner on a disc. So a density map held across a commit must be keyed
geometrically -- rounded edge midpoints, or one remembered ``(gkey(u), gkey(v))``
-- exactly as §4.2 says.

The polylines keep their job as the bake/export format and as the seam for
TOPOLOGICAL edits -- deleting a patch, splitting one -- which are awkward on a
live mesh and trivial on polylines. Both routes end in the same call.

**The edit is not finished until :meth:`CoarseLayoutObject.commit` runs.**
Moving vertices in place produces a layout whose edges no longer match any
traced separatrix by geometric key. Densifying THAT directly gives straight
chords and throws away the field alignment the whole front end exists to
produce. ``edit_coarse`` re-matches each edge to the separatrix it came from and
warps that curve onto the moved corners (``edit.warp_polyline``), so the nudge
costs the nudge and not the curvature. It also snaps boundary corners onto the
domain walls, which is load-bearing rather than cosmetic -- see ``edit.py``.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


try:
    import rhinoscriptsyntax as rs
    import scriptcontext as sc
    import Rhino
    from Rhino.Geometry import Line
    from Rhino.Geometry import Point3d
    # ``from System.Drawing.Color import FromArgb`` is an IronPython idiom and
    # raises under Rhino 8's CPython, where Color is a class and not a module.
    # It raised at the END of this block, so everything above it was bound and
    # the failure only surfaced later as ``NameError: FromArgb``.
    from System.Drawing import Color
except ImportError:
    pass


__all__ = ['CoarseLayoutObject']


def _closest_on_loop():
    """``compas_singular.framefield.edit.closest_on_loop``, imported late.

    The import is deliberately kept inside the function. Rhino's script engine
    re-initialises ``sys.path`` between runs while keeping ``sys.modules`` warm,
    so a late import cannot assume whatever put the path there earlier is still
    in effect -- see ``ff_common.ensure_paths``.
    """
    try:
        from compas_singular.framefield.edit import closest_on_loop
    except ImportError:
        raise ImportError(
            "cannot import compas_singular.framefield.edit -- is compas_singular/src "
            "on sys.path? See ff_common.ensure_paths / CMD_start.")
    return closest_on_loop


class CoarseLayoutObject(object):
    """Scene object for a coarse quad layout being hand-edited in Rhino.

    Parameters
    ----------
    decomposition : :class:`FieldDecomposition`
        The one that produced the layout. Kept because :meth:`commit` hands the
        edited mesh back to its :meth:`~FieldDecomposition.edit_coarse`, which
        needs the separatrices and the domain loops that only it has.
    coarse : :class:`CoarsePseudoQuadMesh`, optional
        The layout to edit. Defaults to ``decomposition.mesh``.

    Examples
    --------
    >>> decomposition = FieldDecomposition.from_boundary(outer, inners)
    >>> decomposition.decomposition_mesh()
    >>> obj = CoarseLayoutObject(decomposition)
    >>> obj.update()                    # move corners, then 'commit'
    >>> mesh = decomposition.quad_mesh()
    """

    settings = {
        'layer': 'ff::coarse',
        'color.vertices': [0, 0, 0],
        'color.vertices.boundary': [0, 120, 200],
        'color.vertices.pole': [200, 0, 120],
        'color.edges': [80, 80, 80],
        'color.preview': [255, 120, 0],
    }

    def __init__(self, decomposition, coarse=None, loops=None, snap_tol=None):
        self.decomposition = decomposition
        self.mesh = coarse if coarse is not None else decomposition.mesh
        if self.mesh is None:
            raise ValueError('no coarse layout -- call decomposition_mesh() first')
        # ``_loops`` is private but it is the only accessor for outer + inners,
        # and ``edit_coarse`` uses the same one internally.
        self.loops = loops if loops is not None else decomposition._loops()
        self.snap_tol = snap_tol
        self._guid_vertices = {}
        self._guid_edges = {}
        self._original = {vkey: self.mesh.vertex_coordinates(vkey)
                          for vkey in self.mesh.vertices()}

    # ------------------------------------------------------------------
    # layer state
    # ------------------------------------------------------------------

    def _layer_path(self):
        """The layer and every parent of it -- a locked parent locks the child."""
        parts = self.settings['layer'].split('::')
        return ['::'.join(parts[:i + 1]) for i in range(len(parts))]

    def unlock(self):
        """Unlock the layout's layer, its parents, and the objects on it.

        Returns what was locked, to hand back to :meth:`relock`. Layer lock and
        object lock are separate in Rhino and either one alone is enough to make
        the corners unpickable, so both are cleared.
        """
        state = {'layers': {}, 'objects': []}
        for layer in self._layer_path():
            if rs.IsLayer(layer):
                state['layers'][layer] = rs.LayerLocked(layer)
                if state['layers'][layer]:
                    rs.LayerLocked(layer, False)
        for guid in list(self._guid_vertices) + list(self._guid_edges):
            if rs.IsObjectLocked(guid):
                state['objects'].append(guid)
                rs.UnlockObject(guid)
        return state

    def relock(self, state):
        """Put back exactly what :meth:`unlock` took off, and nothing else."""
        for guid in state.get('objects', []):
            if rs.IsObject(guid):
                rs.LockObject(guid)
        for layer, was_locked in state.get('layers', {}).items():
            if was_locked and rs.IsLayer(layer):
                rs.LayerLocked(layer, True)

    # ------------------------------------------------------------------
    # the edit
    # ------------------------------------------------------------------

    def update(self):
        """Command loop. ``move``, ``reset``, ``commit``, ``finish``.

        Everything runs with the layer unlocked and re-locked on the way out,
        including on an exception -- a command that leaves the user's layers in
        a state they did not set them to is worse than one that fails.
        """
        self.draw()
        state = self.unlock()
        try:
            while True:
                operation = rs.GetString(
                    'next', 'finish', ['move', 'reset', 'commit', 'finish'])
                if operation == 'move':
                    self.move_vertices()
                elif operation == 'reset':
                    self.reset()
                    self.draw()
                elif operation == 'commit':
                    if self.commit():
                        break
                else:
                    break
        finally:
            self.relock(state)

    def move_vertices(self):
        """Pick a corner, drag it, repeat until Enter."""
        while True:
            guid = rs.GetObject('Select a coarse corner', preselect=True,
                                filter=rs.filter.point)
            if not guid:
                return
            vkey = self._guid_vertices.get(guid)
            if vkey is None:
                print('Not a corner of the coarse layout.')
                continue
            if self._move(vkey):
                self.draw()

    def _move(self, vkey):
        """Drag one corner, previewing the layout as it deforms.

        The preview shows the SNAPPED position, not the cursor: a boundary
        corner is projected onto its domain wall as it moves, so what is drawn
        is what ``edit_coarse`` will get. The point is constrained to the world
        XY plane because the pipeline is planar -- ``edit.mesh_from_faces``
        zeroes z anyway, and letting the cursor leave the plane only makes the
        preview lie.
        """
        boundary = vkey in set(self.mesh.vertices_on_boundary())
        neighbours = [self.mesh.vertex_coordinates(nbr)
                      for nbr in self.mesh.vertex_neighbors(vkey)]
        start = self.mesh.vertex_coordinates(vkey)
        color = Color.FromArgb(*self.settings['color.preview'])

        gp = Rhino.Input.Custom.GetPoint()
        gp.SetCommandPrompt('New position for corner {}'.format(vkey))
        gp.SetBasePoint(Point3d(*start), True)
        gp.Constrain(Rhino.Geometry.Plane.WorldXY, False)

        def OnDynamicDraw(sender, e):
            cp = e.CurrentPoint
            xyz = self._project([cp.X, cp.Y, 0.0], boundary)
            point = Point3d(*xyz)
            for nbr in neighbours:
                e.Display.DrawLine(Line(point, Point3d(*nbr)), color, 2)
            e.Display.DrawPoint(point, color)

        gp.DynamicDraw += OnDynamicDraw
        gp.Get()
        if gp.CommandResult() != Rhino.Commands.Result.Success:
            return False

        ep = gp.Point()
        self.mesh.vertex_attributes(
            vkey, 'xyz', self._project([ep.X, ep.Y, 0.0], boundary))
        return True

    def _project(self, xyz, boundary):
        """A boundary corner goes onto the nearest wall; an interior one does not.

        Eligibility is TOPOLOGICAL, exactly as in ``edit.snap_to_loops``, and for
        the same reason: projecting anything merely near a wall would drag an
        interior corner that legitimately sits in a narrow slot onto it and
        collapse a patch nobody touched. Doing it here as well as in
        ``edit_coarse`` is not redundant -- it is what makes the preview honest.
        """
        if not boundary or not self.loops:
            return xyz
        closest_on_loop = _closest_on_loop()
        best_d, best_p = float('inf'), None
        for loop in self.loops:
            d, q = closest_on_loop(xyz, loop)
            if d < best_d:
                best_d, best_p = d, q
        return best_p if best_p is not None else xyz

    def reset(self):
        """Put every corner back where the generated layout had it."""
        for vkey, xyz in self._original.items():
            if vkey in list(self.mesh.vertices()):
                self.mesh.vertex_attributes(vkey, 'xyz', xyz)

    # ------------------------------------------------------------------
    # handing it back
    # ------------------------------------------------------------------

    def commit(self):
        """**Hand the edited layout to** :meth:`FieldDecomposition.edit_coarse`.

        This is the step that makes the edit worth anything: the moved edges are
        re-matched to the separatrices they were traced from and those curves are
        warped onto the new corners, so the layout keeps its field alignment.
        Skipping it and densifying the mutated mesh directly gives straight
        chords between the moved corners.

        ``edit_coarse`` deliberately RAISES on a layout that will not densify
        rather than falling back, so the failure is caught here and reported
        while the user is still in Rhino and can fix it. Returns whether the
        layout was adopted; on failure the mesh is left as edited so the offending
        corner can be moved again.
        """
        kwargs = {} if self.snap_tol is None else {'snap_tol': self.snap_tol}
        try:
            mesh = self.decomposition.edit_coarse(self.mesh, **kwargs)
        except ValueError as e:
            print('The edited layout was NOT adopted:')
            print('  {}'.format(e))
            return False

        # ``edit_coarse`` welds and repairs, so it returns a DIFFERENT mesh
        # object whose vertex keys are RENUMBERED -- measured 2 of 8 unchanged on
        # a disc. Rebind, or the next move edits an orphan the decomposition no
        # longer knows about, and re-key anything held by index (see the module
        # docstring).
        self.mesh = mesh
        self._original = {vkey: mesh.vertex_coordinates(vkey)
                          for vkey in mesh.vertices()}
        notes = self.decomposition.edit_notes
        print('Layout adopted: {} patch(es) in, {} out, {} corner(s) snapped'.format(
            notes.get('faces_in'), notes.get('faces_out'), notes.get('snapped')))
        if notes.get('off_wall'):
            print('  {} corner(s) sit on no wall -- expected only where a patch '
                  'was deleted'.format(notes['off_wall']))
        # 3 is a pseudo-quad pole and is correct; anything else means a patch
        # went out with more points than it has corners, which the repair then
        # fanned into patches nobody drew.
        odd = {n: count for n, count in notes.get('sides', {}).items()
               if n not in (3, 4)}
        if odd:
            print('  patches with {} side(s) -- one point per CORNER, not per '
                  'sample'.format(sorted(odd)))
        self.draw()
        return True

    # ------------------------------------------------------------------
    # display
    # ------------------------------------------------------------------

    def draw(self):
        """Corners as points, edges as lines, on ``ff::coarse``.

        Redraws everything: a coarse layout is tens of patches, so there is
        nothing to gain from updating in place and a stale guid map is a real
        way to move the wrong vertex.
        """
        self.clear()
        layer = self.settings['layer']
        if not rs.IsLayer(layer):
            rs.AddLayer(layer)

        boundary = set(self.mesh.vertices_on_boundary())
        poles = set(self.mesh.poles()) if hasattr(self.mesh, 'poles') else set()

        rs.EnableRedraw(False)
        for vkey in self.mesh.vertices():
            guid = rs.AddPoint(Point3d(*self.mesh.vertex_coordinates(vkey)))
            rs.ObjectLayer(guid, layer)
            if vkey in poles:
                key = 'color.vertices.pole'
            elif vkey in boundary:
                key = 'color.vertices.boundary'
            else:
                key = 'color.vertices'
            rs.ObjectColor(guid, self.settings[key])
            self._guid_vertices[guid] = vkey

        color = self.settings['color.edges']
        for u, v in self.mesh.edges():
            guid = rs.AddLine(Point3d(*self.mesh.vertex_coordinates(u)),
                              Point3d(*self.mesh.vertex_coordinates(v)))
            rs.ObjectLayer(guid, layer)
            rs.ObjectColor(guid, color)
            self._guid_edges[guid] = (u, v)
        rs.EnableRedraw(True)
        sc.doc.Views.Redraw()

    def clear(self):
        guids = list(self._guid_vertices) + list(self._guid_edges)
        if guids:
            rs.DeleteObjects(guids)
        self._guid_vertices = {}
        self._guid_edges = {}
