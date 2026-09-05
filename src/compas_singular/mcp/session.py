"""**One mesh being improved, and the record of how it got that way.**

Much smaller than ``agent/core/session.py``, and for a reason worth stating:
that class holds a three-state machine because it can generate a layout from a
field and a dense mesh from a layout, and every backward move destroys work. This
one improves a mesh that already exists. There is no stage to be in, so there is
no stage to gate, no ``StageError``, and no transition that costs anything.

**Undo is a map of positions, not a pickle of the world.** Smoothing moves
vertices and changes nothing else, so a snapshot is ``{vertex: [x, y, z]}`` --
kilobytes, instant, and no dependence on pickle being able to reach the object
graph. ``agent``'s session pickles the decomposition and both editors together
because it has to; this one does not have to, and the difference is one of the
things the two approaches are meant to be compared on.

**That cheapness is bought with an assumption, so it is checked.** A position map
only restores a mesh whose topology has not changed underneath it. Version 1 has
no tool that changes topology, so the assumption holds -- but :meth:`undo`
verifies the vertex and face counts anyway and refuses rather than silently
writing coordinates into the wrong mesh. When a topology tool is added, that
refusal is the thing that will say so.

**The history is the deliverable, not a side effect.** Every step records what
was asked, what changed, and the quality before and after. The model reads it
back through ``inspect`` and ``history``, a person reads it in the transcript,
and ``save_example`` turns it into a worked example other sessions can learn
from. Remarks -- prose the model or a person attaches -- ride along on it,
because a number without the reason it was accepted is not much of a record.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import time

from ..framefield.quality import mesh_quality


__all__ = ['MeshSession', 'UNDO_DEPTH']


#: How many snapshots to keep. Position maps are small, but a long session with
#: an undo per step would still grow without bound.
UNDO_DEPTH = 20


class MeshSession(object):
    """The mesh in hand, the geometry constraining it, and what has been done."""

    def __init__(self):
        #: The authoritative mesh, in double precision. Rhino only gets a copy.
        self.mesh = None
        #: Outer and inner boundaries as compas polylines. What the boundary
        #: vertices are allowed to slide along.
        self.walls = []
        #: Drawn guide curves, kept for reporting and for a later guide tool.
        self.guides = []
        #: Input point features, as [x, y, z]. Each one is a place the mesh is
        #: supposed to carry a pole; whether it does is a thing you see rather
        #: than measure, which is why they are drawn.
        self.points = []
        #: Where the mesh came from: ``{'kind': 'rhino'|'file', ...}``.
        self.source = None
        self.history = []
        self.remarks = []
        self.started = time.time()
        self._undo = []
        #: Monotonic counters behind :attr:`visually_current`. Numbers rather
        #: than a bool because "seen since the last change" is the question, and
        #: a bool would need resetting in every tool that moves a vertex.
        self._changed_at = 0
        self._seen_at = 0

    # --------------------------------------------------------------------
    # loading
    # --------------------------------------------------------------------

    def adopt(self, mesh, walls=None, guides=None, points=None, source=None):
        """Take a new mesh, discarding whatever was held before.

        Clears the undo stack: a snapshot of the previous mesh's positions
        cannot be applied to this one, and keeping it would only offer an undo
        that must refuse.
        """
        self.mesh = mesh
        if walls is not None:
            self.walls = _as_polylines(walls)
        if guides is not None:
            self.guides = _as_polylines(guides)
        if points is not None:
            self.points = [list(p) for p in points]
        self.source = source
        self._undo = []
        # A newly adopted mesh has never been looked at.
        self._changed_at = 1
        self._seen_at = 0
        return self

    @property
    def loaded(self):
        return self.mesh is not None

    # --------------------------------------------------------------------
    # measuring
    # --------------------------------------------------------------------

    def quality(self, low_angle=None):
        """``mesh_quality`` of the mesh in hand, or ``None`` if there is none."""
        if self.mesh is None:
            return None
        if low_angle is None:
            return mesh_quality(self.mesh)
        return mesh_quality(self.mesh, low_angle=low_angle)

    def triple(self):
        """``(min angle, max angle, max aspect)`` -- the three a step is judged on.

        The same three ``framefield.relax`` gates on, deliberately: a step that
        improves one and wrecks another has not improved the mesh.
        """
        metrics = self.quality()
        if metrics is None:
            return None
        return (metrics.get('min_angle'), metrics.get('max_angle'),
                metrics.get('aspect_max'))

    # --------------------------------------------------------------------
    # has anybody actually looked at it
    # --------------------------------------------------------------------

    def mark_changed(self):
        """Record that the mesh moved. Every tool that moves a vertex calls this."""
        self._changed_at += 1

    def mark_seen(self):
        """Record that the mesh has been rendered and looked at."""
        self._seen_at = self._changed_at

    @property
    def visually_current(self):
        """Whether the mesh has been LOOKED at since it last changed.

        The numbers can say a mesh improved while it is visibly wrong -- off its
        boundary, ignoring a guide, a point feature that never became a pole.
        None of that is in ``mesh_quality``, and all of it is obvious in a
        picture. So delivering a result is gated on somebody having seen one:
        :func:`rhino_push` and :func:`save_mesh` refuse until this is true.
        """
        return self._seen_at >= self._changed_at

    # --------------------------------------------------------------------
    # undo
    # --------------------------------------------------------------------

    def positions(self):
        return dict((key, list(self.mesh.vertex_coordinates(key)))
                    for key in self.mesh.vertices())

    def snapshot(self, label=''):
        """Remember where every vertex is, so a step can be taken back."""
        if self.mesh is None:
            return None
        self._undo.append({
            'label': label,
            'positions': self.positions(),
            'vertices': self.mesh.number_of_vertices(),
            'faces': self.mesh.number_of_faces(),
            'step': len(self.history),
        })
        while len(self._undo) > UNDO_DEPTH:
            self._undo.pop(0)
        return label

    def can_undo(self):
        return bool(self._undo)

    def undo(self):
        """Restore the most recent snapshot.

        Returns
        -------
        tuple
            ``(True, label)``, or ``(False, reason)``.
        """
        if self.mesh is None:
            return False, 'nothing is loaded'
        if not self._undo:
            return False, 'nothing to undo -- no snapshot has been taken'
        entry = self._undo.pop()
        if (entry['vertices'] != self.mesh.number_of_vertices()
                or entry['faces'] != self.mesh.number_of_faces()):
            # Not reachable in version 1, which has no topology tool. Kept
            # because the day one is added, this is the check that says the
            # position-map undo no longer covers it.
            return False, ('the mesh has {} vertices and {} faces, but the '
                           'snapshot was taken at {} and {} -- a position undo '
                           'cannot restore a changed topology'.format(
                               self.mesh.number_of_vertices(),
                               self.mesh.number_of_faces(),
                               entry['vertices'], entry['faces']))
        for key, xyz in entry['positions'].items():
            self.mesh.vertex_attributes(key, 'xyz', xyz)
        del self.history[entry['step']:]
        # An undo moves every vertex, so the last picture is out of date too.
        self.mark_changed()
        return True, entry['label'] or 'the last snapshot'

    # --------------------------------------------------------------------
    # the record
    # --------------------------------------------------------------------

    def record(self, action, before=None, after=None, **detail):
        """Append a step. Returns the entry, so a tool can report it back."""
        entry = {'step': len(self.history), 'action': action,
                 'at': time.time(), 'before': before, 'after': after}
        entry.update(detail)
        self.history.append(entry)
        return entry

    def add_remark(self, text, about=None):
        remark = {'step': len(self.history), 'text': text, 'about': about,
                  'at': time.time()}
        self.remarks.append(remark)
        return remark

    def remarks_for(self, step):
        return [r['text'] for r in self.remarks if r['step'] == step]

    def state(self):
        """A small summary. The reading in prose is ``describe``'s job."""
        mesh = self.mesh
        return {
            'loaded': self.loaded,
            'source': self.source,
            'vertices': mesh.number_of_vertices() if mesh else None,
            'faces': mesh.number_of_faces() if mesh else None,
            'walls': len(self.walls),
            'guides': len(self.guides),
            'points': len(self.points),
            'steps': len(self.history),
            'remarks': len(self.remarks),
            'undo_depth': len(self._undo),
            'visually_current': self.visually_current,
        }

    def __repr__(self):
        state = self.state()
        return '<MeshSession {} faces, {} step(s)>'.format(
            state['faces'], state['steps'])


def _as_polylines(curves):
    """Point lists as compas polylines, which is what the smoothers want.

    ``automated_boundary_constraints`` accepts either, but a ``Polyline`` gives
    ``closest_point_on_constraint`` a fast path and makes the type obvious at
    every later use.
    """
    from compas.geometry import Polyline
    out = []
    for curve in curves or []:
        if isinstance(curve, Polyline):
            out.append(curve)
            continue
        points = getattr(curve, 'points', curve)
        points = [list(p) for p in points]
        if len(points) >= 2:
            out.append(Polyline(points))
    return out
