"""**The pipeline as one session, with a state machine over it.**

Six stages, three states. The stages are things a caller DOES -- solve a field,
edit the layout, densify, set densities, smooth, edit the dense mesh -- but the
session is only ever in one of three states, because that is how many meshes
there are:

``'field'``
    A field is solved and its separatrices traced. There is no layout yet.
``'coarse'``
    A coarse layout exists and :class:`CoarseLayoutEditor` is holding it.
``'dense'``
    A dense mesh exists and :class:`DenseMeshEditor` is holding it.

**The states are not commutative, and every backward move destroys work.** That
is not a limitation to design around, it is what the pipeline is: each stage is
generated from the one above it, so going back means regenerating and the edit
below has nowhere to live.

* ``dense -> coarse`` throws away every dense edit. A hand-edited dense mesh has
  no layout that reproduces it -- see ``editing/dense_mesh.py``.
* ``coarse -> field`` throws away the layout, its edits, the densities and the
  dense mesh. A re-solve returns a NEW ``FieldDecomposition``; it does not
  mutate the one in hand.

So a backward move is never implicit. :meth:`~MeshEditSession.revert_to_coarse`
and :meth:`~MeshEditSession.adopt_decomposition` are the only two, both explicit,
and :meth:`~MeshEditSession.discarded_by` says what either would cost BEFORE it
is called -- which is what a confirmation prompt, human or otherwise, needs.

**The session never calls** :meth:`FieldDecomposition.quad_mesh`. It densifies
through :func:`~compas_singular.agent.core.rebuild.densify_preserving`, which
neither resets the densities nor falls back to the triangulation backstop. Those
two behaviours are correct for a person pressing a button and catastrophic for a
session holding twenty steps of state; ``rebuild.py`` says why at length.

**Undo is a whole-state snapshot, and it is not free.** The decomposition, both
editors and the density map are pickled together in ONE call, so the references
they share -- the editor's ``decomposition``, the layout it points at -- come
back shared rather than duplicated. Pickle is the mechanism ``framefield.cache``
already relies on for exactly this object, so it is proven on the type; it is
still a copy of the field, so snapshot deliberately rather than on every step.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import pickle

from ...editing.coarse_layout import CoarseLayoutEditor
from ...editing.dense_mesh import DenseMeshEditor
from .address import AddressBook
# ``DensifyRefused`` is deliberately NOT caught here: ``rebuild_dense`` lets it
# propagate, because the caller has to be told a density was rejected rather
# than discover it from an unchanged face count.
from .rebuild import densify_preserving


__all__ = ['STAGES', 'StageError', 'MeshEditSession']


#: In order. A stage is only ever entered from the one before it.
STAGES = ('field', 'coarse', 'dense')


class StageError(Exception):
    """An operation was asked for in a state where it is not defined."""


class MeshEditSession(object):
    """One meshing job being driven by something that is not a person.

    Parameters
    ----------
    decomposition : FieldDecomposition
        The solved field. From ``cache.solve(...)``, Rhino's
        ``get_decomposition()``, or ``FieldDecomposition.from_boundary(...)``.
    outer : list[[x, y, z]], optional
        The outer boundary the field was solved from. Kept so a re-solve has
        something to vary -- the decomposition does not carry its own inputs
        back. Without it :meth:`inputs` is incomplete and stage 0 cannot run.
    inners : list, optional
    guides : list, optional
        The guide curves -- cables, force lines -- the field was solved with.
        Defaults to ``decomposition.guides``, which IS carried.
    solve_params : dict, optional
        The keyword arguments the field was solved with, so a re-solve can vary
        one and keep the rest.
    target_length : float, optional
        Base element size for densification. Defaults to the background spacing,
        which is what ``quad_mesh`` defaults to. **This is the quad size, not
        the background triangulation spacing** -- they share a default and are
        not the same knob.

    Attributes
    ----------
    stage : str
        One of :data:`STAGES`.
    densities : dict
        ``{strip address: density}``. Keyed geometrically because strip keys
        renumber -- see ``address.py``.
    history : list[dict]
        Every operation recorded through :meth:`record`, in order.
    """

    def __init__(self, decomposition, outer=None, inners=None, guides=None,
                 solve_params=None, target_length=None):
        if decomposition is None:
            raise ValueError('a session needs a solved FieldDecomposition')
        self.decomposition = decomposition
        self.outer = outer
        self.inners = list(inners or [])
        self.guides = list(guides if guides is not None
                           else (decomposition.guides or []))
        self.solve_params = dict(solve_params or {})
        self.target_length = (target_length
                              or decomposition.background.target_length)

        self.stage = 'field'
        self.coarse_editor = None
        self.dense_editor = None
        self.densities = {}
        self.history = []

        self._undo = []
        self._book = None
        self._last_carry = {}
        #: ``(chain, guide)`` from the last ``select_guide_chain``, waiting for
        #: an ``attach_guide_chain``. Selection and attachment are two steps on
        #: purpose: attaching MOVES vertices, so a second guide chosen after the
        #: first had already pulled the mesh about would depend on the order the
        #: guides were processed in. Select, measure, then attach.
        self._pending_chain = None
        #: Set by :meth:`adopt_dense`. A dense mesh handed in from outside has
        #: no layout that reproduces it, so the backward transition is refused.
        self._dense_adopted = False
        #: The most recent dense mesh produced or adopted, kept across a revert.
        #: See :meth:`rebuild_dense`.
        self.last_dense = None

    # ------------------------------------------------------------------
    # what the session is looking at
    # ------------------------------------------------------------------

    @property
    def mesh(self):
        """The mesh of the current stage, or ``None`` at ``'field'``."""
        if self.stage == 'dense' and self.dense_editor is not None:
            return self.dense_editor.mesh
        if self.stage == 'coarse' and self.coarse_editor is not None:
            return self.coarse_editor.mesh
        return None

    @property
    def editor(self):
        """The editor of the current stage, or ``None`` at ``'field'``."""
        if self.stage == 'dense':
            return self.dense_editor
        if self.stage == 'coarse':
            return self.coarse_editor
        return None

    def inputs(self):
        """What a re-solve would need. ``complete`` is False when it cannot run."""
        return {
            'outer': self.outer,
            'inners': self.inners,
            'guides': self.guides,
            'params': dict(self.solve_params),
            'complete': self.outer is not None,
        }

    # ------------------------------------------------------------------
    # addressing
    # ------------------------------------------------------------------

    def book(self):
        """The :class:`AddressBook` for the current mesh, built on demand."""
        mesh = self.mesh
        if mesh is None:
            return None
        if self._book is None:
            self._book = AddressBook(mesh)
        return self._book

    def refresh_book(self, tol=None):
        """Rebuild the labels for the current mesh. Returns the carry report.

        Call after anything that renumbers -- a cut, a strip deletion, a commit.
        The report says which of the OLD labels survived, which moved and which
        are gone, so a caller holding a plan across the change is told rather
        than left pointing an old name at a new corner. See
        :meth:`AddressBook.carry`.
        """
        mesh = self.mesh
        if mesh is None:
            self._book = None
            self._last_carry = {}
            return {}
        if self._book is None:
            self._book = AddressBook(mesh)
            self._last_carry = {}
            return {}
        self._book, report = self._book.carry(mesh, tol=tol)
        self._last_carry = report
        return report

    @property
    def last_carry(self):
        """The report from the most recent :meth:`refresh_book`."""
        return dict(self._last_carry)

    # ------------------------------------------------------------------
    # forward transitions
    # ------------------------------------------------------------------

    def enter_coarse(self, loops=None, snap_tol=None):
        """Move to ``'coarse'``. Builds the layout if there is not one yet.

        Idempotent: calling it again while already at ``'coarse'`` returns the
        editor that is already there rather than rebuilding, so an edit in
        progress is not silently discarded.
        """
        if self.stage == 'dense':
            raise StageError(
                'already at the dense stage -- use revert_to_coarse(), which '
                'says what it discards, rather than entering coarse again')
        if self.coarse_editor is not None:
            self.stage = 'coarse'
            return self.coarse_editor

        if self.decomposition.mesh is None:
            # Note this clears ``repair_notes`` -- decomposition_mesh() does, and
            # that is the documented behaviour, not an accident here.
            self.decomposition.decomposition_mesh()
        self.coarse_editor = CoarseLayoutEditor(
            self.decomposition, loops=loops, snap_tol=snap_tol)
        self.stage = 'coarse'
        self._book = None
        self.record('enter_coarse',
                    faces=self.coarse_editor.mesh.number_of_faces())
        return self.coarse_editor

    def rebuild_dense(self, densities=None, target_length=None, density=None,
                      walls=None):
        """**Densify the layout and move to** ``'dense'``.

        Goes through :func:`densify_preserving`, so per-strip densities survive
        and a failure raises :class:`DensifyRefused` **without replacing the
        layout**. On a refusal the session stays exactly where it was and the
        coarse edit is intact -- which is the whole reason this method is not
        ``quad_mesh``.

        Parameters
        ----------
        densities : dict, optional
            ``{strip address: density}``. Defaults to :attr:`densities`, so
            per-strip work set earlier is re-applied automatically.
        target_length : float, optional
            Base element size. Defaults to :attr:`target_length`.
        density : int, optional
            Uniform base, takes precedence over ``target_length``.
        walls : list, optional
            Domain walls for the dense editor, so a moved boundary vertex is
            held on the outline. Defaults to the decomposition's own loops.

        Returns
        -------
        (QuadMesh, dict)
            The dense mesh and :func:`densify_preserving`'s report.
        """
        if self.stage == 'field':
            self.enter_coarse()

        if densities is None:
            densities = self.densities
        if target_length is None and density is None:
            target_length = self.target_length

        dense, report = densify_preserving(
            self.decomposition, densities=densities,
            target_length=target_length, density=density,
            coarse=self.coarse_editor.mesh if self.coarse_editor else None)

        # Only now, once it worked. A refusal must leave the map alone: the
        # densities the caller set are still what they want, and overwriting
        # them with a failed attempt's would lose the work twice.
        self.densities = dict(report['densities'])
        if walls is None:
            walls = self.decomposition._loops()
        self.dense_editor = DenseMeshEditor(dense, walls=walls)
        #: The most recent dense mesh this session produced, kept even after a
        #: revert. ``revert_to_coarse`` clears ``decomposition.dense``, so
        #: without this a run that densifies well and then steps back to the
        #: layout at the end leaves its caller with nothing to show -- measured
        #: on a real run that spent 37 tool calls on densities and then reverted
        #: on its last turn. What to DELIVER and what stage the session is IN
        #: are different questions.
        self.last_dense = dense
        self.stage = 'dense'
        self._book = None
        self.record('rebuild_dense', faces=report['faces'],
                    route=report['route'],
                    unmatched=len(report['density_report']['unmatched']))
        return dense, report

    def adopt_dense(self, mesh, walls=None):
        """**Edit a dense mesh that already exists**, rather than densifying one.

        For a front end holding a mesh from somewhere else -- baked in a
        document, loaded from a file, produced by an earlier run. The session
        goes straight to ``'dense'`` and the layout is not consulted.

        **The mesh has to be a QuadMesh**, not a plain ``Mesh``: the line
        operations select from its polyedges and strips.

        **There is no going back from here.** A mesh adopted this way has no
        coarse layout that reproduces it, so :meth:`revert_to_coarse` refuses
        unless a layout was separately built -- see :meth:`discarded_by`. That
        is the same rule ``CMD_edit_quad_mesh`` states: re-densifying regenerates
        from the layout and throws the edit away, and no correspondence exists
        once the layout changes.
        """
        if walls is None:
            walls = self.decomposition._loops()
        self.dense_editor = DenseMeshEditor(mesh, walls=walls)
        self.decomposition.dense = mesh
        self.last_dense = mesh
        self.stage = 'dense'
        self._book = None
        #: Whether the dense mesh was handed in rather than densified from the
        #: layout in hand. Read by :meth:`discarded_by`.
        self._dense_adopted = True
        self.record('adopt_dense', faces=mesh.number_of_faces(),
                    vertices=mesh.number_of_vertices())
        return self.dense_editor

    # ------------------------------------------------------------------
    # backward transitions -- both explicit, both priced first
    # ------------------------------------------------------------------

    def discarded_by(self, transition):
        """What ``transition`` would destroy. Call this before calling that.

        Parameters
        ----------
        transition : {'revert_to_coarse', 'adopt_decomposition'}

        Returns
        -------
        dict
            ``allowed`` and, when it is not, ``reason``. Otherwise the counts
            that make the cost concrete: ``dense_faces``, ``dense_edits``,
            ``coarse_faces``, ``coarse_edits``, ``densities``, and ``final``
            -- True when nothing can bring the work back.
        """
        if transition == 'revert_to_coarse':
            if self.stage != 'dense':
                return {'allowed': False,
                        'reason': 'not at the dense stage -- nothing to revert'}
            if self.coarse_editor is None:
                # Do NOT suggest enter_coarse as a way out: it refuses at the
                # dense stage too, so the advice would send a caller round a
                # loop. A dense mesh adopted from outside is genuinely the end
                # of the line, and saying so plainly is the useful answer.
                return {'allowed': False,
                        'reason': 'this dense mesh was adopted rather than '
                                  'densified from a layout, so there is no '
                                  'layout to go back to and no way to make one '
                                  'that would reproduce it. Edit the mesh where '
                                  'it is.'}
            return {
                'allowed': True,
                'final': True,
                'dense_faces': (self.dense_editor.mesh.number_of_faces()
                                if self.dense_editor else 0),
                'dense_edits': self._count('dense'),
                'coarse_faces': None, 'coarse_edits': None,
                'densities': 0,
            }

        if transition == 'adopt_decomposition':
            return {
                'allowed': True,
                'final': True,
                'dense_faces': (self.dense_editor.mesh.number_of_faces()
                                if self.dense_editor else 0),
                'dense_edits': self._count('dense'),
                'coarse_faces': (self.coarse_editor.mesh.number_of_faces()
                                 if self.coarse_editor else 0),
                'coarse_edits': self._count('coarse'),
                'densities': len(self.densities),
            }

        return {'allowed': False,
                'reason': 'unknown transition {!r}'.format(transition)}

    def revert_to_coarse(self):
        """Drop the dense mesh and go back to the layout. **Destroys dense edits.**

        Returns the ``discarded_by`` report for what it just did, so a caller
        can say what was lost after the fact as well as before it.
        """
        cost = self.discarded_by('revert_to_coarse')
        if not cost['allowed']:
            raise StageError(cost['reason'])
        self.dense_editor = None
        self.decomposition.dense = None
        self.stage = 'coarse'
        self._book = None
        self.record('revert_to_coarse', **{k: v for k, v in cost.items()
                                           if k != 'allowed'})
        return cost

    def adopt_decomposition(self, decomposition, guides=None,
                            solve_params=None):
        """**Take a newly solved field.** Destroys the layout and everything under it.

        The hook stage 0 uses. A re-solve returns a NEW ``FieldDecomposition``
        rather than mutating this one, so there is nothing to merge: the layout,
        its edits, the densities and the dense mesh all belonged to the old
        field and none of them describe the new one.

        Returns the ``discarded_by`` report.
        """
        cost = self.discarded_by('adopt_decomposition')
        self.decomposition = decomposition
        if guides is not None:
            self.guides = list(guides)
        elif decomposition.guides:
            self.guides = list(decomposition.guides)
        if solve_params is not None:
            self.solve_params = dict(solve_params)
        self.coarse_editor = None
        self.dense_editor = None
        self.densities = {}
        self.stage = 'field'
        self._book = None
        self._last_carry = {}
        self.record('adopt_decomposition', **{k: v for k, v in cost.items()
                                              if k != 'allowed'})
        return cost

    # ------------------------------------------------------------------
    # undo
    # ------------------------------------------------------------------

    def snapshot(self, label=''):
        """Push the whole session state onto the undo stack.

        One ``pickle.dumps`` over the decomposition, both editors, the density
        map and the stage, so shared references survive -- see the module
        docstring. Not free; take one before something risky, not every step.
        """
        blob = pickle.dumps((
            self.decomposition, self.coarse_editor, self.dense_editor,
            self.densities, self.stage, self.target_length, self.guides,
        ), protocol=pickle.HIGHEST_PROTOCOL)
        self._undo.append((label, blob))
        return len(self._undo)

    def can_undo(self):
        return bool(self._undo)

    def undo(self):
        """Pop the last snapshot back. ``(True, label)``, or ``(False, reason)``."""
        if not self._undo:
            return False, 'nothing to undo'
        label, blob = self._undo.pop()
        (self.decomposition, self.coarse_editor, self.dense_editor,
         self.densities, self.stage, self.target_length,
         self.guides) = pickle.loads(blob)
        self._book = None
        self._last_carry = {}
        self.record('undo', label=label, depth=len(self._undo))
        return True, label

    # ------------------------------------------------------------------
    # history
    # ------------------------------------------------------------------

    def record(self, action, **detail):
        """Append one step to :attr:`history`. Returns the entry."""
        entry = {'step': len(self.history), 'stage': self.stage,
                 'action': action}
        entry.update(detail)
        self.history.append(entry)
        return entry

    def _count(self, stage):
        """How many recorded steps were edits at ``stage``."""
        skip = ('enter_coarse', 'rebuild_dense', 'revert_to_coarse',
                'adopt_decomposition', 'undo', 'snapshot')
        return sum(1 for e in self.history
                   if e.get('stage') == stage and e.get('action') not in skip)

    # ------------------------------------------------------------------

    def state(self):
        """A small summary. The full per-stage digest is ``digest.py``'s job."""
        mesh = self.mesh
        return {
            'stage': self.stage,
            'faces': mesh.number_of_faces() if mesh is not None else None,
            'vertices': mesh.number_of_vertices() if mesh is not None else None,
            'densities': len(self.densities),
            'target_length': self.target_length,
            'guides': len(self.guides),
            'steps': len(self.history),
            'undo_depth': len(self._undo),
            'inputs_complete': self.outer is not None,
            'last_reason': getattr(self.editor, 'last_reason', '') or '',
        }

    def __repr__(self):
        s = self.state()
        return '<MeshEditSession {} stage, {} faces, {} step(s)>'.format(
            s['stage'], s['faces'], s['steps'])
