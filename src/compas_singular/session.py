"""**A whole project in one object: settings, inputs, coarse layout, field and dense mesh.**

Runs without Rhino, and saves as one JSON file::

    session = SingularSession.load('floor.json')
    editor = CoarseEditor(session.coarse, field=session.field)
    editor.remove_strip(skey=3)
    session.coarse = editor.mesh
    session.record('remove strip 3')
    session.dump('floor.json')

**Each item is stored once, here.** In Rhino, a scene object draws an item and
holds a reference to it, but the scene is never saved as data: compas would write
the item a second time, and it would come back as a different object with the
same guid. The scene is rebuilt from the items instead (see
``compas_singular.rhino.session.RhinoSession``).

**Undo is a list of snapshots.** :meth:`~SingularSession.record` stores the whole
state as JSON after a change; :meth:`~SingularSession.undo` and
:meth:`~SingularSession.redo` step through them. A snapshot of a 3 200-face
project is 371 KB and 9 ms. In Rhino, ``RhinoSession`` records into the document
instead, so Rhino's own Ctrl+Z restores it.

The method names follow ``compas_session.Session`` (``dump``, ``load``,
``record``, ``undo``, ``redo``), but this does not subclass it: its ``dump``
writes the scene with its items, and its instances are shared per Python
interpreter rather than per document.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import annotations

from typing import Any
from typing import TYPE_CHECKING

import compas
from compas.data import Data

from compas_singular.settings import Settings

if TYPE_CHECKING:
    from compas_singular.datastructures import CoarsePseudoQuadMesh
    from compas_singular.datastructures import QuadMesh
    from compas_singular.framefield import CrossField
    from compas_singular.symmetry import Domain


__all__ = ['SingularSession', 'UNDO_DEPTH']


#: How many recorded states :meth:`SingularSession.undo` can step back through.
UNDO_DEPTH = 50


class SingularSession(Data):
    """Everything a project is.

    Parameters
    ----------
    settings : :class:`compas_singular.settings.Settings`, optional
    domain : :class:`compas_singular.symmetry.Domain`, optional
        The inputs: outer boundary, holes, guides and poles.
    coarse : :class:`compas_singular.datastructures.CoarsePseudoQuadMesh`, optional
    field : :class:`compas_singular.framefield.CrossField`, optional
        ``None`` on the skeleton route.
    dense : :class:`compas_singular.datastructures.QuadMesh`, optional

    Notes
    -----
    Call :meth:`record` after every change. Undo stops at the first recorded
    state, so :meth:`load` records one, and a new session should record one
    before its first change.
    """

    #: The attributes that hold items. Everything else is settings or history.
    ITEMS = ('domain', 'coarse', 'field', 'dense')

    def __init__(
        self,
        settings: Settings | None = None,
        domain: Domain | None = None,
        coarse: CoarsePseudoQuadMesh | None = None,
        field: CrossField | None = None,
        dense: QuadMesh | None = None,
        name: str | None = None,
    ) -> None:
        super(SingularSession, self).__init__(name=name)
        self.settings = settings if settings is not None else Settings()
        self.domain = domain
        self.coarse = coarse
        self.field = field
        self.dense = dense
        self._history = []      # [(step name, JSON snapshot)], oldest first
        self._current = -1      # index into _history of the state we are in

    # --------------------------------------------------------------------------
    # compas Data
    # --------------------------------------------------------------------------

    @property
    def __dtype__(self) -> str:
        # A project file is a SingularSession wherever it was written, so a file
        # dumped by RhinoSession loads in a script with no Rhino.
        return SingularSession.__clstype__()

    @property
    def __data__(self) -> dict[str, Any]:
        data = {'settings': self.settings.model_dump()}
        for item in self.ITEMS:
            data[item] = getattr(self, item)
        return data

    @classmethod
    def __from_data__(cls, data: dict[str, Any]) -> SingularSession:
        settings = Settings.model_validate(data.get('settings') or {})
        return cls(settings=settings, **{item: data.get(item) for item in cls.ITEMS})

    def clear(self, *items: str) -> None:
        """Empty the named items (``'coarse'``, ``'dense'``, ...). Settings stay.

        Only the ones named: ``clear(*session.ITEMS)`` empties them all.
        """
        for item in items:
            if item not in self.ITEMS:
                raise ValueError('not a session item: {!r}; the items are {}'.format(item, self.ITEMS))
            setattr(self, item, None)

    # --------------------------------------------------------------------------
    # files
    # --------------------------------------------------------------------------

    def dump(self, filepath: str) -> str:
        """Write the whole project to one JSON file. Returns ``filepath``."""
        compas.json_dump(self, filepath)
        return filepath

    @classmethod
    def load(cls, filepath: str) -> SingularSession:
        """Read a project written by :meth:`dump`, recorded as the first undo state."""
        session = cls()
        session.take(compas.json_load(filepath))
        session.record('load')
        return session

    # --------------------------------------------------------------------------
    # undo
    # --------------------------------------------------------------------------

    def record(self, name: str) -> None:
        """Remember the current state as the step called ``name``. Drops any redo."""
        del self._history[self._current + 1:]
        self._history.append((name, compas.json_dumps(self)))
        del self._history[:-UNDO_DEPTH]
        self._current = len(self._history) - 1

    def undo(self) -> bool:
        """Go back to the state before the last recorded step. ``False`` if there is none."""
        if self._current < 1:
            return False
        self._current -= 1
        self.take(compas.json_loads(self._history[self._current][1]))
        return True

    def redo(self) -> bool:
        """Go forward again after :meth:`undo`. ``False`` if there is nothing to redo."""
        if self._current >= len(self._history) - 1:
            return False
        self._current += 1
        self.take(compas.json_loads(self._history[self._current][1]))
        return True

    @property
    def history(self) -> list[str]:
        """The names of the recorded steps, oldest first."""
        return [name for name, _ in self._history]

    def take(self, other: SingularSession) -> None:
        """Make this session hold ``other``'s settings and items -- a session
        loaded from a file, say. Keeps ``self``, so whatever refers to this
        session (a Rhino document) stays valid."""
        if not isinstance(other, SingularSession):
            raise TypeError('not a compas_singular session: {}'.format(type(other).__name__))
        self.settings = other.settings
        for item in self.ITEMS:
            setattr(self, item, getattr(other, item))
