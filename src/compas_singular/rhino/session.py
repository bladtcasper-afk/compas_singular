"""**The session of one Rhino document, kept inside that document.**

Every command starts with :meth:`RhinoSession.current` and, if it changed
something, ends with :meth:`RhinoSession.record`::

    session = RhinoSession.current()
    coarse = session.coarse.copy()          # edit a copy ...
    coarse.set_strip_density(3, 8)
    session.coarse = coarse                 # ... and hand it back when done
    session.record('Densities')

``record`` brings the permanent display of the layout and the dense mesh up to
date (:meth:`~RhinoSession.draw`) and writes a snapshot of the whole session into
the document (see :mod:`compas_singular.rhino.document`), both inside the
command's undo record. Rhino's own Ctrl+Z and Ctrl+Y then restore the two
together, and the next :meth:`~RhinoSession.current` notices the document holds
a different revision and reloads. So :meth:`undo` and :meth:`redo` are not used
here; they are Rhino's.

**Edit a copy, assign it back, then record.** The session in memory is reused
from one command to the next as long as the document holds the revision it last
wrote. A command that changed an item in place and was then cancelled would
leave that change in memory, and the next command would build on it. Assigning
back only on the way to ``record`` keeps what is in memory equal to what is in
the document -- and it is how :meth:`~RhinoSession.draw` knows what changed.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import uuid

import compas

from compas_singular.session import SingularSession

from . import document

try:
    import scriptcontext as sc
except ImportError:
    sc = None

try:
    import rhinoscriptsyntax as rs
except ImportError:
    rs = None


__all__ = ['RhinoSession']


def display_options():
    """How the permanent display draws each item, by item name.

    The same layers the commands always baked to. The layout's parts -- poles,
    edge curves, the polylines those come from -- are drawn by its scene object
    too. The dense mesh leaves the sublayers of ``QuadMesh`` alone: ``Edited``,
    ``Smoothened`` and ``Dual`` belong to other commands.
    """
    from .project import layer_path
    return {
        'coarse': dict(layer=layer_path('Mesh'), show_faces=True, joined=True,
                       poles_layer=layer_path('Poles'),
                       curves_layer=layer_path('EdgeCurves'),
                       polylines_layer=layer_path('Polylines')),
        'dense': dict(layer=layer_path('QuadMesh'), show_faces=True, joined=True,
                      show_vertices=False, show_edges=False),
    }


class RhinoSession(SingularSession):
    """A :class:`~compas_singular.session.SingularSession` stored in a Rhino document.

    Attributes
    ----------
    doc : Rhino.RhinoDoc
        The document this session belongs to.
    revision : str or None
        The revision of the snapshot this session was loaded from or last wrote.
        ``None`` for a document that holds no snapshot yet.
    """

    #: ``doc.RuntimeSerialNumber`` -> the session of that document.
    _sessions = {}

    #: The items that have a permanent display.
    DISPLAYED = ('coarse', 'dense')

    def __init__(self, *args, **kwargs):
        super(RhinoSession, self).__init__(*args, **kwargs)
        self.doc = None
        self.revision = None
        self._scene = None
        self._drawn = {}                # item name -> the object the display shows

    @property
    def scene(self):
        """The compas ``Scene`` that draws this session's items. Never saved: it
        is rebuilt from the items (see :mod:`compas_singular.session`)."""
        if self._scene is None:
            from compas.scene import Scene

            from .scene import ensure_registered

            ensure_registered()
            self._scene = Scene(context='Rhino')
        return self._scene

    @classmethod
    def current(cls, doc=None):
        """The session of ``doc`` (default: the active document).

        Reloaded from the document when it holds a different revision than the
        one in memory: after Ctrl+Z or Ctrl+Y, after opening a ``.3dm``, and on
        the first command of a Rhino session.
        """
        doc = doc or sc.doc
        revision = document.read_revision(doc)
        session = cls._sessions.get(doc.RuntimeSerialNumber)
        if session is None or session.revision != revision:
            session = cls._load(doc, revision)
            cls._sessions[doc.RuntimeSerialNumber] = session
        return session

    @classmethod
    def _load(cls, doc, revision):
        session = cls()
        session.doc = doc
        session.revision = revision
        text = document.read_snapshot(doc) if revision else None
        if text:
            session.take(compas.json_loads(text))
        else:
            # A document from before sessions existed, or a new one.
            from .legacy import read_legacy
            read_legacy(session)
        # The document already shows what it holds: Rhino restored the display
        # together with the snapshot, or it was drawn when the snapshot was.
        session._drawn = {name: getattr(session, name) for name in cls.DISPLAYED}
        return session

    def record(self, name):
        """Bring the display up to date and write this session into its document,
        both as part of the running command.

        ``name`` is for the reader of the command: Rhino's undo list names the
        step after the command itself.
        """
        previous = rs.EnableRedraw(False)
        try:
            self.draw()
        finally:
            rs.EnableRedraw(previous)
        self.doc.Views.Redraw()
        self.revision = str(uuid.uuid4())
        document.write_snapshot(self.doc, self.revision, compas.json_dumps(self))

    # --------------------------------------------------------------------------
    # the permanent display
    # --------------------------------------------------------------------------

    def draw(self):
        """Redraw every displayed item that is not the object drawn last time.

        By the copy rule a changed item is always a NEW object, so an identity
        check is the whole test, and a record that only changed settings redraws
        nothing.
        """
        options = display_options()
        for name in self.DISPLAYED:
            item = getattr(self, name)
            if name in self._drawn and self._drawn[name] is item:
                continue
            self._draw_item(item, options[name])
            self._drawn[name] = item

    def _draw_item(self, item, options):
        """Clear the item's layers and draw it there; with ``None``, only clear."""
        from .project import ensure_layers

        ensure_layers()
        if item is None:
            import compas_rhino.layers

            for key in ('layer', 'poles_layer', 'curves_layer', 'polylines_layer'):
                if options.get(key):
                    compas_rhino.layers.clear_layer(options[key], include_children=False, purge=False)
            return
        display = self.scene.add(item, **options)
        display.clear()
        display.draw()
        self.scene.remove(display)
