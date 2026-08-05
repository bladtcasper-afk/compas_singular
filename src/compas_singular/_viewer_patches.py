"""Runtime patches for known compas_viewer bugs.

compas_viewer (2.0.2) exposes per-object visibility through a parent-index
chain baked into a GPU settings texture: a leaf object's *effective* show
value is computed in the vertex shader by walking ``object -> parent ->
parent's parent -> ...`` and multiplying their ``show`` flags together
(``compas_viewer/renderer/shaders/model.vert:getEffectiveShow``). For this to
work, every ancestor -- including :class:`compas_viewer.scene.Group` -- needs
a row in the buffer manager's settings/transform textures so children can
find their parent's index.

``Renderer.init()`` (runs once, the first time the viewer renders) builds
that buffer correctly and includes ``Group`` objects. ``Renderer.
rebuild_buffers()`` (runs on every later ``viewer.scene.add()``/``.remove()``
call, e.g. adding new elements after the window has been shown once) instead
excludes ``Group`` from the buffer entirely. The first ``scene.add()`` call
made after the viewer has ever been shown -- which is exactly what "add more
elements, then reopen the viewer" requires -- drops every ``Group`` row from
the buffer, so every object's parent index resolves to -1 from then on. The
group's "Show" checkbox in the sidebar keeps setting ``group.show``, but
nothing reads it anymore: group-based show/hide silently and permanently
stops working for the rest of the process.

Call :func:`patch_group_visibility` once after importing
``compas_viewer.Viewer`` to make ``rebuild_buffers`` match ``init()``.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

__all__ = ["patch_group_visibility"]


def patch_group_visibility():
    """Make ``Renderer.rebuild_buffers`` keep Groups in the visibility buffer.

    Idempotent -- safe to call more than once (e.g. from multiple example
    scripts in the same process).
    """
    from compas_viewer.renderer.renderer import Renderer
    from compas_viewer.scene import Group
    from compas_viewer.scene import TagObject

    if getattr(Renderer.rebuild_buffers, "_compas_singular_patched", False):
        return

    import OpenGL.GL as GL

    def rebuild_buffers(self):
        """Rebuild the buffers (patched: keeps Group rows, matching Renderer.init())."""
        GL.glBindVertexArray(self._vao)
        self.buffer_manager.clear()

        for obj in self.viewer.scene.objects:
            if not isinstance(obj, (Group, TagObject)) and not obj._inited:
                obj.init()

        for obj in self.viewer.scene.objects:
            if not isinstance(obj, TagObject):
                self.buffer_manager.add_object(obj)
        self.buffer_manager.create_buffers()
        GL.glBindVertexArray(0)

    rebuild_buffers._compas_singular_patched = True
    Renderer.rebuild_buffers = rebuild_buffers
