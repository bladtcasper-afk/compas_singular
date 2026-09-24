"""**How the session's items are drawn and picked in Rhino: compas scene objects.**

``scene.add(coarse)`` gives a ``RhinoCoarseObject`` for a
``CoarsePseudoQuadMesh`` and ``scene.add(dense)`` a ``RhinoDenseObject`` for a
``QuadMesh``, the way ``scene.add(mesh)`` gives compas_rhino's
``RhinoMeshObject`` for a plain mesh. That mapping is compas's scene registry,
filled by :func:`register_scene_objects_rhino`.

It has to be filled twice over, and :func:`ensure_registered` does both:

* compas fills the registry only while it is EMPTY -- once, with every
  installed package's plugins, ours included (``__all_plugins__``). Registering
  ours first would leave compas_rhino's own scene objects out;
* ``CMD_dev_reload`` re-imports compas_singular, and the registry still maps the
  PREVIOUS import's classes. Registering again replaces them.

Nothing here imports Rhino at module level, so compas can discover the plugin
without Rhino present; ``requires=["Rhino"]`` keeps it from running there.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import annotations

from compas.plugins import plugin


__all__ = ['register_scene_objects_rhino', 'ensure_registered']


@plugin(category='factories', pluggable_name='register_scene_objects', requires=['Rhino'])
def register_scene_objects_rhino() -> None:
    from compas.scene.context import register

    from compas_singular.datastructures import CoarsePseudoQuadMesh
    from compas_singular.datastructures import QuadMesh

    from compas_singular.rhino.scene.coarseobject import RhinoCoarseObject
    from compas_singular.rhino.scene.denseobject import RhinoDenseObject

    # The registry follows the item's class hierarchy, nearest first: a layout
    # is a QuadMesh too, and still gets RhinoCoarseObject.
    register(CoarsePseudoQuadMesh, RhinoCoarseObject, context='Rhino')
    register(QuadMesh, RhinoDenseObject, context='Rhino')


def ensure_registered() -> None:
    """compas's own scene objects, then ours from THIS import. Safe to call every time."""
    from compas.scene.context import ITEM_SCENEOBJECT
    from compas.scene.context import register_scene_objects

    if not ITEM_SCENEOBJECT:
        register_scene_objects()
    register_scene_objects_rhino()
