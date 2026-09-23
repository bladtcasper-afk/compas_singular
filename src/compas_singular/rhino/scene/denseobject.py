"""The dense quad mesh (``QuadMesh``) in Rhino."""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from compas_singular.rhino.scene.meshobject import RhinoSingularMeshObject


__all__ = ['RhinoDenseObject']


class RhinoDenseObject(RhinoSingularMeshObject):
    """A dense mesh as pickable points and lines, kept up to date with :meth:`sync`.

    Thousands of edges, so an editor redraws after a click with ``sync`` (only
    what changed) rather than ``redraw`` (everything). No faces by default: the
    points and lines are what is picked.
    """

    def __init__(self, **kwargs):
        kwargs.setdefault('show_faces', False)
        kwargs.setdefault('show_vertices', True)
        kwargs.setdefault('show_edges', True)
        super(RhinoDenseObject, self).__init__(**kwargs)
