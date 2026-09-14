"""Moved to :mod:`compas_singular.editing.rebuild`.

Nothing in it ever touched the field, the tracer or the background: it welds a
hand-edited layout from plain geometry, snaps its boundary corners to the domain
walls and repairs the result to all-quad. That is what ``CoarseEditor.commit``
does, so it belongs beside the editor -- and the editor must not have to import
the field solver to reach it.

This module stays so existing callers keep working. Import from
``compas_singular.editing`` in anything new.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ..editing.rebuild import PRECISION  # noqa: F401
from ..editing.rebuild import closest_on_loop  # noqa: F401
from ..editing.rebuild import coarse_from_skeleton  # noqa: F401
from ..editing.rebuild import face_polylines  # noqa: F401
from ..editing.rebuild import faces_from_geometry  # noqa: F401
from ..editing.rebuild import mesh_from_faces  # noqa: F401
from ..editing.rebuild import snap_to_loops  # noqa: F401
from ..editing.rebuild import warp_polyline  # noqa: F401


__all__ = ['coarse_from_skeleton', 'warp_polyline', 'face_polylines',
           'faces_from_geometry', 'mesh_from_faces', 'snap_to_loops',
           'closest_on_loop', 'PRECISION']
