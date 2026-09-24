"""
********************************************************************************
editing
********************************************************************************

.. currentmodule:: compas_singular.editing

**Hand-editing a layout, with no CAD in it.**

Every operation a user performs on a coarse layout in Rhino -- dragging a
corner, cutting patches with a drawn line, deleting a strip, committing the
result back -- is defined here, on the mesh, and can be run from a plain script
with no Rhino open. ``compas_singular.rhino`` and the ``CMD_`` commands supply
the picks, the previews and the prompts; they supply nothing else.

The split is not tidiness. This is the half that can be wrong in ways a user
cannot see -- a cut that leaves a five-sided patch, a curve map that no longer
matches its edges -- so it has to be testable without a CAD session.

**This package imports nothing from** ``framefield``, and that is enforced by
where things live rather than by convention: :mod:`~compas_singular.editing.rebuild`
(the weld/snap/repair a commit performs) and :mod:`~compas_singular.editing.repair`
(making a mesh all-quad again) moved here out of ``framefield`` precisely because
neither ever needed a field. ``framefield`` re-exports both for existing callers.

Classes
=======

.. autosummary::
    :toctree: generated/
    :nosignatures:

    MeshEditor
    CoarseEditor
    DenseMeshEditor
    GuideCurve

Functions
=========

.. autosummary::
    :toctree: generated/
    :nosignatures:

    guide_chain
    attach_chain
    chain_quality
    collect_polyedges
    mean_edge_length
    coarse_from_skeleton
    snap_to_loops
    warp_polyline
    densifiable
    solve_non_quad_faces

"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import annotations

# The shared base. Imported first: the editors below are its subclasses.
from compas_singular.editing.editor import MeshEditor  # noqa: F401

from compas_singular.editing.coarseeditor import *  # noqa: F401 F403
from compas_singular.editing.denseeditor import *  # noqa: F401 F403
from compas_singular.editing.guide_chain import *  # noqa: F401 F403

# NOT star-imported: ``rebuild`` also defines ``PRECISION``, which would rebind
# ``coarse_layout``'s. Same value today, and that is exactly why a clash here
# would go unnoticed.
from compas_singular.editing.rebuild import coarse_from_skeleton  # noqa: F401
from compas_singular.editing.rebuild import face_polylines  # noqa: F401
from compas_singular.editing.rebuild import faces_from_geometry  # noqa: F401
from compas_singular.editing.rebuild import mesh_from_faces  # noqa: F401
from compas_singular.editing.rebuild import snap_to_loops  # noqa: F401
from compas_singular.editing.rebuild import warp_polyline  # noqa: F401
from compas_singular.editing.repair import densifiable  # noqa: F401
from compas_singular.editing.repair import solve_non_quad_faces  # noqa: F401
from compas_singular.editing.repair import topological_quad_split  # noqa: F401

import types  # noqa: E402

# Never re-export the submodules themselves: ``from .foo import *`` binds ``foo``
# here too, and a module object of the same name shadows a function further up.
__all__ = [name for name, obj in list(globals().items())
           if not name.startswith('_') and not isinstance(obj, types.ModuleType)]
