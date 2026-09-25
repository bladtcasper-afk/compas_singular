"""
********************************************************************************
compas_singular.datastructures
********************************************************************************

.. currentmodule:: compas_singular.datastructures


Mesh
=====

Mesh class.

Core
----

.. autosummary::
    :toctree: generated/
    :nosignatures:

    Mesh

Operations
----------

.. autosummary::
    :toctree: generated/
    :nosignatures:

    mesh_move_by
    mesh_move_vertices_by
    mesh_move_vertex_to
    is_face_degenerate
    trimesh_face_circle
    mesh_weld
    meshes_join
    meshes_join_and_weld

Coloring
--------

.. autosummary::
    :toctree: generated/
    :nosignatures:

    mesh_vertex_2_coloring
    mesh_vertex_n_coloring
    mesh_face_2_coloring
    mesh_face_n_coloring

Smoothing
---------

Constrained smoothing, without Rhino: the projections are computed with
:mod:`compas.geometry`.

.. autosummary::
    :toctree: generated/
    :nosignatures:

    closest_point_on_constraint
    mesh_boundary_loops
    mesh_boundary_polylines
    mesh_boundary_corners
    constrained_smoothing
    automated_boundary_constraints
    boundary_constrained_smoothing
    boundary_smoothing
    region_smoothing
    relaxation

On a surface: interior vertices on the surface, boundary vertices on its borders.

.. autosummary::
    :toctree: generated/
    :nosignatures:

    automated_smoothing_surface_constraints
    automated_smoothing_constraints
    surface_constrained_smoothing



Quad Mesh
=========

QuadMesh class.

Core
----

.. autosummary::
    :toctree: generated/
    :nosignatures:

    QuadMesh

Morphing
--------

.. autosummary::
    :toctree: generated/
    :nosignatures:

    fold
    fold_vertex_group

Coloring
--------

.. autosummary::
    :toctree: generated/
    :nosignatures:

    quad_mesh_strip_2_coloring
    quad_mesh_strip_n_coloring
    quad_mesh_polyedge_2_coloring
    quad_mesh_polyedge_n_coloring

Shape grammar
-------------

.. autosummary::
    :toctree: generated/
    :nosignatures:

    add_opening
    add_handle

Pattern grammar
---------------

.. autosummary::
    :toctree: generated/
    :nosignatures:

    add_strip
    add_strips
    delete_strip
    delete_strips
    split_strip
    split_strips
    strip_polyedge_update
    strips_to_split_to_prevent_boundary_collapse
    collateral_strip_deletions
    total_boundary_deletions

Coarse Quad Mesh
================

CoarseQuadMesh class.

Core
----

.. autosummary::
    :toctree: generated/
    :nosignatures:

    CoarseQuadMesh

Coloring
--------

.. autosummary::
    :toctree: generated/
    :nosignatures:

    dense_quad_mesh_polyedge_2_coloring

Pseudo Quad Mesh
================

PseudoQuadMesh class.

Core
----

.. autosummary::
    :toctree: generated/
    :nosignatures:

    PseudoQuadMesh

Pole grammar
------------

.. autosummary::
    :toctree: generated/
    :nosignatures:

    split_quad_in_pseudo_quads
    merge_pseudo_quads_in_quad

Coarse Pseudo Quad Mesh
=======================

CoarsePseudoQuadMesh class.

.. autosummary::
    :toctree: generated/
    :nosignatures:

    CoarsePseudoQuadMesh

"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from compas_singular.datastructures.mesh import *  # noqa: F401 F403
from compas_singular.datastructures.mesh_quad import *  # noqa: F401 F403
from compas_singular.datastructures.mesh_quad_coarse import *  # noqa: F401 F403
from compas_singular.datastructures.mesh_quad_pseudo import *  # noqa: F401 F403
from compas_singular.datastructures.mesh_quad_pseudo_coarse import *  # noqa: F401 F403
from compas_singular.datastructures.skeleton import *  # noqa: F401 F403

import types  # noqa: E402

# Only re-export names bound in this namespace, never the submodules themselves:
# ``from .foo import *`` also binds ``foo`` as an attribute of this package, and
# re-exporting that module object shadows any function or subpackage of the same
# name further up the chain.
__all__ = [name for name, obj in list(globals().items())
           if not name.startswith('_') and not isinstance(obj, types.ModuleType)]
