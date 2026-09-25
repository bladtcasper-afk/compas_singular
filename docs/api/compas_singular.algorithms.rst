********************************************************************************
compas_singular.algorithms
********************************************************************************

.. currentmodule:: compas_singular.algorithms

.. rst-class:: lead

Decomposition of a domain into a coarse quad layout along its topological skeleton, and the propagation and triangulation helpers it uses.


Classes
=======

.. autosummary::
    :toctree: generated/
    :nosignatures:

    SkeletonDecomposition


Triangulation
=============

.. autosummary::
    :toctree: generated/
    :nosignatures:

    boundary_triangulation
    weld_polyline_features
    arrange_polyline_features
    as_points
    as_curves


Propagation
===========

.. autosummary::
    :toctree: generated/
    :nosignatures:

    quadrangulate_mesh
    quadrangulate_faces
    quadrangulate_face
    discrete_coons_patch_mesh
    update_adjacent_face


Dual
====

.. currentmodule:: compas_singular.algorithms.dual_mesh

.. autosummary::
    :toctree: generated/
    :nosignatures:

    dual_mesh
