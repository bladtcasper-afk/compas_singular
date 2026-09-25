********************************************************************************
compas_singular.datastructures
********************************************************************************

.. currentmodule:: compas_singular.datastructures

.. rst-class:: lead

Mesh data structures for quad meshes, coarse layouts and pseudo-quads, with the strip grammar, smoothing and colouring that act on them.


Classes
=======

.. autosummary::
    :toctree: generated/
    :nosignatures:

    Mesh
    QuadMesh
    CoarseQuadMesh
    PseudoQuadMesh
    CoarsePseudoQuadMesh
    Skeleton
    BoundaryLoop


Strip grammar
=============

.. autosummary::
    :toctree: generated/
    :nosignatures:

    add_strip
    add_strips
    split_strip
    split_strips
    open_added_strip
    strip_polyedge_update
    is_polyedge_valid_for_strip_addition
    polyedge_sides
    delete_strip
    delete_strips
    collateral_strip_deletions
    total_boundary_deletions
    strips_to_split_to_prevent_boundary_collapse


Pole grammar
============

.. autosummary::
    :toctree: generated/
    :nosignatures:

    split_quad_in_pseudo_quads
    merge_pseudo_quads_in_quad


Shape grammar
=============

.. autosummary::
    :toctree: generated/
    :nosignatures:

    add_opening
    add_handle


Smoothing
=========

.. autosummary::
    :toctree: generated/
    :nosignatures:

    constrained_smoothing
    boundary_constrained_smoothing
    boundary_smoothing
    region_smoothing
    relaxation
    automated_boundary_constraints
    closest_point_on_constraint
    mesh_boundary_loops
    mesh_boundary_polylines
    mesh_boundary_corners
    automated_smoothing_surface_constraints
    automated_smoothing_constraints
    surface_constrained_smoothing


Coarse edge curves
==================

.. autosummary::
    :toctree: generated/
    :nosignatures:

    coarse_edges_to_curves
    snap_corners_to_walls
    mean_edge_length


Coarse mesh from polylines
==========================

.. autosummary::
    :toctree: generated/
    :nosignatures:

    split_at_corners
    split_at_junctions
    weld_network
    check_network
    faces_from_network
    check_faces


Colouring
=========

.. autosummary::
    :toctree: generated/
    :nosignatures:

    mesh_vertex_2_coloring
    mesh_vertex_n_coloring
    mesh_face_2_coloring
    mesh_face_n_coloring
    quad_mesh_strip_2_coloring
    quad_mesh_strip_n_coloring
    quad_mesh_polyedge_2_coloring
    quad_mesh_polyedge_n_coloring
    dense_quad_mesh_polyedge_2_coloring


Morphing
========

.. autosummary::
    :toctree: generated/
    :nosignatures:

    fold
    fold_vertex_group


Operations
==========

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
