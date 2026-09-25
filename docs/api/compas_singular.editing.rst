********************************************************************************
compas_singular.editing
********************************************************************************

.. currentmodule:: compas_singular.editing

.. rst-class:: lead

Hand-editing a coarse layout or a dense mesh, with no CAD in it. Rhino supplies only picks and prompts.


Classes
=======

.. autosummary::
    :toctree: generated/
    :nosignatures:

    MeshEditor
    CoarseEditor
    DenseMeshEditor
    GuideCurve


Guide chains
============

.. autosummary::
    :toctree: generated/
    :nosignatures:

    guide_chain
    attach_chain
    chain_quality
    collect_polyedges
    mean_edge_length


Rebuilding
==========

.. autosummary::
    :toctree: generated/
    :nosignatures:

    coarse_from_skeleton
    face_polylines
    faces_from_geometry
    mesh_from_faces
    snap_to_loops
    warp_polyline


Repair
======

.. autosummary::
    :toctree: generated/
    :nosignatures:

    densifiable
    solve_non_quad_faces
    topological_quad_split
