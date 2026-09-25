********************************************************************************
compas_singular.framefield
********************************************************************************

.. currentmodule:: compas_singular.framefield

.. rst-class:: lead

The frame-field route: solve a cross field, trace its separatrices into a coarse layout, and densify along the field.


Classes
=======

.. autosummary::
    :toctree: generated/
    :nosignatures:

    FieldDecomposition
    CrossField
    FieldInputs
    BackgroundMesh
    Constraint
    PointLocator
    Symmetry


Constraints
===========

.. autosummary::
    :toctree: generated/
    :nosignatures:

    from_boundary
    from_curves


Densification and relaxation
============================

.. autosummary::
    :toctree: generated/
    :nosignatures:

    field_densification
    relax_mesh


Quality
=======

.. autosummary::
    :toctree: generated/
    :nosignatures:

    mesh_quality
    curve_alignment
    hard_floor
    guide_metrics
