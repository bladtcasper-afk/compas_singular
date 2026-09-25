********************************************************************************
compas_singular.geometry
********************************************************************************

.. currentmodule:: compas_singular.geometry

.. rst-class:: lead

Geometry functions on plain point lists: arrays of points, and projection onto and discretisation of polylines.


Arrays
======

.. autosummary::
    :toctree: generated/
    :nosignatures:

    circle_evaluate
    archimedean_spiral_evaluate
    line_array
    rectangular_array
    circular_array
    spiral_array


Polylines
=========

.. autosummary::
    :toctree: generated/
    :nosignatures:

    polyline_length
    point_at_length
    closest_on_polyline
    project_on_polyline
    distance_to_polyline
    distance_to_loop
    near_loop
    loop_arc_lengths
    loop_parameter
    resample_loop
    discretise_line
    discretise_boundary
    bounding_box_diagonal
    points_in_polygon_xy
    signed_area
