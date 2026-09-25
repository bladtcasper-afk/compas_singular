********************************************************************************
Smoothing the quad mesh
********************************************************************************

.. rst-class:: lead

A dense mesh straight from ``densify()`` is made of patches that were each filled on their own, so the quads kink where patches meet. Smoothing moves the vertices to even out the quads, without changing the topology.

.. figure:: ../../examples/GitHub/images/06_quad_mesh_smoothing.png
    :figclass: figure
    :class: figure-img img-fluid

.. literalinclude:: ../../examples/GitHub/06_quad_mesh_smoothing.py
    :language: python
