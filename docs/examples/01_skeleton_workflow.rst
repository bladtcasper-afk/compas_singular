********************************************************************************
The basic skeleton workflow
********************************************************************************

.. rst-class:: lead

The skeleton route from boundary to dense quad mesh: the medial axis of the domain cuts it into a few large quad patches, and each patch is filled with a grid of quads.

.. figure:: ../../examples/GitHub/images/01_skeleton_workflow.png
    :figclass: figure
    :class: figure-img img-fluid

.. literalinclude:: ../../examples/GitHub/01_skeleton_workflow.py
    :language: python
