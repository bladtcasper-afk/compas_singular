********************************************************************************
COMPAS Singular
********************************************************************************

.. figure:: /_images/01_front.jpg
     :figclass: figure
     :class: figure-img img-fluid

.. rst-class:: lead

:mod:`compas_singular` turns a 2D domain -- a boundary, optional holes, and optional
point features, curve features and guides -- into a quad mesh with controlled topology.
You choose where the irregular vertices (the *singularities*) sit,
instead of getting them wherever a generic quadrangulation happens to put them.

A coarse quad layout is built from the domain along one of two routes:
the topological skeleton of the domain, or a cross field solved over it.
The layout can then be edited strip by strip, given densities and patterns,
densified into a dense quad mesh, smoothed, and turned into its dual.
The same pipeline runs as a Python API, as a set of Rhino 8 commands, and as an MCP server.

:mod:`compas_singular` stems from the doctoral research of Robin Oval
on topology finding of patterns for structural design,
and is built on the `COMPAS <https://compas.dev>`_ framework.


Table of Contents
=================

.. toctree::
   :maxdepth: 3
   :titlesonly:

   Introduction <self>
   installation
   tutorial
   examples
   api
   publications
   citing
   license
   acknowledgements


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
