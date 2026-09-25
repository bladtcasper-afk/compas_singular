********************************************************************************
Python workflow
********************************************************************************

.. rst-class:: lead

The pipeline as a plain Python script, step by step.
Every step also has a worked example with a picture in :doc:`/examples`.


Input
=====

A domain is given as lists of ``[x, y, z]`` points.
A curved boundary is given as a sampled point list: the points *are* the boundary from then on.

.. code-block:: python

    import math

    outer = [[-6.0, -4.0, 0.0], [6.0, -4.0, 0.0], [6.0, 4.0, 0.0], [-6.0, 4.0, 0.0]]
    hole = [[1.6 * math.cos(2 * math.pi * i / 64), 1.6 * math.sin(2 * math.pi * i / 64), 0.0] for i in range(64)]

Besides the outer boundary, a domain can have

* ``inner_boundaries``: closed loops, the holes;
* ``point_features``: interior points that become poles, for a column head or the apex of a dome;
* ``polyline_features`` *(skeleton route)*: curves the layout is cut along, for a crease or a rib;
* ``guides`` *(frame-field route)*: curves the quads should run along, for a cable or a force line.

Sample a curved boundary finer than ``target_length``.
Sampled too coarsely, a circle is a polygon, and every point is a corner (:doc:`/examples/12_curved_boundaries`).


Coarse layout
=============

The skeleton route:

.. code-block:: python

    from compas_singular.algorithms import SkeletonDecomposition

    decomposition = SkeletonDecomposition.from_boundary(outer, inner_boundaries=[hole], target_length=0.5)
    coarse = decomposition.coarse_mesh()

The frame-field route:

.. code-block:: python

    from compas_singular.framefield import FieldDecomposition

    decomposition = FieldDecomposition.from_boundary(outer, inner_boundaries=[hole], target_length=0.5, relax=True)
    coarse = decomposition.coarse_mesh()

    for warning in decomposition.warnings():
        print(warning)

``target_length`` is the spacing of the triangulation the layout is computed on, not the size of the quads.
Finer is more accurate and slower.
On the frame-field route, ``relax=True`` gives more accurate singularity positions than the default solve,
and ``warnings()`` says whether any part of the layout was repaired: a layout is returned either way.


Editing
=======

:class:`~compas_singular.editing.CoarseEditor` edits a copy of the coarse mesh, strip by strip.
Corners and edges are picked by a point near them, so the points below fit the skeleton layout of this domain;
``locate`` returns ``None`` when nothing is near the point.
An edit that cannot be made returns ``ok=False`` with the reason and changes nothing.

.. code-block:: python

    from compas_singular.editing import CoarseEditor

    editor = CoarseEditor(coarse, loops=[outer, hole], polylines=decomposition.polylines)

    ok, notes = editor.divide(points=[[4.75, -4.0, 0.0], [4.75, 4.0, 0.0]])
    kind, edge, point = editor.locate([2.55, -1.12, 0.0])
    ok, notes = editor.remove_strip(edge=edge)

    coarse, notes = editor.commit()
    edges_to_curves, tally = editor.edge_curves()

An edit that changes the strips renumbers them, so densities are set after editing.
See :doc:`/examples/03_coarse_mesh_editing`.


Densities and patterns
======================

.. code-block:: python

    coarse.collect_strips()
    coarse.set_strips_density_target(0.4)

``set_strips_density_target`` derives the density of every strip from a target edge length.
``set_strip_density`` sets an exact number for one strip (:doc:`/examples/04_densities`).
``set_face_pattern`` fills a patch with ``'ortho'``, ``'diagonal'`` or ``'fan'`` instead of a plain grid (:doc:`/examples/05_patterns`).


Dense quad mesh
===============

.. code-block:: python

    dense = coarse.densify()

Each patch is blended from its own four sides.
The patch edges keep the curved shape of the boundary and skeleton they came from;
``boundary_curvature=False`` or ``skeleton_curvature=False`` use straight edges instead.
After an edit, pass ``overwrite_edges_to_curves=edges_to_curves`` from the editor.

On the frame-field route, passing the field fills each patch interior along the field,
so the quads follow it everywhere, including near a guide:

.. code-block:: python

    dense = coarse.densify(overwrite_edges_to_curves=decomposition.edges_to_curves(), field=decomposition.get_field())
    print(decomposition.quality(dense))


Smoothing and the dual
======================

:doc:`smoothing` compares the smoothing functions.
The dual turns every vertex of the quad mesh into a block (:doc:`/examples/13_dual_blocks`):

.. code-block:: python

    from compas.datastructures import Mesh
    from compas_singular.algorithms.dual_mesh import dual_mesh

    plain = Mesh.from_vertices_and_faces(*dense.to_vertices_and_faces(keep_keys=False))
    blocks = dual_mesh(plain, redistribute=True)

The dual wants a plain COMPAS mesh, renumbered from zero.
