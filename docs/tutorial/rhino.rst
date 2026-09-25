********************************************************************************
Rhino 8
********************************************************************************

.. rst-class:: lead

In Rhino, the pipeline is a sequence of commands, numbered in workflow order.
Each reads and writes named layers under ``TopologyProblem::``,
and the session state is kept with the document, so it survives between commands and a save and reopen.

The commands are the scripts in ``rhino_plugin/commands/``.
They run on Rhino 8's CPython 3.9; see :doc:`/installation` to install the package there.


Commands
========

.. list-table::
   :header-rows: 1
   :widths: 30 8 62

   * - Command
     - Step
     - What it does
   * - ``CMD00_start``
     - 1
     - Reset the project: empty the ``TopologyProblem`` layers and recreate them. The settings stay.
   * - ``CMD99_settings``
     - 1
     - Show and edit the settings of this document.
   * - ``CMD01_boundary_selection``
     - 2
     - Pick the outer boundary, holes, point features and guides.
   * - ``CMD02_coarse_mesh``
     - 3
     - Choose the frame-field or the skeleton route, and build the coarse layout.
   * - ``CMD02_read_coarse_mesh``
     - 3
     - Read a coarse layout drawn by hand over the domain boundaries, instead of generating one.
   * - ``CMD03_edit_coarse_mesh``
     - 4
     - Edit the coarse layout by hand: move a corner or a pole, cut a new line along a drawn polyline or arc,
       unzip a run of corners into a strip, remove a strip, undo. Nothing is written until you commit.
   * - ``CMD04_densities``
     - 5
     - Pick a strip and set its density, or set a target length for every strip.
   * - ``CMD05_dense_pattern``
     - 5
     - Set the pattern of each coarse patch: ortho, diagonal or fan.
   * - ``CMD06_quad_mesh``
     - 6
     - Densify the coarse layout into the quad mesh.
   * - ``CMD07_smooth``
     - 7
     - Smooth the whole mesh or a region, optionally with guide curves, from one option line.
   * - ``CMD08_dual``
     - 8
     - Take the dual of a picked quad mesh.
   * - ``CMD09_edit_quad_mesh``
     -
     - Edit the final mesh by hand: move or remove vertices, edges and faces, draw new edges,
       add or remove a line of quads. The result no longer has to be all quads.
       Densifying the layout again discards the edit.
   * - ``CMD99_session_save`` / ``CMD99_session_open``
     -
     - Save the whole project to one JSON file, or open one into this document.
   * - ``CMD99_mesh_export``
     -
     - Export one picked mesh to a JSON file a script can load again.
   * - ``CMD99_clear``
     -
     - Clear chosen layers to start over.

Densities are matched to strips by geometry, never by strip number,
because a layout read back from Rhino can come back with its strips renumbered.
