********************************************************************************
Overview
********************************************************************************

:mod:`compas_singular` turns a 2D domain — a boundary, optional holes, and
optional internal points/curves/cables — into a **quad mesh with controlled
topology**: you choose where the irregular vertices (*singularities*) sit,
instead of getting them wherever a generic quadrangulation algorithm happens to
put them.

There is one pipeline, with **two ways to build the first, coarse layout**, and
everything downstream of that layout — editing, setting densities, densifying,
smoothing, taking the dual — is shared. The whole pipeline is available three
ways: as a **Python API** (headless, scriptable), as a **Rhino plugin**
(a sequence of commands, ``CMD_*``, driven by picks and prompts), and as an
**AI-assisted session** (an MCP server driving the same underlying tools from
natural language, with the Rhino document attached through ``CMD_mcp_link``).

This page is the map of that workflow: what each stage does, what it can be
told, and where to find it in each of the three front ends. It is deliberately
not a function-by-function reference — that lives in the API docs.

|

----


The two routes to a coarse layout
==================================

Everything starts from a **coarse mesh**: one quad face per topological patch,
each patch edge carrying a *strip density*. It is cheap to edit and carries no
geometry detail yet. Two independent algorithms build one from a domain:

.. list-table::
   :header-rows: 1
   :widths: 18 41 41

   * -
     - **Skeleton route**
     - **Field route**
   * - class
     - ``SkeletonDecomposition``
     - ``FieldDecomposition``
   * - idea
     - the topological skeleton (medial axis) of a Delaunay triangulation of
       the domain partitions it into quad patches
     - a cross field is solved over the domain and traced into separatrices,
       which partition it into quad patches
   * - can be steered by guides (cables, force lines)?
     - no — a medial axis is equidistant from the walls by definition
     - **yes** — a guide is the reason to pick this route
   * - respects input symmetry?
     - only incidentally
     - yes, with ``symmetry='auto'``
   * - line ("curve") features
     - yes — a curve becomes a topological cut, patches lie along it
     - not applicable (use guides instead)
   * - point features → poles
     - yes
     - yes
   * - measured quality, 18-domain baseline
     - wins 1 of 18
     - wins 15 of 18
   * - weak spot
     - a free curve-feature tip still isn't fully resolved (see the curve
       feature notes in the developer docs)
     - a domain of *uniform curvature* loses quality at patch **corners**; a
       sampled arc that turns more than 45° between two points reads as a
       corner to the field

**Rule of thumb:** start with the field route. Use the skeleton route only when
you need a line feature the mesh must be cut along (not just steered by), or
when you deliberately want the medial-axis layout.

|

----


The eight-step pipeline
========================

Both routes share the same eight stages. This is the spine of the whole
workflow — the same eight steps exist as eight Rhino commands and as eight
sections of a plain Python script.

.. list-table::
   :header-rows: 1
   :widths: 6 22 36 36

   * - #
     - Stage
     - What happens
     - Where
   * - 1
     - **Settings**
     - background spacing, target edge length, guide mode, symmetry, caching
     - ``CMD_start`` / a settings dict at the top of the script
   * - 2
     - **Input**
     - outer boundary, holes, point features, line features and/or guides
     - ``CMD_boundary_selection`` / plain point lists or curves
   * - 3
     - **Coarse layout**
     - build the patch layout — this is where Skeleton vs Field is chosen
     - ``CMD_coarse_mesh`` / ``SkeletonDecomposition`` or ``FieldDecomposition``
   * - 4
     - **Edit** *(optional)*
     - move a corner, cut a new division, unzip a strip, delete a strip
     - ``CMD_edit_coarse_mesh`` / ``CoarseEditor``
   * - 5
     - **Densities**
     - how many elements across each strip — per strip, or a global target
     - ``CMD_densities`` (+ ``CMD_dense_pattern``) / ``set_strip(s)_density*``
   * - 6
     - **Quad mesh**
     - fill every patch with an actual grid of quads
     - ``CMD_quad_mesh`` / ``coarse.densification(...)`` or ``decomposition.quad_mesh(...)``
   * - 7
     - **Smoothing** *(optional)*
     - relax the dense mesh to improve element shape
     - ``CMD_smoothen`` / ``CMD_smoothen_guide`` / ``relax`` and friends
   * - 8
     - **Dual mesh** *(optional)*
     - turn the quad mesh into its dual (e.g. for a block layout)
     - ``CMD_dual`` / ``dual_mesh(...)``

The two canonical, runnable references for this pipeline are
``examples/11_skeleton_workflow.py`` and ``examples/12_framefield_workflow.py``:
the same eight sections, in the same order, with the route as the only
structural difference. Read them side by side to see exactly what each route
needs and returns.

|

----


The Python workflow
====================

A minimal field-route script is, at its core, four lines::

    from compas_singular.framefield.field_decomposition import FieldDecomposition

    decomposition = FieldDecomposition.from_boundary(outer, inner_boundaries=holes,
                                                       guides=guides, target_length=0.5)
    coarse = decomposition.coarse_mesh(point_features)
    coarse.set_strips_density_target(t=0.5)
    dense = decomposition.densify()            # or coarse.densification(...)

The skeleton route mirrors it, one import and one call to build the layout::

    from compas_singular.algorithms import SkeletonDecomposition

    decomposition = SkeletonDecomposition.from_boundary(outer, inner_boundaries=holes,
                                                          point_features=points,
                                                          target_length=0.5)
    coarse = decomposition.coarse_mesh()
    coarse.set_strips_density_target(t=0.5)
    dense = coarse.densification()

Everything past step 3 is genuinely shared code. What varies at each step:

**Step 2 — input.** A list of ``[x, y, z]`` points (or ``compas.geometry``
``Point``/``Polyline`` objects — both are accepted and coerced) for:

- ``outer_boundary`` — required, closed;
- ``inner_boundaries`` — a list of closed loops (holes); faces inside one are removed;
- ``point_features`` — interior points that become **poles** (an odd number of
  strips fanning into one vertex — e.g. a column head);
- ``polyline_features`` *(skeleton route only)* — curves the layout is **cut
  along**: the curve survives as a continuous course of mesh edges, not by
  snapping;
- ``guides`` *(field route only)* — cables / force lines the field is steered
  by, with ``mode='tangent'`` (mesh runs *along* the guide) or
  ``'perpendicular'`` (runs *across* it).

A curved boundary or hole should be sampled *finer* than the background
spacing before it is handed in — the wall is the boundary from that point on,
so under-sampling costs coverage (a densified circle can come out tens of
percent short of its true radius if sampled too coarsely, see
``05_boundary_curvature`` in the developer notes).

**Step 3 — coarse layout, key options.**

- ``target_length`` — the *background* triangulation/field spacing, **not**
  the eventual quad size; finer is a better layout and a slower solve;
- ``symmetry='auto'`` *(field route)* — detects the symmetry group of the
  input and holds it through the solve;
- ``relax=True`` *(field route)* — diffusion + normalisation instead of a
  single Dirichlet solve; use it whenever guides are present (the default
  solver is measurably wrong about where singularities sit without it) —
  except with a **closed** guide, where relaxation is measurably worse;
- always read ``decomposition.route()`` and ``decomposition.warnings()``
  afterwards — the front end always returns *a* mesh, including when the
  separatrix network did not close and it silently fell back to a coarser
  covering.

**Step 4 — editing the coarse layout, via ``CoarseEditor``.** Anchored on the
layout itself (not on a particular route), it takes ``loops`` (the domain
walls) and ``polylines`` (the traced branches/separatrices, so an edited edge
keeps its curve instead of falling back to a straight chord), and offers:

- ``move_vertex`` — drag a corner; a boundary corner is projected back onto
  its own wall, an interior one moves freely;
- ``divide`` — split a band of patches in two, either by strip key and a
  parameter ``t``, or along a drawn curve (``points=``); must enter/leave
  every patch through opposite sides;
- ``add_strip`` / ``remove_strip`` — unzip a new band along a chosen run of
  corners, or delete one entirely;
- ``commit()`` — writes the edit back **in place** and returns the mesh (or
  ``None`` on refusal — test that, not truthiness). A moves-only edit keeps
  every strip key and density; a cut/divide/strip edit re-welds and
  **renumbers everything**.

**Step 5 — densities.** ``coarse.set_strip_density(skey, n)`` for one strip,
``set_strips_density(...)`` for several, or ``set_strips_density_target(t)``
to derive a density per strip from a target edge length. Optionally,
``coarse.set_pattern(...)`` / per-patch pattern assignment chooses ``ortho``
(default grid), ``diagonal`` (both diagonals cut through the grid) or ``fan``
(four polar fans meeting at the centre) — diagonal and fan need an even
density on both strips crossing the patch, and cost measured quality on
purpose (they are triangulated/fanned by construction), so set patterns
*before* final densities.

**Step 6 — quad mesh.** ``coarse.densification()`` blends each patch from its
own four (subdivided) sides — a discrete Coons patch. By default
(``boundary_curvature=True, skeleton_curvature=True``) it uses whatever shape
``coarse.edges_to_curves()`` has stored for each edge — the two toggles let a
caller keep curvature on the layout's own boundary while chording the interior
separatrices, or the reverse. Without a stored (or overwritten, via
``overwrite_edges_to_curves=...``) shape every edge densifies as a straight
chord and a curved boundary loses real area. On the field route,
``decomposition.densify()`` (or
``coarse.densification(..., field=decomposition.get_field())``) additionally
integrates the **patch interior** from the field — the only way a guide
reaches a patch it forced no topology in. ``decomposition.quad_mesh(...)`` is
steps 5+6 combined, with two fallbacks if the separatrix network did not
close; read ``decomposition.route()`` after calling it.

**Step 7 — smoothing.** Not a single call but a small decision ladder — see
*Smoothing and quality* below.

**Step 8 — dual mesh.** ``compas_singular.algorithms.dual_mesh.dual_mesh(mesh,
redistribute=True)`` (despite the module path, it needs no Rhino) turns the
quad mesh into its dual, redistributing boundary faces so they come out full
size rather than half.

|

----


Smoothing and quality
======================

Quality here means *element* quality — the shape of the individual quads —
not whether the layout/topology is right. Four numbers describe it:
``min_angle`` / ``max_angle`` (worst corner angle anywhere), ``aspect_max``
(longest/shortest edge, worst face) and ``share_below`` (the fraction of
corners below a low-angle threshold — the number that says whether a problem
is *one face* or *the whole mesh*). A high aspect ratio and an irregular
vertex are not automatically defects: a graded mesh is supposed to have
stretched elements, and a quad mesh over anything but a rectangle must have
some non-4-valent vertices.

Five smoothing tools, not interchangeable, tried in this order:

1. **relax** — the safe first move. Tries a schedule of settings
   gentlest-first and **keeps a result only if min angle, max angle and max
   aspect all improve together**; otherwise it leaves the mesh untouched. Not
   a failure if nothing changes — it means the layout, not the smoothing, is
   the limit.
2. **smooth_boundary_constrained** — the workhorse for a mesh that is
   generally slack. Interior vertices relax freely, boundary vertices slide
   along the walls, kinks stay pinned. No gate — it can make things worse, so
   measure before/after.
3. **smooth_region** — for *one* bad spot (``share_below`` near zero). Area-
   weighted and tapered over a few rings so it blends rather than creasing at
   its edge.
4. **smooth_guides** — when the mesh has guide curves it should hug. Attaches
   the longest run of one polyedge that already follows each guide and
   smooths with that chain (plus the boundary) held.
5. **relax_fdm** — force-density relaxation, a specialised alternative to
   ``relax`` rather than a general constrained solve (its ``constraints`` /
   ``algorithm`` arguments are not wired through).

**Never smooth a mesh whose whole purpose is interface perpendicularity to a
boundary or a force line** (e.g. a block-coursing mesh) — global smoothing has
been measured to degrade that from ~90° down to ~45–65° while barely improving
regularity.

When the layout itself is the limit rather than the element shapes,
``dense_add_line`` / ``dense_remove_line`` (adding or removing a whole strip
on the *dense* mesh) change topology directly — but any such edit is lost the
next time the coarse layout is re-densified, and is refused outright on a
mesh with a pole.

|

----


The Rhino plugin workflow
===========================

The Rhino side is the same eight-step pipeline as a sequence of commands
(``rhino_plugin/commands/CMD_*.py``), each one reading and writing named
layers under ``TopologyProblem::`` plus a small JSON side-car cached next to
the document, so state survives between commands (and a save/reopen) without
re-deriving anything that was computed and is expensive to recompute.

.. list-table::
   :header-rows: 1
   :widths: 22 14 64

   * - Command
     - Step
     - What it does
   * - ``CMD_start``
     - 1
     - not run directly — holds the settings dict every other command reads
       (background spacing, target length, guide mode, symmetry, relax,
       field caching) and the ``import_compas_singular`` bootstrap
   * - ``CMD_boundary_selection``
     - 2
     - pick/draw the outer boundary, holes, point features and guide curves
       onto their layers
   * - ``CMD_coarse_mesh``
     - 3
     - asks **FrameField** or **Skeleton**, then builds the coarse layout and
       bakes it (mesh, poles, separatrices/branches, edge curves)
   * - ``CMD_read_coarse_mesh``
     - 3 *(alt.)*
     - the mirror of the above: reads a coarse layout **you drew by hand**
       (walls on ``Outer``/``Inner``, division lines on ``EdgeCurves``)
       instead of generating one
   * - ``CMD_edit_coarse_mesh``
     - 4
     - move a corner; draw a polyline/arc across the layout to cut a new
       division; pick a run of corners to unzip into a new strip; pick a
       strip's ribbon to delete it; ``commit`` or ``reset``
   * - ``CMD_densities``
     - 5
     - pick a strip (shown as a ribbon shaded by density) and set its count,
       or set a global target for every strip at once
   * - ``CMD_dense_pattern``
     - 5 *(opt.)*
     - assign ``ortho`` / ``diagonal`` / ``fan`` to picked patches or all of
       them, before final densities
   * - ``CMD_quad_mesh``
     - 6
     - densifies the coarse layout into the dense quad mesh, using the
       document's own boundary curves for curved edges and, on the field
       route, the cached field for patch interiors
   * - ``CMD_edit_quad_mesh``
     - 6 *(alt.)*
     - hand-edit the **dense** mesh directly: drag a vertex, or pick an edge
       to add/remove the whole strip through it. Ends the pipeline for that
       mesh — re-running ``CMD_quad_mesh`` regenerates from the layout and
       discards a hand edit
   * - ``CMD_smoothen``
     - 7
     - relax / boundary-constrained / region smoothing on a picked dense mesh
   * - ``CMD_smoothen_guide``
     - 7 *(opt.)*
     - attach chains of mesh vertices to guide curves, automatically or by
       hand, then smooth with them held
   * - ``CMD_dual``
     - 8
     - takes the dual of a picked dense mesh, redistributing boundary faces
   * - ``CMD_symmetry_report``
     - diagnostic
     - read-only: reports the worst symmetry-mismatch distance per input
       element, to tell drawing imprecision from a genuinely asymmetric input
   * - ``CMD_clear``
     - utility
     - clears chosen layers (inputs, field, skeleton, quad mesh, trash) to
       start over

Settings persist in the document (``CMD_start``'s settings dict, stored in
document user text) and the field/layout are cached to disk, keyed on the
inputs, so re-opening a file does not force a re-solve. Strip **keys** are not
stable across a bake — a layout read back from Rhino geometry can come back
with its strips renumbered — so densities and patterns are always matched
geometrically (by edge midpoint), never by raw index.

|

----


AI-assisted editing
=====================

The same tools underneath ``CMD_edit_coarse_mesh``, ``CMD_densities``,
``CMD_quad_mesh`` and the smoothers are also exposed over the **Model Context
Protocol**, so they can be driven from natural language by an MCP client
(``python -m compas_singular.mcp``; the Rhino document is attached with
``CMD_mcp_link``).

They arrive as typed tools — ``create_coarse_mesh``, ``coarse_set_density``,
``coarse_densify``, ``relax``, ``smooth_boundary_constrained``,
``smooth_guides``, ``smooth_region``, ``dense_add_line`` /
``dense_remove_line``, ``inspect``, ``rhino_pull`` / ``rhino_push`` and their
``_coarse`` counterparts, among others — plus ``check_inputs``, which is run
first, because the skeleton route never raises on a bad input, it just builds a
bad layout, and a picture-based ``inspect(image=True)`` for judging what the
numbers can't say (boundary deviation, ignored point features, guides the mesh
ran across instead of along, kinks and propagated distortion).

The MCP server does not solve a frame field: it builds a layout from the
topological skeleton, and a field-route layout reaches it from Rhino, through
``CMD_coarse_mesh`` and ``rhino_pull_coarse``. (An earlier in-process agent,
``CMD_ai_edit`` and ``compas_singular.agent``, could solve one; it was removed
on 2026-09-16 in favour of this server.)

In every case, the underlying decision ladder is the same one described above
under *Smoothing and quality*: try ``relax`` first, escalate to a global or
regional pass only if it did not help, and treat two passes that change
nothing as proof the *layout* is the limit rather than something to keep
retrying.

|

----


What this page does not cover yet
====================================

This is the workflow foundation: what exists, how the pieces fit together,
and what can be configured at each stage. It deliberately stops short of a
function-by-function reference (parameter lists, return types, edge cases) —
that belongs in the generated API pages (:doc:`04_api`) and in individual
docstrings, and is the natural next layer to build on top of this one.
