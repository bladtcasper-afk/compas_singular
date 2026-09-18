"""Step 9 -- **the hand-edited coarse layout**. The one seam in the pipeline.

``FieldDecomposition.quad_mesh`` owns the whole chain in one call -- solve the
field, trace the separatrices, cut the coarse layout, densify it -- and that is
what makes it usable and also what makes it closed. The coarse layout is the
one artefact a designer actually wants to touch: merge two patches, drag a
corner off a wall, delete a patch the tracer put somewhere silly. Until now
there was nowhere to reach in.

**Moved here from ``framefield/edit.py``.** Nothing in it ever touched the field,
the tracer or the background -- it welds, snaps and repairs a layout from plain
geometry -- so it belongs beside the editor that calls it rather than beside the
solver.

This module opens exactly one seam, and no more. The generated layout goes out
to Rhino, comes back edited, and densification proceeds from THAT layout. See
:meth:`FieldDecomposition.edit_coarse`, which is the entry point; everything
here is the machinery it needs.

THE FIELD IS NOT RE-SOLVED, AND THE SEPARATRICES ARE KEPT
---------------------------------------------------------

The edit enters at ``decomposition.mesh`` and nowhere else. ``field``,
``tracer``, ``background`` and ``separatrices`` are untouched, which is the
requirement -- an edit to the layout is not evidence about the field, and
letting a dragged corner feed back into the solve would mean the user could
never move anything without changing what they were moving it relative to.

Keeping the separatrices is the other half of that, and it is what makes the
edit worth anything. A coarse edge whose endpoint moved no longer matches any
traced polyline by geometric key, so the naive answer is to densify it as a
straight chord -- which throws away precisely the field alignment the front end
exists to produce. Instead the edge is re-matched to the separatrix it came
from and that curve is WARPED onto its new endpoints
(:func:`warp_polyline`). A nudged corner costs the nudge, not the curvature.

WHAT COMES BACK FROM RHINO
--------------------------

Two forms, because both are natural to produce there and neither is more
correct than the other:

* **a mesh** -- ``compas_rhino.conversions.mesh_to_compas`` of an edited Rhino
  mesh, or a ``CoarsePseudoQuadMesh`` handed straight back;
* **closed polylines, one per patch** -- easier to edit, and the only one of the
  two in which deleting a patch or splitting one in half is a two-second
  operation.

Both reduce to the same thing: a list of faces, each a list of CORNER points,
welded by rounded coordinate. That welding is the whole topology recovery --
two patches share a corner because their corner points coincide, at the same
3-decimal precision ``from_polylines`` matches endpoints at. Nothing carries an
index across the round trip, which is what lets the face count change.

ONE POINT PER CORNER
--------------------

A patch outline must have one point per corner, not one point per sample of a
curved edge. The round trip cannot tell the difference -- a 4-sided patch whose
edges were baked as their traced separatrices comes back as a 40-gon, and
:func:`repair.solve_non_quad_faces` will dutifully fan it into 38 patches. That
is why :func:`face_polylines` bakes CORNERS only and the curved separatrices go
out separately as reference geometry. The side-count histogram in the returned
notes says immediately when this has gone wrong.

SNAPPING IS LOAD-BEARING, NOT COSMETIC
--------------------------------------

``decomposition._on_loop`` tests membership of a boundary loop at ``1e-6``. A
corner dragged in Rhino to what looks like the wall is not on the wall, so
``edges_to_curves`` finds no boundary arc for it, densifies the edge as a
chord, and the mesh quietly loses the area between chord and arc -- the same
failure that took an ellipse down to 65% coverage before ``_arc_between``
existed. So incoming corners are PROJECTED onto the wall before anything is
measured or validated.

Which corners, though, is decided by TOPOLOGY and not by distance, and
:func:`snap_to_loops` sets out why: "everything within ``snap_tol`` of a wall"
would drag interior corners that are legitimately near one -- the comb plate's
teeth are 2 units wide -- onto it, silently collapsing patches nobody touched.
A vertex on the coarse mesh's own boundary belongs on a domain loop by
definition; an interior one is never moved however close it is.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from compas.geometry import Polyline
from compas.geometry import distance_point_point
from compas.itertools import pairwise

from ..datastructures import CoarsePseudoQuadMesh

from .repair import solve_non_quad_faces


__all__ = ['coarse_from_skeleton', 'warp_polyline', 'face_polylines',
           'faces_from_geometry', 'mesh_from_faces', 'snap_to_loops',
           'closest_on_loop', 'PRECISION']


#: Decimals the corner weld rounds to. Matches ``_mesh_from_faces`` and is the
#: resolution ``from_polylines`` matches polyline endpoints at, so a layout that
#: went out through one route and came back through the other welds the same.
PRECISION = 3


# ------------------------------------------------------------------
# reading what Rhino handed back
# ------------------------------------------------------------------

def _is_mesh(thing):
    """Duck-typed: a compas mesh of any class, including the pseudo-quad ones."""
    return all(hasattr(thing, name)
               for name in ('vertices', 'faces', 'vertex_coordinates',
                            'face_vertices'))


def _as_points(thing):
    """A polyline, a compas ``Polyline``, or a bare list of points -> point list.

    The closing point of a closed polyline is dropped: Rhino writes it, a face
    must not have it, and a face that keeps it has a repeated vertex and is
    thrown away by :func:`mesh_from_faces` as degenerate.
    """
    points = [list(p)[:3] for p in getattr(thing, 'points', thing)]
    points = [p + [0.0] * (3 - len(p)) for p in points]
    if len(points) > 1 and distance_point_point(points[0], points[-1]) < 1e-9:
        points = points[:-1]
    return points


def faces_from_geometry(geometry):
    """``[[corner, corner, ...], ...]`` from a mesh or from closed polylines.

    Accepts, in order of how it is checked:

    * a mesh -- anything with ``vertices``/``faces``/``vertex_coordinates``;
    * an iterable of closed polylines, ``Polyline`` or bare point lists;
    * an iterable of faces already given as point lists, which is the same
      thing and is what makes this idempotent.
    """
    if _is_mesh(geometry):
        return [[list(geometry.vertex_coordinates(v))
                 for v in geometry.face_vertices(f)]
                for f in geometry.faces()]

    faces = []
    for item in geometry:
        points = _as_points(item)
        if len(points) >= 3:
            faces.append(points)
    return faces


# ------------------------------------------------------------------
# putting an edited corner back on the wall
# ------------------------------------------------------------------

def closest_on_loop(point, loop):
    """``(distance, projected point)`` for the closest point of a closed loop.

    The loop is treated as closed whether or not its last point repeats its
    first, matching every other loop consumer in this package.
    """
    best_d, best_p = float('inf'), None
    ring = list(loop) + list(loop[:1])
    for a, b in pairwise(ring):
        abx, aby = b[0] - a[0], b[1] - a[1]
        length2 = abx * abx + aby * aby
        if length2 == 0.0:
            continue
        t = ((point[0] - a[0]) * abx + (point[1] - a[1]) * aby) / length2
        t = max(0.0, min(1.0, t))
        q = [a[0] + abx * t, a[1] + aby * t, 0.0]
        d = distance_point_point(point, q)
        if d < best_d:
            best_d, best_p = d, q
    return best_d, best_p


def snap_to_loops(mesh, loops, tol):
    """Project the mesh's own BOUNDARY corners onto the domain walls.

    **Which corners are eligible is decided topologically, not by distance.**
    That distinction is the whole safety of this step. Snapping anything within
    ``tol`` of a wall is destructive on the domains where it matters: a patch
    corner that legitimately sits a fifth of a background spacing inside a
    narrow slot -- the comb plate has several -- would be yanked onto the wall
    it was near, silently collapsing a patch nobody touched. Whereas a vertex
    on the coarse mesh's own boundary belongs on a domain loop by definition,
    whatever its distance, so moving it there cannot destroy information; it can
    only restore it.

    So: a vertex with a naked edge gets projected, an interior one never does,
    however close it is.

    Returns
    -------
    (int, int)
        How many boundary corners were moved, and how many were further from
        every loop than ``tol`` and so left where they were. The second number
        should be 0; anything else means a patch is hanging off the domain.
    """
    if not loops or tol <= 0.0:
        return 0, 0

    moved = off = 0
    # ``vertices_on_boundarIES`` -- PLURAL. The singular form returns only the
    # LONGEST boundary, so every corner of every HOLE was reported as interior and
    # silently skipped: measured on an annulus, 12 of 16 boundary vertices seen and
    # all four of the hole's missed. It went unnoticed because the editor's own
    # preview used the same singular form, so the two agreed with each other while
    # both being wrong -- and because a generated layout already has its corners on
    # the walls, which is why no baseline row moves either way. It bites exactly
    # when a user DRAGS a hole corner, which is the case this seam exists for.
    for vkey in set(v for ring in mesh.vertices_on_boundaries() for v in ring):
        point = mesh.vertex_coordinates(vkey)
        best_d, best_p = float('inf'), None
        for loop in loops:
            d, q = closest_on_loop(point, loop)
            if d < best_d:
                best_d, best_p = d, q
        if best_p is None or best_d <= 1e-9:
            continue
        if best_d > tol:
            off += 1
            continue
        mesh.vertex_attributes(vkey, 'xyz', best_p)
        moved += 1
    return moved, off


# ------------------------------------------------------------------
# faces -> mesh
# ------------------------------------------------------------------

def _key(point):
    return (round(point[0], PRECISION), round(point[1], PRECISION))


def mesh_from_faces(faces, cls=CoarsePseudoQuadMesh):
    """A coarse mesh from faces given as lists of corner points.

    Corners are welded by rounded coordinate at :data:`PRECISION`, the same
    resolution ``from_polylines`` matches endpoints at, so two patches meeting
    at a cut share that vertex instead of each carrying its own copy.

    A face that visits one corner twice is degenerate -- a patch with no
    interior -- and is dropped rather than welded into a mesh that will not
    densify. ``None`` when nothing is left.
    """
    if not faces:
        return None, 0
    index = {}
    vertices = []
    cells = []
    dropped = 0
    for face in faces:
        cell = []
        for point in face:
            key = _key(point)
            if key not in index:
                index[key] = len(vertices)
                vertices.append([point[0], point[1], 0.0])
            cell.append(index[key])
        if len(set(cell)) == len(cell) and len(cell) >= 3:
            cells.append(cell)
        else:
            dropped += 1
    if not cells:
        return None, dropped
    return cls.from_vertices_and_faces(vertices, cells), dropped


# ------------------------------------------------------------------
# the entry point
# ------------------------------------------------------------------

def coarse_from_skeleton(geometry, loops=None, poles=(), snap_tol=0.0,
                         cls=CoarsePseudoQuadMesh):
    """**A coarse quad layout from an edited patch skeleton.**

    ``skeleton`` here means THE PATCH SKELETON HANDED BACK FROM RHINO -- the
    corner-and-edge wireframe of a coarse layout. It has nothing to do with
    ``SkeletonDecomposition`` or the medial-axis front end this package
    replaces, which is what the word means everywhere else in this repository.

    Parameters
    ----------
    geometry : mesh or list
        A compas mesh, or closed polylines with ONE POINT PER PATCH CORNER.
        See the module docstring on why the distinction matters.
    loops : list[list[[x, y, z]]], optional
        Boundary loops -- outer first, then holes. Used for snapping and handed
        to :func:`repair.solve_non_quad_faces` so a split vertex lands on the
        wall rather than on a chord.
    poles : list[[x, y, z]], optional
        Preferred pole positions. Pass the previous layout's poles so an
        untouched pseudo-quad keeps its collapsed corner where it was.
    snap_tol : float, optional
        How far a BOUNDARY corner may be projected onto a wall -- not which
        corners are eligible, which is topological. ``0`` disables snapping.
        A boundary corner further than this from every loop is left alone and
        counted in ``off_wall``, because at that distance the layout is wrong
        rather than imprecise and moving it would hide the fact.
    cls : type, optional

    Returns
    -------
    (mesh or None, dict)
        The layout, and what happened to it: ``faces_in``, ``faces_out``,
        ``vertices``, ``snapped``, ``off_wall``, ``degenerate``, ``sides`` (the
        side-count histogram BEFORE repair -- anything other than ``{4: n}``
        means the skeleton was not four-sided patches) and ``repair`` (the note
        :func:`repair.solve_non_quad_faces` returned).
    """
    faces = faces_from_geometry(geometry)
    notes = {'faces_in': len(faces), 'faces_out': 0, 'vertices': 0,
             'snapped': 0, 'off_wall': 0, 'degenerate': 0, 'sides': {},
             'repair': ''}
    if not faces:
        return None, notes

    mesh, dropped = mesh_from_faces(faces, cls)
    notes['degenerate'] = dropped
    if mesh is None:
        return None, notes

    # Snapping needs the topology -- see :func:`snap_to_loops` -- so it happens
    # after the weld and before anything is measured or validated.
    notes['snapped'], notes['off_wall'] = snap_to_loops(mesh, loops or [], snap_tol)

    sides = {}
    for fkey in mesh.faces():
        n = len(mesh.face_vertices(fkey))
        sides[n] = sides.get(n, 0) + 1
    notes['sides'] = sides

    # Whatever came back, it leaves here all-quad -- or all-quad plus registered
    # pseudo-quad poles, which is the same thing to ``densification``. Reused
    # wholesale: this is the identical repair the generated layout goes through.
    mesh, note = solve_non_quad_faces(mesh, cls, loops, poles=poles)
    notes['repair'] = note
    notes['faces_out'] = mesh.number_of_faces()
    notes['vertices'] = mesh.number_of_vertices()
    return mesh, notes


# ------------------------------------------------------------------
# moving a traced curve onto new endpoints
# ------------------------------------------------------------------

def warp_polyline(points, pa, pb, limit=1.0):
    """**End-anchored warp of a traced separatrix onto moved endpoints.**

    The reason an edited layout keeps its field alignment. ``points`` is the
    curve as traced; ``pa`` and ``pb`` are where its two ends have to be now.
    Every sample is displaced by a blend of the two end displacements, weighted
    by normalised arc length::

        d0 = pa - c[0]
        d1 = pb - c[-1]
        c'[k] = c[k] + (1 - t_k) * d0 + t_k * d1

    so the ends land exactly on ``pa`` and ``pb``, and everything between keeps
    the shape it was traced with, rigidly translated where both ends moved the
    same way and sheared where they did not. This is the standard curve warp
    (Sederberg's, in its one-dimensional case); nothing here is novel and it is
    the cheapest thing that preserves curvature exactly under translation.

    Parameters
    ----------
    points : list[[x, y, z]]
    pa, pb : [x, y, z]
    limit : float, optional
        Reject the warp when either end has to move further than ``limit``
        times the curve's own length. A displacement of that size is not this
        curve moved, it is a different curve, and warping it produces a
        self-crossing polyline that densifies into folded quads.

    Returns
    -------
    list[[x, y, z]] or None
        ``None`` when the warp was rejected, so the caller can fall back to the
        straight chord.
    """
    curve = [list(p) for p in points]
    if len(curve) < 2:
        return None

    cum = [0.0]
    for a, b in pairwise(curve):
        cum.append(cum[-1] + distance_point_point(a, b))
    total = cum[-1]
    if total <= 1e-12:
        return None

    d0 = [pa[i] - curve[0][i] for i in range(3)]
    d1 = [pb[i] - curve[-1][i] for i in range(3)]
    reach = max((d0[0] ** 2 + d0[1] ** 2) ** 0.5, (d1[0] ** 2 + d1[1] ** 2) ** 0.5)
    if reach > limit * total:
        return None

    out = []
    for point, s in zip(curve, cum):
        t = s / total
        out.append([point[i] + (1.0 - t) * d0[i] + t * d1[i] for i in range(3)])
    out[0] = [pa[0], pa[1], pa[2] if len(pa) > 2 else 0.0]
    out[-1] = [pb[0], pb[1], pb[2] if len(pb) > 2 else 0.0]
    return out


# ------------------------------------------------------------------
# what to bake into Rhino
# ------------------------------------------------------------------

def face_polylines(coarse):
    """**One closed polyline per patch: the thing to bake and edit.**

    CORNERS ONLY -- four points and the closing repeat, or three for a
    pseudo-quad. Deliberately not the curved separatrix geometry, however much
    better that looks: the round trip recovers a patch's corners from the
    polyline's points, so a polyline carrying forty samples of a traced curve
    comes back as a forty-sided patch. Bake
    ``decomposition.decomposition_polylines()`` on a separate reference layer
    if the curvature needs to be visible while editing -- it is drawn, not
    read.

    Returns
    -------
    list[:class:`compas.geometry.Polyline`]
    """
    out = []
    for fkey in coarse.faces():
        corners = [list(coarse.vertex_coordinates(v))
                   for v in coarse.face_vertices(fkey)]
        out.append(Polyline(corners + corners[:1]))
    return out
