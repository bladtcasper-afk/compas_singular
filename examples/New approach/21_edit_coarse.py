"""EXAMPLE 9 -- **editing the coarse layout by hand, between the two halves.**

``quad_mesh`` owns the whole chain in one call, which is what makes it usable
and also what makes it closed. The coarse layout is the one artefact a designer
actually wants to touch -- drag a corner the tracer put somewhere awkward,
delete a patch, split one in two -- and until now there was nowhere to reach in.

This opens exactly one seam::

    d = FieldDecomposition.from_boundary(outer, holes, guides=cables, target_length=0.5)
    coarse = d.decomposition_mesh()        # <- generated layout, bake it to Rhino
    ...                                    # <- edit it there
    d.edit_coarse(edited)                  # <- hand it back
    dense = d.quad_mesh(target_length=1.0) # <- densified from YOUR layout

**The field is not re-solved and not touched.** The edit enters at ``d.mesh``
and nowhere else. That is checked here rather than asserted: panel 1 compares
``field.theta`` vertex by vertex across the edit and prints the largest
difference, which is 0.

THE EDITS ARE SCRIPTED, AND THAT IS THE POINT
---------------------------------------------

Nothing below opens Rhino. Each edit is applied in-process to the baked
polylines by :func:`drag` / :func:`slide` / ``del``, which is exactly what
Rhino would hand back -- a list of closed polylines with moved points or a
missing patch. That is what makes this a CHECK you can run rather than a
screenshot of somebody's session, and it is the only reason the numbers below
mean anything.

Four edits, chosen because between them they exercise every branch:

    identity     bake it and hand it straight back, unchanged
    drag         move the interior singularity corner 1.8 units
    slide        move a boundary corner along the wall, 0.05 sloppily off it
    delete       remove a patch entirely

WHAT A DRAGGED CORNER COSTS -- THE NUDGE, NOT THE CURVATURE
-----------------------------------------------------------

The interesting failure this avoids is silent. A coarse edge whose endpoint
moved no longer matches any traced separatrix by geometric key, so the obvious
implementation densifies it as a straight chord -- and a straight chord is
precisely the alignment the front end exists to produce, thrown away. Instead
the edge is re-matched to the separatrix it came from, anchored on the end that
did NOT move, and that curve is warped onto the new corner
(``edit.warp_polyline``).

Panel 3 measures it. ``sagitta`` is a curve's largest deviation from its own
chord -- 0 for a straight line -- so comparing the warped edge's sagitta with
the separatrix it came from says directly whether the curvature survived.

THE RHINO RECIPE
----------------

The code runs inside Rhino, so the layout does not travel through a file. Two
runs against the same document::

    # run 1 -- generate and bake
    d = FieldDecomposition.from_boundary(outer, holes, target_length=0.5)
    coarse = d.decomposition_mesh()
    for pl in face_polylines(coarse):          # ONE POINT PER CORNER
        sc.doc.Objects.AddPolyline([Rhino.Geometry.Point3d(*p) for p in pl])
    for pl in d.decomposition_polylines():     # reference layer, drawn not read
        ...

    # ... edit the patch polylines in Rhino ...

    # run 2 -- hand them back
    d = FieldDecomposition.from_boundary(outer, holes, target_length=0.5)
    d.edit_coarse([[list(p) for p in pl.ToArray()] for pl in edited])
    dense = d.quad_mesh(target_length=1.0)

Re-solving in run 2 is safe: the pipeline is deterministic by construction --
``background.py``'s nudge is a *fixed* pseudo-random one, precisely so this
holds -- so run 2's field is run 1's field and the edited layout is measured
against the separatrices it was generated from. ``edit_coarse`` also accepts a
mesh (``compas_rhino.conversions.mesh_to_compas``) if that is easier.

BAKE CORNERS, NOT CURVES
------------------------

``face_polylines`` gives four points and the closing repeat, not the traced
separatrix geometry, however much better the latter looks. The round trip
recovers a patch's corners from its polyline's points, so a patch baked as its
curved edges comes back as a 40-gon and ``solve_non_quad_faces`` fans it into
38 patches. Bake ``decomposition_polylines()`` on a separate layer if the
curvature needs to be visible while editing -- it is drawn, not read.

Run:
    python 21_edit_coarse.py
    python 21_edit_coarse.py --no-view
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from compas.geometry import distance_point_point                   # noqa: E402
from compas.itertools import pairwise                              # noqa: E402

from compas_singular.framefield.decomposition import FieldDecomposition  # noqa: E402
from compas_singular.framefield.edit import face_polylines               # noqa: E402


_spec = importlib.util.spec_from_file_location(
    'baseline_suite', os.path.join(HERE, '15_baseline.py'))
suite = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(suite)


#: Background spacing -- the resolution the FIELD is solved at.
SPACING = 0.6

#: Target quad edge length, applied to the coarse layout.
QUAD = 1.0

#: Panel pitch for the viewer grid.
PITCH = 22.0


# ------------------------------------------------------------------
# standing in for Rhino
# ------------------------------------------------------------------

def bake(coarse):
    """What run 1 would put in the Rhino document: one closed polyline per patch."""
    return [[list(p) for p in polyline] for polyline in face_polylines(coarse)]


def drag(polylines, frm, to, tol=1e-6):
    """Move every occurrence of a corner. What dragging one vertex in Rhino does.

    Every occurrence, because a corner shared by four patches appears in four
    polylines and Rhino's own point editing moves them together. Moving one copy
    and not the others would tear the layout apart at that corner, which the
    weld would then report as four separate vertices.
    """
    moved = 0
    for polyline in polylines:
        for point in polyline:
            if (abs(point[0] - frm[0]) < tol and abs(point[1] - frm[1]) < tol):
                point[0], point[1] = to[0], to[1]
                moved += 1
    return moved


def interior_corner(coarse):
    """A patch corner that is not on a wall -- a singularity, in these domains."""
    for vkey in coarse.vertices():
        if not coarse.is_vertex_on_boundary(vkey):
            return coarse.vertex_coordinates(vkey)
    return None


def along_loop(loop, point, distance):
    """Walk ``distance`` along a closed loop from the point nearest ``point``.

    What sliding a corner in Rhino with the boundary curve snapped actually
    does, and NOT the same as adding an offset in x. On a curved wall the two
    diverge fast: sliding 0.8 in +x from a point on the suite's disc lands 0.58
    OUTSIDE it, which is past any sane snap tolerance and is correctly refused
    as a corner hanging off the domain. Measuring the snap needs an edit that a
    user would plausibly make.
    """
    ring = list(loop) + list(loop[:1])
    cum = [0.0]
    for a, b in pairwise(ring):
        cum.append(cum[-1] + distance_point_point(a, b))
    total = cum[-1]

    # arc-length parameter of the point on the loop nearest ``point``
    best_d, best_s = float('inf'), 0.0
    for i, (a, b) in enumerate(pairwise(ring)):
        abx, aby = b[0] - a[0], b[1] - a[1]
        length2 = abx * abx + aby * aby
        if length2 == 0.0:
            continue
        t = max(0.0, min(1.0, ((point[0] - a[0]) * abx
                               + (point[1] - a[1]) * aby) / length2))
        q = [a[0] + abx * t, a[1] + aby * t, 0.0]
        d = distance_point_point(point, q)
        if d < best_d:
            best_d, best_s = d, cum[i] + t * (cum[i + 1] - cum[i])

    s = (best_s + distance) % total
    for i, (c0, c1) in enumerate(zip(cum, cum[1:])):
        if c0 <= s <= c1:
            a, b = ring[i], ring[i + 1]
            t = 0.0 if c1 == c0 else (s - c0) / (c1 - c0)
            return [a[k] + (b[k] - a[k]) * t for k in range(3)]
    return list(ring[-1])


def wall_corner(coarse, outline):
    """A corner ON a wall but not at a corner OF the domain.

    A separatrix landing point, in other words -- the kind that can be slid
    along the wall without changing the shape of the domain.
    """
    for vkey in coarse.vertices():
        if not coarse.is_vertex_on_boundary(vkey):
            continue
        p = coarse.vertex_coordinates(vkey)
        if not any(distance_point_point(p, c) < 1e-6 for c in outline):
            return p
    return None


# ------------------------------------------------------------------
# measuring
# ------------------------------------------------------------------

def sagitta(points):
    """A polyline's largest deviation from its own chord. 0 for a straight line.

    The measure of "did the curvature survive", because it is exactly what a
    chord throws away and exactly what a warp is supposed to keep.
    """
    a, b = points[0], points[-1]
    abx, aby = b[0] - a[0], b[1] - a[1]
    length = (abx * abx + aby * aby) ** 0.5
    if length < 1e-12:
        return 0.0
    worst = 0.0
    for p in points:
        # perpendicular distance from the chord
        cross = abs(abx * (p[1] - a[1]) - aby * (p[0] - a[0])) / length
        worst = max(worst, cross)
    return worst


def polyline_length(points):
    return sum(distance_point_point(a, b) for a, b in pairwise(points))


def nearest_vertex_distance(mesh, point):
    """How close the dense mesh comes to a point -- 0 if it has a vertex there."""
    return min(distance_point_point(mesh.vertex_coordinates(v), point)
               for v in mesh.vertices())


def theta_snapshot(field):
    """The field itself, as a plain dict, for a before/after comparison."""
    return dict(field.theta)


def theta_difference(before, after):
    """Largest per-vertex change in the field. Must be exactly 0."""
    if set(before) != set(after):
        return float('inf')
    return max((abs(after[k] - before[k]) for k in before), default=0.0)


# ------------------------------------------------------------------
# the four edits
# ------------------------------------------------------------------

def fresh(outer, holes=None):
    d = FieldDecomposition.from_boundary(
        outer, inner_boundaries=holes, target_length=SPACING)
    d.decomposition_mesh()
    return d


def edit_identity(d, outer):
    """Bake it and hand it straight back. Nothing may change."""
    return bake(d.mesh), 'baked and returned unchanged'


def edit_drag(d, outer):
    """Move the interior singularity corner. The case the warp exists for."""
    polylines = bake(d.mesh)
    corner = interior_corner(d.mesh)
    if corner is None:
        return None, 'no interior corner to drag'
    target = [corner[0] + 1.5, corner[1] - 1.0, 0.0]
    n = drag(polylines, corner, target)
    return polylines, 'interior corner {} -> {} ({} patch corner(s))'.format(
        _fmt(corner), _fmt(target), n)


def edit_slide(d, outer):
    """Slide a boundary corner 0.8 ALONG the wall, landing 0.05 off it.

    The sloppiness is deliberate: it is what an eyeballed drag in Rhino
    produces, and it is what the snap exists to absorb. Without the snap the
    corner is not on the wall, ``edges_to_curves`` finds no boundary arc for it,
    and the domain quietly loses the area between chord and arc.

    Along the wall rather than along x, because on a curved domain those are
    different edits -- see :func:`along_loop`.
    """
    polylines = bake(d.mesh)
    corner = wall_corner(d.mesh, outer)
    if corner is None:
        return None, 'no slidable wall corner'
    slid = along_loop(d.background.outer, corner, 0.8)
    target = [slid[0], slid[1] + 0.05, 0.0]
    n = drag(polylines, corner, target)
    return polylines, ('wall corner {} -> {}, slid 0.8 along the wall and left '
                       '0.05 off it ({} corner(s))'.format(
                           _fmt(corner), _fmt(target), n))


def edit_delete(d, outer):
    """Remove a patch. A topology change, and a deliberate void."""
    polylines = bake(d.mesh)
    if len(polylines) < 3:
        return None, 'too few patches to delete one'
    index = len(polylines) // 2
    del polylines[index]
    return polylines, 'patch {} of {} deleted'.format(index, len(polylines) + 1)


EDITS = [
    ('identity', edit_identity),
    ('drag', edit_drag),
    ('slide', edit_slide),
    ('delete', edit_delete),
]

#: Three domains, for three different reasons. ``pentagon`` has one interior
#: singularity and straight walls, so a drag is easy to read. ``disc`` has four
#: interior singularities and a CURVED wall, which is what puts the boundary-arc
#: branch under load -- slide a corner there and the edge either follows the arc
#: or the mesh loses area. ``plus-plate`` has no interior corner at all and its
#: layout is five patches meeting in a cross, so deleting the middle one leaves
#: four patches touching at corners only: the natural way to make a
#: non-manifold layout, and the natural place to see the refusal.
DOMAINS = [
    ('pentagon', suite.PENTAGON, None),
    ('disc', suite.DISC, None),
    ('plus-plate', suite.PLUS, None),
]


def _fmt(point):
    return '({:.2f}, {:.2f})'.format(point[0], point[1])


# ------------------------------------------------------------------

def run(label, outer, holes):
    """Every edit against one domain, each from its own fresh decomposition."""
    base = fresh(outer, holes)
    reference = base.quad_mesh(target_length=QUAD)

    print()
    print('=== {} '.format(label) + '=' * max(0, 66 - len(label)))
    print('  generated   {} patches, {} dense faces, route {!r}'.format(
        base.mesh.number_of_faces(), reference.number_of_faces(), base.route()))

    results = [('generated', base, reference, '')]
    for name, apply_edit in EDITS:
        d = fresh(outer, holes)
        before = theta_snapshot(d.field)
        tracer, separatrices = d.tracer, d.separatrices

        polylines, what = apply_edit(d, outer)
        if polylines is None:
            print('  {:11s} skipped -- {}'.format(name, what))
            continue

        try:
            d.edit_coarse(polylines)
        except ValueError as exc:
            print('  {:11s} REFUSED -- {}'.format(name, exc))
            continue
        dense = d.quad_mesh(target_length=QUAD)

        moved = theta_difference(before, theta_snapshot(d.field))
        q = d.quality()
        edges = d.edit_notes.get('edges', {})

        print('  {:11s} {}'.format(name, what))
        print('              {} patches -> {} dense faces, route {!r}, '
              'coverage {:.0%}'.format(
                  d.mesh.number_of_faces(), dense.number_of_faces(),
                  d.route(), q['coverage']))
        print('              angles {:.2f} to {:.2f} deg, aspect {:.2f}'.format(
            q['min_angle'], q['max_angle'], q['aspect_max']))
        print('              edges: {} exact, {} boundary arc, {} warped, '
              '{} chord'.format(edges.get('exact', 0), edges.get('boundary', 0),
                                edges.get('warped', 0), edges.get('chord', 0)))
        print('              FIELD UNCHANGED: max |dtheta| = {:.3e}, '
              'tracer {}, separatrices {}'.format(
                  moved,
                  'same object' if d.tracer is tracer else 'REPLACED',
                  'same object' if d.separatrices is separatrices else 'REPLACED'))
        if d.edit_notes.get('snapped'):
            print('              {} corner(s) snapped onto a wall'.format(
                d.edit_notes['snapped']))

        if name == 'identity':
            same = (d.mesh.number_of_faces() == base.mesh.number_of_faces()
                    and dense.number_of_faces() == reference.number_of_faces())
            print('              IDENTITY: {} ({} vs {} dense faces)'.format(
                'unchanged' if same else 'CHANGED',
                dense.number_of_faces(), reference.number_of_faces()))

        if name == 'drag':
            corner = interior_corner(d.mesh)
            if corner is not None:
                gap = nearest_vertex_distance(dense, corner)
                print('              the dense mesh reaches the moved corner: '
                      '{:.2e} away'.format(gap))
            report_curvature(d, base)

        for note in d.warnings():
            print('              ! {}'.format(note))

        results.append((name, d, dense, what))
    return results


def incident_sagittas(d, corner):
    """Sagitta of every coarse edge meeting ``corner``, as it will densify.

    THE SAME EDGES before and after the drag, which is the only comparison that
    means anything. Comparing the warped edges against the average of all the
    traced polylines does not: half of those are boundary arcs of a different
    length and curvature, and the mean says more about the domain than about
    the warp.
    """
    out = []
    for (u, v), curve in d.edges_to_curves().items():
        pa = d.mesh.vertex_coordinates(u)
        pb = d.mesh.vertex_coordinates(v)
        if (distance_point_point(pa, corner) < 1e-6
                or distance_point_point(pb, corner) < 1e-6):
            out.append(sagitta(curve))
    return out


def report_curvature(d, base):
    """**Did the drag cost the nudge, or the curvature?**

    The five edges around the interior corner, measured before the drag on the
    generated layout and after it on the edited one. If the warp is doing its
    job the two are of the same order; if those edges had fallen back to straight
    chords the second number would be 0.000 exactly, and that is the failure
    this whole branch exists to prevent.
    """
    was = incident_sagittas(base, interior_corner(base.mesh))
    now = incident_sagittas(d, interior_corner(d.mesh))
    if not was or not now:
        print('              (no edge at the corner to measure)')
        return
    print('              curvature at that corner -- mean sagitta of its {} '
          'edges: {:.3f} generated -> {:.3f} after the drag ({} straight). '
          'A chord would be 0.000.'.format(
              len(now), sum(was) / len(was), sum(now) / len(now),
              sum(1 for s in now if s < 1e-9)))


def refusal_check():
    """A layout that will not densify must RAISE, not fall back silently.

    ``quad_mesh`` promises to always return a mesh, and honouring that promise
    on a hand edit would mean discarding the user's work and handing back a
    triangulation that looks nothing like what they drew -- at the one moment
    they are still in Rhino and could fix it. So this is the deliberate
    exception, and it is checked.
    """
    print()
    print('=== a bad edit is refused, not quietly replaced ' + '=' * 21)

    # A BOWTIE: one patch's corners reordered so its outline crosses itself.
    # Chosen because it is the failure a hand edit actually produces -- dragging
    # a corner past its diagonal neighbour -- and because it is a failure of the
    # PATCH rather than of the collection, so nothing about it can be repaired
    # away. Note that adding a stray extra patch is NOT such a case: it fans and
    # is absorbed, which is the repair working correctly, not a hole in the
    # check.
    d = fresh(suite.PENTAGON)
    polylines = bake(d.mesh)
    quad = polylines[0]
    polylines[0] = [quad[0], quad[2], quad[1], quad[3], quad[0]]
    try:
        d.edit_coarse(polylines)
    except ValueError as exc:
        print('  bowtied one patch -> RAISED as it should:')
        print('    {}'.format(exc))
        print('  The generated layout is untouched and the field was never '
              'consulted; fix the patch in Rhino and hand it back.')
        return
    print('  bowtied one patch -> accepted, {} patches, route {!r}. That is a '
          'HOLE IN THE CHECK, not a success.'.format(
              d.mesh.number_of_faces(), d.route()))


def main(argv):
    show = '--no-view' not in argv
    wanted = [a for a in argv[1:] if not a.startswith('-')]

    print('Hand-edited coarse layouts. Background spacing {}, quad target {}.'
          .format(SPACING, QUAD))
    print('Each edit starts from its OWN fresh decomposition, so the field it '
          'is measured against is the field that generated the layout.')

    scenes = []
    for label, outer, holes in DOMAINS:
        if wanted and label not in wanted:
            continue
        scenes.append((label, outer, run(label, outer, holes)))

    refusal_check()

    if show and scenes:
        view(scenes)
    return 0


def view(scenes):
    """One row per domain: the generated result, then one panel per edit.

    Each panel carries the dense mesh with its coarse layout drawn over it, so
    what the edit did to the LAYOUT and what it did to the ELEMENTS are visible
    together -- which is the whole question.
    """
    from compas_viewer import Viewer

    from compas_singular.framefield.viz import Grid
    from compas_singular.framefield.viz import add_boundary
    from compas_singular.framefield.viz import add_dense
    from compas_singular.framefield.viz import add_layout

    viewer = Viewer()
    cols = 1 + len(EDITS)
    grid = Grid(pitch=PITCH, cols=cols)

    for row, (label, outer, results) in enumerate(scenes):
        group = viewer.scene.add_group(name=label)
        xs = [p[0] for p in outer]
        ys = [p[1] for p in outer]
        cx, cy = (min(xs) + max(xs)) * 0.5, (min(ys) + max(ys)) * 0.5

        for col, (name, d, dense, what) in enumerate(results):
            dx, dy = grid.cell(row * cols + col)
            dx += 0.35 * PITCH - cx
            dy += 0.35 * PITCH - cy
            sub = viewer.scene.add_group(
                name='{} {} -- {} faces'.format(col + 1, name,
                                                dense.number_of_faces()),
                parent=group)
            add_dense(sub, dense, dx, dy)
            add_layout(sub, d, dx, dy, faces=False)
            add_boundary(sub, d, dx, dy)

    grid.frame(viewer)
    viewer.show()


if __name__ == '__main__':
    sys.exit(main(sys.argv))
