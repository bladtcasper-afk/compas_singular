"""EXAMPLE 8 -- the method itself, one panel per step.

Every other example shows RESULTS: what decomposes, how it compares, what the
numbers are. None of them shows the METHOD. ``quad_mesh`` is deliberately one
call that owns the whole chain, which is what makes it usable and also what
makes it opaque -- an outline goes in and a quad mesh comes out, with the six
things that happened in between visible only as attributes on the object.

This lays those six out side by side, left to right, for eight domains:

    1  domain          the outline as handed in -- nothing computed yet
    2  background      the triangulation the field will live on
    3  cross field     one cross per vertex, plus the singularities it forces
    4  separatrices    field lines traced out of the singularities and corners
    5  coarse layout   the patches those lines cut the domain into
    6  quad mesh       each patch densified WITH the field steering its interior

Read a row left to right and it is the whole front end. Read a column top to
bottom and it is one step across eight shapes, which is where the differences
live -- the square has no singularity at all and one patch; the cross-shaped
plate's four reflex corners each launch a separatrix; the stadium's curved end
has to be approximated by a field that only knows about the triangles under it.

The eight domains, chosen because between them they cover every mechanism:

    square       the trivial case -- no singularity, no reflex corner, ONE patch
    cross        four reflex corners, so four separatrices and a real layout
    comb-plate   narrow teeth: the case where a tolerance scaled to the whole
                 domain used to swallow a 2-unit feature (see repair.SHARP_TURN)
    pentagon     odd corner count -- the field MUST put a singularity inside
    stadium      a curved wall, which is what needed arrangement.py to work
    annulus      the only domain here with a HOLE -- a second wall the field has
                 to close around, and a layout that has to wrap a ring rather
                 than fill a disc
    square+circle  and
    circle+square  the two MIXED-curvature annuli: a round hole in a straight
                 domain, and a square hole in a round one. The annulus above has
                 both walls curved, so its field is smooth everywhere and every
                 singularity it has comes from the topology. These two put a
                 corner on one wall and none on the other, so the field has to
                 reconcile the two -- the four corners of the square hole in
                 ``circle+square`` each want a separatrix in a domain whose
                 outer wall offers nothing to end it on

Outlines are IMPORTED from ``15_baseline.py`` rather than retyped, so these are
the same polygons the locked baseline measures -- see ``17_vs_skeleton.py`` for
why that matters (the obvious guess at a sampled outline is a different polygon
that measures differently and convincingly).

Nothing here is a check: no assertions, no pass/fail. It is a picture of the
pipeline, with the counts at each step printed beside it.

Run:
    python 18_pipeline.py
    python 18_pipeline.py --no-view
    python 18_pipeline.py pentagon stadium      # a subset, same layout
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from compas_singular.framefield.decomposition import FieldDecomposition  # noqa: E402


#: The outlines, imported so there is exactly one definition of each shape.
#: ``15_baseline`` is not a legal module name and is loaded by path; it guards
#: its CLI behind ``__main__``, so importing it runs nothing.
_spec = importlib.util.spec_from_file_location(
    'baseline_suite', os.path.join(HERE, '15_baseline.py'))
suite = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(suite)


#: Background triangulation spacing -- NOT the quad size. This is the resolution
#: the FIELD is solved at, and it is what steps 2 to 4 are pictures of.
SPACING = 0.5

#: Target quad edge length, applied to the coarse layout in step 6.
QUAD = 1

#: ``(label, outer, holes, guides, extra kwargs)``. ``cross`` is the suite's
#: ``plus-plate`` and ``annulus`` its ``disc+round hole`` -- same polygons, the
#: names they are asked for by; the last two keep the suite's own names.
#: ``ROUND_HOLE`` is wound the OTHER way round from the outer wall and
#: ``SQUARE_HOLE`` the SAME way, and both work: the background triangulation
#: classifies points by ``is_point_in_polygon_xy``, which does not care.
DOMAINS = [
    ('square',        suite.SQUARE,   None, None, {}),
    ('cross',         suite.PLUS,     None, None, {}),
    ('comb-plate',    suite.COMB,     None, None, {}),
    ('pentagon',      suite.PENTAGON, None, None, {}),
    ('stadium',       suite.STADIUM,  None, None, {}),
    ('annulus',       suite.DISC,     [suite.ROUND_HOLE],  None, {}),
    ('square+circle', suite.SQUARE,   [suite.ROUND_HOLE],  None, {}),
    ('circle+square', suite.DISC,     [suite.SQUARE_HOLE], None, {}),
    ('ellipse',       suite.ELLIPSE,  None, None, {}),
    ('ellipse+hole',       suite.ELLIPSE,  [suite.ROUND_HOLE], None, {}),
    ('decomposition',  suite.DECOMPOSITION_OUTER, [suite.DECOMPOSITION_HOLES], None, {}),
    ('RV',              suite.RV_OUTER, suite.RV_HOLES, None, {}),
    ('RV+hole',         suite.RV_WITH_HOLE_OUTER, suite.RV_WITH_HOLE_HOLES, None, {}),
]

#: One column per step, in the order they happen. The scene tree names match.
STAGES = [
    ('domain', 'the outline as handed in'),
    ('background', 'triangulation the field lives on'),
    ('cross field', 'one cross per vertex + singularities'),
    ('separatrices', 'field lines traced out of the singularities'),
    ('coarse layout', 'the patches those lines cut'),
    ('quad mesh', 'densified, field steering patch interiors'),
]

#: Panel pitch. The widest domain here is 15 units (pentagon, stadium), so this
#: leaves a clear gap between panels at every step.
PITCH = 40.0


def bbox_centre(loop):
    """Centre of a loop's bounding box, used to sit each shape in its panel.

    The domains are different sizes and none of them is centred on its own
    origin -- the pentagon starts at x = -2 -- so placing them all at the raw
    cell offset puts a lopsided grid in front of the camera.
    """
    xs = [p[0] for p in loop]
    ys = [p[1] for p in loop]
    return (min(xs) + max(xs)) * 0.5, (min(ys) + max(ys)) * 0.5


def run(label, outer, holes, guides, kwargs):
    """Everything the six panels need, in one pass over the pipeline.

    Steps 1 to 4 all come out of ``from_boundary`` -- it solves the field and
    traces the separatrices -- and steps 5 and 6 out of ``quad_mesh``. Nothing
    is recomputed for the picture: the network is cached on the decomposition
    and the viz helpers read the same objects the measurements below do.
    """
    d = FieldDecomposition.from_boundary(
        outer, inner_boundaries=holes, guides=guides,
        target_length=SPACING, **(kwargs or {}))

    # steps 5 and 6. ``quad_mesh`` leaves the coarse layout on ``d.mesh`` and
    # the densification on ``d.dense``; on a fallback route ``d.mesh`` IS the
    # fallback, which is the honest thing to draw, and ``route`` says so.
    dense = d.quad_mesh(target_length=QUAD)
    return d, dense


def report(label, outer, d, dense):
    """The counts at each step, printed in the same order as the panels."""
    _, others, _ = d._build()
    singular = list(d.field.singularities())
    positive = sum(1 for _f, k in singular if k > 0)
    negative = len(singular) - positive
    q = d.quality()

    print()
    print('=== {} '.format(label) + '=' * max(0, 66 - len(label)))
    # the INPUT count, not ``background.outer`` -- that loop has already been
    # resampled to the background spacing, which is step 2's doing, not step 1's
    print('  1 domain        {} outline points'.format(len(outer)))
    print('  2 background    {} vertices, {} triangles at spacing {} '
          '(wall resampled to {} points)'.format(
              d.background.mesh.number_of_vertices(),
              d.background.mesh.number_of_faces(), SPACING,
              len(d.background.outer)))
    print('  3 cross field   {} singularit{} ({} wanting valence 3, {} valence 5)'
          .format(len(singular), 'y' if len(singular) == 1 else 'ies',
                  positive, negative))
    print('  4 separatrices  {} interior line{}'.format(
        len(others), '' if len(others) == 1 else 's'))
    print('  5 coarse layout {} patch{} (route {!r})'.format(
        d.mesh.number_of_faces(),
        '' if d.mesh.number_of_faces() == 1 else 'es', q['route']))
    print('  6 quad mesh     {} faces, angles {:.2f} to {:.2f} deg, aspect {:.2f}, '
          'coverage {:.0%}'.format(q['faces'], q['min_angle'], q['max_angle'],
                                   q['aspect_max'], q['coverage']))
    stats = d.densify_stats or {}
    if stats:
        print('                  field-integrated interiors: {} of {} patches'
              .format(stats.get('patches', 0) - stats.get('guarded', 0),
                      stats.get('patches', 0)))
    for warning in d.warnings():
        print('  ! {}'.format(warning))


def main(argv):
    show = '--no-view' not in argv
    wanted = [a for a in argv[1:] if not a.startswith('-')]

    print('The frame-field front end, step by step. Background spacing {}, '
          'quad target {}.'.format(SPACING, QUAD))
    for i, (name, what) in enumerate(STAGES):
        print('  column {}  {:14s} {}'.format(i + 1, name, what))

    scenes = []
    for label, outer, holes, guides, kwargs in DOMAINS:
        if wanted and label not in wanted:
            continue
        try:
            d, dense = run(label, outer, holes, guides, kwargs)
        except Exception as exc:                              # noqa: BLE001
            print()
            print('=== {} -- RAISED {}: {}'.format(label, type(exc).__name__, exc))
            continue
        report(label, outer, d, dense)
        scenes.append((label, outer, d, dense))

    if show and scenes:
        view(scenes)
    return 0


def view(scenes):
    """One row per domain, one panel per step, in the order they happen.

    Each domain is a group in the scene tree with the six steps nested inside
    it, so a step can be switched off across the whole row -- turning the dense
    mesh off everywhere leaves the five layouts, which is the comparison worth
    making.
    """
    from compas.geometry import Polyline
    from compas_viewer import Viewer                                # noqa: E402

    from compas_singular.framefield.viz import Grid               # noqa: E402
    from compas_singular.framefield.viz import WALL               # noqa: E402
    from compas_singular.framefield.viz import add_background     # noqa: E402
    from compas_singular.framefield.viz import add_boundary       # noqa: E402
    from compas_singular.framefield.viz import add_dense          # noqa: E402
    from compas_singular.framefield.viz import add_field          # noqa: E402
    from compas_singular.framefield.viz import add_guides         # noqa: E402
    from compas_singular.framefield.viz import add_layout         # noqa: E402
    from compas_singular.framefield.viz import add_separatrices   # noqa: E402
    from compas_singular.framefield.viz import add_singularities  # noqa: E402

    viewer = Viewer()
    cols = len(STAGES)
    grid = Grid(pitch=PITCH, cols=cols)

    for row, (label, outer, d, dense) in enumerate(scenes):
        group = viewer.scene.add_group(name=label)
        cx, cy = bbox_centre(outer)

        def panel(step):
            """Offset of this domain's ``step``-th panel, shape centred in it.

            ``Grid.frame`` aims the camera assuming a shape sits at roughly
            ``0.35 * pitch`` from its cell origin, so centre on that point and
            the framing stays correct for any set of shapes.
            """
            dx, dy = grid.cell(row * cols + step)
            return dx + 0.35 * PITCH - cx, dy + 0.35 * PITCH - cy

        # 1 -- the outline, before anything has been computed from it. Drawn
        # from the INPUT points, not background.outer, which is already the
        # resampled version step 2 needs.
        dx, dy = panel(0)
        sub = viewer.scene.add_group(name='1 domain', parent=group)
        loop = [[p[0] + dx, p[1] + dy, 0.0] for p in outer]
        sub.add(Polyline(loop + loop[:1]), linecolor=WALL, linewidth=3,
                show_points=True, pointcolor=WALL, pointsize=12,
                name='outline ({} points)'.format(len(outer)))
        add_guides(sub, d.guides, dx, dy)

        # 2 -- the background triangulation
        dx, dy = panel(1)
        sub = viewer.scene.add_group(name='2 background', parent=group)
        add_background(sub, d, dx, dy)
        add_boundary(sub, d, dx, dy)

        # 3 -- the field, and the singularities its topology forces
        dx, dy = panel(2)
        sub = viewer.scene.add_group(name='3 cross field', parent=group)
        # a lower tick budget than the default 140: this panel is a picture of a
        # DIRECTION field, and at the zoom that fits thirty panels on screen the
        # default thins to something that reads as hatching. ``add_field`` grows
        # the tick length as it thins, so fewer ticks stay legible zoomed in too
        add_field(sub, d, dx, dy, budget=300)
        add_singularities(sub, d, dx, dy)
        add_boundary(sub, d, dx, dy)

        # 4 -- the separatrices traced through it
        dx, dy = panel(3)
        sub = viewer.scene.add_group(name='4 separatrices', parent=group)
        add_separatrices(sub, d, dx, dy)
        add_singularities(sub, d, dx, dy)
        add_boundary(sub, d, dx, dy)

        # 5 -- the coarse layout they cut
        dx, dy = panel(4)
        sub = viewer.scene.add_group(name='5 coarse layout', parent=group)
        add_layout(sub, d, dx, dy, faces=True)

        # 6 -- the deliverable
        dx, dy = panel(5)
        sub = viewer.scene.add_group(name='6 quad mesh', parent=group)
        add_dense(sub, dense, dx, dy)
        add_boundary(sub, d, dx, dy)
        add_guides(sub, d.guides, dx, dy)

    grid.frame(viewer)
    viewer.show()


if __name__ == '__main__':
    sys.exit(main(sys.argv))
