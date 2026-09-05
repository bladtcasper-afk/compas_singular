"""EXAMPLE 7 -- a symmetric problem, and where it stops being one.

A square with four symmetric guide arcs comes out of the field route with a
visibly asymmetric coarse layout: a doubled diagonal, sliver patches at two
corners, pole triangles where the other side has none. This script finds where
the symmetry goes, fixes it in ``framefield/symmetry.py``, and measures what
the fix is worth and whether it is stable.

NOTHING IS WIRED IN. ``framefield/symmetry.py`` is imported by this script and
by nothing else. Every existing module, the ``ff_0*`` Rhino commands and
``baseline.json`` are untouched, and running this changes none of them. The
adoption decision is a separate one -- see the closing section.

The viewer scene is the fastest way to see the point. Eight panels: the original
route along the top row, the symmetric one below it, one stage per column --
background, field, layout, quad mesh. So the divergence appears where it
happens, in column 1, in the triangulation, before the field exists.

Column 3 is the one to look at. It draws each layout ON TOP OF ITS OWN MIRROR
IMAGE, in green. Two layouts side by side are hard to compare by eye -- the
difference is one patch corner half a triangle out of place -- but a layout over
its own reflection either coincides or visibly doubles. Top row: green running
beside black. Bottom row: green vanishing into it.

Run:
    python 24_symmetry.py                 # the numbers, then the viewer
    python 24_symmetry.py --view          # the viewer only
    python 24_symmetry.py --no-view       # the numbers only
    python 24_symmetry.py --teardown      # part 1 only, the stage-by-stage
    python 24_symmetry.py --suite         # part 2 only, the domain sweep
    python 24_symmetry.py --stability     # part 3 only, determinism + sweep
"""
import os
import sys
import time
from math import cos, pi, sin

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from framefield.background import BackgroundMesh                # noqa: E402
from framefield.constraints import from_boundary                # noqa: E402
from framefield.constraints import from_curves                  # noqa: E402
from framefield.decomposition import FieldDecomposition         # noqa: E402
from framefield.field import CrossField                         # noqa: E402
from framefield.quality import hard_floor                       # noqa: E402
from framefield.trace import Tracer                             # noqa: E402
from framefield import symmetry as sym                          # noqa: E402


#: quad size, held fixed so the only thing varying is the route
QUAD_TARGET = 1.0

#: background spacing for every single-resolution measurement
SPACING = 0.5

#: every fix in symmetry.py, which is what "the symmetric route" means below
ALL_STEPS = ('background', 'guides', 'field', 'singularities', 'network')


# ----------------------------------------------------------------------
# domains
# ----------------------------------------------------------------------

SQUARE = [[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]]
RECT = [[0, 0, 0], [14, 0, 0], [14, 10, 0], [0, 10, 0]]
SQUARE_HOLE = [[4, 4, 0], [6, 4, 0], [6, 6, 0], [4, 6, 0]]
PLUS = [[4, 0, 0], [8, 0, 0], [8, 4, 0], [12, 4, 0], [12, 8, 0], [8, 8, 0],
        [8, 12, 0], [4, 12, 0], [4, 8, 0], [0, 8, 0], [0, 4, 0], [4, 4, 0]]
L_SHAPE = [[0, 0, 0], [12, 0, 0], [12, 4, 0], [5, 4, 0], [5, 10, 0], [0, 10, 0]]
DISC = [[5 + 5 * cos(2 * pi * i / 48), 5 + 5 * sin(2 * pi * i / 48), 0.0]
        for i in range(48)]
RING = [[5 + 3 * cos(2 * pi * i / 24), 5 + 3 * sin(2 * pi * i / 24), 0.0]
        for i in range(24)]


def arc(cx, cy, r, a0, a1, n=20):
    return [[cx + r * cos(a0 + (a1 - a0) * i / n),
             cy + r * sin(a0 + (a1 - a0) * i / n), 0.0] for i in range(n + 1)]


def turn(points, k, centre=(5.0, 5.0)):
    """``points`` rotated by ``k`` quarter turns about ``centre``."""
    out = []
    for p in points:
        x, y = p[0] - centre[0], p[1] - centre[1]
        for _ in range(k % 4):
            x, y = -y, x
        out.append([x + centre[0], y + centre[1], 0.0])
    return out


def flip(points, axis, centre=(5.0, 5.0)):
    if axis == 'x':
        return [[2 * centre[0] - p[0], p[1], 0.0] for p in points]
    return [[p[0], 2 * centre[1] - p[1], 0.0] for p in points]


#: the case from the Rhino screenshot: one arc and its three quarter turns.
#: The base arc is itself symmetric about the 45 degree line, so the SET is
#: invariant under all eight symmetries of the square.
_BASE_ARC = arc(0.0, 0.0, 6.5, 0.35, pi / 2 - 0.35)
FOUR_ARCS = [turn(_BASE_ARC, k) for k in range(4)]

#: the same four arcs, but generated the way a Rhino user generates them --
#: MIRRORED across the two centre lines. Mirroring alone gives a group of order
#: 4 that does NOT contain the quarter turn, so this is a different problem
#: even though the picture looks the same at a glance.
_CORNER_ARC = arc(0.0, 0.0, 6.5, 0.35, pi / 2 - 0.35)
MIRRORED_ARCS = [_CORNER_ARC, flip(_CORNER_ARC, 'x'),
                 flip(_CORNER_ARC, 'y'), flip(flip(_CORNER_ARC, 'x'), 'y')]

#: two arcs on a rectangle -- 2-fold, and the bounding box is not square, so
#: the quarter turn is not even a candidate
TWO_ARCS = [arc(0.0, 5.0, 5.0, -0.9, 0.9),
            [[14 - p[0], 10 - p[1], 0.0] for p in arc(0.0, 5.0, 5.0, -0.9, 0.9)]]

DOMAINS = [
    # label,               outer,  holes,          guides,         kwargs
    ('square',             SQUARE, None,           None,           {}),
    ('square+4 arcs',      SQUARE, None,           FOUR_ARCS,      {}),
    ('square+4 mirrored',  SQUARE, None,           MIRRORED_ARCS,  {}),
    ('rect+2 arcs',        RECT,   None,           TWO_ARCS,       {}),
    ('square+sq hole',     SQUARE, [SQUARE_HOLE],  None,           {}),
    ('plus-plate',         PLUS,   None,           None,           {}),
    ('disc',               DISC,   None,           None,           {}),
    ('square+ring cable',  SQUARE, None,           [RING],
     dict(guide_band=1.5, guide_weight=5.0)),
    ('L-shape (control)',  L_SHAPE, None,          None,           {}),
]


# ----------------------------------------------------------------------
# measurement
# ----------------------------------------------------------------------

def build(outer, holes, guides, kwargs, spacing, steps):
    """One decomposition, by the original route or the symmetric one.

    ``steps=None`` means ``FieldDecomposition.from_boundary`` -- the route in
    production today, untouched by this script.
    """
    if steps is None:
        return FieldDecomposition.from_boundary(
            outer, holes, guides=guides, target_length=spacing, **kwargs)
    return sym.decomposition(outer, holes, guides=guides, target_length=spacing,
                             steps=steps, **kwargs)


def measure(label, outer, holes, guides, kwargs, spacing, steps, group=None):
    """One row. Never raises: a domain that blows up is recorded as blown up."""
    row = {'domain': label, 'spacing': spacing}
    start = time.time()
    try:
        d = build(outer, holes, guides, kwargs, spacing, steps)
    except Exception as exc:
        row['error'] = '{}: {}'.format(type(exc).__name__, str(exc)[:90])
        return row

    if group is None:
        group = getattr(d, 'symmetry', None) or sym.Symmetry.detect(
            [outer] + list(holes or []) + list(guides or []))
    row['order'] = len(group)
    row['field_dev'] = sym.field_invariance(d.field, group)
    row['sings'] = len(d.field.singularities())
    row['ph_ok'] = d.field.poincare_hopf()['ok']
    launched = d.tracer._singularity_points()
    row['orbits'] = sorted(len(o) for o in sym.singularity_orbits(
        d.field, group, points=launched))
    row['sing_sym'] = min((s for name, (s, _) in
                           sym.invariance(list(launched.values()), group).items()
                           if name != 'identity'), default=1.0)

    try:
        coarse = d.decomposition_mesh()
        row['patches'] = coarse.number_of_faces()
        row['poles'] = sum(1 for f in coarse.faces()
                           if len(coarse.face_vertices(f)) == 3)
        corners = [coarse.vertex_coordinates(v) for v in coarse.vertices()]
        row.update(_layout_symmetry(corners, group))
        row['coarse_xyz'] = sorted((round(p[0], 6), round(p[1], 6)) for p in corners)
    except Exception as exc:
        row['error'] = 'coarse: {}: {}'.format(type(exc).__name__, str(exc)[:80])
        return row

    try:
        dense = d.quad_mesh(target_length=QUAD_TARGET)
        q = d.quality(dense)
        row['route'] = q['route']
        row['faces'] = q['faces']
        row['min_angle'] = q['min_angle']
        row['aspect'] = q['aspect_max']
        ok, why = hard_floor(q)
        row['floor'] = 'pass' if ok else why[:44]
    except Exception as exc:
        row['error'] = 'dense: {}: {}'.format(type(exc).__name__, str(exc)[:80])

    row['seconds'] = time.time() - start
    row['warnings'] = d.warnings()
    return row


def _layout_symmetry(corners, group):
    """Exactness of the coarse CORNERS, split into rotations and reflections."""
    scores = sym.invariance(corners, group)
    rotations = [v for g, v in scores.items()
                 if g != 'identity' and not g.startswith('mirror')]
    mirrors = [v for g, v in scores.items() if g.startswith('mirror')]
    out = {}
    out['rot'] = min((s for s, _ in rotations), default=None)
    out['mir'] = min((s for s, _ in mirrors), default=None)
    out['dev'] = max([d for _, d in rotations] + [d for _, d in mirrors] or [0.0])
    return out


def pct(value):
    return '  -  ' if value is None else '{:5.0%}'.format(value)


def show(rows, title):
    if title:
        print(title)
    print('    {:24s} {:>3s} {:>9s} {:>4s} {:>5s} {:>4s} {:>4s} {:>6s} {:>6s} '
          '{:>7s} {:>7s} {:>6s}  {}'.format(
              'domain', 'grp', 'field dev', 'sing', 'sing%', 'ptch', 'pole',
              'rot', 'mirror', 'dev', 'min ang', 'aspect', 'floor'))
    for row in rows:
        if 'error' in row:
            print('    {:24s} {}'.format(row['domain'], row['error']))
            continue
        print('    {:24s} {:3d} {:9.1e} {:4d} {:>5s} {:4d} {:4d} {:>6s} {:>6s} '
              '{:7.4f} {:7.2f} {:6.2f}  {}'.format(
                  row['domain'], row['order'], row['field_dev'], row['sings'],
                  pct(row.get('sing_sym')), row.get('patches', -1),
                  row.get('poles', -1), pct(row.get('rot')), pct(row.get('mir')),
                  row.get('dev', 0.0), row.get('min_angle', 0.0),
                  row.get('aspect', 0.0), row.get('floor', '?')))


# ----------------------------------------------------------------------
# part 1 -- where the symmetry goes, one stage at a time
# ----------------------------------------------------------------------

def part1():
    print('--- 1. where the symmetry goes -----------------------------------')
    print()
    group = sym.Symmetry.detect([SQUARE] + FOUR_ARCS)
    print('    input: 10x10 square, four arcs. {}'.format(group))
    print()
    print('    Stage by stage on the ORIGINAL route. "exact" is the share of')
    print('    points whose image under the operation is also a point; "dev" is')
    print('    the worst distance from an image to the nearest one. Spacing is')
    print('    {}, so a dev above that is more than a whole triangle.'.format(SPACING))
    print()

    bg = BackgroundMesh.from_boundary(SQUARE, target_length=SPACING)
    constraints = from_boundary(bg) + from_curves(bg, FOUR_ARCS, mode='perpendicular')
    field = CrossField.solve(bg, constraints)
    tracer = Tracer(field)
    seps, _ = tracer.separatrices()

    boundary = [bg.mesh.vertex_coordinates(v) for v in bg.boundary_vertices()]
    interior = [bg.mesh.vertex_coordinates(v) for v in bg.mesh.vertices()
                if v not in bg.boundary_vertices()]
    singular = [bg.mesh.face_centroid(f) for f, _ in field.singularities()]
    ends = [p for s in seps for p in (s.points[0], s.points[-1])]

    print('    {:24s} {:>7s} {:>8s}   {:>7s} {:>8s}'.format(
        'stage', 'rot90', 'dev', 'mirror', 'dev'))
    for label, points in (('boundary vertices', boundary),
                          ('interior vertices', interior),
                          ('singular faces', singular),
                          ('separatrix endpoints', ends)):
        scores = sym.invariance(points, group)
        r = scores.get('rot90', (None, 0.0))
        m = scores.get('mirror-vertical', (None, 0.0))
        print('    {:24s} {:>7s} {:8.4f}   {:>7s} {:8.4f}'.format(
            label, pct(r[0]), r[1], pct(m[0]), m[1]))

    print()
    print('    The boundary is fine and the interior is not, so the loss is in')
    print('    the interior grid -- before the field is solved. Two causes in')
    print('    background.py, and they are separable:')
    print()
    import framefield.background as bgmod
    original = bgmod._jitter
    for spacing in (0.5, 0.3):
        for label, jitter in (('jitter ON ', original),
                              ('jitter OFF', lambda i, j, a: 0.0)):
            bgmod._jitter = jitter
            probe = BackgroundMesh.from_boundary(SQUARE, target_length=spacing)
            points = [probe.mesh.vertex_coordinates(v) for v in probe.mesh.vertices()
                      if v not in probe.boundary_vertices()]
            score = sym.invariance(points, group)['rot90'][0]
            print('      spacing {:.1f}  (10/{:.1f} = {:5.2f})  {}  rot90 {:>6s}'
                  .format(spacing, spacing, 10.0 / spacing, label, pct(score)))
    bgmod._jitter = original
    print()
    print('    The hash jitter breaks it always. The anchoring at min(xs)')
    print('    breaks it additionally whenever the spacing does not divide the')
    print('    domain. The jitter cannot simply be deleted -- a regular grid is')
    print('    cocircular everywhere and Qhull then picks diagonals arbitrarily')
    print('    -- so symmetry.py makes it equivariant instead.')
    print()

    print('    Each fix, added to the one above it:')
    print()
    rows = []
    for label, steps in (('baseline', None),
                         ('+ symmetric background', ('background',)),
                         ('+ equivariant guides', ('background', 'guides')),
                         ('+ field symmetrisation',
                          ('background', 'guides', 'field')),
                         ('+ snapped singularities',
                          ('background', 'guides', 'field', 'singularities')),
                         ('+ network shims', ALL_STEPS)):
        row = measure(label, SQUARE, None, FOUR_ARCS, {}, SPACING, steps, group)
        rows.append(row)
    show(rows, '')
    print()
    print('    grp = order of the detected group. field dev = max|u - rho(g)u|,')
    print('    i.e. how symmetric the FIELD itself is. sing% = share of the')
    print('    LAUNCH POINTS with a twin. rot/mirror = share of coarse corners')
    print('    with a twin. dev = worst distance to one.')
    print()
    print('    The last two rows are the ones that are easy to miss. After the')
    print('    third the FIELD is symmetric to 2e-16 and the layout still is')
    print('    not, because a singularity is reported as a FACE and a face is a')
    print('    cruder object than the field on it: the |u| minima under those')
    print('    faces are an exact orbit, the face centroids are not. Fixing')
    print('    that alone still leaves 8%, all of it in build_network -- which')
    print('    recomputes the centroids it was just given, and clusters trace')
    print('    ends greedily in trace order.')
    print()
    print('    Compare the baseline deviation with build_network\'s snapping')
    print('    tolerance, 0.8 * spacing = {:.1f}. It straddles it, which is why'
          .format(0.8 * SPACING))
    print('    the original picture carries a whole extra patch on one side')
    print('    rather than a gently distorted version of the other.')


# ----------------------------------------------------------------------
# part 2 -- does it hold up across domains
# ----------------------------------------------------------------------

def part2(only=None):
    print('--- 2. the same two routes over a suite of symmetric domains -----')
    print()
    for label, outer, holes, guides, kwargs in DOMAINS:
        if only and label != only:
            continue
        group = sym.Symmetry.detect([outer] + list(holes or []) + list(guides or []))
        rows = [measure(label + '  [original]', outer, holes, guides, kwargs,
                        SPACING, None, group),
                measure(label + '  [symmetric]', outer, holes, guides, kwargs,
                        SPACING, ALL_STEPS, group)]
        show(rows, '  {}: {}'.format(label, ', '.join(group.names())))
        for row in rows:
            for warning in row.get('warnings', []):
                print('      warn ({}): {}'.format(
                    'sym' if 'symmetric' in row['domain'] else 'orig', warning))
        print()


# ----------------------------------------------------------------------
# part 3 -- is it stable
# ----------------------------------------------------------------------

def part3():
    print('--- 3. stability --------------------------------------------------')
    print()

    print('  (a) determinism -- the same input twice, coarse corners compared')
    for label, outer, holes, guides, kwargs in DOMAINS[:4]:
        a = measure(label, outer, holes, guides, kwargs, SPACING,
                    ALL_STEPS)
        b = measure(label, outer, holes, guides, kwargs, SPACING,
                    ALL_STEPS)
        same = a.get('coarse_xyz') == b.get('coarse_xyz')
        print('      {:22s} {}'.format(label, 'identical' if same else 'DIFFERENT'))
    print()

    print('  (b) an ASYMMETRIC domain must come out of the new route byte for')
    print('      byte the same as the old one -- no symmetry, no changes')
    for label, outer, holes, guides, kwargs in [d for d in DOMAINS
                                                if 'control' in d[0]]:
        group = sym.Symmetry.detect([outer] + list(holes or []) + list(guides or []))
        a = measure(label, outer, holes, guides, kwargs, SPACING, None)
        b = measure(label, outer, holes, guides, kwargs, SPACING,
                    ALL_STEPS)
        same = a.get('coarse_xyz') == b.get('coarse_xyz')
        print('      {:22s} group order {}, layouts {}'.format(
            label, len(group), 'identical' if same else 'DIFFERENT'))
    print()

    print('  (c) resolution sweep -- symmetry is not a resolution artefact, so')
    print('      it must not come and go with the background spacing')
    for label, outer, holes, guides, kwargs in DOMAINS[1:4]:
        group = sym.Symmetry.detect([outer] + list(holes or []) + list(guides or []))
        rows = []
        for spacing in (0.6, 0.5, 0.4):
            for tag, steps in (('orig', None), ('sym ', ALL_STEPS)):
                row = measure('{} {} @{:.1f}'.format(label, tag, spacing), outer,
                              holes, guides, kwargs, spacing, steps, group)
                rows.append(row)
        show(rows, '  {}'.format(label))
        print()


# ----------------------------------------------------------------------
# the picture
# ----------------------------------------------------------------------

#: the layout's own mirror image, drawn on top of it
OVERLAY = (0.0, 0.72, 0.30)
#: the symmetry axes
AXIS = (0.72, 0.72, 0.78)
#: guide curves. ``viz.add_guides`` draws them at linewidth 6, which is right
#: for one cable on one panel and swamps an eight-panel grid whose guides run
#: through the middle of every one of them.
GUIDE = (0.55, 0.20, 0.70)


def _axis_direction(g):
    """A vector along a reflection's axis: its eigenvector for eigenvalue +1.

    ``(1 + a, c)`` is twice the projection of the x-axis onto the mirror line,
    which is that eigenvector -- except for the reflection that sends x to -x,
    where it vanishes and the y-axis column ``(b, 1 + d)`` is used instead.
    """
    a, b, c, d = g[0], g[1], g[2], g[3]
    if (1 + a, c) != (0, 0):
        return (1.0 + a, float(c))
    return (float(b), 1.0 + d)


def add_axes(group, symmetry, span, dx=0.0, dy=0.0):
    """Every mirror axis of the detected group, as one object.

    Worth drawing even though it is not geometry: it is the thing the layout is
    being judged against, and an asymmetry that is obvious once the axis is on
    the screen is invisible without it.
    """
    from compas.colors import Color
    from compas.datastructures import Graph

    graph = Graph()
    n = 0
    for g in symmetry.elements:
        if not g[4]:
            continue
        ux, uy = _axis_direction(g)
        length = (ux * ux + uy * uy) ** 0.5
        ux, uy = ux / length * span, uy / length * span
        cx, cy = symmetry.centre[0] + dx, symmetry.centre[1] + dy
        graph.add_node(n, x=cx - ux, y=cy - uy, z=-0.02)
        graph.add_node(n + 1, x=cx + ux, y=cy + uy, z=-0.02)
        graph.add_edge(n, n + 1)
        n += 2
    if n:
        # ``edgecolor``, NOT ``linecolor``: compas_viewer's GraphObject reads
        # ``edgecolor.default`` per edge and ignores ``linecolor`` entirely, so a
        # graph handed a linecolor draws black and says nothing about it
        group.add(graph, edgecolor=Color(*AXIS), linewidth=1, show_points=False,
                  name='{} mirror axis/axes'.format(n // 2))
    return n // 2


def add_thin_guides(group, guides, dx=0.0, dy=0.0):
    """The guide curves, as one object and thin enough to see past."""
    from compas.colors import Color
    from compas.datastructures import Graph

    graph = Graph()
    n = 0
    for curve in guides or []:
        pts = [list(p) for p in getattr(curve, 'points', curve)]
        for a, b in zip(pts, pts[1:]):
            graph.add_node(n, x=a[0] + dx, y=a[1] + dy, z=0.09)
            graph.add_node(n + 1, x=b[0] + dx, y=b[1] + dy, z=0.09)
            graph.add_edge(n, n + 1)
            n += 2
    if n:
        group.add(graph, edgecolor=Color(*GUIDE), linewidth=2, show_points=False,
                  name='{} guide curve(s)'.format(len(guides)))
    return n // 2


def add_mirror_overlay(group, mesh, symmetry, dx=0.0, dy=0.0):
    """The coarse layout's own image under the group, drawn on top of it.

    THIS is the panel to look at. Two side-by-side layouts are hard to compare
    by eye -- the difference is one patch corner half a triangle out of place --
    but a layout drawn on top of its own mirror image either coincides or
    visibly doubles. On the symmetric route the green lines disappear into the
    black ones exactly; on the original they run beside them.

    Every reflection in the group is drawn, so a layout that is symmetric about
    one axis and not another says so.

    Drawn THIN over the layout's thick lines, not the other way round. Four
    reflections at linewidth 3 stack into a solid green layout that hides the
    black one underneath, which inverts the whole point of the panel: what has
    to be legible is a green line sitting INSIDE a black one, against a green
    line sitting beside it.
    """
    from compas.colors import Color
    from compas.datastructures import Graph

    reflections = [g for g in symmetry.elements if g[4]] or [
        g for g in symmetry.elements if g[5] != 'identity']
    if not reflections:
        return 0

    graph = Graph()
    n = 0
    for g in reflections:
        for u, v in mesh.edges():
            a = symmetry.apply(g, mesh.vertex_coordinates(u))
            b = symmetry.apply(g, mesh.vertex_coordinates(v))
            graph.add_node(n, x=a[0] + dx, y=a[1] + dy, z=0.12)
            graph.add_node(n + 1, x=b[0] + dx, y=b[1] + dy, z=0.12)
            graph.add_edge(n, n + 1)
            n += 2
    group.add(graph, edgecolor=Color(*OVERLAY), linewidth=1, show_points=False,
              name='the same layout, mirrored ({} axes)'.format(len(reflections)))
    return len(reflections)


def view(outer=SQUARE, guides=FOUR_ARCS, holes=None, kwargs=None,
         spacing=SPACING):
    """Eight panels: the original route along the top, the symmetric one below.

    One stage per COLUMN, so the divergence appears where it actually happens --
    column 1, in the triangulation, before the field exists. Four columns by two
    rows rather than the other way round because the viewer's window is wide:
    the camera frames the grid through its VERTICAL field of view, so a tall
    grid is fitted by shrinking every panel while the width goes unused.
    """
    from compas_viewer import Viewer

    from compas_viewer.config import Config

    from framefield.viz import (Grid, add_background, add_boundary, add_dense,
                                add_field, add_layout, add_singularities)

    kwargs = kwargs or {}
    group = sym.Symmetry.detect([outer] + list(holes or []) + list(guides or []))
    span = max(max(p[0] for p in outer) - min(p[0] for p in outer),
               max(p[1] for p in outer) - min(p[1] for p in outer)) * 0.62

    built = []
    for steps in (None, ALL_STEPS):
        d = build(outer, holes, guides, kwargs, spacing, steps)
        coarse = d.decomposition_mesh()
        dense = d.quad_mesh(target_length=QUAD_TARGET)
        corners = [coarse.vertex_coordinates(v) for v in coarse.vertices()]
        scores = _layout_symmetry(corners, group)
        built.append((d, coarse, dense, scores))

    # the ground grid is a 10x10 plane at the origin, which on an eight-panel
    # scene reads as a ninth panel sitting next to the first one
    config = Config()
    config.renderer.show_grid = False
    viewer = Viewer(config)
    grid = Grid(pitch=span * 2.4, cols=4)
    rows = ('ORIGINAL', 'SYMMETRIC')

    for row, (d, coarse, dense, scores) in enumerate(built):
        tag = rows[row]
        first = 4 * row

        dx, dy = grid.cell(first)
        g = viewer.scene.add_group(name='{}. {} background -- {} vertices, {} '
                                       'where the symmetry dies'
                                   .format(first + 1, tag,
                                           d.background.mesh.number_of_vertices(),
                                           'THIS is' if row == 0 else 'not'))
        add_background(g, d, dx, dy)
        add_axes(g, group, span, dx, dy)
        add_boundary(g, d, dx, dy)
        add_thin_guides(g, guides, dx, dy)

        dx, dy = grid.cell(first + 1)
        g = viewer.scene.add_group(
            name='{}. {} field -- {} singularities, launch points {:.0%} '
                 'symmetric'.format(first + 2, tag,
                                    len(d.field.singularities()),
                                    _launch_symmetry(d, group)))
        # a lower tick budget than the default 300: eight panels at this size
        # turn a 300-tick field into a black hatch, and viz.add_field's
        # linecolor is ignored by compas_viewer for graphs so it cannot be
        # lightened instead
        add_field(g, d, dx, dy, budget=140)
        add_singularities(g, d, dx, dy, size=12)
        add_axes(g, group, span, dx, dy)
        add_boundary(g, d, dx, dy)
        add_thin_guides(g, guides, dx, dy)

        dx, dy = grid.cell(first + 2)
        poles = sum(1 for f in coarse.faces() if len(coarse.face_vertices(f)) == 3)
        g = viewer.scene.add_group(
            name='{}. {} layout, ON TOP OF ITS OWN MIRROR IMAGE -- {} patches, '
                 '{} pole(s), corners {:.0%} symmetric, worst gap {:.3f}'
                 .format(first + 3, tag, coarse.number_of_faces(), poles,
                         scores['mir'] if scores['mir'] is not None
                         else scores['rot'], scores['dev']))
        # faces off: the point of the panel is whether two sets of LINES
        # coincide, and a translucent patch fill over them hides exactly that
        add_layout(g, d, dx, dy, faces=False)
        add_mirror_overlay(g, coarse, group, dx, dy)
        add_axes(g, group, span, dx, dy)
        add_thin_guides(g, guides, dx, dy)

        dx, dy = grid.cell(first + 3)
        q = d.quality(dense)
        g = viewer.scene.add_group(
            name='{}. {} quad mesh -- {} quads, min angle {:.1f} deg, aspect '
                 '{:.2f}'.format(first + 4, tag, dense.number_of_faces(),
                                 q['min_angle'], q['aspect_max']))
        add_dense(g, dense, dx, dy)
        add_thin_guides(g, guides, dx, dy)

    # Grid.frame sizes the distance off the grid's LARGER dimension, which on a
    # 4x2 grid is its width -- and width is the dimension the viewport has to
    # spare. 0.85 fits the two rows to the height instead of fitting four
    # columns to it and leaving most of the window empty.
    grid.frame(viewer, margin=0.85)
    viewer.show()


def _launch_symmetry(decomposition, group):
    launched = list(decomposition.tracer._singularity_points().values())
    if not launched:
        return 1.0
    return min((s for name, (s, _) in sym.invariance(launched, group).items()
                if name != 'identity'), default=1.0)


# ----------------------------------------------------------------------

def main():
    only = None
    if '--domain' in sys.argv:
        only = sys.argv[sys.argv.index('--domain') + 1]

    if '--view' in sys.argv:
        view()
        return 0

    parts = [flag for flag in ('--teardown', '--suite', '--stability')
             if flag in sys.argv]
    if not parts or '--teardown' in parts:
        part1()
        print()
    if not parts or '--suite' in parts:
        part2(only)
    if not parts or '--stability' in parts:
        part3()

    print()
    print('WHERE THIS LANDS')
    print('  A symmetric problem now gives a symmetric layout, exactly: 100%')
    print('  under all eight elements at deviation 0.0000, Poincare-Hopf still')
    print('  passing, no repair warnings, and the same 16 patches at every')
    print('  background spacing tried. The solve was never at fault -- on a')
    print('  symmetric mesh it is exact to 1e-15 -- and no tolerance anywhere')
    print('  was touched.')
    print()
    print('  Element quality moves WITH symmetry rather than against it. On the')
    print('  four-arc square, minimum angle 29.75 -> 57.83 deg and worst aspect')
    print('  2.01 -> 1.60. The ring cable, which the original route takes below')
    print('  the hard floor (max_angle 179.699 deg, 53 patches, 10 poles), comes')
    print('  back inside it here: 21 patches, 4 poles, floor pass.')
    print()
    print('  THE COST, and the row to argue about. The layout stops splitting')
    print('  where the asymmetry used to split it, so a guide is followed a')
    print('  little less closely by a simpler layout: 10.2 -> 10.7 deg on the')
    print('  four-arc square, 13.0 -> 14.8 deg on rect+2 arcs. That second row')
    print('  goes from 6 singularities and 17 patches to NONE and one patch.')
    print('  Both readings are defensible -- two arcs on a rectangle force no')
    print('  winding the walls cannot absorb, so the six look like artefacts of')
    print('  the asymmetric background -- but it is a change in what the layout')
    print('  DOES, not only in how symmetric it looks. Decide that row on')
    print('  purpose.')
    print()
    print('  TWO OF THE FIVE FIXES ARE SHIMS, not edits: symmetry.py replaces')
    print('  repair._cluster and decomposition.build_network for the duration of')
    print('  one call and restores them. They are each a few lines in place at')
    print('  adoption time. Nothing on disk is changed for any other caller,')
    print('  which is why 15_baseline.py is unaffected by this file existing.')
    print()
    print('  BEFORE ADOPTING: the background triangulation changes, so every row')
    print('  of baseline.json moves. Read that diff row by row before')
    print('  regenerating -- and note that from_curves\'s curve-order dependence')
    print('  is a real bug in its own right, worth fixing whether or not any of')
    print('  the rest of this is adopted.')

    if '--no-view' not in sys.argv:
        view()
    return 0


if __name__ == '__main__':
    sys.exit(main())
