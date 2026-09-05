"""EXAMPLE 4 -- the constraints, each with a reproducer.

Read this before pointing the field front end at your own outline. Every entry
is a MEASURED failure of the current implementation, not a caveat in principle.

The shape of the limits changed once the plate case was fixed. It is no longer
"polygons work, plates do not" -- rectilinear plates are now the ROBUST case.
What remains fragile is CURVED boundaries.

Opens a compas_viewer scene putting each failure next to the equivalent success.

Run:
    python 13_limits.py
    python 13_limits.py --no-view          # text only
"""
import os
import sys
from math import atan2, cos, degrees, pi, sin

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from compas_singular.framefield.decomposition import FieldDecomposition  # noqa: E402


SQUARE = [[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]]
T_SHAPE = [[0, 0, 0], [12, 0, 0], [12, 4, 0], [8, 4, 0],
           [8, 10, 0], [4, 10, 0], [4, 4, 0], [0, 4, 0]]


def circle(cx, cy, r, n=48):
    return [[cx + r * cos(2 * pi * i / n), cy + r * sin(2 * pi * i / n), 0.0]
            for i in range(n)]


def ellipse(cx, cy, a, b, n=48):
    """A circle's two radii pulled apart -- LIMIT 4's knob.

    Same winding as ``circle``, so between the round hole and the ovals below it
    the only things that change are the aspect ratio and, at fixed ``n``, how
    hard the tip is to sample. ``n`` is the second knob and the decisive one.
    """
    return [[cx + a * cos(2 * pi * i / n), cy + b * sin(2 * pi * i / n), 0.0]
            for i in range(n)]


def max_turn(loop):
    """Sharpest turn between consecutive edges of a closed loop, in degrees.

    LIMIT 4's whole measurement. Above 45 degrees -- the half-period a cross
    field folds every angle difference into -- a sampled arc and a corner become
    the same object, and the arc loses.
    """
    n = len(loop)
    worst = 0.0
    for i in range(n):
        p, q, r = loop[i - 1], loop[i], loop[(i + 1) % n]
        a = atan2(q[1] - p[1], q[0] - p[0])
        b = atan2(r[1] - q[1], r[0] - q[0])
        worst = max(worst, abs(degrees((b - a + pi) % (2 * pi) - pi)))
    return worst


DISC = circle(5, 5, 5)
ELLIPSE = [[5 + 7 * cos(2 * pi * i / 60), 5 + 3 * sin(2 * pi * i / 60), 0.0]
           for i in range(60)]
STADIUM = ([[0, 0, 0], [12, 0, 0]]
           + [[12 + 3 * cos(a), 3 + 3 * sin(a), 0.0]
              for a in [-pi / 2 + pi * i / 10 for i in range(1, 10)]]
           + [[12, 6, 0], [0, 6, 0]])
ROUND_HOLE = circle(5, 5, 2)[::-1]


def attempt(boundary, holes=None, target_length=0.5, guides=None, **kw):
    """Run the pipeline; return (decomposition, one-line outcome, ok)."""
    try:
        d = FieldDecomposition.from_boundary(boundary, holes, guides=guides,
                                             target_length=target_length, **kw)
    except Exception as exc:
        return None, 'field/trace raised {}: {}'.format(type(exc).__name__, str(exc)[:55]), False

    warn = '  [{}]'.format('; '.join(d.warnings())) if d.warnings() else ''
    try:
        coarse = d.decomposition_mesh()
    except Exception as exc:
        return d, '{}: {}'.format(type(exc).__name__, str(exc)[:70]), False

    sides = {}
    for f in coarse.faces():
        n = len(coarse.face_vertices(f))
        sides[n] = sides.get(n, 0) + 1
    if set(sides) != {4}:
        return d, 'NOT all quads: sides {}{}'.format(sides, warn), False

    try:
        coarse.collect_strips()
        coarse.set_strips_density_target(0.8)
        coarse.densification(edges_to_curves=d.edges_to_curves())
        dense = coarse.get_quad_mesh()
    except Exception as exc:
        return d, 'densification raised {}'.format(type(exc).__name__), False

    return d, 'OK -- {} patches, {} quads, manifold {}{}'.format(
        coarse.number_of_faces(), dense.number_of_faces(), dense.is_manifold(), warn), True


def section(n, title, why):
    print()
    print('=' * 74)
    print('LIMIT {} -- {}'.format(n, title))
    print('=' * 74)
    print(why)
    print()


def main():
    failures = []

    section(1, 'Curved boundaries -- FIXED, and what it cost',
            "This used to read 'the one real geometric gap', and the diagnosis was\n"
            "wrong. The separatrix network DOES close on a disc; it was rejected\n"
            "because from_polylines never computes an intersection, so two\n"
            "separatrices crossing in their interiors got no node at the crossing\n"
            "and the recovered faces repeated vertices. See framefield/arrangement.py.\n"
            "The tell was that the hexagon failed at tl=0.5 and worked at 0.35 --\n"
            "a geometric limit does not come and go with background resolution.\n"
            "\n"
            "Every case below now reaches route 'field', and the layouts are coarse\n"
            "again. They were not: one stray non-quad patch used to send the WHOLE\n"
            "layout through the global quad split, so the hexagon's 10 quads and one\n"
            "triangle at tl=0.5 came out as 43 patches, the ellipse's 30, and the\n"
            "stadium's 73. The non-monotonicity was the tell -- 43 patches at 0.5,\n"
            "11 at 0.4, 8 at 0.35 -- and, as with the arrangement, it was three\n"
            "upstream defects rather than a geometric limit: a landing-point\n"
            "tolerance that merged two distinct wall nodes 0.35 apart, corner\n"
            "detection that read every vertex of a discretised arc as a corner, and\n"
            "a duplicate test that merged two different arms of one singularity.\n"
            "See repair.boundary_corners, repair._cluster and build_network.\n"
            "Now: hexagon 11 / 11 / 8, ellipse 10, stadium 4.")
    for label, bnd, holes in [('disc', DISC, None),
                              ('ellipse', ELLIPSE, None),
                              ('stadium', STADIUM, None),
                              ('disc with round hole', DISC, [ROUND_HOLE])]:
        for tl in (0.6, 0.5, 0.4):
            d, msg, ok = attempt(bnd, holes, target_length=tl)
            print('  {:22s} tl={:.1f}  {}  route {}'.format(
                label, tl, msg, d.route() if d else '-'))
            if not ok and label == 'disc' and tl == 0.4:
                failures.append((label + ' @0.4', d))
        print()
    print('  All four are on route "field" now, the disc with a round hole')
    print('  included -- it was the last one on a fallback, and for a different')
    print('  reason: an annulus has no corner and no singularity, so nothing')
    print('  volunteered a cut. LIMIT 4 is where that got its own source, and')
    print('  where the part of it that is still a limit is measured. Note the')
    print('  4 patches here hold across all three resolutions, which the disc')
    print('  itself (5/5/9) does not.')

    section(2, 'A cable steers the field but not the mesh',
            "The headline limitation against the stated objective, and the reason\n"
            "12_cables.py exists. Densification is a discrete Coons patch between\n"
            "coarse edges and never consults the field, so a cable reaches the\n"
            "output only if it moves a singularity or bends a separatrix.")
    _, base, _ = attempt(SQUARE)
    _, guided, _ = attempt(SQUARE, guides=[[[1, 2, 0], [5, 5, 0], [9, 8, 0]]],
                           guide_weight=None, guide_band=3.0)
    print('  square, no cable        : {}'.format(base))
    print('  square, hard wide cable : {}'.format(guided))
    print('  -> byte-identical. See 12_cables.py for the full measurement.')

    section(3, 'Non-orthogonal cable families are not supported',
            "Cases A and B of the evaluation doc -- an orthogonal stress field, or a\n"
            "SINGLE cable family -- are exactly what a cross field represents. Two\n"
            "families meeting at anything but 90 degrees (case C) need a frame field\n"
            "and the warp step: milestone 3.")
    try:
        FieldDecomposition.from_boundary(SQUARE, orthogonal=False)
        print('  orthogonal=False: unexpectedly accepted')
    except NotImplementedError as exc:
        print('  orthogonal=False raises, as it should:')
        print('    {}'.format(str(exc)))

    section(4, 'A sampled arc whose turn exceeds 45 deg is read as a corner',
            "This used to read 'regions needing a cut from neither a singularity nor\n"
            "a corner', with the round-holed disc as the case that had no source for\n"
            "one and fell to the triangulation backstop. That is now false: there is\n"
            "a THIRD source. trace.hole_launches fires cuts off any inner loop with\n"
            "no corners of its own -- an annulus has index 0, so the field is smooth\n"
            "and singularity-free and will never volunteer one -- and the disc with a\n"
            "round hole gives 4 patches on route 'field' at every resolution below.\n"
            "\n"
            "Flatten the hole and it breaks, but NOT because an oval is hard. The\n"
            "limit is 45 degrees, and it is a limit on the INPUT POLYGON:\n"
            "\n"
            "A cross is invariant under 90 deg rotation, so field.wrap_to_period folds\n"
            "every angle difference into +/-45 deg -- past that a direction is better\n"
            "matched by the NEXT arm round, and the two readings are indistinguishable.\n"
            "That is not an approximation, it is what a cross IS. So a boundary turn\n"
            "above 45 deg cannot be read as an arc continuing; it can only be read as\n"
            "a corner, where the turning is absorbed and the arm carries straight on.\n"
            "For a real corner that is correct -- a square's 90 deg turns fold to 0,\n"
            "which is why its field is constant and it has no singularity at all. For\n"
            "a SAMPLED ARC it is a lie the discretisation told, and nothing downstream\n"
            "can catch it: repair.SHARP_TURN is 45 deg too, so the corner detector\n"
            "agrees with the field, consistently and wrongly.\n"
            "\n"
            "The 4.2x0.5 ellipse has tip curvature radius b^2/a = 0.0595. At 48 points\n"
            "its two tips turn 57.7 deg, and three things flip at once:\n"
            "  - each aliased tip loses one full period of boundary winding, so the\n"
            "    measured winding is 2 instead of 0 -- exactly one per tip.\n"
            "  - the field MATCHES that winding (poincare_hopf ok=True, residual 0).\n"
            "    It is not violating anything; it is solving the wrong problem.\n"
            "  - boundary_corners finds those tips, so hole_launches skips the loop\n"
            "    as already-cornered -- and they are CONVEX, so corner_launches will\n"
            "    not fire on them either. The hole gets no cut from anywhere.\n"
            "\n"
            "Sample the SAME ellipse at 72 points instead and the worst turn is 40.3\n"
            "deg. Nothing else changes. Everything else does.\n"
            "\n"
            "Note this is not the free-magnitude solver of 19_field_accuracy.py -- the\n"
            "rows below are relax=True, and the aliasing is identical under both. That\n"
            "solver is a SECOND and independent cause of excess singularities here:\n"
            "on the 3.5x0.8 oval, relax=False wanders 8/12/18/8 across tl 0.6/0.5/0.4/\n"
            "0.3 with |u| collapsing to 0.025, while relax=True holds 8 at every\n"
            "resolution. Eight is that shape's real answer; the rest was noise.")
    for label, hole in [
            ('round r=2.0       ', ROUND_HOLE),
            ('oval  3.5x0.8  n=48', ellipse(5, 5, 3.5, 0.8)[::-1]),
            ('oval  4.2x0.5  n=48', ellipse(5, 5, 4.2, 0.5)[::-1]),
            ('oval  4.2x0.5  n=72', ellipse(5, 5, 4.2, 0.5, n=72)[::-1])]:
        d, msg, _ = attempt(circle(5, 5, 5), [hole], target_length=0.5, relax=True)
        if d is None:
            print('  {} {}'.format(label, msg))
            continue
        ph = d.field.poincare_hopf()
        _, rep = d.tracer.separatrices()
        print('  {}  max turn {:5.1f} deg{}'.format(
            label, max_turn(hole), '  <-- ALIASES' if max_turn(hole) > 45.0 else ''))
        print('  {:19s}  {}'.format('', msg))
        print('  {:19s}  {} singularities, winding {}, {} of 4 cuts landed'.format(
            '', ph['singular_faces'], ph['boundary_winding'], rep.get('hole_cuts', 0)))

    section(5, 'Poles appear only as repair, never from the field itself',
            "This used to read 'poles are not produced', and the poles argument was\n"
            "accepted for signature compatibility and ignored. Both are now false.\n"
            "solve_non_quad_faces resolves a triangular patch it cannot otherwise\n"
            "repair as a pseudo-quad pole -- the same fix INTERIOR_TRIANGLE_POLE_FIX\n"
            "documents for the medial-axis front end -- and a position passed in\n"
            "poles decides which of that triangle's corners the collapsed side sits\n"
            "at. What is still true is that a POLE IS NOT A FIELD SINGULARITY: the\n"
            "field never asks for one, so none of the domains below produces one,\n"
            "and a pole here always means the layout needed patching. 12_cables.py's\n"
            "ring cable is the live case -- 8 poles, and 43 patches instead of 160.")
    d, _, _ = attempt(T_SHAPE)
    coarse = d.decomposition_mesh(poles=[[6, 2, 0]])
    print('  T-plate, all quads already: {} faces, {} pseudo-quads (poles unused)'.format(
        coarse.number_of_faces(),
        sum(1 for f in coarse.faces() if coarse.is_face_pseudo_quad(f))))
    d2, _, _ = attempt(SQUARE, guides=[circle(5, 5, 3, n=24)])
    c2 = d2.decomposition_mesh()
    print('  ring cable, needs repair : {} faces, {} pseudo-quads'.format(
        c2.number_of_faces(),
        sum(1 for f in c2.faces() if c2.is_face_pseudo_quad(f))))

    section(6, 'Planar domains only',
            "background.py builds a planar Delaunay triangulation and every tangent\n"
            "space is world XY, so there is no parallel transport and no per-edge\n"
            "connection. face_basis() returns the identity and exists so the surface\n"
            "case will not need an API change. Z coordinates are dropped silently.")
    print('  (no reproducer -- a 3D outline is flattened without warning)')

    print()
    print('=' * 74)
    print('WHAT IS ROBUST: rectilinear plates -- square, L, T, U, plus, comb,')
    print('rectangle with a rectangular hole or slot -- AND curved boundaries:')
    print('disc, ellipse, stadium, hexagon, and a disc with a ROUND hole. All on')
    print('route "field". 14_arrangement.py measures crossings, dangling ends and')
    print('route for each.')
    print('WHAT IS NOT   : cables reaching the mesh (LIMIT 2); and any boundary')
    print('sampled so coarsely that one turn exceeds 45 deg, which a cross field')
    print('cannot tell from a corner (LIMIT 4). Neither the patch count nor the')
    print('round-holed disc is among them any more -- see LIMIT 1.')
    print()
    print('LIMIT 4 is a condition on how you GENERATE an outline, not a test you')
    print('can run on one. "No turn above 45 degrees" flags 15 of the 19 domains')
    print('here, including every rectilinear plate: a square turns 90 degrees at')
    print('each corner and is the best-behaved shape in the suite. A cross field')
    print('WANTS that -- 90 folds to 0, the arm continues into the next wall, and')
    print('the field comes out constant.')
    print()
    print('The distinction is not in the polygon. A 48-gon approximating a sharp')
    print('ellipse and a genuine 48-sided polygon ARE THE SAME OBJECT, and only')
    print('the author knows which was meant. Three guards were measured and all')
    print('three fail: turn magnitude flags the square; neighbour support cannot')
    print('separate the stadium corner (55% of peak) from the ellipse tip (52%);')
    print('and "the inner loop received no separatrix" never fires, because the')
    print('aliased loop is not starved but FLOODED -- 26 endpoints land on it.')
    print()
    print('So: if your outline is SAMPLED from a curve, sample it so no turn')
    print('exceeds 45 degrees -- refine by curvature, not by chord length, which')
    print('is what background._densify_loop does and why it subdivides everywhere')
    print('except the tip, where chords are already 6x below target and the turn')
    print('is worst. If you hand in a polygon with a turn above 45 degrees, it is')
    print('a corner. That is correct, and it is the only thing it can be.')
    print('=' * 74)

    if '--no-view' not in sys.argv:
        show()
    return 0


def show():
    """Each limit next to the equivalent success, so the gap is legible.

    Panels 3 to 5 are the ones that changed: they used to be drawn WITHOUT patch
    faces because a many-sided patch renders as a meaningless polygon. They now
    have a real layout to draw, so they are drawn like the working cases.
    """
    from compas_viewer import Viewer

    from compas_singular.framefield.viz import (Grid, add_dense, add_layout,
                                                add_separatrices, add_singularities)

    viewer = Viewer()
    grid = Grid(pitch=17.0, cols=3)

    # 1. a rectilinear plate, working
    d, _, _ = attempt(T_SHAPE, target_length=0.5)
    coarse = d.decomposition_mesh()
    coarse.collect_strips()
    coarse.set_strips_density_target(0.8)
    coarse.densification(edges_to_curves=d.edges_to_curves())
    dx, dy = grid.cell(0)
    g = viewer.scene.add_group(name='1. T-plate WORKS -- separatrices land on corners')
    add_dense(g, coarse.get_quad_mesh(), dx, dy)
    add_layout(g, d, dx, dy, faces=False)

    # 2. the disc at the background it always survived
    d6, _, _ = attempt(DISC, target_length=0.6)
    dx, dy = grid.cell(1)
    g = viewer.scene.add_group(name='2. disc tl=0.6 -- 5 patches, always worked')
    add_layout(g, d6, dx, dy)
    add_singularities(g, d6, dx, dy)

    # 3. the same disc, finer background: 3 separatrix crossings, and the case
    # that used to be discarded for them
    d4, _, _ = attempt(DISC, target_length=0.4)
    dx, dy = grid.cell(2)
    g = viewer.scene.add_group(
        name='3. disc tl=0.4 -- 3 crossings resolved, 9 patches, route field')
    add_layout(g, d4, dx, dy)
    add_singularities(g, d4, dx, dy)

    # 4, 5. the curved boundaries that used to fail at every resolution
    for i, (label, bnd) in enumerate(
            [('4. ellipse -- was triangular patches, now route field', ELLIPSE),
             ('5. stadium -- round end was drawn twice, now route field',
              STADIUM)], start=3):
        dd, _, _ = attempt(bnd, target_length=0.5)
        dx, dy = grid.cell(i)
        g = viewer.scene.add_group(name=label)
        add_layout(g, dd, dx, dy)
        add_separatrices(g, dd, dx, dy)
        add_singularities(g, dd, dx, dy)

    # 6. the cable that changes nothing
    dc, _, _ = attempt(SQUARE, guides=[[[1, 2, 0], [5, 5, 0], [9, 8, 0]]],
                       guide_weight=None, guide_band=3.0)
    cc = dc.decomposition_mesh()
    cc.collect_strips()
    cc.set_strips_density_target(0.8)
    cc.densification(edges_to_curves=dc.edges_to_curves())
    dx, dy = grid.cell(5)
    g = viewer.scene.add_group(
        name='6. hard cable, and the mesh ignores it -- see 12_cables.py')
    add_dense(g, cc.get_quad_mesh(), dx, dy)
    from compas_singular.framefield.viz import add_guides
    add_guides(g, [[[1, 2, 0], [5, 5, 0], [9, 8, 0]]], dx, dy)

    grid.frame(viewer)
    viewer.show()


if __name__ == '__main__':
    sys.exit(main())
