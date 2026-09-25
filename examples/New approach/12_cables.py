"""EXAMPLE 3 -- cables and force lines as INPUT to the field.

This is the objective of the whole exercise, so be precise about what is being
asked for:

    WANTED -- a quad mesh GENERATED FROM the frame field, where cables and force
    vectors are inputs that shape the field, and the mesh comes out already
    aligned to them.

    NOT -- an existing quad mesh whose lines are afterwards pulled onto a curve.
    That is post-processing, it is what ``src/compas_singular/guide_lines.py``
    does, and it is a different problem. guide_lines snaps polyedges the
    decomposition already contains, so it can only ever reuse directions the
    layout already has; its own docstring says "refinement subdivides the lines
    that exist, it does not invent lines in new directions". Useful for final
    adjustment. Not a substitute for generating the right layout.

Where the implementation stands against that objective:

    the FIELD  takes the constraints exactly -- 0.0 deg off a hard-constrained
               cable, tangent or perpendicular
    the MESH   follows -- ``framefield/densify.py`` integrates the field inside
               each patch instead of blending the patch's four edges, so a
               cable steers the interior without the layout having to split

This example measures both halves rather than describing them. Part 2 used to
report a median edge angle of 0.0 for every row; that was not a measurement,
it was two byte-identical meshes being compared to each other.

What the mesh does NOT do is reach 0.0 degrees off the cable, and the reason is
geometric. A patch boundary is held fixed so the patches still weld and the
strip grammar still applies, and on a plain square that boundary IS the four
walls -- so the mesh must turn to meet them however well the field is
integrated. Part 2 prints the profile along the cable, which shows the mesh
sitting near 25 degrees off where it is clear of a wall and near 33 at the ends.

Opens a compas_viewer scene contrasting the two halves: the field with and
without the cable, then the meshes they produce. Panels 3 and 4 were identical
for the whole life of this project and are not any more.

Run:
    python 12_cables.py
    python 12_cables.py --no-view          # numbers only
"""
import os
import statistics
import sys
from math import atan2, cos, degrees, pi, sin

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from compas_singular.framefield import BackgroundMesh                    # noqa: E402
from compas_singular.framefield import CrossField                        # noqa: E402
from compas_singular.framefield import from_boundary                     # noqa: E402
from compas_singular.framefield import from_curves                       # noqa: E402
from compas_singular.framefield.decomposition import FieldDecomposition  # noqa: E402
from compas_singular.framefield.quality import curve_alignment           # noqa: E402
from compas_singular.framefield.quality import curve_alignment_profile   # noqa: E402


SQUARE = [[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]]

#: a single cable running corner to corner -- case B of the evaluation doc
CABLE = [[1.0, 2.0, 0.0], [5.0, 5.0, 0.0], [9.0, 8.0, 0.0]]


def arc(cx, cy, r, a0, a1, n=24):
    return [[cx + r * cos(a0 + (a1 - a0) * i / n),
             cy + r * sin(a0 + (a1 - a0) * i / n), 0.0] for i in range(n + 1)]


#: a hoop -- two cable families that between them wind around a point
RING = arc(5, 5, 3, 0, 2 * pi)


def solve(guides=None, mode='perpendicular', **kw):
    bg = BackgroundMesh.from_boundary(SQUARE, target_length=0.5)
    constraints = from_boundary(bg)
    if guides:
        constraints += from_curves(bg, guides, mode=mode, **kw)
    return bg, CrossField.solve(bg, constraints)


# --- 1. the field takes the constraint ---------------------------------------

def alignment(bg, field, curve):
    """Mean angle between the cross field and a curve, sampled along it."""
    offs = []
    for (a, b) in zip(curve, curve[1:]):
        mid = [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2, 0.0]
        best = min(bg.mesh.faces(),
                   key=lambda f: (bg.mesh.face_centroid(f)[0] - mid[0]) ** 2
                   + (bg.mesh.face_centroid(f)[1] - mid[1]) ** 2)
        theta = field.angle_in_face(best, (1 / 3., 1 / 3., 1 / 3.))
        tangent = atan2(b[1] - a[1], b[0] - a[0])
        d = (tangent - theta) % (pi / 2)
        offs.append(degrees(min(d, pi / 2 - d)))
    return statistics.mean(offs)


def part1():
    print('--- 1. the FIELD takes cables and force lines -------------------')
    print('    mean angle between the cross field and the cable, along it:')
    print()
    for label, mode, kw in [
        ('hard, tangent', 'tangent', dict(weight=None, band=3.0)),
        ('hard, perpendicular', 'perpendicular', dict(weight=None, band=3.0)),
        ('soft (weight 1)', 'tangent', dict(weight=1.0)),
        ('soft, wide band', 'tangent', dict(weight=5.0, band=3.0)),
        ('no cable at all', None, None),
    ]:
        if mode is None:
            bg, field = solve()
        else:
            bg, field = solve([CABLE], mode=mode, **kw)
        print('    {:24s} {:5.1f} deg off'.format(label, alignment(bg, field, CABLE)))
    print()
    print('    Tangent and perpendicular give the SAME number. That is correct,')
    print('    not a bug: a cross is 4-fold symmetric, so exp(i*4*theta) is')
    print('    unchanged by a 90 degree rotation and the two constraints are one')
    print('    constraint. They only separate once the frame is non-orthogonal')
    print('    -- two cable families crossing at 60 deg, case C, warp.py.')


# --- 2. and now the mesh does too --------------------------------------------

def part2():
    """The other half of the objective, on the same scale as part 1.

    This used to report the median edge angle of each mesh and print 0.0 for
    every row. That number was not measuring anything: the meshes it compared
    were byte-identical, because densification was a Coons patch between coarse
    edges and never consulted the field. The metric here is the one part 1 uses
    on the field -- angle to the cable, sampled along the cable -- so the two
    halves are directly comparable.
    """
    print('--- 2. and the MESH now does too ---------------------------------')
    print()
    print('    Angle between dense mesh edges and the cable, sampled along it,')
    print('    on the same scale as part 1. Two reference points from there:')
    print('      a hard-constrained cable puts the FIELD    0.0 deg off')
    print('      no cable at all leaves it                 36.9 deg off')
    print()
    print('    {:26s} {:>4s} {:>5s} {:>8s} {:>6s} {:>9s}'.format(
        'cable', 'sing', 'seps', 'patches', 'quads', 'off cable'))
    for label, guides, kw in [
        ('no cable', None, {}),
        ('cable, soft', [CABLE], {}),
        ('cable, hard', [CABLE], dict(guide_weight=None)),
        ('cable, hard + wide band', [CABLE], dict(guide_band=3.0, guide_weight=5.0)),
    ]:
        d = FieldDecomposition.from_boundary(SQUARE, guides=guides,
                                             target_length=0.5, **kw)
        dense = d.quad_mesh(target_length=1.0)
        print('    {:26s} {:4d} {:5d} {:8d} {:6d} {:9.1f}'.format(
            label, len(d.field.singularities()), len(d.separatrices),
            d.mesh.number_of_faces(), dense.number_of_faces(),
            curve_alignment(dense, CABLE)))

    print()
    print('    The chain, end to end:')
    print('      cable -> field        yes, measured in part 1')
    print('      field -> LAYOUT       only through singularities and corners')
    print('      layout -> mesh        yes -- framefield/densify.py integrates')
    print('                            the field INSIDE each patch instead of')
    print('                            blending its four edges')
    print()
    print('    Still one patch: on a square the cable creates no singularity, so')
    print('    the LAYOUT does not change. It no longer has to. The cable now')
    print('    steers the patch interior, which is what option (b) in part 3 was.')
    print()
    print('    It does not reach 0.0, and the reason is geometric rather than a')
    print('    shortfall in the integrator. The patch boundary is fixed -- here')
    print('    it IS the four walls -- so the mesh must turn to meet them however')
    print('    well the field is integrated. Where the cable is clear of a wall')
    print('    the mesh is much closer to it than the mean suggests:')
    print()
    d = FieldDecomposition.from_boundary(SQUARE, guides=[CABLE],
                                         target_length=0.5,
                                         guide_band=3.0, guide_weight=None)
    dense = d.quad_mesh(target_length=1.0)
    plain = FieldDecomposition.from_boundary(SQUARE, target_length=0.5)
    flat = plain.quad_mesh(target_length=1.0)
    guided = curve_alignment_profile(dense, CABLE)
    unguided = curve_alignment_profile(flat, CABLE)
    print('        {:>18s} {:>10s} {:>10s}'.format(
        'point on cable', 'unguided', 'guided'))
    for (pt, a), (_, b) in zip(unguided, guided):
        print('        {:>18s} {:10.1f} {:10.1f}'.format(
            '({:.1f}, {:.1f})'.format(*pt), a, b))


# --- 3. what would have to change --------------------------------------------

def part3():
    print('--- 3. the two ways to close it ----------------------------------')
    d = FieldDecomposition.from_boundary(SQUARE, guides=[RING],
                                         target_length=0.5,
                                         guide_band=1.5, guide_weight=5.0)
    try:
        coarse = d.decomposition_mesh()
        sides = {}
        for f in coarse.faces():
            n = len(coarse.face_vertices(f))
            sides[n] = sides.get(n, 0) + 1
        state = '{} patches, sides {}'.format(coarse.number_of_faces(), sides)
    except Exception as exc:
        state = '{}: {}'.format(type(exc).__name__, str(exc)[:60])
    print('    (a) make the cable change the LAYOUT.')
    print('        A closed hoop forces winding the walls cannot absorb, so')
    print('        singularities appear and the layout does change:')
    print('          ring cable: {} singularities, {} separatrices, {}'.format(
        len(d.field.singularities()), len(d.separatrices), state))
    print('        Real, but a blunt instrument -- it only works when the cable')
    print('        forces topology, and the patches then need repair.')
    print()
    print('    (b) make DENSIFICATION field-aware.  <-- DONE, this is part 2')
    print('        Replace discrete_coons_patch inside a patch with an')
    print('        integration of the field, so a cable steers the interior of a')
    print('        patch without needing to split it. The cable is an input to')
    print('        the mesh, not a post-hoc correction. It lives in')
    print('        framefield/densify.py, so compas_singular is untouched.')
    print()
    print('        (a) was NOT used to get there, and this ring is why. Its')
    print('        layout is the only one in the suite a cable really reaches,')
    print('        and it is in poor shape: the pole faces above carry the')
    print("        suite's three remaining hard-floor breaches, and at the two")
    print('        finer backgrounds the layout does not close at all and falls')
    print('        to the polygon route. densify.py routes AROUND that -- a')
    print('        pseudo-quad patch keeps its Coons interior, because a pole')
    print('        sits at a singularity and the field has no continuous branch')
    print('        there to integrate. The ring cable is therefore essentially')
    print('        unchanged by field-aware densification, deliberately.')


def show():
    """The two halves of the objective, side by side.

    Panels 1 and 2 are the FIELD with and without the cable -- the cable
    clearly bends it. Panels 3 and 4 are the resulting MESHES, and they now
    differ: panel 4's elements swing onto the purple diagonal in the middle of
    the plate and unwind to meet the walls. Those two panels were byte-identical
    for the whole life of this project; that they are not any more is the
    deliverable.
    """
    from compas_viewer import Viewer

    from compas_singular.framefield.viz import (Grid, add_background, add_boundary, add_dense,
                                                add_field, add_guides, add_layout,
                                                add_singularities)

    viewer = Viewer()
    grid = Grid(pitch=14.0, cols=2)

    def built(guides=None, **kw):
        d = FieldDecomposition.from_boundary(SQUARE, guides=guides,
                                             target_length=0.5, **kw)
        return d, d.quad_mesh(target_length=1.0)

    d0, mesh0 = built()
    d1, mesh1 = built([CABLE], guide_band=3.0, guide_weight=None)

    dx, dy = grid.cell(0)
    g = viewer.scene.add_group(name='1. FIELD, no cable -- constant')
    add_background(g, d0, dx, dy)
    add_field(g, d0, dx, dy)
    add_boundary(g, d0, dx, dy)

    dx, dy = grid.cell(1)
    g = viewer.scene.add_group(name='2. FIELD, hard cable -- bends onto it, 0.0 deg off')
    add_background(g, d1, dx, dy)
    add_field(g, d1, dx, dy)
    add_boundary(g, d1, dx, dy)
    add_guides(g, [CABLE], dx, dy)

    dx, dy = grid.cell(2)
    g = viewer.scene.add_group(
        name='3. MESH, no cable -- {} quads, {:.1f} deg off the cable'.format(
            mesh0.number_of_faces(), curve_alignment(mesh0, CABLE)))
    add_dense(g, mesh0, dx, dy)

    dx, dy = grid.cell(3)
    g = viewer.scene.add_group(
        name='4. MESH, hard cable -- {} quads, {:.1f} deg off the cable'.format(
            mesh1.number_of_faces(), curve_alignment(mesh1, CABLE)))
    add_dense(g, mesh1, dx, dy)
    add_guides(g, [CABLE], dx, dy)

    # a cable that DOES change the topology, for contrast
    dr = FieldDecomposition.from_boundary(SQUARE, guides=[RING], target_length=0.5,
                                          guide_band=1.5, guide_weight=5.0)
    dx, dy = grid.cell(4)
    g = viewer.scene.add_group(
        name='5. ring cable -- {} singularities forced by the winding'.format(
            len(dr.field.singularities())))
    add_background(g, dr, dx, dy)
    add_field(g, dr, dx, dy)
    add_boundary(g, dr, dx, dy)
    add_singularities(g, dr, dx, dy)
    add_guides(g, [RING], dx, dy)

    dx, dy = grid.cell(5)
    g = viewer.scene.add_group(name='6. ring cable layout -- topology changes, '
                                    'patches need repair')
    try:
        dr.decomposition_mesh()
        add_layout(g, dr, dx, dy)
    except Exception:
        add_boundary(g, dr, dx, dy)
    add_singularities(g, dr, dx, dy)
    add_guides(g, [RING], dx, dy)

    grid.frame(viewer)
    viewer.show()


def main():
    part1()
    print()
    part2()
    print()
    part3()
    print()
    print('IN SHORT')
    print('  Both halves of the objective now work. The cable shapes the field')
    print('  (part 1) and the field shapes the mesh (part 2), with no step in')
    print('  between that pulls an existing mesh onto a curve. guide_lines.py is')
    print('  still not the answer and was not used: it can only reuse directions')
    print('  the layout already has. framefield/densify.py invents them.')

    if '--no-view' not in sys.argv:
        show()
    return 0


if __name__ == '__main__':
    sys.exit(main())
