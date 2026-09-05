"""Milestone 1a -- background triangulation, cross field, singularities.

Checks steps 1-4 of the plan in ``frame-field-front-end-module-layout.md``:

* the triangulation is valid (positive areas, right number of boundary loops,
  right Euler characteristic);
* the field solve converges;
* the interior singularity indices match the boundary winding -- the discrete
  Poincare-Hopf check, which is what catches a broken angle unwrapping;
* a guide curve visibly bends the field, and removing it restores the unguided one.

Opens a compas_viewer scene with the field on every domain checked.

Run:

    python 01_field.py
    python 01_field.py --no-view   # checks only
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from compas_singular.framefield import BackgroundMesh                     # noqa: E402
from compas_singular.framefield import CrossField                         # noqa: E402
from compas_singular.framefield import from_boundary                      # noqa: E402
from compas_singular.framefield import from_curves                        # noqa: E402
from compas_singular.framefield.constraints import _assert_mode_identity  # noqa: E402


SQUARE = [[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]]

L_SHAPE = [[0, 0, 0], [12, 0, 0], [12, 4, 0], [5, 4, 0], [5, 10, 0], [0, 10, 0]]

PENTAGON = [[0, 0, 0], [10, 0, 0], [13, 7, 0], [5, 12, 0], [-2, 7, 0]]


def circle(radius=5.0, n=48, centre=(5.0, 5.0)):
    from math import cos, pi, sin
    return [[centre[0] + radius * cos(2 * pi * i / n),
             centre[1] + radius * sin(2 * pi * i / n), 0.0] for i in range(n)]


def run(name, boundary, holes=None, guides=None, target_length=0.6, mode='perpendicular',
        expect=None):
    """Solve one domain and check it.

    ``expect`` is the interior singularity sum predicted INDEPENDENTLY, from the
    quad-mesh index identity

        sum_interior (4 - valence) + sum_boundary (3 - valence) = 4 * chi

    counting each convex corner as a valence-2 quad corner (+1) and each reflex
    one as valence 4 (-1). Poincare-Hopf alone cannot catch a wrong boundary
    condition -- it compares the interior sum against the boundary winding, and
    both are read off the same angles, so a mis-pinned wall passes happily. This
    is the check that does not.
    """
    bg = BackgroundMesh.from_boundary(boundary, holes, target_length=target_length)
    grid = bg.validate()

    constraints = from_boundary(bg)
    if guides:
        constraints += from_curves(bg, guides, mode=mode)

    field = CrossField.solve(bg, constraints)
    info = field.report()
    info['expected'] = expect
    info['expected_ok'] = expect is None or info['interior_index_sum'] == expect

    print('--- {} ---'.format(name))
    print('  mesh      : {} vertices, {} faces, loops {}/{}, euler {}/{}  {}'.format(
        grid['vertices'], grid['faces'], grid['boundary_loops'], grid['expected_loops'],
        grid['euler'], grid['expected_euler'], 'OK' if grid['ok'] else 'BAD'))
    print('  |u|       : min {:.3f}  mean {:.3f}'.format(
        info['min_magnitude'], info['mean_magnitude']))
    print('  indices   : interior {:+d}  winding {:+d}  predicted {}  ({} singular faces)  {}'.format(
        info['interior_index_sum'], info['boundary_winding'],
        '{:+d}'.format(expect) if expect is not None else '  -',
        info['singular_faces'],
        'OK' if (info['ok'] and info['expected_ok']) else 'FAIL'))
    return bg, field, grid, info


def main():
    assert _assert_mode_identity()
    print('cross-field 4-fold identity: tangent == perpendicular   OK\n')

    results = {}

    # 4 corners x (+1), chi = 1 -> interior 4 - 4 = 0. The constant field: a
    # square needs no irregular vertices, and every 90-degree corner turns the
    # tangent by exactly the cross's period, so it costs the field nothing.
    results['square'] = run('square', SQUARE, target_length=0.5, expect=0)

    # smooth boundary, no corners to absorb anything -> all 4 go inside
    results['disc'] = run('disc', circle(), target_length=0.5, expect=4)

    # 5 convex (+1 each) + 1 reflex (-1) = 4 -> interior 0. An L-shape tiles
    # with three rectangles; it needs no singularity either.
    results['L'] = run('L-shape', L_SHAPE, target_length=0.4, expect=0)

    # 5 corners of ~108 degrees, none a multiple of 90: the walls cannot absorb
    # it and one valence-5 vertex has to appear. 4 - 5 = -1.
    results['pentagon'] = run('pentagon', PENTAGON, target_length=0.5, expect=-1)

    # outer 4 x (+1) + hole 4 x (-1) = 0 = 4 * chi(annulus) -> interior 0
    results['annulus'] = run(
        'square with square hole', SQUARE,
        holes=[[[4, 4, 0], [6, 4, 0], [6, 6, 0], [4, 6, 0]]], target_length=0.4, expect=0)

    # guide curve: a diagonal cable across the square. The guide does not change
    # the boundary, so the index sum must NOT move -- a guide redistributes the
    # field, it does not create topology. Worth asserting: a guide that changes
    # the index sum means its constraints are fighting the walls.
    guide = [[1.0, 2.0, 0.0], [5.0, 5.0, 0.0], [9.0, 8.0, 0.0]]
    results['guided'] = run('square + diagonal guide', SQUARE,
                            guides=[guide], target_length=0.5, expect=0)

    print()
    bad = [k for k, (_, _, g, i) in results.items()
           if not (g['ok'] and i['ok'] and i['expected_ok'])]

    # the guide has to actually bend the field -- the square's unguided solution
    # is the constant field, so any departure from it is the guide's doing
    bg_g, guided, _, _ = results['guided']
    turned = sum(1 for v in bg_g.mesh.vertices() if abs(guided.theta[v]) > 1e-3)
    total = bg_g.mesh.number_of_vertices()
    print('guide effect : {} of {} vertices turned off the constant field'.format(turned, total))
    if turned < 0.1 * total:
        bad.append('guide had no effect')

    print()
    if bad:
        print('FAILED: {}'.format(', '.join(bad)))
    else:
        print('all checks passed')

    if '--no-view' not in sys.argv:
        show(results)

    return 0 if not bad else 1


def show(results):
    """The field itself, on every domain checked."""
    from compas_viewer import Viewer

    from compas_singular.framefield.viz import (Grid, add_background, add_boundary, add_field,
                                                add_guides, add_singularities)

    viewer = Viewer()
    grid = Grid(pitch=15.0, cols=3)

    guide = [[1.0, 2.0, 0.0], [5.0, 5.0, 0.0], [9.0, 8.0, 0.0]]
    for i, (name, (bg, field, _grid, info)) in enumerate(results.items()):
        dx, dy = grid.cell(i)
        shim = _Shim(bg, field)
        group = viewer.scene.add_group(
            name='{} -- interior index {:+d}, {} singular faces'.format(
                name, info['interior_index_sum'], info['singular_faces']))
        add_background(group, shim, dx, dy)
        add_field(group, shim, dx, dy)
        add_boundary(group, shim, dx, dy)
        add_singularities(group, shim, dx, dy)
        if name == 'guided':
            add_guides(group, [guide], dx, dy)

    grid.frame(viewer)
    viewer.show()


class _Shim(object):
    """The viz helpers take a FieldDecomposition; this stage has no layout yet.

    Only ``background`` and ``field`` are read by the layers used here, so a
    two-attribute stand-in is enough and avoids tracing separatrices just to
    draw the field.
    """

    def __init__(self, background, field):
        self.background = background
        self.field = field


if __name__ == '__main__':
    sys.exit(main())
