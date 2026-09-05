"""EXAMPLE 1 -- what the field front end can decompose today.

The whole pipeline in four lines:

    d = FieldDecomposition.from_boundary(outline)
    coarse = d.decomposition_mesh()
    coarse.collect_strips()
    coarse.densification(edges_to_curves=d.edges_to_curves())

Everything after ``decomposition_mesh`` is stock compas_singular. The only thing
replaced is where the coarse layout comes from: a cross field's separatrices
instead of a medial-axis skeleton.

The domains below are the ones that work end to end. For the ones that do not,
see ``13_limits.py`` -- run that before trusting this on your own outline.

Opens a compas_viewer scene: one group per domain on a grid, each with its dense
mesh, coarse patch layout and singularities. Toggle groups in the scene tree.

Run:
    python 10_domains.py
    python 10_domains.py pentagon disc     # only these domains
    python 10_domains.py --no-view         # numbers only
"""
import os
import sys
from math import cos, pi, sin

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from compas_singular.framefield.decomposition import FieldDecomposition  # noqa: E402


# --- domains -----------------------------------------------------------------

SQUARE = [[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]]

L_SHAPE = [[0, 0, 0], [12, 0, 0], [12, 4, 0], [5, 4, 0], [5, 10, 0], [0, 10, 0]]

PENTAGON = [[0, 0, 0], [10, 0, 0], [13, 7, 0], [5, 12, 0], [-2, 7, 0]]

SQUARE_HOLE = [[4, 4, 0], [6, 4, 0], [6, 6, 0], [4, 6, 0]]

T_SHAPE = [[0, 0, 0], [12, 0, 0], [12, 4, 0], [8, 4, 0],
           [8, 10, 0], [4, 10, 0], [4, 4, 0], [0, 4, 0]]

U_SHAPE = [[0, 0, 0], [12, 0, 0], [12, 10, 0], [9, 10, 0],
           [9, 4, 0], [3, 4, 0], [3, 10, 0], [0, 10, 0]]

PLUS = [[4, 0, 0], [8, 0, 0], [8, 4, 0], [12, 4, 0], [12, 8, 0], [8, 8, 0],
        [8, 12, 0], [4, 12, 0], [4, 8, 0], [0, 8, 0], [0, 4, 0], [4, 4, 0]]

DISC = [[5 + 5 * cos(2 * pi * i / 48), 5 + 5 * sin(2 * pi * i / 48), 0.0]
        for i in range(48)]

ROUND_HOLE = [[5 + 2 * cos(-2 * pi * i / 48), 5 + 2 * sin(-2 * pi * i / 48), 0.0]
              for i in range(48)]

HEXAGON = [[5 + 5 * cos(2 * pi * i / 6), 5 + 5 * sin(2 * pi * i / 6), 0.0]
           for i in range(6)]

ELLIPSE = [[5 + 7 * cos(2 * pi * i / 60), 5 + 3 * sin(2 * pi * i / 60), 0.0]
           for i in range(60)]


#: ``target_length`` is the BACKGROUND TRIANGULATION spacing, not the quad size.
#: The quad size is set later, on the coarse mesh, by
#: ``set_strips_density_target``. They are independent, and the pentagon is
#: pinned to 0.5 deliberately -- see ``13_limits.py``, resolution sensitivity.
DOMAINS = [
    ('square', SQUARE, None, 0.5,
     'No singularities, no reflex corners: one quad. The correct answer, but '
     'from_polylines cannot produce it -- see decomposition._single_patch.'),
    ('L-shape', L_SHAPE, None, 0.5,
     'No singularities either. The two cuts come from the REFLEX CORNER, which '
     'emits k-1 = 2 separatrices. Gives the 3-rectangle decomposition.'),
    ('T-plate', T_SHAPE, None, 0.5,
     'Two reflex corners. They launch the same line from both ends, so it is '
     'traced twice and deduplicated; 4 patches.'),
    ('U-plate', U_SHAPE, None, 0.5,
     'Two reflex corners again, further apart. 5 patches.'),
    ('plus-plate', PLUS, None, 0.5,
     'Four reflex corners, four duplicate pairs, 5 patches -- the cross of the '
     'plus plus its four arms.'),
    ('pentagon', PENTAGON, None, 0.5,
     'Corners of ~108 deg cannot absorb the boundary turning, so one valence-5 '
     'singularity appears in the interior and fans 5 separatrices to the wall.'),
    ('disc', DISC, None, 0.5,
     'A smooth boundary has no corners to absorb anything, so the whole +4 goes '
     'into the interior as four +1/4 singularities.'),
    ('square with hole', SQUARE, [SQUARE_HOLE], 0.4,
     'Each hole corner is REFLEX seen from the domain (270 deg) and emits 2 '
     'separatrices. 8 patches, no interior singularity needed.'),
    ('hexagon', HEXAGON, None, 0.5,
     'Two separatrices cross in their interiors. from_polylines computes no '
     'intersection, so that used to become a face repeating a vertex and the '
     'layout was discarded; arrangement.py inserts the node and it is now the '
     'FIELD route. 11 patches -- it was 43, because the two singularities '
     'landed 0.35 apart on the top wall and the landing tolerance merged them, '
     'collapsing one quad to a triangle that then quad-split the whole layout. '
     'A landing is pinned to the wall, so repair._cluster no longer merges '
     'two of them. 11 / 11 / 8 at tl 0.5 / 0.4 / 0.35, and monotone at last.'),
    ('ellipse', ELLIPSE, None, 0.5,
     'Also the FIELD route now. 10 patches -- it was 30, because the two arms '
     'leaving each tip singularity run to the same far tip and the duplicate '
     'test merged them, leaving a valence-3 singularity with two arms and a '
     'pentagon either side. Arms of one singularity leave in different '
     'directions; two traces of one connection do not. Boundary-aware '
     'splitting keeps new vertices ON the curve, so coverage stays exact.'),
    ('disc with round hole', DISC, [ROUND_HOLE], 0.5,
     'Multiply connected with no corners anywhere: no route but the '
     'TRIANGULATION backstop. Ugly and high-valence -- and a mesh.'),
]


def build(name, boundary, holes, target_length, density=0.8, guides=None):
    d = FieldDecomposition.from_boundary(
        boundary, holes, guides=guides, target_length=target_length)
    # quad_mesh owns the whole chain and every fallback in it: it ALWAYS returns
    # an all-quad mesh of the domain. d.route() says how it got there.
    dense = d.quad_mesh(target_length=density)
    return d, d.mesh, dense


def main():
    wanted = [a for a in sys.argv[1:] if not a.startswith('--')]

    scenes = []
    print('{:18s} {:>5s} {:>5s} {:>7s} {:>7s}  {}'.format(
        'domain', 'sing', 'seps', 'patches', 'quads', 'route'))
    print('-' * 60)

    for name, boundary, holes, target_length, _note in DOMAINS:
        if wanted and name not in wanted:
            continue
        try:
            d, coarse, dense = build(name, boundary, holes, target_length)
        except Exception as exc:
            print('{:18s} FAILED {}: {}'.format(name, type(exc).__name__, exc))
            continue

        for w in d.warnings():
            print('    warning: {}'.format(w))

        print('{:18s} {:5d} {:5d} {:7d} {:7d}  {}'.format(
            name, len(d.field.singularities()), len(d.separatrices),
            coarse.number_of_faces(), dense.number_of_faces(), d.route()))
        scenes.append((name, d, coarse, dense))

    print()
    for name, boundary, holes, tl, note in DOMAINS:
        if wanted and name not in wanted:
            continue
        print('{}:\n    {}'.format(name, note))

    # density is a property of the COARSE mesh, not of the field -- the same
    # layout drives any resolution
    print()
    d, coarse, _ = build('pentagon', PENTAGON, None, 0.5)
    for target in (2.0, 1.0, 0.5):
        print('pentagon at quad size {:.1f}: {} quads from the same {} patches'.format(
            target, d.quad_mesh(target_length=target).number_of_faces(),
            coarse.number_of_faces()))

    if '--no-view' not in sys.argv and scenes:
        show(scenes)
    return 0


def show(scenes):
    """One group per domain, laid out on a grid.

    Each panel carries the dense mesh with its coarse layout drawn over it, so
    the patches the field produced and the quads they densify into are visible
    together. Toggle groups in the viewer's scene tree to isolate a stage.
    """
    from compas_viewer import Viewer

    from compas_singular.framefield.viz import Grid, add_dense, add_layout, add_singularities

    viewer = Viewer()
    grid = Grid(pitch=17.0, cols=3)

    for i, (name, d, coarse, dense) in enumerate(scenes):
        dx, dy = grid.cell(i)
        group = viewer.scene.add_group(name='{} -- {} patches, {} quads'.format(
            name, coarse.number_of_faces(), dense.number_of_faces()))
        add_dense(group, dense, dx, dy)
        add_layout(group, d, dx, dy, faces=False)
        add_singularities(group, d, dx, dy)

    grid.frame(viewer)
    viewer.show()


if __name__ == '__main__':
    sys.exit(main())
