"""Milestone 1c -- separatrix network to coarse quad mesh.

Checks step 6 of the plan: the network reaches
``CoarsePseudoQuadMesh.from_polylines``, the result is a quad mesh, and the
EXISTING downstream pipeline accepts it -- ``collect_strips`` then
``densification`` -- which is the whole point of the seam.

Opens a compas_viewer scene with each coarse layout over its dense mesh.

Run:
    python 03_layout.py
    python 03_layout.py --no-view          # checks only
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from compas_singular.framefield.decomposition import FieldDecomposition  # noqa: E402

SQUARE = [[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]]
PENTAGON = [[0, 0, 0], [10, 0, 0], [13, 7, 0], [5, 12, 0], [-2, 7, 0]]
L_SHAPE = [[0, 0, 0], [12, 0, 0], [12, 4, 0], [5, 4, 0], [5, 10, 0], [0, 10, 0]]


def circle(radius=5.0, n=48, centre=(5.0, 5.0)):
    from math import cos, pi, sin
    return [[centre[0] + radius * cos(2 * pi * i / n),
             centre[1] + radius * sin(2 * pi * i / n), 0.0] for i in range(n)]


def run(name, boundary, holes=None, guides=None, target_length=0.5):
    print('--- {} ---'.format(name))
    d = FieldDecomposition.from_boundary(boundary, holes, guides=guides,
                                         target_length=target_length)
    b, o, info = d._build()
    print('  network   : {} boundary arcs, {} separatrices'.format(len(b), len(o)))

    try:
        coarse = d.decomposition_mesh()
    except Exception as exc:
        print('  from_polylines FAILED: {}: {}'.format(type(exc).__name__, exc))
        return name, False, d, None, None

    faces = [len(coarse.face_vertices(f)) for f in coarse.faces()]
    counts = {}
    for n in faces:
        counts[n] = counts.get(n, 0) + 1
    print('  coarse    : {} vertices, {} faces, sides {}'.format(
        coarse.number_of_vertices(), coarse.number_of_faces(), counts))

    if not faces:
        print('  EMPTY MESH')
        return name, False, d, coarse, None

    ok = set(faces) == {4}
    dense = None
    try:
        coarse.collect_strips()
        coarse.set_strips_density_target(1.0)
        coarse.densification(edges_to_curves=d.edges_to_curves())
        dense = coarse.get_quad_mesh()
        print('  dense     : {} vertices, {} faces, manifold {}'.format(
            dense.number_of_vertices(), dense.number_of_faces(), dense.is_manifold()))
    except Exception as exc:
        print('  downstream FAILED: {}: {}'.format(type(exc).__name__, exc))
        ok = False

    print('  {}'.format('OK' if ok else 'NOT ALL QUADS'))
    return name, ok, d, coarse, dense


def main():
    results = []
    scenes = []
    for label, bnd, holes in [('square', SQUARE, None),
                              ('pentagon', PENTAGON, None),
                              ('disc', circle(), None),
                              ('L-shape', L_SHAPE, None),
                              ('square + guide', SQUARE, None)]:
        guides = [[[1.0, 2.0, 0.0], [5.0, 5.0, 0.0], [9.0, 8.0, 0.0]]]             if label == 'square + guide' else None
        out = run(label, bnd, holes, guides=guides,
                  target_length=0.4 if label == 'L-shape' else 0.5)
        results.append(out[:2])
        if len(out) > 2 and out[2] is not None:
            scenes.append((label,) + out[2:])

    print()
    bad = [n for n, ok in results if not ok]
    print('FAILED: {}'.format(', '.join(bad)) if bad else 'all checks passed')

    if '--no-view' not in sys.argv and scenes:
        show(scenes)
    return 1 if bad else 0


def show(scenes):
    """The coarse layout each network produced, over its dense mesh."""
    from compas_viewer import Viewer

    from compas_singular.framefield.viz import (Grid, add_dense, add_layout, add_singularities)

    viewer = Viewer()
    grid = Grid(pitch=16.0, cols=3)

    for i, (label, d, coarse, dense) in enumerate(scenes):
        dx, dy = grid.cell(i)
        group = viewer.scene.add_group(name='{} -- {} patches'.format(
            label, coarse.number_of_faces()))
        if dense is not None:
            add_dense(group, dense, dx, dy)
        add_layout(group, d, dx, dy, faces=dense is None)
        add_singularities(group, d, dx, dy)

    grid.frame(viewer)
    viewer.show()


if __name__ == '__main__':
    sys.exit(main())
