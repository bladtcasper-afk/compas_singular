"""EXAMPLE 5 -- the planar arrangement, and the route it buys.

``CoarsePseudoQuadMesh.from_polylines`` does not compute intersections. It hands
every polyline SEGMENT to ``Mesh.from_lines``, which keys vertices off segment
endpoints, so two separatrices crossing in their interiors contribute no node at
the crossing and the graph the face search walks is not a planar embedding. The
faces it recovers then VISIT THE SAME VERTEX TWICE, ``topological_quad_split``
turns each of those into a fan of inverted quads, and ``densifiable`` rejects the
layout at ``area2 <= 0.0``. The field layout is discarded and the domain is
meshed as a plain polygon -- which is why every curved domain used to stop
following its own field.

``framefield/arrangement.py`` computes the arrangement before ``from_polylines``
sees it. This script measures, per domain:

* **crossings** -- polyline pairs meeting anywhere but at a shared node. Must be 0.
* **dangling** -- network nodes only one polyline reaches. Must be 0.
* **repeat**   -- recovered faces visiting a vertex twice. Must be 0.
* **route**    -- ``field`` / ``polygon`` / ``triangulation``. Must be ``field``.
* **cover**    -- the dense mesh's area against the domain's. Must be ~100%.

``route`` is the number that matters. A mesh always came out before this change
too; it just was not a mesh of the field.

Run:
    python 14_arrangement.py
"""
import os
import sys
from math import cos, pi, sin

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from compas_singular.framefield.arrangement import (count_dangling_ends,  # noqa: E402
                                                    count_interior_crossings,
                                                    faces_with_repeated_vertices)
from compas_singular.framefield.decomposition import FieldDecomposition  # noqa: E402


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
HEXAGON = [[5 + 5 * cos(2 * pi * i / 6), 5 + 5 * sin(2 * pi * i / 6), 0.0]
           for i in range(6)]
ELLIPSE = [[5 + 7 * cos(2 * pi * i / 60), 5 + 3 * sin(2 * pi * i / 60), 0.0]
           for i in range(60)]
STADIUM = ([[0, 0, 0], [12, 0, 0]]
           + [[12 + 3 * cos(a), 3 + 3 * sin(a), 0.0]
              for a in [-pi / 2 + pi * i / 10 for i in range(1, 10)]]
           + [[12, 6, 0], [0, 6, 0]])


#: The curved domains that used to fall back, at several resolutions each --
#: the hexagon failed at 0.5 and 0.4 and worked at 0.35, and a geometric limit
#: does not come and go with background resolution. Then the rectilinear plates
#: and the pentagon, which always worked and must not regress.
CASES = [
    ('disc', DISC, None, 0.6),
    ('disc', DISC, None, 0.5),
    ('disc', DISC, None, 0.4),
    ('disc', DISC, None, 0.3),
    ('hexagon', HEXAGON, None, 0.5),
    ('hexagon', HEXAGON, None, 0.4),
    ('hexagon', HEXAGON, None, 0.35),
    ('ellipse', ELLIPSE, None, 0.5),
    ('stadium', STADIUM, None, 0.5),
    ('L-shape', L_SHAPE, None, 0.5),
    ('T-plate', T_SHAPE, None, 0.5),
    ('U-plate', U_SHAPE, None, 0.5),
    ('plus-plate', PLUS, None, 0.5),
    ('pentagon', PENTAGON, None, 0.5),
    ('square', SQUARE, None, 0.5),
    ('square with hole', SQUARE, [SQUARE_HOLE], 0.4),
]


def polygon_area(loop):
    pts = list(loop)
    return abs(0.5 * sum(p[0] * q[1] - q[0] * p[1]
                         for p, q in zip(pts, pts[1:] + pts[:1])))


def coverage(d, mesh):
    want = polygon_area(d.background.outer)
    want -= sum(polygon_area(loop) for loop in d.background.inners)
    got = 0.0
    for fkey in mesh.faces():
        pts = [mesh.vertex_coordinates(v) for v in mesh.face_vertices(fkey)]
        got += 0.5 * sum(p[0] * q[1] - q[0] * p[1]
                         for p, q in zip(pts, pts[1:] + pts[:1]))
    return abs(got) / want if want else 1.0


def main():
    print('{:18s} {:>4s} {:>5s} {:>4s} {:>4s} {:>7s} {:>6s} {:>6s}  {}'.format(
        'domain', 'tl', 'cross', 'dang', 'rep', 'patches', 'quads', 'cover',
        'route'))
    print('-' * 84)

    failed = []
    for label, bnd, holes, tl in CASES:
        try:
            d = FieldDecomposition.from_boundary(bnd, holes, target_length=tl)
            boundary, others, _info = d._build()
            # the node tolerance the arrangement itself worked to; anything
            # closer than this it was entitled to merge
            tol = d.background.target_length * 0.2
            cross = count_interior_crossings(boundary + others, tol)
            dang = count_dangling_ends(boundary, others, tol)

            coarse = d.decomposition_mesh()
            rep = faces_with_repeated_vertices(coarse)
            dense = d.quad_mesh(target_length=0.8)
            route = d.route()
            cover = coverage(d, dense)
        except Exception as exc:
            print('{:18s} {:4.2f}  EXCEPTION {}: {}'.format(
                label, tl, type(exc).__name__, str(exc)[:40]))
            failed.append('{} @{}'.format(label, tl))
            continue

        bad = (route != 'field' or cross or dang or rep
               or not 0.97 <= cover <= 1.03)
        print('{:18s} {:4.2f} {:5d} {:4d} {:4d} {:7d} {:6d} {:5.0f}%  {}{}'.format(
            label, tl, len(cross), len(dang), len(rep),
            coarse.number_of_faces(), dense.number_of_faces(), cover * 100,
            route, '   <-- FAIL' if bad else ''))
        if bad:
            failed.append('{} @{}'.format(label, tl))

    print()
    if failed:
        print('NOT on route field / not clean: {}'.format(', '.join(failed)))
        return 1
    print('all {} domains reach route "field" with a clean arrangement'.format(
        len(CASES)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
