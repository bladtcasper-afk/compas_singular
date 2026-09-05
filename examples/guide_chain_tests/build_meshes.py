"""The dense quad meshes and guide curves the guide-chain measurement runs on.

Built headless with the recipe the Rhino commands use -- boundary triangulation,
skeleton decomposition, strip densification -- so the thing being measured is a mesh of
the kind a guide is actually attached to: singularities where the decomposition put them,
edge lengths that vary, and a boundary that a chain can be tempted onto.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from math import cos
from math import pi
from math import sin

from compas.geometry import Line

from compas_singular.algorithms import SkeletonDecomposition
from compas_singular.algorithms import boundary_triangulation


def densify(points, n=5):
    """A closed polygon as a densified point loop, corners kept."""
    out = []
    loop = list(points) + [points[0]]
    for a, b in zip(loop[:-1], loop[1:]):
        out.extend([list(point)[:3] for point in Line(a, b).to_polyline(n=n).points[:-1]])
    return out


def build(outer, holes=(), density=4):
    """Outer boundary (+ holes) -> coarse decomposition -> dense quad mesh."""
    trimesh = boundary_triangulation(
        outer_boundary=densify(outer),
        inner_boundaries=[densify(hole) for hole in holes],
        point_features=[])
    decomposition = SkeletonDecomposition.from_mesh(trimesh)
    coarse = decomposition.decomposition_mesh([])
    coarse.collect_strips()
    coarse.set_strips_density(d=density)
    coarse.densification()
    return coarse.get_quad_mesh()


def arc(centre, radius, start, end, count=24):
    """A circular arc as points. ``start``/``end`` in degrees."""
    points = []
    for index in range(count + 1):
        angle = (start + (end - start) * index / float(count)) * pi / 180.0
        points.append([centre[0] + radius * cos(angle), centre[1] + radius * sin(angle), 0.0])
    return points


SQUARE = [[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]]
HOLE = [[4, 4, 0], [6, 4, 0], [6, 6, 0], [4, 6, 0]]
L_PLATE = [[0, 0, 0], [10, 0, 0], [10, 4, 0], [4, 4, 0], [4, 10, 0], [0, 10, 0]]


def cases():
    """[(mesh name, mesh, [(guide name, what it is for, points)])].

    Every guide has a known right answer, so the table can be read without opening a
    viewer:

    * ``axis``     -- runs the way a course does: must come out one whole polyedge.
    * ``short``    -- half the length: the chain must be TRIMMED, not extended to the wall.
    * ``diagonal`` -- 45 degrees to both families: no polyedge follows it, and the honest
      answer is a short, poorly aligned stub rather than a staircase.
    * ``arc``      -- curved: a polyedge cannot bend, so coverage is the question.
    * ``hug``      -- one row in from the wall: the boundary rule.
    * ``wall``     -- drawn ON the wall: the pathological case, and what saves the outline.
    * ``ring``     -- closed, around the hole: the arc-length wrap.
    * ``through``  -- crosses the hole: half the guide has no mesh under it at all.
    * ``elbow``    -- a polyline with a corner in it, following the reflex corner.
    """
    square = build(SQUARE)
    holed = build(SQUARE, [HOLE])
    plate = build(L_PLATE)

    return [
        ('square', square, [
            ('axis', [[0.2, 5.0, 0.0], [9.8, 5.0, 0.0]]),
            ('short', [[3.0, 5.0, 0.0], [7.0, 5.0, 0.0]]),
            ('diagonal', [[0.5, 0.5, 0.0], [9.5, 9.5, 0.0]]),
            ('arc', arc([5.0, -3.0], 8.0, 60.0, 120.0)),
            ('hug', [[0.5, 1.1, 0.0], [9.5, 1.1, 0.0]]),
            ('wall', [[0.5, 0.0, 0.0], [9.5, 0.0, 0.0]]),
        ]),
        ('square+hole', holed, [
            ('axis', [[0.2, 2.0, 0.0], [9.8, 2.0, 0.0]]),
            ('through', [[0.2, 5.0, 0.0], [9.8, 5.0, 0.0]]),
            ('ring', arc([5.0, 5.0], 3.5, 0.0, 360.0, count=48)),
            ('hug', [[0.5, 0.7, 0.0], [9.5, 0.7, 0.0]]),
        ]),
        ('L-plate', plate, [
            ('axis', [[0.2, 2.0, 0.0], [9.8, 2.0, 0.0]]),
            ('elbow', [[0.5, 2.0, 0.0], [2.0, 2.0, 0.0], [2.0, 9.5, 0.0]]),
            ('diagonal', [[0.5, 0.5, 0.0], [3.5, 3.5, 0.0]]),
            ('hug', [[0.5, 0.6, 0.0], [9.5, 0.6, 0.0]]),
        ]),
    ]
