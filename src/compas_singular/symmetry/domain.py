"""The input every meshing route takes: walls, holes, guides and poles."""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ._geometry import area_centroid
from ._geometry import as_points
from ._geometry import bbox_diagonal
from ._geometry import open_loop
from ._geometry import point_in_polygon


__all__ = ['Domain']


class Domain(object):
    """``outer``, ``inners``, ``guides`` and ``poles`` as plain point lists.

    Loops are stored OPEN (the closing point is not repeated), guides as open
    polylines, poles as points. A compas ``Polyline`` or ``Point`` is accepted
    anywhere a point list or point is.

    ``guides`` means any curve that is not a wall: a cable for the field route, a
    polyline feature for the skeleton route. Symmetry does not care which.
    """

    def __init__(self, outer=None, inners=None, guides=None, poles=None):
        self.outer = open_loop(as_points(outer)) if outer is not None else []
        self.inners = [open_loop(as_points(loop)) for loop in (inners or [])]
        self.guides = [as_points(curve) for curve in (guides or [])]
        self.guides = [curve for curve in self.guides if len(curve) >= 2]
        self.poles = [[float(p[0]), float(p[1]), 0.0] for p in (poles or [])]

    @property
    def diagonal(self):
        return bbox_diagonal([self.outer] + self.inners + self.guides + [self.poles])

    def area_centroid(self):
        """The centroid of the REGION: outer minus holes. ``None`` with no outer."""
        if len(self.outer) < 3:
            return None
        area, centre = area_centroid(self.outer)
        cx, cy = centre[0] * area, centre[1] * area
        total = area
        for loop in self.inners:
            if len(loop) < 3:
                continue
            a, c = area_centroid(loop)
            cx -= c[0] * a
            cy -= c[1] * a
            total -= a
        if total <= 1e-300:
            return [centre[0], centre[1], 0.0]
        return [cx / total, cy / total, 0.0]

    def contains(self, point):
        if len(self.outer) < 3 or not point_in_polygon(point[0], point[1], self.outer):
            return False
        return not any(point_in_polygon(point[0], point[1], loop) for loop in self.inners if len(loop) >= 3)

    def to_data(self):
        return {'outer': self.outer, 'inners': self.inners,
                'guides': self.guides, 'poles': self.poles}

    @classmethod
    def from_data(cls, data):
        return cls(data.get('outer'), data.get('inners'), data.get('guides'), data.get('poles'))
