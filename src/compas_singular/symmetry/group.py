"""**A finite symmetry group of the plane: rotations and mirrors about a centre.**

Every finite group of isometries that maps a bounded domain onto itself fixes a
point -- the centre -- and is one of two families:

* ``Cn`` -- the ``n`` rotations by multiples of ``2 pi / n``;
* ``Dn`` -- those rotations plus ``n`` mirrors, whose axes pass through the
  centre ``pi / n`` apart.

``C1`` is the trivial group and ``D1`` is a single mirror.

**Elements are named, and the names are geometry.** A rotation is ``R<degrees>``
and a mirror is ``M<axis angle in degrees>``, measured anticlockwise from +x
about the centre and rounded to 0.1 degree -- ``R90``, ``M45``, ``M7.5``. The
names are how a script picks a subgroup (:meth:`SymmetryGroup.subgroup`), so
they must not depend on the order anything was found in; for ``n <= 12`` two
axes are at least 15 degrees apart, so rounding never merges two of them.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from collections import namedtuple
from math import atan2
from math import cos
from math import degrees
from math import pi
from math import radians
from math import sin


__all__ = ['Element', 'SymmetryGroup', 'rotation_key', 'mirror_key']


#: One group element.
#:
#: ``kind`` is ``'identity'``, ``'rotation'`` or ``'mirror'``; ``angle`` is the
#: rotation angle or the mirror's AXIS angle, in radians; ``matrix`` is the linear
#: part ``(a, b, c, d)`` acting as ``x' = a x + b y``, ``y' = c x + d y`` about the
#: centre.
Element = namedtuple('Element', ['key', 'kind', 'angle', 'matrix', 'reflects'])


def _fmt(value):
    value = round(value, 1)
    if value == 0.0:
        value = 0.0          # no '-0'
    return '{:g}'.format(value)


def rotation_key(angle):
    """``R<degrees>`` for a rotation by ``angle`` radians, in [0, 360)."""
    deg = degrees(angle) % 360.0
    if abs(deg - 360.0) < 0.05:
        deg = 0.0
    return 'R' + _fmt(deg)


def mirror_key(axis):
    """``M<degrees>`` for a mirror whose axis is at ``axis`` radians, in [0, 180)."""
    deg = degrees(axis) % 180.0
    if abs(deg - 180.0) < 0.05:
        deg = 0.0
    return 'M' + _fmt(deg)


def _rotation(angle):
    c, s = cos(angle), sin(angle)
    kind = 'identity' if rotation_key(angle) == 'R0' else 'rotation'
    key = 'I' if kind == 'identity' else rotation_key(angle)
    return Element(key, kind, angle % (2 * pi), (c, -s, s, c), False)


def _mirror(axis):
    axis = axis % pi
    c, s = cos(2 * axis), sin(2 * axis)
    return Element(mirror_key(axis), 'mirror', axis, (c, s, s, -c), True)


class SymmetryGroup(object):
    """``Cn`` or ``Dn`` about a centre.

    Parameters
    ----------
    centre : [x, y, z]
    n : int
        Number of rotations, the identity included. ``1`` for no rotation.
    mirrors : bool
        ``True`` for ``Dn``.
    axis : float
        Angle of one mirror axis, radians. Only meaningful for ``Dn``; the
        others are at ``axis + k pi / n``. Stored as the smallest of them.
    """

    def __init__(self, centre=(0.0, 0.0, 0.0), n=1, mirrors=False, axis=0.0):
        self.centre = [float(centre[0]), float(centre[1]), 0.0]
        self.n = max(1, int(n))
        self.mirrors = bool(mirrors)
        self.axis = (axis % (pi / self.n)) if self.mirrors else 0.0

    # ------------------------------------------------------------------
    # identity
    # ------------------------------------------------------------------

    @classmethod
    def trivial(cls, centre=(0.0, 0.0, 0.0)):
        return cls(centre, 1, False)

    @property
    def name(self):
        return '{}{}'.format('D' if self.mirrors else 'C', self.n)

    @property
    def order(self):
        return self.n * (2 if self.mirrors else 1)

    @property
    def is_trivial(self):
        return self.order == 1

    def __len__(self):
        return self.order

    def __repr__(self):
        return '<SymmetryGroup {} about ({:.3f}, {:.3f}): {}>'.format(
            self.name, self.centre[0], self.centre[1], ' '.join(self.keys()) or '-')

    def __eq__(self, other):
        return (isinstance(other, SymmetryGroup) and self.name == other.name
                and sorted(self.keys()) == sorted(other.keys())
                and abs(self.centre[0] - other.centre[0]) < 1e-9
                and abs(self.centre[1] - other.centre[1]) < 1e-9)

    def __ne__(self, other):
        return not self.__eq__(other)

    # ------------------------------------------------------------------
    # elements
    # ------------------------------------------------------------------

    @property
    def elements(self):
        """Identity first, then the rotations by increasing angle, then the mirrors
        by increasing axis angle."""
        out = [_rotation(2 * pi * k / self.n) for k in range(self.n)]
        if self.mirrors:
            axes = sorted((self.axis + pi * k / self.n) % pi for k in range(self.n))
            out += [_mirror(a) for a in axes]
        return out

    def keys(self):
        """Every non-identity element's key."""
        return [e.key for e in self.elements if e.kind != 'identity']

    def element(self, key):
        for e in self.elements:
            if e.key == key:
                return e
        raise KeyError('{} is not an element of {} -- it has {}'.format(
            key, self.name, ', '.join(['I'] + self.keys())))

    @property
    def identity(self):
        return self.elements[0]

    def rotations(self):
        return [e for e in self.elements if not e.reflects]

    def mirror_elements(self):
        return [e for e in self.elements if e.reflects]

    # ------------------------------------------------------------------
    # action
    # ------------------------------------------------------------------

    def apply(self, element, point):
        """``element`` applied to a point, about :attr:`centre`."""
        a, b, c, d = element.matrix
        x, y = point[0] - self.centre[0], point[1] - self.centre[1]
        z = point[2] if len(point) > 2 else 0.0
        return [a * x + b * y + self.centre[0], c * x + d * y + self.centre[1], z]

    def apply_points(self, element, points):
        return [self.apply(element, p) for p in points]

    @staticmethod
    def linear(element, vector):
        a, b, c, d = element.matrix
        return [a * vector[0] + b * vector[1], c * vector[0] + d * vector[1], 0.0]

    def orbit(self, point):
        return [self.apply(e, point) for e in self.elements]

    def compose(self, first, second):
        """The element ``first o second`` -- ``second`` applied first."""
        a1, b1, c1, d1 = first.matrix
        a2, b2, c2, d2 = second.matrix
        matrix = (a1 * a2 + b1 * c2, a1 * b2 + b1 * d2,
                  c1 * a2 + d1 * c2, c1 * b2 + d1 * d2)
        return self._match(matrix)

    def inverse(self, element):
        if element.reflects:
            return element
        return self._match(_rotation(-element.angle).matrix)

    def _match(self, matrix, tol=1e-6):
        best, gap = None, None
        for e in self.elements:
            d = max(abs(x - y) for x, y in zip(e.matrix, matrix))
            if gap is None or d < gap:
                best, gap = e, d
        if gap > tol:
            raise ValueError('matrix {} is not an element of {}'.format(matrix, self.name))
        return best

    # ------------------------------------------------------------------
    # subgroups
    # ------------------------------------------------------------------

    def subgroup(self, keys=None):
        """The subgroup GENERATED by ``keys``.

        Arbitrary keys are generators, not a group: two mirrors 45 degrees apart
        generate the quarter turns as well. The result is closed under
        composition, so ``subgroup(['M0', 'M45'])`` of ``D4`` is all of ``D4``.
        ``None`` is the whole group; ``[]`` is the trivial group.
        """
        if keys is None:
            return SymmetryGroup(self.centre, self.n, self.mirrors, self.axis)
        if isinstance(keys, str):
            keys = [keys]
        generators = [self.element(k) for k in keys if k != 'I']
        closure = {self.identity.key: self.identity}
        frontier = list(closure.values())
        while frontier:
            fresh = []
            for e in frontier:
                for g in generators:
                    h = self.compose(g, e)
                    if h.key not in closure:
                        closure[h.key] = h
                        fresh.append(h)
            frontier = fresh
        rotations = [e for e in closure.values() if not e.reflects]
        mirrors = [e for e in closure.values() if e.reflects]
        axis = min(m.angle for m in mirrors) if mirrors else 0.0
        return SymmetryGroup(self.centre, len(rotations), bool(mirrors), axis)

    def subgroups(self):
        """Every subgroup, largest first -- what a user can choose between."""
        seen, out = set(), []
        rot_orders = [m for m in range(self.n, 0, -1) if self.n % m == 0]
        for m in rot_orders:
            if self.mirrors:
                for k in range(self.n // m):
                    axis = (self.axis + pi * k / self.n)
                    g = SymmetryGroup(self.centre, m, True, axis)
                    sig = (g.name, tuple(sorted(g.keys())))
                    if sig not in seen:
                        seen.add(sig)
                        out.append(g)
            g = SymmetryGroup(self.centre, m, False)
            sig = (g.name, tuple(sorted(g.keys())))
            if sig not in seen:
                seen.add(sig)
                out.append(g)
        out.sort(key=lambda g: -g.order)
        return out

    # ------------------------------------------------------------------
    # data
    # ------------------------------------------------------------------

    def to_data(self):
        return {'centre': [self.centre[0], self.centre[1]], 'n': self.n,
                'mirrors': self.mirrors, 'axis': degrees(self.axis),
                'name': self.name, 'keys': self.keys()}

    @classmethod
    def from_data(cls, data):
        return cls(data['centre'] + [0.0], data['n'], data['mirrors'],
                   radians(data.get('axis', 0.0)))


def angle_of(vector):
    return atan2(vector[1], vector[0])
