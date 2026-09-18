"""**Symmetry detection: which rotations and mirrors map the input onto itself.**

HOW
---

1. **Centre = area centroid of the region** (outer minus holes). Every isometry
   that maps a region onto itself fixes its centroid, so there is exactly one
   candidate, for any order. The bounding-box centre the old detector used is
   right for square-like shapes and wrong for a triangle.

2. **The matching test.** An element ``g`` is applied to every sample -- each
   vertex and each segment midpoint of every included curve -- and the distance
   from the image to the nearest SEGMENT of the curves of the same kind (walls
   to walls, holes to holes, guides to guides, poles to poles) is measured.
   ``g`` is accepted when the worst distance is within ``tol``. A segment that
   samples a smooth curve (no vertex turning 22.5 degrees or more) is allowed
   its chord sag on top: an image landing between the chord and the arc is on
   the curve as far as the samples can say. Corners get no slack, so a polygon
   is matched exactly.

   These are real distances from :class:`~._geometry.SegmentHash`, which
   searches every cell a query can reach. The old detector rounded coordinates
   into buckets, so two points 1e-9 apart on either side of a bucket edge
   disagreed: it lost symmetries at +-1e-7 noise and once rejected the identity.
   Measuring to a segment rather than to a sample also makes the sampling
   irrelevant -- a circle divided into 26 points, or rotated 5 degrees, still
   matches its own images.

3. **Rotations** by ``2 pi / n`` for ``n = 2 .. max_order``; the order is the
   largest ``n`` whose divisors all pass.

4. **Mirrors.** A mirror axis passes through the centre and, for a polygon,
   through a vertex or an edge midpoint -- where the distance from the centre is
   locally extreme. Those directions and the bisectors between them are the
   candidates, reduced modulo ``pi / n``, each tested with step 2 (on a small
   subset first, then in full).

5. **The group** is the largest subgroup of the result whose every element
   passes: elements tested one at a time near the tolerance need not compose.

Nothing is enforced here. See :mod:`.unit` for that.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from math import atan2
from math import hypot
from math import pi

from ._geometry import SegmentHash
from ._geometry import chord_sags
from .domain import Domain
from .group import SymmetryGroup
from .group import _mirror
from .group import _rotation
from .report import SymmetryReport


__all__ = ['find_symmetry', 'Matcher']


INCLUDE = ('walls', 'holes', 'guides', 'poles')

_KIND = {'walls': 'outer', 'holes': 'hole', 'guides': 'guide', 'poles': 'pole'}


class Matcher(object):
    """Samples and per-kind segment hashes for one domain, about one centre."""

    def __init__(self, domain, centre, include, near):
        self.centre = centre
        self.near = near
        self.samples = []
        include = set(include)
        cell = max(near, 1e-9)
        self.hashes = {}

        def add_curve(points, closed, kind, tag):
            h = self.hashes.setdefault(kind, SegmentHash(cell))
            sags = chord_sags(points, closed)
            h.add_polyline(points, closed, tag, sags)
            n = len(points)
            last = n if closed else n - 1
            for i, p in enumerate(points):
                self.samples.append((p[0], p[1], kind, tag, 0.0))
            for i in range(last):
                a, b = points[i], points[(i + 1) % n]
                # A midpoint is itself off the curve by its own segment's sag.
                self.samples.append((0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]), kind, tag, sags[i]))

        if 'walls' in include and len(domain.outer) >= 2:
            add_curve(domain.outer, True, 'outer', 'outer')
        if 'holes' in include:
            for i, loop in enumerate(domain.inners):
                add_curve(loop, True, 'hole', 'hole {}'.format(i + 1))
        if 'guides' in include:
            for i, curve in enumerate(domain.guides):
                add_curve(curve, False, 'guide', 'guide {}'.format(i + 1))
        if 'poles' in include:
            h = self.hashes.setdefault('pole', SegmentHash(cell))
            for i, p in enumerate(domain.poles):
                tag = 'pole {}'.format(i + 1)
                h.add_point(p, tag)
                self.samples.append((p[0], p[1], 'pole', tag, 0.0))

    def deviation(self, element, limit=None, subset=None):
        """``(worst distance, (tag, point))`` for one element.

        Stops as soon as the worst exceeds ``limit`` (default: the near-miss
        limit) and returns ``inf`` then -- the exact size of a large miss is
        not worth the time.
        """
        limit = self.near if limit is None else limit
        a, b, c, d = element.matrix
        cx, cy = self.centre[0], self.centre[1]
        worst, where = 0.0, None
        samples = self.samples if subset is None else self.samples[::subset]
        for x, y, kind, tag, own in samples:
            dx, dy = x - cx, y - cy
            qx, qy = a * dx + b * dy + cx, c * dx + d * dy + cy
            dist, _ = self.hashes[kind].nearest(qx, qy, limit + own)
            dist = max(0.0, dist - own)
            if dist > worst:
                worst, where = dist, (tag, [x, y, 0.0])
                if worst > limit:
                    return float('inf'), where
        return worst, where


def find_symmetry(outer=None, inners=None, guides=None, poles=None, tol=None,
                  include=INCLUDE, max_order=12, near=None, centre=None, domain=None):
    """Detect the rotations and mirrors that map the input onto itself.

    Parameters
    ----------
    outer, inners, guides, poles
        The domain. Or pass ``domain`` (a :class:`Domain`) instead.
    tol : float, optional
        Drawing precision: the largest distance at which an image still counts
        as landing on the input. Default ``1e-4`` times the bounding-box
        diagonal.
    include : sequence of {'walls', 'holes', 'guides', 'poles'}
        What must be symmetric. Leave out ``'guides'`` to get the symmetry of
        the plate regardless of an off-centre cable.
    max_order : int
        Highest rotation order tested.
    near : float, optional
        Deviations up to this are reported as near-misses. Default ``1e-2``
        times the diagonal.
    centre : [x, y, z], optional
        Override the area centroid.

    Returns
    -------
    :class:`SymmetryReport`
    """
    if domain is None:
        domain = Domain(outer, inners, guides, poles)
    include = tuple(k for k in include if k in INCLUDE)
    diagonal = domain.diagonal or 1.0
    tol = 1e-4 * diagonal if tol is None else float(tol)
    near = max(1e-2 * diagonal, 10 * tol) if near is None else float(near)

    if centre is None:
        centre = domain.area_centroid() if 'walls' in include else None
    if centre is None:
        pts = [p for p in _included_points(domain, include)]
        if not pts:
            raise ValueError('nothing to detect a symmetry from -- the domain is empty')
        centre = [sum(p[0] for p in pts) / len(pts), sum(p[1] for p in pts) / len(pts), 0.0]
    centre = [float(centre[0]), float(centre[1]), 0.0]

    matcher = Matcher(domain, centre, include, near)
    if not matcher.samples:
        raise ValueError('nothing to detect a symmetry from -- no included geometry')
    subset = max(1, len(matcher.samples) // 64)
    tested = {}

    def test(element):
        if element.key in tested:
            return tested[element.key]
        quick, where = matcher.deviation(element, limit=near, subset=subset)
        if quick > tol:
            result = (quick, where)
        else:
            result = matcher.deviation(element, limit=near)
        tested[element.key] = result
        return result

    # --- rotations ---------------------------------------------------------
    passing = set()
    for n in range(2, max_order + 1):
        fwd = test(_rotation(2 * pi / n))
        back = test(_rotation(-2 * pi / n))
        if max(fwd[0], back[0]) <= tol:
            passing.add(n)
    order = 1
    for n in sorted(passing):
        if all((n % m != 0) or (m in passing) for m in range(2, n)):
            order = n
    circle_like = len(passing) == max_order - 1

    # --- mirrors -----------------------------------------------------------
    period = pi / order
    best_axis, best_dev = None, None
    for axis in _mirror_candidates(domain, include, centre, tol, period):
        dev, _ = test(_mirror(axis))
        if dev <= tol and (best_dev is None or dev < best_dev):
            best_axis, best_dev = axis, dev

    group = SymmetryGroup(centre, order, best_axis is not None, best_axis or 0.0)

    # --- the largest subgroup whose every element passes ------------------
    for candidate in group.subgroups():
        if all(test(e)[0] <= tol for e in candidate.elements if e.kind != 'identity'):
            group = candidate
            break

    deviations = {}
    for e in group.elements:
        if e.kind == 'identity':
            continue
        dev, where = test(e)
        deviations[e.key] = (dev, where)

    near_misses = []
    for key, (dev, where) in sorted(tested.items(), key=lambda item: item[1][0]):
        if key in deviations or dev <= tol or dev == float('inf'):
            continue
        near_misses.append((key, dev, where))

    return SymmetryReport(domain=domain, group=group, tol=tol, near=near,
                          include=include, deviations=deviations,
                          near_misses=near_misses, circle_like=circle_like,
                          max_order=max_order)


def _included_points(domain, include):
    if 'walls' in include:
        for p in domain.outer:
            yield p
    if 'holes' in include:
        for loop in domain.inners:
            for p in loop:
                yield p
    if 'guides' in include:
        for curve in domain.guides:
            for p in curve:
                yield p
    if 'poles' in include:
        for p in domain.poles:
            yield p


def _mirror_candidates(domain, include, centre, tol, period):
    """Axis angles worth testing, reduced modulo ``period`` and de-duplicated."""
    cx, cy = centre[0], centre[1]
    angles = []

    def extremes(loop, closed):
        n = len(loop)
        if n == 0:
            return
        r = [hypot(p[0] - cx, p[1] - cy) for p in loop]
        for i in range(n):
            if closed:
                prev, nxt = r[i - 1], r[(i + 1) % n]
            else:
                prev = r[i - 1] if i > 0 else r[i]
                nxt = r[i + 1] if i < n - 1 else r[i]
            is_max = r[i] >= prev - tol and r[i] >= nxt - tol
            is_min = r[i] <= prev + tol and r[i] <= nxt + tol
            if (is_max or is_min) and r[i] > tol:
                angles.append(atan2(loop[i][1] - cy, loop[i][0] - cx))

    if 'walls' in include and len(domain.outer) >= 3:
        extremes(domain.outer, True)
    else:
        for loop in domain.inners if 'holes' in include else []:
            extremes(loop, True)
        for curve in domain.guides if 'guides' in include else []:
            extremes(curve, False)
        for p in domain.poles if 'poles' in include else []:
            if hypot(p[0] - cx, p[1] - cy) > tol:
                angles.append(atan2(p[1] - cy, p[0] - cx))
    if not angles:
        return []

    ordered = sorted(a % (2 * pi) for a in angles)
    bisectors = [0.5 * (ordered[i] + ordered[i + 1]) for i in range(len(ordered) - 1)]
    bisectors.append(0.5 * (ordered[-1] + ordered[0] + 2 * pi))
    reduced = sorted((a % pi) % period for a in ordered + bisectors)

    radius = max(hypot(p[0] - cx, p[1] - cy) for p in _included_points(domain, include)) or 1.0
    atol = max(1e-9, 0.5 * tol / radius)
    out = []
    for a in reduced:
        if out and a - out[-1] <= atol:
            continue
        out.append(a)
    if len(out) > 1 and (out[0] + period) - out[-1] <= atol:
        out.pop()
    return out
