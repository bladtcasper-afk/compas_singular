"""Step 5 -- separatrix tracing.

Separatrices are the field lines that leave a singularity. They partition the
domain into four-sided patches, which is what lets this front end skip the
mixed-integer parametrization entirely (see the evaluation doc): the patches go
straight to ``CoarsePseudoQuadMesh.from_polylines``.

This is the module the design notes flag as highest-risk, and the risk is not in
the integration -- it is that a field line has no obligation to terminate. Three
failure modes, all guarded here rather than left for later:

* LIMIT CYCLES -- a separatrix that spirals forever without reaching a boundary
  or a singularity. Caught by an arc-length cap plus self-proximity detection.
* NEAR-MISSES -- two separatrices passing within epsilon, which produces a sliver
  patch instead of a clean junction. Caught by proximity snapping in ``repair``.
* GRAZING BOUNDARY HITS -- a separatrix approaching the wall tangentially, where
  the exit point is numerically ill-conditioned. Caught by clipping against the
  boundary segment rather than extrapolating the last step.
"""
from collections import namedtuple
from math import atan2
from math import cos
from math import pi
from math import sin

from compas.geometry import distance_point_point
from compas.geometry import intersection_segment_segment_xy
from compas.itertools import pairwise

from .constraints import PERIOD
from .field import wrap_to_period


__all__ = ['Tracer', 'Separatrix']


#: ``reason`` is one of 'boundary', 'singularity', 'cycle', 'length', 'lost'.
Separatrix = namedtuple('Separatrix', 'points source sink reason')

#: Two cut sites whose field alignment differs by less than this count as
#: equally aligned, and ``_cut_launches_from`` picks between them on arc length
#: instead. Sized to sit above discretisation noise and below any real
#: preference: on a round hole every site is normal to the field to within about
#: 1e-3, while on an irregular hole the alignment that matters differs by tenths.
ALIGN_TIE = 0.02


def _bary(p, a, b, c):
    """Barycentric coordinates of ``p`` in triangle ``abc``, in XY."""
    v0x, v0y = b[0] - a[0], b[1] - a[1]
    v1x, v1y = c[0] - a[0], c[1] - a[1]
    v2x, v2y = p[0] - a[0], p[1] - a[1]
    den = v0x * v1y - v1x * v0y
    if abs(den) < 1e-18:
        return None
    v = (v2x * v1y - v1x * v2y) / den
    w = (v0x * v2y - v2x * v0y) / den
    return (1.0 - v - w, v, w)


class Tracer(object):
    """Traces the separatrices of a :class:`CrossField`.

    Parameters
    ----------
    field : CrossField
    step : float, optional
        Integration step. Defaults to a third of the background target length --
        small enough that a step never crosses more than one triangle, which is
        what keeps the branch matching stable.
    max_length : float, optional
        Arc-length cap. Defaults to 40x the domain diagonal.
    """

    def __init__(self, field, step=None, max_length=None):
        self.field = field
        self.background = field.background
        self.mesh = field.background.mesh
        self.step = step or self.background.target_length / 3.0

        xs = [p[0] for p in self.background.outer]
        ys = [p[1] for p in self.background.outer]
        self.diagonal = ((max(xs) - min(xs)) ** 2 + (max(ys) - min(ys)) ** 2) ** 0.5
        self.max_length = max_length or 40.0 * self.diagonal

        self._grid = None
        self._cell = self.background.target_length * 1.5
        self._origin = (min(xs), min(ys))

    # ------------------------------------------------------------------
    # point location
    # ------------------------------------------------------------------

    def _build_grid(self):
        grid = {}
        for fkey in self.mesh.faces():
            pts = [self.mesh.vertex_coordinates(v) for v in self.mesh.face_vertices(fkey)]
            i0 = int((min(p[0] for p in pts) - self._origin[0]) // self._cell)
            i1 = int((max(p[0] for p in pts) - self._origin[0]) // self._cell)
            j0 = int((min(p[1] for p in pts) - self._origin[1]) // self._cell)
            j1 = int((max(p[1] for p in pts) - self._origin[1]) // self._cell)
            for i in range(i0, i1 + 1):
                for j in range(j0, j1 + 1):
                    grid.setdefault((i, j), []).append(fkey)
        self._grid = grid

    def locate(self, point, hint=None):
        """Containing face and barycentric coordinates, or ``None`` if outside.

        ``hint`` is checked first -- consecutive integration steps almost always
        land in the same or an adjacent triangle, so this turns point location
        from the tracer's bottleneck into a couple of dot products.
        """
        if hint is not None:
            for fkey in [hint] + list(self.mesh.face_neighbors(hint)):
                found = self._test(fkey, point)
                if found:
                    return found

        if self._grid is None:
            self._build_grid()
        i = int((point[0] - self._origin[0]) // self._cell)
        j = int((point[1] - self._origin[1]) // self._cell)
        for fkey in self._grid.get((i, j), ()):
            found = self._test(fkey, point)
            if found:
                return found
        return None

    def _test(self, fkey, point, tol=-1e-9):
        vkeys = self.mesh.face_vertices(fkey)
        a, b, c = [self.mesh.vertex_coordinates(v) for v in vkeys]
        bary = _bary(point, a, b, c)
        if bary is None:
            return None
        if bary[0] >= tol and bary[1] >= tol and bary[2] >= tol:
            return fkey, bary
        return None

    # ------------------------------------------------------------------
    # field sampling
    # ------------------------------------------------------------------

    def direction_at(self, point, reference, hint=None):
        """Unit field direction at ``point``, on the arm closest to ``reference``.

        Returns ``(direction, fkey)`` or ``(None, None)`` outside the domain.

        Picking the arm by proximity to where the tracer came from is the whole
        trick. A cross has four arms and no intrinsic orientation; matching to the
        incoming direction is what makes a streamline continue straight instead of
        turning 90 degrees at an arbitrary triangle.
        """
        found = self.locate(point, hint)
        if found is None:
            return None, None
        fkey, bary = found
        ref_angle = atan2(reference[1], reference[0])
        theta = self.field.angle_in_face(fkey, bary, reference=ref_angle)
        theta += PERIOD * round((ref_angle - theta) / PERIOD)
        return [cos(theta), sin(theta), 0.0], fkey

    # ------------------------------------------------------------------
    # launch directions
    # ------------------------------------------------------------------

    def _singularity_points(self):
        return {fkey: self.mesh.face_centroid(fkey)
                for fkey, _ in self.field.singularities()}

    def launch_directions(self, centre, index):
        """Directions in which separatrices leave a singularity.

        Walk a small circle around the singularity; a separatrix leaves wherever
        the field points RADIALLY, i.e. wherever the residual ``phi - theta``
        hits a multiple of the cross's period.

        The residual is UNWRAPPED along the circle and the crossings counted as
        level crossings, rather than hunting sign changes in the wrapped residual.
        That distinction matters: the wrapped residual is a sawtooth, so a
        sign-change search cannot tell a genuine root from the -45/+45 jump, and
        any noise wiggle near a root registers as an extra pair. Measured on a
        disc, the sign-change version produced two launches 1.5 degrees apart and
        four arms where a +1 singularity has three.

        Unwrapped, the residual is monotone by construction: over one lap it
        advances by exactly ``(4 - index)`` periods, so counting level crossings
        gives the right number of arms and cannot double-count.
        """
        target = 4 - index

        # Several radii AND several sample counts. A single probe circle is not
        # reliable: too small and it sits in the handful of triangles where
        # |u| -> 0 and the angle is noise; too large and it leaves the domain or
        # swallows a neighbouring singularity. Under-sampling aliases -- if the
        # field turns more than half a period between adjacent samples the
        # unwrapping loses a period and an arm goes missing.
        #
        # Measured before this: the same pentagon gave 5 arms at target_length
        # 0.50 and 0.35 but only 4 at 0.60, 0.45 and 0.40, with no indication in
        # the output that anything was wrong.
        best = []
        attempts = [(s, n) for n in (360, 720, 1440)
                    for s in (1.2, 0.9, 1.6, 0.6, 2.0, 0.45)]
        for scale, samples in attempts:
            radius = self.background.target_length * scale
            phis, residual = [], []
            theta = None
            complete = True
            for i in range(samples + 1):                  # +1 closes the lap at 2*pi
                phi = 2.0 * pi * i / samples
                p = [centre[0] + cos(phi) * radius, centre[1] + sin(phi) * radius, 0.0]
                found = self.locate(p)
                if found is None:
                    complete = False
                    break
                fkey, bary = found
                raw = self.field.angle_in_face(fkey, bary, reference=phi)
                theta = raw if theta is None else theta + wrap_to_period(raw - theta)
                phis.append(phi)
                residual.append(phi - theta)
            if not complete:
                continue                    # circle left the domain; try another radius

            found_dirs = self._level_crossings(phis, residual, target)
            if len(found_dirs) == target:
                return [[cos(a), sin(a), 0.0] for a in found_dirs]
            if len(found_dirs) > len(best):
                best = found_dirs

        return [[cos(a), sin(a), 0.0] for a in best]

    @staticmethod
    def _level_crossings(phis, residual, target):
        """Where the unwrapped residual first reaches each multiple of the period.

        Over one lap the residual advances by exactly ``target`` periods, so there
        are exactly ``target`` separatrices -- the count is known in advance and
        does not need to be discovered. What has to be found is only WHERE.

        The search runs against the RUNNING MAXIMUM of the residual rather than
        the residual itself. Close to a singularity ``|u| -> 0``, so
        ``theta = arg(u)/4`` is ill-conditioned and the residual acquires small
        oscillations; against the raw signal every wiggle across a level counts as
        another crossing. Measured on a pentagon, that turned 5 arms into 11, some
        of them 0.2 degrees apart. The running maximum is monotone by
        construction, so each level is consumed exactly once.
        """
        out = []
        start = residual[0]
        levels = []
        k = int(start // PERIOD) + 1
        while len(levels) < target:
            levels.append(k * PERIOD)
            k += 1

        running = start
        li = 0
        for i in range(1, len(residual)):
            previous = running
            running = max(running, residual[i])
            while li < len(levels) and levels[li] <= running:
                level = levels[li]
                # interpolate within the segment that actually crossed it
                span = running - previous
                t = 0.0 if span <= 0 else max(0.0, min(1.0, (level - previous) / span))
                out.append(phis[i - 1] + t * (phis[i] - phis[i - 1]))
                li += 1
        return out

    # ------------------------------------------------------------------
    # integration
    # ------------------------------------------------------------------

    def _inward_normal(self, loop, i):
        """Unit normal at ``loop[i]`` pointing INTO the domain, or ``None``.

        Taken from the edge tangent rather than from the loop's centroid, so it
        is still correct on a hole that is not convex. The sign is decided by
        probing: whichever side locates inside a background face is the domain.
        """
        n = len(loop)
        a, b = loop[i - 1], loop[(i + 1) % n]
        tx, ty = b[0] - a[0], b[1] - a[1]
        length = (tx * tx + ty * ty) ** 0.5
        if length < 1e-12:
            return None
        for sign in (1.0, -1.0):
            normal = [-ty / length * sign, tx / length * sign, 0.0]
            probe = [loop[i][0] + normal[0] * self.background.target_length * 0.6,
                     loop[i][1] + normal[1] * self.background.target_length * 0.6, 0.0]
            if self.locate(probe) is not None:
                return normal
        return None

    def _cut_launches_from(self, loop, count):
        """``count`` launch sites spread around ``loop``, each on a field line.

        Two independent choices, and both matter:

        **Spacing is by ARC LENGTH, not vertex index.** Indexing assumes uniform
        sampling; on any loop that is not uniformly sampled it bunches the cuts
        together and they fail to separate the domain. This is the same mistake
        that made ``_split_loop`` emit overlapping boundary arcs.

        **Within each slot the site is chosen by FIELD ALIGNMENT** -- the vertex
        whose best cross arm is most nearly normal to the loop. A cut launched
        across the field turns immediately and stops being the straight radial
        line the layout needs; launched along it, the cut IS a field line. On a
        round hole the field is polar, so every site scores alike and the arc
        spacing decides; on an irregular hole the score does the work.
        """
        from .repair import _arc_lengths

        _seg, cum = _arc_lengths(loop)
        total = cum[-1]
        if total <= 0.0:
            return []

        step = self.background.target_length * 0.8
        window = total / (2.0 * count)
        launches = []
        for k in range(count):
            target = total * k / count
            candidates = []
            for i in range(len(loop)):
                # circular distance from this vertex to the slot centre
                offset = abs(cum[i] - target)
                offset = min(offset, total - offset)
                if offset > window:
                    continue
                normal = self._inward_normal(loop, i)
                if normal is None:
                    continue
                probe = [loop[i][0] + normal[0] * step,
                         loop[i][1] + normal[1] * step, 0.0]
                found = self.locate(probe)
                if found is None:
                    continue
                fkey, bary = found
                theta = self.field.angle_in_face(fkey, bary)
                arm, score = None, -2.0
                for a in range(4):
                    angle = theta + a * PERIOD
                    d = [cos(angle), sin(angle), 0.0]
                    dot = d[0] * normal[0] + d[1] * normal[1]
                    if dot > score:
                        arm, score = d, dot
                if arm is not None:
                    # Start ON the loop vertex, not at the probe. The probe is
                    # inset only to sample the field somewhere safely inside the
                    # domain; tracing from it displaces the whole cut sideways by
                    # that inset, so its landing snaps to a different point of the
                    # loop than the launch vertex and the loop ends up marked
                    # TWICE per cut -- which turns each quarter-annulus into a
                    # 5-gon and every fan of one into slivers. This is the same
                    # mistake, and the same fix, as corner launches.
                    candidates.append((score, offset, list(loop[i]), arm))
            if not candidates:
                continue
            # Best alignment wins, but only when the difference MEANS something.
            # A plain argmax over ``score`` is an argmax over discretisation
            # noise wherever the field is already normal to the loop everywhere
            # -- which is exactly the round hole this method exists for. It put
            # the annulus's four cuts at 45, 97.5, 225 and 255 degrees: sectors
            # of 30 and 150 degrees side by side, an aspect ratio of 5 baked into
            # the layout before densification ever ran. Within ``ALIGN_TIE`` the
            # sites are equally good field lines, so the arc-length offset
            # decides, which is what the slots were computed for in the first
            # place.
            top = max(c[0] for c in candidates)
            _, _, point, arm = min(
                (c for c in candidates if c[0] > top - ALIGN_TIE),
                key=lambda c: c[1])
            launches.append((point, arm))
        return launches

    def hole_launches(self, corner_limit=pi / 12.0, count=4):
        """Cuts opening a boundary loop that emits no separatrix of its own.

        Separatrices otherwise come from two places only: interior singularities
        and reflex corners. A SMOOTH hole has neither, so nothing cuts the domain
        and an annulus -- which is not simply connected and cannot be one patch
        however much the layout is repaired -- gets no layout at all. That is
        ``13_limits`` LIMIT 4, and it is why a round-holed disc falls to the
        triangulation backstop.

        **Why two and not one.** One cut is a SLIT: both of its sides belong to
        the same face, so the planar face walk runs up it and back down and
        emits a face with a repeated vertex -- exactly the pathology
        ``arrangement.py`` exists to remove. Two cuts leave two genuine regions,
        each simply connected and each bounded by distinct edges.

        **Why this is now worth trying again.** An earlier version of this method
        was written and left disabled: two radial cuts did not give
        ``from_polylines`` a network it could resolve, and enabling it broke
        annulus cases that had worked. That failure was ``from_polylines``
        recovering faces from an embedding whose crossings had no node and whose
        loose ends were never landed -- which is precisely what
        ``planar_arrangement`` was later built to fix. The blocker was removed
        by other work; this method was not retried until now.

        Parameters
        ----------
        corner_limit : float, optional
            A loop with corners above this is left alone -- its corners already
            supply the cuts. A square hole in a square plate gives 8 and
            decomposes cleanly without any of this.
        count : int, optional
            Cuts per loop. Two is the minimum that separates; more only adds
            patches.

        Returns
        -------
        list[(point, direction)]
        """
        from .repair import boundary_corners

        launches = []
        for loop in self.background.inners:
            if boundary_corners(loop, corner_limit):
                continue                    # its corners already supply the cuts
            launches.extend(self._cut_launches_from(loop, count))
        return launches

    def corner_points(self, corner_limit=pi / 12.0):
        """Every boundary corner, as a flat list of points."""
        from .repair import boundary_corners
        out = []
        for loop in [self.background.outer] + list(self.background.inners):
            for i in boundary_corners(loop, corner_limit):
                out.append(list(loop[i]))
        return out

    def trace(self, start, direction, source=None, singular_points=None, corners=None):
        """Integrate one separatrix from ``start`` along ``direction``.

        Returns
        -------
        Separatrix
        """
        singular_points = singular_points or {}
        corners = corners or []
        stop_radius = self.background.target_length * 1.2
        h = self.step

        points = [list(start)]
        d = list(direction)
        fkey = None
        length = 0.0
        reason = 'length'
        sink = None

        while length < self.max_length:
            # RK2 (midpoint): a plain Euler step drifts off a curving field line
            # badly enough to matter over the long traces this produces
            d1, f1 = self.direction_at(points[-1], d, fkey)
            if d1 is None:
                reason = 'lost'
                break
            mid = [points[-1][k] + d1[k] * h * 0.5 for k in range(3)]
            d2, _ = self.direction_at(mid, d1, f1)
            if d2 is None:
                d2 = d1
            nxt = [points[-1][k] + d2[k] * h for k in range(3)]

            found = self.locate(nxt, f1)
            if found is None:
                clipped = self._clip_to_boundary(points[-1], nxt)
                if clipped is not None:
                    points.append(clipped)
                    length += distance_point_point(points[-2], clipped)
                reason = 'boundary'
                break

            length += distance_point_point(points[-1], nxt)
            points.append(nxt)
            d, fkey = d2, found[0]

            # reached another singularity
            hit = None
            for skey, centre in singular_points.items():
                if skey == source and length < stop_radius * 3.0:
                    continue                       # our own launch point
                if distance_point_point(nxt, centre) < stop_radius:
                    hit = skey
                    break
            if hit is not None:
                points.append(list(singular_points[hit]))
                reason, sink = 'singularity', hit
                break

            # Arrived at the wall without crossing it. Leaving this to the
            # "stepped outside the domain" test alone is not enough: a separatrix
            # can run ALONG the boundary and never leave. On a T-shaped plate the
            # line from one reflex corner reaches the other and then continues
            # down the boundary edge, overlapping it exactly -- and
            # from_polylines' planar face search cannot separate a region from a
            # line lying on its own edge, so the whole plate comes back as one
            # 26-sided face.
            #
            # Proximity to the wall, not proximity to a CORNER: terminating near
            # corners instead stops separatrices that merely pass one, which
            # collapses patches on a pentagon.
            if length > stop_radius * 2.0:
                near = self._project_to_boundary(nxt)
                if near is not None and distance_point_point(nxt, near) < h * 0.75:
                    points.append(near)
                    reason = 'boundary'
                    break

            # limit cycle: back near a point we already passed, long ago
            if len(points) > 40 and self._revisits(points, stop_radius * 0.5):
                reason = 'cycle'
                break

        return Separatrix(points, source, sink, reason)

    def _revisits(self, points, tol):
        """Whether the head has returned to an early part of its own trail."""
        head = points[-1]
        for p in points[:-20]:
            if distance_point_point(head, p) < tol:
                return True
        return False

    def _project_to_boundary(self, point):
        """Closest point on any boundary loop."""
        best, best_d = None, float('inf')
        for loop in [self.background.outer] + list(self.background.inners):
            for a, b in pairwise(loop + loop[:1]):
                abx, aby = b[0] - a[0], b[1] - a[1]
                length2 = abx * abx + aby * aby
                if length2 == 0.0:
                    continue
                t = ((point[0] - a[0]) * abx + (point[1] - a[1]) * aby) / length2
                t = max(0.0, min(1.0, t))
                q = [a[0] + abx * t, a[1] + aby * t, 0.0]
                d = distance_point_point(point, q)
                if d < best_d:
                    best, best_d = q, d
        return best

    def _clip_to_boundary(self, inside, outside):
        """Where the segment leaves the domain.

        Clipping against the actual boundary segment rather than extrapolating
        the last integration step is what keeps a grazing exit well-conditioned:
        near-tangential approach makes the step direction useless but leaves the
        intersection itself perfectly determined.
        """
        best, best_d = None, float('inf')
        for loop in [self.background.outer] + list(self.background.inners):
            for a, b in pairwise(loop + loop[:1]):
                x = intersection_segment_segment_xy((inside, outside), (a, b))
                if x is None:
                    continue
                x = [x[0], x[1], 0.0]
                d = distance_point_point(inside, x)
                if d < best_d:
                    best, best_d = x, d
        return best

    # ------------------------------------------------------------------
    # the whole set
    # ------------------------------------------------------------------

    def corner_launches(self, corner_limit=pi / 12.0):
        """Separatrices emitted by the domain's own corners.

        Interior singularities are not the only sources. A boundary corner with
        interior angle ``alpha`` is a BOUNDARY singularity: the quad layout puts
        ``k = round(alpha / 90)`` patch corners there, and it emits ``k - 1``
        separatrices into the domain.

        * 90 degrees (k=1) -- a plain convex corner, emits nothing.
        * 180 degrees (k=2) -- a straight wall, not a corner at all, never reaches
          here because the turn is below ``corner_limit``.
        * 270 degrees (k=3) -- a reflex corner, emits two.

        Without this an L-shape has no singularities, hence no separatrices, hence
        no layout -- yet it obviously decomposes into three rectangles, and those
        cuts start exactly at the reflex corner. Tracing only from interior
        singularities misses every domain whose structure lives on the wall.

        Returns
        -------
        list[(point, direction)]
        """
        from .repair import boundary_corners

        launches = []
        for loop_index, loop in enumerate([self.background.outer] + list(self.background.inners)):
            n = len(loop)
            is_hole = loop_index > 0
            for i in boundary_corners(loop, corner_limit):
                here = loop[i]
                prev, nxt = loop[(i - 1) % n], loop[(i + 1) % n]
                back = [prev[0] - here[0], prev[1] - here[1], 0.0]
                fwd = [nxt[0] - here[0], nxt[1] - here[1], 0.0]
                interior_angle, ccw = self._interior_angle(back, fwd, loop)
                if is_hole:
                    # a hole's own interior is OUTSIDE the domain, so the angle the
                    # loop reports is the complement of the one that matters. A
                    # square hole's corners are 90 degrees from the hole's side and
                    # 270 from the domain's -- read the wrong one and every hole
                    # corner looks convex and emits nothing.
                    interior_angle = 2.0 * pi - interior_angle
                    ccw = not ccw
                k = int(round(interior_angle / PERIOD))
                if k < 2:
                    continue

                # The interior direction is the outgoing wall rotated by HALF the
                # interior angle, not the bisector of the two wall vectors. At a
                # reflex corner those are opposite: back + fwd points into the
                # 90-degree wedge the domain does NOT occupy, so a probe placed
                # there lands outside, the field sample fails, and the corner
                # silently emits nothing.
                a_fwd = atan2(fwd[1], fwd[0])
                turn = 1.0 if ccw else -1.0
                a_in = a_fwd + turn * interior_angle * 0.5
                inset = self.background.target_length * 0.8
                probe = [here[0] + cos(a_in) * inset, here[1] + sin(a_in) * inset, 0.0]
                found = self.locate(probe)
                if found is None:
                    continue
                fkey, bary = found
                theta = self.field.angle_in_face(fkey, bary)

                # keep the arms lying strictly inside the interior sector --
                # a dot product against the interior direction is too blunt for a
                # 270-degree corner, where a valid arm can be 135 degrees off it
                margin = 0.15
                candidates = []
                for arm in range(4):
                    a = theta + arm * PERIOD
                    rel = turn * (a - a_fwd) % (2.0 * pi)
                    if margin < rel < interior_angle - margin:
                        candidates.append(([cos(a), sin(a), 0.0], rel))
                candidates.sort(key=lambda c: c[1])

                # Launch from the CORNER, stepped along the launch direction --
                # not from the probe. The probe is inset along the interior
                # bisector and exists only to sample the field somewhere safely
                # inside; starting the trace there displaces the whole separatrix
                # sideways by that inset. On a T-shaped plate the line leaving one
                # reflex corner then runs 0.28 BELOW the wall it should lie on,
                # passes under the opposite corner instead of terminating at it,
                # and carries on to the far edge -- which is what turned the plate
                # into a single 18-sided face.
                step = self.background.target_length * 0.5
                launches.extend(
                    ([here[0] + d[0] * step, here[1] + d[1] * step, 0.0], d)
                    for d, _ in candidates[:k - 1])
        return launches

    def _interior_angle(self, back, fwd, loop):
        """Interior angle at a corner, measured on the domain side.

        On a counter-clockwise loop the interior is to the left, so the interior
        angle is the sweep from the OUTGOING edge round to the INCOMING one --
        that order, not the reverse. Getting it backwards returns the exterior
        angle, which reads a 270-degree reflex corner as a 90-degree convex one
        and emits no separatrices from it at all: an L-shape then has no interior
        lines and no layout.
        """
        a_back = atan2(back[1], back[0])
        a_fwd = atan2(fwd[1], fwd[0])
        pts = loop
        area2 = sum(p[0] * q[1] - q[0] * p[1] for p, q in zip(pts, pts[1:] + pts[:1]))
        ccw = area2 > 0
        sweep = (a_back - a_fwd) % (2.0 * pi)
        return (sweep if ccw else 2.0 * pi - sweep), ccw

    def separatrices(self):
        """Every separatrix of the field -- from singularities and from corners.

        Returns
        -------
        (list[Separatrix], dict)
            The traces, and a report with the per-reason counts and the
            arm-count check.
        """
        singular_points = self._singularity_points()
        indices = dict(self.field.singularities())
        corners = self.corner_points()

        out = []
        arm_mismatch = []
        for fkey, centre in singular_points.items():
            dirs = self.launch_directions(centre, indices[fkey])
            expected = 4 - indices[fkey]
            if len(dirs) != expected:
                arm_mismatch.append((fkey, len(dirs), expected))
            for d in dirs:
                start = [centre[0] + d[0] * self.background.target_length * 0.5,
                         centre[1] + d[1] * self.background.target_length * 0.5, 0.0]
                out.append(self.trace(start, d, source=fkey,
                                      singular_points=singular_points, corners=corners))

        for start, d in self.corner_launches():
            out.append(self.trace(start, d, source=None,
                                  singular_points=singular_points, corners=corners))

        # A loop with neither a corner nor a singularity to cut it emits nothing
        # above, and a domain it bounds cannot be decomposed at all. Only reached
        # for smooth holes; a hole with corners is skipped inside.
        cuts = 0
        for start, d in self.hole_launches():
            # ``source='hole'`` rather than None. It is never a face key, so the
            # self-hit guard in ``trace`` ignores it, but it distinguishes a cut
            # from a CORNER launch -- which is the other thing that carries
            # ``source=None`` -- and ``build_network`` has to tell them apart:
            # a corner launch starts one inset inside a corner it should be
            # snapped back onto, a cut starts on a smooth wall and must be left
            # exactly where it is AND recorded as a split point of that wall.
            trace = self.trace(start, d, source='hole',
                               singular_points=singular_points, corners=corners)
            # A cut that does not reach another boundary leaves the loop open and
            # is worse than none: it is a dangling edge, and the face walk runs up
            # it and back. Keep only the ones that landed.
            if trace.reason == 'boundary':
                out.append(trace)
                cuts += 1

        reasons = {}
        for s in out:
            reasons[s.reason] = reasons.get(s.reason, 0) + 1

        report = {
            'count': len(out),
            'reasons': reasons,
            'arm_mismatch': arm_mismatch,
            'hole_cuts': cuts,
            'ok': not arm_mismatch and not reasons.get('lost') and not reasons.get('cycle'),
        }
        return out, report
