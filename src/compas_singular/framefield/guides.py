"""Directional guide metrics: does the mesh run ALONG the guide, and mid-block?

``quality.curve_alignment`` cannot answer either question, deliberately. It folds
the angle into the cross's period, so *"an edge crossing the curve at 90 degrees
counts as aligned as one running along it, because a quad mesh has two families and
a cable may be either."* That is the right convention for asking whether the FIELD
is aligned -- it is what makes the mesh number and the field number comparable, which
is what the whole field -> mesh objective is stated on.

It is the wrong convention for asking whether the guide got what it was for. A cable
is a cable and the quads are blocks, so the cable must run down the **middle of a row
of faces**; a cable sitting on the shared edges between two rows is a cable in the
joint between two blocks. Those two meshes score **identically** under a
cross-symmetric metric. So does a mesh whose courses run across the cable instead of
along it.

WHY THIS IS NOT ``guide_lines.py``'s ``measure``
------------------------------------------------

That one (in the ``__main__`` block of ``src/compas_singular/guide_lines.py``) selects
its edges by INCIDENCE: vertices within ``on_tol=0.05`` of the guide, i.e. essentially
exactly on it. ``guide_lines`` may assume that, because it *creates* the incidence by
snapping a polyedge onto the guide. **The frame-field route never snaps.** The guide is
a field constraint; no vertex is placed on it on purpose. At ``target_length=1.0`` on a
10-unit domain, ``on_tol`` is half a percent of an edge, so the incident set is expected
to be empty and ``along``/``across`` undefined. That is very likely why
``curve_alignment`` was written radius-based in the first place.

Loosening ``on_tol`` into a radius is the wrong repair: a radius around a curve collects
a ragged band of vertices from BOTH sides, which is the selection that folded faces when
it was tried as a smoothing constraint.

So select by FACE instead, the way ``guide_lines``' ``off_centre`` already does -- it
has no incidence assumption and transfers unchanged. A face counts if the guide passes
through it. Then:

* the face's two edge families give **along** (the family nearer the tangent, want 0
  degrees) and **across** (the other one, want 90);
* the guide's offset from the face centre, in units of half the face width across the
  guide, gives **off-centre** -- 0 is dead centre (mid-block), 1 is on the face edge
  (in the joint).

All three from one pass, no incidence, no chord tracking, no radius averaging.

READING THE NUMBERS
-------------------

**Always report the no-guide control.** Where a guide runs parallel to what the mesh
would have done anyway, ``across`` scores ~90 whether or not the guide did anything, so
an impressive number on a flat guide means nothing on its own. ``across_tilted`` bins to
the stations where the guide is more than ``steep`` degrees off horizontal, which on an
axis-aligned domain is where a course actually has to turn. Horizontal is a proxy, and
only a meaningful one when the domain's natural grid is axis-aligned -- the control is
the real defence, not the binning.

**No boundary margin.** Excluding a fan of vertices near the boundary looks reasonable
and, on a bowed guide, deletes exactly the tilted part -- the only part the metric is
about. Measured on one such mesh: margin 0 -> 76.4 deg / 62% within 15; 1.5 -> 78.7 / 79;
2.0 -> 80.6 / 88. None of those three numbers is about the boundary.

**Poles are skipped, not reconstituted.** A pseudo-quad has three corners and no two
edge families, so there is nothing to call along or across. They are counted separately
rather than silently dropped, for the same reason ``quality.py`` never reconstitutes the
phantom fourth corner.
"""
from math import acos
from math import degrees

from compas_singular.geometry.polyline import closest_on_polyline


__all__ = ['guide_metrics', 'format_guide_metrics']

#: A guide station is "tilted" past this many degrees off horizontal. See the
#: module docstring on why this is a proxy and the control is the real defence.
STEEP = 12.0

#: ``across`` counts as on target within this many degrees of 90.
ACROSS_TOLERANCE = 15.0


def _median(values):
    s = sorted(values)
    k = len(s)
    if not k:
        return None
    return s[k // 2] if k % 2 else 0.5 * (s[k // 2 - 1] + s[k // 2])


def _acute(u, v):
    """Angle between two directions ignoring sense, in degrees: 0 to 90."""
    d = abs(u[0] * v[0] + u[1] * v[1])
    return degrees(acos(max(-1.0, min(1.0, d))))


def _nearest_on_polyline(point, polyline):
    """``(distance, unit tangent)`` of the closest point of a polyline, or ``(inf, None)``.

    The tangent is the direction of the segment the point landed on, which is
    why this needs the index :func:`closest_on_polyline` returns and not just
    the closest point.
    """
    index, _t, _q, distance = closest_on_polyline(point, polyline)
    if distance == float('inf'):
        return distance, None
    a, b = polyline[index], polyline[index + 1]
    abx, aby = b[0] - a[0], b[1] - a[1]
    length = (abx * abx + aby * aby) ** 0.5
    return distance, (abx / length, aby / length)


def _families(points):
    """Mean direction of a quad's two opposite edge families, or ``None``.

    ``(v0 v1, v2 v3)`` is one family and ``(v1 v2, v3 v0)`` the other. Directions are
    averaged without sense -- the two edges of a family run opposite ways round the
    face, so the second is flipped onto the first before averaging.
    """
    if len(points) != 4:
        return None

    out = []
    for start in (0, 1):
        vectors = []
        for k in (start, start + 2):
            a, b = points[k], points[(k + 1) % 4]
            dx, dy = b[0] - a[0], b[1] - a[1]
            length = (dx * dx + dy * dy) ** 0.5
            if length <= 1e-12:
                continue
            vectors.append((dx / length, dy / length))
        if not vectors:
            return None
        ux, uy = vectors[0]
        for vx, vy in vectors[1:]:
            if ux * vx + uy * vy < 0.0:      # opposite sense round the face
                vx, vy = -vx, -vy
            ux, uy = ux + vx, uy + vy
        length = (ux * ux + uy * uy) ** 0.5
        if length <= 1e-12:
            return None
        out.append((ux / length, uy / length))
    return out


def guide_metrics(mesh, guides, steep=STEEP):
    """**How well a mesh carries its guides**, as numbers rather than a verdict.

    Parameters
    ----------
    mesh : Mesh
        A quad or pseudo-quad mesh.
    guides : sequence[sequence[point]]
        The guide curves, as polylines -- the same objects fed to the field.
    steep : float, optional
        Degrees off horizontal past which a guide station counts as tilted.

    Returns
    -------
    dict
        ``faces`` -- faces the guides pass through, the population everything else is
        measured over. **Zero means the metric could not see the guide at all**, and
        every other value is None.

        ``along`` -- median angle between the guide and the face family nearer to it,
        in degrees. 0 is a course running along the guide.

        ``across`` / ``across_within`` / ``across_tilted`` -- median angle of the other
        family, the share of it within :data:`ACROSS_TOLERANCE` of 90 degrees, and the
        median over tilted stations only. 90 is an interface meeting the guide square.

        ``off_centre`` / ``off_centre_worst`` -- the guide's offset from the face centre
        in units of half the face width across it. 0 is mid-block, 1 is on the joint.

        ``poles`` -- pseudo-quad faces skipped, having no two families to compare.
    """
    empty = {'faces': 0, 'along': None, 'across': None, 'across_within': None,
             'across_tilted': None, 'off_centre': None, 'off_centre_worst': None,
             'poles': 0}
    if not guides:
        return empty

    along, across, ratios, poles = [], [], [], 0

    for guide in guides:
        points = [p for p in guide]
        if len(points) < 2:
            continue

        for face in mesh.faces():
            corners = [mesh.vertex_coordinates(v) for v in mesh.face_vertices(face)]
            if len(corners) < 3:
                continue
            centre = [sum(p[i] for p in corners) / len(corners) for i in range(3)]

            distance, tangent = _nearest_on_polyline(centre, points)
            if tangent is None:
                continue
            # only the faces the guide actually passes through
            if distance >= max(((centre[0] - p[0]) ** 2 + (centre[1] - p[1]) ** 2) ** 0.5
                               for p in corners):
                continue

            # off-centre: the guide's offset in units of half the face width across it
            normal = (-tangent[1], tangent[0])
            projections = [(p[0] - centre[0]) * normal[0] + (p[1] - centre[1]) * normal[1]
                           for p in corners]
            half = (max(projections) - min(projections)) / 2.0
            if half > 1e-9:
                ratios.append(distance / half)

            families = _families(corners)
            if families is None:
                poles += 1
                continue

            angles = sorted(_acute(f, tangent) for f in families)
            along.append(angles[0])
            tilt = _acute(tangent, (1.0, 0.0))
            across.append((tilt, angles[1]))

    if not along and not ratios:
        return dict(empty, poles=poles)

    angles = [a for _, a in across]
    tilted = [a for t, a in across if t >= steep]
    return {
        'faces': len(along),
        'along': _median(along),
        'across': _median(angles),
        'across_within': (100.0 * sum(1 for a in angles if a >= 90.0 - ACROSS_TOLERANCE)
                          / len(angles)) if angles else None,
        'across_tilted': _median(tilted),
        'off_centre': _median(ratios),
        'off_centre_worst': max(ratios) if ratios else None,
        'poles': poles,
    }


def format_guide_metrics(metrics):
    """One line of the numbers from :func:`guide_metrics`, for a harness to print."""
    if not metrics['faces']:
        return 'guide unseen (0 faces{})'.format(
            ', {} poles'.format(metrics['poles']) if metrics['poles'] else '')

    def num(value, fmt='{:5.1f}'):
        return '  -- ' if value is None else fmt.format(value)

    out = 'faces {:3d} | along {} deg | across {} deg, {}% within {:.0f}'.format(
        metrics['faces'], num(metrics['along']), num(metrics['across']),
        num(metrics['across_within'], '{:3.0f}'), ACROSS_TOLERANCE)
    out += ' | tilted {} deg'.format(num(metrics['across_tilted']))
    out += ' | off-centre {} med, {} worst'.format(
        num(metrics['off_centre'], '{:4.2f}'), num(metrics['off_centre_worst'], '{:4.2f}'))
    if metrics['poles']:
        out += ' | {} poles skipped'.format(metrics['poles'])
    return out
