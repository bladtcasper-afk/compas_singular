"""Element quality, in two tiers -- a hard floor and a regression baseline.

Until this module existed there was no automated quality gate at all.
``FieldDecomposition._acceptable`` checked all-quad, manifold, and 97-103% area
coverage, and ``route()`` was the de-facto gate on top of that. Both stopped
discriminating: every test domain now reaches route ``'field'``, and those three
structural properties certified a mesh carrying a literal **180 degree angle**
as fine. That mesh is live in the suite -- ``12_cables``' ring cable, face 128,
a four-sided non-pole quad whose corners are
``(7.2455, 3.0752) (7.5500, 3.2143) (7.5341, 2.8200) (7.4960, 2.8537)``: three
of them collinear. Nothing structural is wrong with it. It is not a usable
element.

TWO TIERS, and one threshold cannot do both jobs
------------------------------------------------

* The **hard floor** (:func:`hard_floor`) is absolute and permanent: an angle at
  or near 0 or 180 degrees, or a non-finite aspect ratio, is a degeneracy rather
  than a taste, and no baseline may bless one.
* The **regression check** lives in ``15_baseline.py``, against committed
  numbers. It is what catches "the change made the ellipse worse", which no
  absolute threshold can.

The reason for the split is measured, not stylistic. Minimum angle across the
field route legitimately spans 24 to 90 degrees over the domain suite, so any
absolute threshold loose enough to pass the ellipse is far too loose to notice
the square dropping off 90.

PSEUDO-QUADS
------------

A pole is the quad ``(p, a, b, p)`` with one side collapsed, stored as the
three-vertex face ``(p, a, b)`` plus an entry in ``mesh.attributes['face_pole']``.
Every metric here measures the corners a face ACTUALLY has and never
reconstitutes the phantom fourth one -- reconstituting it yields a zero-length
edge, hence a 0 degree angle and a division by zero, and would fail every domain
carrying a pole the moment it saw one. Measured on the ring cable's 14 dense
pole faces: read as triangles their angles run 7.4 to 117.5 degrees and their
shortest edge is 0.051, all finite; read as collapsed quads every one of them
reports 0 degrees and an infinite aspect ratio.

A 7.4 degree pole is bad, and it is the regression check's business, not the
floor's.
"""
from math import acos
from math import atan2
from math import degrees
from math import pi


__all__ = ['mesh_quality', 'hard_floor', 'face_angles', 'curve_alignment',
           'curve_alignment_profile',
           'HARD_MIN_ANGLE', 'HARD_MAX_ANGLE', 'LOW_ANGLE']


#: Below this, an angle is a degeneracy rather than a poor element. The worst
#: legitimate value in the suite is the ring cable's 7.4 degree pole, so this
#: sits an order of magnitude clear of anything real.
HARD_MIN_ANGLE = 0.5

#: Likewise at the top. The worst legitimate value is the triangulation
#: backstop's 160.2 degrees on the round-holed disc; the ring cable's 180.000
#: is the case this catches.
HARD_MAX_ANGLE = 179.5

#: Default "low angle" for the share-of-angles metric. Reported, never gated on.
LOW_ANGLE = 20.0


def face_angles(points):
    """Interior angles of a polygon, in degrees, one per corner given.

    ``points`` are the corners the face ACTUALLY has -- for a pseudo-quad that
    is three, not four. A repeated or coincident corner yields 0.0 rather than
    raising, so a degenerate face is reported by :func:`hard_floor` instead of
    crashing the harness.
    """
    n = len(points)
    out = []
    for i in range(n):
        a, b, c = points[i - 1], points[i], points[(i + 1) % n]
        ux, uy = a[0] - b[0], a[1] - b[1]
        vx, vy = c[0] - b[0], c[1] - b[1]
        lu = (ux * ux + uy * uy) ** 0.5
        lv = (vx * vx + vy * vy) ** 0.5
        if lu < 1e-12 or lv < 1e-12:
            out.append(0.0)
            continue
        out.append(degrees(acos(max(-1.0, min(1.0, (ux * vx + uy * vy) / (lu * lv))))))
    return out


def _face_edges(points):
    """Edge lengths of a polygon, in the order its corners were given."""
    n = len(points)
    return [((points[i][0] - points[(i + 1) % n][0]) ** 2
             + (points[i][1] - points[(i + 1) % n][1]) ** 2) ** 0.5
            for i in range(n)]


def mesh_quality(mesh, low_angle=LOW_ANGLE):
    """Element quality of a quad (or pseudo-quad) mesh, as plain numbers.

    Parameters
    ----------
    mesh : Mesh
    low_angle : float, optional
        Threshold for ``share_below``, in degrees.

    Returns
    -------
    dict
        ``faces``, ``poles``, ``min_angle``, ``max_angle``, ``aspect_max``,
        ``share_below``, ``low_angle``, ``irregular_interior``, and
        ``worst_face`` -- the key of the face carrying ``min_angle``, so a
        failure can be looked at rather than merely counted.

    Aspect ratio is **longest edge over shortest edge**, per face, worst over
    the mesh. A perfect grid reads 1.00. It is ``inf`` when a face has a
    zero-length edge, which is what :func:`hard_floor` rejects on.
    """
    face_pole = mesh.attributes.get('face_pole') or {}

    lo, hi = 180.0, 0.0
    aspect = 0.0
    total = 0
    below = 0
    worst_face = None

    for fkey in mesh.faces():
        # The corners the face HAS. Never (p, a, b, p) -- see the module
        # docstring; that is the whole pseudo-quad trap.
        points = [mesh.vertex_coordinates(v) for v in mesh.face_vertices(fkey)]
        if len(points) < 3:
            lo, worst_face = 0.0, fkey
            aspect = float('inf')
            continue

        angles = face_angles(points)
        total += len(angles)
        below += sum(1 for a in angles if a < low_angle)
        if min(angles) < lo:
            lo, worst_face = min(angles), fkey
        hi = max(hi, max(angles))

        edges = _face_edges(points)
        shortest = min(edges)
        aspect = float('inf') if shortest <= 0.0 else max(aspect, max(edges) / shortest)

    irregular = 0
    for vkey in mesh.vertices():
        if mesh.is_vertex_on_boundary(vkey):
            continue
        if len(mesh.vertex_neighbors(vkey)) != 4:
            irregular += 1

    return {
        'faces': mesh.number_of_faces(),
        'poles': len(face_pole),
        'min_angle': lo,
        'max_angle': hi,
        'aspect_max': aspect,
        'share_below': (below / float(total)) if total else 0.0,
        'low_angle': low_angle,
        'irregular_interior': irregular,
        'worst_face': worst_face,
    }


#: Radius, as a multiple of the sampling spacing, within which a dense mesh
#: edge counts as being "at" a point on the curve.
ALIGNMENT_RADIUS = 1.0


def _curve_samples(curve, spacing):
    """``(point, tangent)`` every ``spacing`` units along a polyline."""
    out = []
    for a, b in zip(curve, curve[1:]):
        dx, dy = b[0] - a[0], b[1] - a[1]
        length = (dx * dx + dy * dy) ** 0.5
        if length < 1e-12:
            continue
        steps = max(1, int(round(length / spacing)))
        tangent = atan2(dy, dx)
        for k in range(steps):
            t = (k + 0.5) / steps
            out.append(((a[0] + dx * t, a[1] + dy * t), tangent))
    return out


def curve_alignment_profile(mesh, curve, radius=ALIGNMENT_RADIUS, spacing=0.5):
    """**How far a MESH is from a curve it was supposed to follow.**

    Per sample along the curve, the mean angle between the curve and every mesh
    edge whose midpoint is within ``radius``, folded into the cross's period.
    An edge crossing the curve at 90 degrees counts as aligned as one running
    along it, because a quad mesh has two families and a cable may be either.
    0 means the mesh runs on the curve; 45 is the worst a quad mesh can do.

    This is the same convention ``12_cables`` part 1 measures the FIELD with, so
    the mesh number and the field number are directly comparable -- which is the
    entire point. The field reaching 0.0 degrees off a hard-constrained cable
    while the mesh sat at 36.9 was how the gap between the two halves of the
    objective was stated; closing it has to be measured on the same scale.

    A plain MEAN over the contributing edges, deliberately. Taking the closest
    few would flatter any mesh dense enough to have one edge pointing the right
    way somewhere.

    Returns
    -------
    list[((float, float), float)]
        Sample point and its mean offset in degrees, in order along the curve.
    """
    edges = []
    for u, v in mesh.edges():
        a = mesh.vertex_coordinates(u)
        b = mesh.vertex_coordinates(v)
        if abs(a[0] - b[0]) < 1e-9 and abs(a[1] - b[1]) < 1e-9:
            continue                      # a pole's collapsed edge
        edges.append((((a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0),
                      atan2(b[1] - a[1], b[0] - a[0])))

    out = []
    for mid, tangent in _curve_samples(curve, spacing):
        vals = []
        for centre, ang in edges:
            if ((centre[0] - mid[0]) ** 2 + (centre[1] - mid[1]) ** 2
                    <= radius * radius):
                d = (tangent - ang) % (pi / 2)
                vals.append(degrees(min(d, pi / 2 - d)))
        if vals:
            out.append((mid, sum(vals) / len(vals)))
    return out


def curve_alignment(mesh, curve, radius=ALIGNMENT_RADIUS, spacing=0.5):
    """Mean of :func:`curve_alignment_profile`, in degrees.

    The mean over the whole curve understates what a guided mesh achieves,
    because a patch boundary is fixed and the mesh cannot turn where the curve
    runs up against one. Look at the profile too before concluding a cable only
    half arrived.
    """
    profile = curve_alignment_profile(mesh, curve, radius, spacing)
    if not profile:
        return float('nan')
    return sum(v for _, v in profile) / len(profile)


def hard_floor(metrics):
    """TIER 1. Is anything in this mesh degenerate rather than merely poor?

    Absolute and permanent: no baseline, and no future retuning, may bless an
    angle at or near 0 or 180 degrees or a non-finite aspect ratio. Everything
    softer than that is the regression check's business.

    Returns
    -------
    (bool, str)
        Verdict, and when false the offending metric NAMED with its value --
        the harness's whole job is to stop a number being described by hand.
    """
    if not metrics.get('faces'):
        return False, 'no faces'

    aspect = metrics['aspect_max']
    if aspect != aspect or aspect in (float('inf'), float('-inf')):
        return False, ('aspect_max is {} -- a face has a zero-length edge '
                       '(face {})'.format(aspect, metrics.get('worst_face')))

    if metrics['min_angle'] < HARD_MIN_ANGLE:
        return False, ('min_angle {:.3f} deg is below the hard floor of {} deg '
                       '(face {})'.format(metrics['min_angle'], HARD_MIN_ANGLE,
                                          metrics.get('worst_face')))

    if metrics['max_angle'] > HARD_MAX_ANGLE:
        return False, ('max_angle {:.3f} deg is above the hard floor of {} deg'.format(
            metrics['max_angle'], HARD_MAX_ANGLE))

    return True, ''
