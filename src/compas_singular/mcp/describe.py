"""**Turning the metric dict into a reading somebody can act on.**

``min_angle: 12.3`` is a number. It does not say where the bad element is,
whether one face is dragging the figure down or the whole mesh is like that, or
whether 12.3 is a defect worth fixing on this mesh at all. A model handed only
the number will either ignore it or over-react to it, and both are wasted steps.

So every report carries a reading as well as the numbers: what the worst element
is, where it is (as a handle -- ``worst_face`` is a face KEY, which must never
cross the protocol), how concentrated the problem is, and a one-word verdict.

**The thresholds are data, not code.** They live in
``library/thresholds.json`` and are read through
:func:`~compas_singular.mcp.library.thresholds`, so what counts as "good" can be
retuned for a project without touching this file or cutting a release. The
defaults below are the fallback when that file is missing or unreadable, so the
server still starts and still says something sensible.

**It reports, it does not decide.** Nothing here refuses a mesh or blocks a
step. ``framefield.relax`` has a real gate with real measurements behind it; a
verdict from this module is a summary for a reader, and calling it a pass or a
fail would give it an authority it has not earned.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from compas_singular.mcp.handle import vertex_handle


__all__ = ['DEFAULT_THRESHOLDS', 'band', 'describe', 'format_triple']


#: Fallback bands. ``min_angle`` and the rest are read as: at least ``good`` is
#: good, at least ``usable`` is usable, at least ``poor`` is poor, below that
#: unusable. ``max_angle`` and ``aspect_max`` run the other way.
DEFAULT_THRESHOLDS = {
    'min_angle': {'good': 45.0, 'usable': 25.0, 'poor': 15.0, 'higher_is_better': True},
    'max_angle': {'good': 135.0, 'usable': 155.0, 'poor': 168.0, 'higher_is_better': False},
    'aspect_max': {'good': 2.0, 'usable': 4.0, 'poor': 8.0, 'higher_is_better': False},
    'share_below': {'good': 0.0, 'usable': 0.02, 'poor': 0.10, 'higher_is_better': False},
}

BANDS = ('good', 'usable', 'poor', 'unusable')


def band(value, limits):
    """Which band a value falls in. ``'unknown'`` when it cannot be judged."""
    if value is None:
        return 'unknown'
    try:
        value = float(value)
    except (TypeError, ValueError):
        return 'unknown'
    if value != value:                       # NaN
        return 'unknown'
    higher_is_better = limits.get('higher_is_better', True)
    for name in ('good', 'usable', 'poor'):
        edge = limits.get(name)
        if edge is None:
            continue
        if higher_is_better and value >= edge:
            return name
        if not higher_is_better and value <= edge:
            return name
    return 'unusable'


def _worst(bands):
    """The least good of the bands seen."""
    for name in reversed(BANDS):
        if name in bands:
            return name
    return 'unknown'


def format_triple(triple):
    """``(min, max, aspect)`` as a short string, or ``'-'``."""
    if not triple:
        return '-'
    lo, hi, aspect = triple
    return '{:.1f} / {:.1f} deg, aspect {:.2f}'.format(
        lo if lo is not None else float('nan'),
        hi if hi is not None else float('nan'),
        aspect if aspect is not None else float('nan'))


def _face_location(mesh, fkey):
    """A handle for where a face is, and whether it touches a boundary.

    The face's own key is useless to a caller -- it renumbers, and it is not a
    place. Its centroid, named the way a vertex is named, is both.
    """
    if fkey is None or mesh is None:
        return None, False
    try:
        vertices = mesh.face_vertices(fkey)
    except Exception:
        return None, False
    if not vertices:
        return None, False
    points = [mesh.vertex_coordinates(key) for key in vertices]
    centroid = [sum(p[i] for p in points) / float(len(points)) for i in range(3)]
    on_boundary = any(mesh.is_vertex_on_boundary(key) for key in vertices)
    return vertex_handle(centroid), on_boundary


def describe(mesh, metrics, limits=None):
    """A prose reading of a quality dict.

    Parameters
    ----------
    mesh : Mesh or None
        Used only to locate the worst face. The reading degrades gracefully
        without it.
    metrics : dict
        From ``framefield.quality.mesh_quality``.
    limits : dict, optional
        Thresholds. Defaults to :data:`DEFAULT_THRESHOLDS`.

    Returns
    -------
    dict
        ``reading`` (a paragraph), ``verdict`` (one of :data:`BANDS`),
        ``bands`` (per metric), and ``worst_face_at`` (a handle, or None).
    """
    if not metrics:
        return {'reading': 'Nothing is loaded.', 'verdict': 'unknown',
                'bands': {}, 'worst_face_at': None}

    limits = limits or DEFAULT_THRESHOLDS
    bands = dict((name, band(metrics.get(name), limits[name]))
                 for name in limits if name in limits)

    where, on_boundary = _face_location(mesh, metrics.get('worst_face'))
    sentences = []

    min_angle = metrics.get('min_angle')
    if min_angle is not None:
        place = ' at {}'.format(where) if where else ''
        edge = ' on the boundary' if on_boundary else ' in the interior'
        sentences.append(
            'Worst element is a {:.1f} degree corner{}{}.'.format(
                min_angle, edge, place))

    max_angle = metrics.get('max_angle')
    if max_angle is not None:
        sentences.append('Widest corner is {:.1f} degrees.'.format(max_angle))

    aspect = metrics.get('aspect_max')
    if aspect is not None:
        if aspect == float('inf'):
            sentences.append(
                'Aspect ratio is infinite: a face has a zero-length edge, which '
                'is degenerate rather than merely poor.')
        else:
            sentences.append('Worst aspect ratio is {:.2f} (a perfect grid '
                             'reads 1.00).'.format(aspect))

    # How WIDESPREAD the low angles are is the part that decides whether to
    # smooth locally or globally, so it is stated in those terms.
    share = metrics.get('share_below')
    low = metrics.get('low_angle')
    if share is not None and low is not None:
        if share > 0.0:
            sentences.append(
                '{:.1f}% of all corners are below {:.0f} degrees, so the problem '
                'is {} rather than a single face.'.format(
                    100.0 * share, low,
                    'widespread' if share > 0.05 else 'in more than one place'))
        elif bands.get('min_angle') in ('poor', 'unusable', 'usable'):
            # Nothing is below the low-angle line, yet the minimum is still not
            # good: the figure is one face, not a general condition.
            sentences.append(
                'No corner is below {:.0f} degrees, so that minimum is one '
                'local defect rather than a general one.'.format(low))
        else:
            sentences.append(
                'No corner is below {:.0f} degrees.'.format(low))

    faces = metrics.get('faces')
    irregular = metrics.get('irregular_interior')
    poles = metrics.get('poles')
    if faces is not None:
        tail = '{} faces'.format(faces)
        if irregular:
            tail += ', {} irregular interior vertices'.format(irregular)
        if poles:
            tail += ', {} poles'.format(poles)
        # NOT ``capitalize()``: it lowercases everything after the first
        # character, which would turn "12 faces, 4 Poles" into mush.
        sentences.append('Mesh has ' + tail + '.')

    verdict = _worst([b for b in bands.values() if b != 'unknown'])
    return {'reading': ' '.join(sentences),
            'verdict': verdict,
            'bands': bands,
            'worst_face_at': where}
