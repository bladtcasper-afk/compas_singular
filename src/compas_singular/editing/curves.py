"""**Giving a moved coarse edge back the curve it was traced with.**

A coarse edge is a straight chord. Everything a layout knows about curvature
lives in a separate mapping -- one polyline per edge -- that is handed to
``densification``. On a GENERATED layout that mapping is easy: every edge either
runs along a domain wall, and gets the wall's arc, or it was traced as a
separatrix, and matches one by the geometric key of its two ends.

**An edit breaks the match.** Drag a corner and the edges around it no longer end
where any traced separatrix ends, so the exact lookup fails and they fall back to
the chord between the moved corners -- which is precisely the field alignment the
whole front end exists to produce, thrown away by a nudge.

This module is the branch in between. It finds the separatrix an edge came FROM
and moves that curve onto the edge's new endpoints, so a nudge costs the nudge and
not the curvature.

**It needs no field.** The warp itself is
:func:`~compas_singular.editing.rebuild.warp_polyline`, which is arithmetic on a
point list, and the search below reads only the traced polylines and one length
scale. This used to live on ``FieldDecomposition`` and looked field-dependent
because it read ``self.background.target_length``; that is a number, not a field.

**The match is ANCHORED on an end that did not move**, and that rule is what makes
it well posed rather than a nearest-neighbour guess. It follows from what an edit
is: dragging one corner moves one end of each edge around it and leaves the other
exactly where the tracer put it, so the right separatrix is the one still ending
at the still point -- and among those, the one whose free end was nearest. Nothing
has to cap the end that MOVED, which is the one the user may have dragged as far
as they liked. An earlier version capped both ends at one background spacing and
rejected every real edit: a 1.8-unit drag on a 0.6 background left five edges as
straight chords, which is the exact loss this branch exists to prevent.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from compas.geometry import distance_point_point
from compas.itertools import pairwise

from compas_singular.editing.rebuild import warp_polyline


__all__ = ['warp_edge_curve', 'warp_chorded_edges']


def warp_edge_curve(polylines, pa, pb, scale, claimed=None):
    """**The traced separatrix this edge came from, moved onto its endpoints.**

    Parameters
    ----------
    polylines : list[list[[x, y, z]]]
        The traced separatrix network, plus any curve the user drew.
    pa, pb : [x, y, z]
        Where the edge's two ends are NOW.
    scale : float
        The length the two tolerances are expressed in -- the background spacing
        the separatrices were traced at, where that is known. ``0.1 * scale`` is
        how close an end has to be to count as not having moved; ``scale`` is how
        far both ends may be when neither is anchored.
    claimed : set, optional
        Indices already taken by another edge. A polyline is only reused when
        nothing unclaimed fits -- the case where one patch was split in two and
        both halves genuinely lie along one separatrix. Mutated in place.

    Returns
    -------
    list[[x, y, z]] or None
        ``None`` when nothing fits, so the caller keeps the straight chord.

    Notes
    -----
    Two guards stay, because the WRONG curve warped is worse than no curve at
    all -- it densifies into a patch that bulges through its neighbour rather
    than one that is merely straight:

    * ``warp_polyline`` rejects a warp reaching further than the candidate's own
      length, which would fold it;
    * a candidate whose length is wildly out of proportion to the new chord is a
      different curve, not this one moved.
    """
    if not polylines:
        return None
    if claimed is None:
        claimed = set()

    anchor = scale * 0.1
    loose = scale
    chord = distance_point_point(pa, pb)

    best = None
    for i, polyline in enumerate(polylines):
        points = [list(p)[:3] for p in polyline]
        if len(points) < 2:
            continue
        for candidate in (points, list(reversed(points))):
            da = distance_point_point(candidate[0], pa)
            db = distance_point_point(candidate[-1], pb)

            if da <= anchor or db <= anchor:
                # anchored: one end never moved, so trust the other freely
                tier, free = 0, (db if da <= anchor else da)
            elif da <= loose and db <= loose:
                # neither end is anchored -- both corners of this edge moved. The
                # correspondence really is a guess here, so it is a timid one.
                tier, free = 1, da + db
            else:
                continue

            length = sum(distance_point_point(a, b) for a, b in pairwise(candidate))
            if chord > 1e-9 and not (0.25 * chord <= length <= 4.0 * chord):
                continue

            score = (i in claimed, tier, free)
            if best is None or score < best[0]:
                best = (score, i, candidate)

    if best is None:
        return None
    _score, index, points = best
    curve = warp_polyline(points, pa, pb)
    if curve is None:
        return None
    claimed.add(index)
    return curve


def warp_chorded_edges(mapping, coarse, polylines, scale):
    """Give every straight-chord edge a warped curve where one fits.

    ``mapping`` is what ``coarse_edges_to_curves`` produced: one polyline per
    coarse edge, complete. An entry of exactly two points IS the chord -- that is
    how an edge that matched nothing comes back -- so those are the entries to
    try, and every other entry is left alone because it already has the shape it
    should.

    Returns ``(mapping, warped)``. The mapping is modified in place and returned
    for convenience.
    """
    claimed = set()
    warped = 0
    for edge, curve in list(mapping.items()):
        if len(curve) > 2:
            continue                      # already a wall arc or an exact match
        u, v = edge
        pa = coarse.vertex_coordinates(u)
        pb = coarse.vertex_coordinates(v)
        out = warp_edge_curve(polylines, pa, pb, scale, claimed)
        if out is not None:
            mapping[edge] = out
            warped += 1
    return mapping, warped
