"""**Reading a mesh, and the three ways of improving one.**

Each tool takes the session first and keyword arguments after, and returns a
plain dict. The descriptions carry the measurements, because an automated caller
with no context picks the worst setting available and every one of those is
cheap to prevent in a sentence and expensive to discover in a loop.

**The ungated passes snapshot themselves.** ``relax`` has a real gate and cannot
make a mesh worse; the other two apply what they are told and can. Rather than
rely on a caller remembering to snapshot first, they take one, and report
``all_improved`` so the decision to undo is a reading rather than a judgement.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import math

from ..datastructures.mesh.smoothing import automated_boundary_constraints
from ..datastructures.mesh.smoothing import boundary_constrained_smoothing
from ..datastructures.mesh.smoothing import constrained_smoothing
from ..datastructures.mesh.smoothing import relaxation
from ..datastructures.mesh.smoothing import smoothing_region
from ..editing.guide_chain import GuideCurve
from ..editing.guide_chain import attach_chain
from ..editing.guide_chain import collect_polyedges
from ..editing.guide_chain import guide_chain
from ..editing.guide_chain import mean_edge_length
from ..framefield.relax import relax_mesh
from . import render
from .describe import describe
from .handle import SelectorError
from .handle import describe_selector
from .handle import resolve
from .handle import select
from .handle import vertex_handle
from .library import thresholds
from .registry import tool


__all__ = []


#: What to actually look at, sent WITH every image. A picture handed over
#: without this gets glanced at; the numbers are already in the payload, so the
#: only reason to render one is the four things below, and none of them are
#: measured anywhere.
VISUAL_CHECKLIST = """Look at the image and judge these four, none of which any
number in this report covers:

1. BOUNDARIES. The input walls are drawn in orange UNDERNEATH the mesh. Where
   the mesh sits on its boundary the orange is hidden. Anywhere orange shows
   along an edge, the mesh has left the wall it was supposed to follow -- say
   where.
2. POINT FEATURES. Each input point feature is a magenta RING. A respected one
   has a magenta DISC (a pole) inside it. A ring with no disc in it is an input
   the mesh ignored.
3. GUIDES. Input guide curves are green, also drawn underneath. Say whether the
   mesh runs along them or crosses them.
4. SMOOTHNESS AND POLYEDGE CONTINUITY. Follow the grey grid lines across the
   mesh. They should run as continuous, gently curving families. Report kinks,
   abrupt direction changes, strips that pinch or fan, and any lumpiness or
   waviness in a region that should be regular. Around a red disc (an irregular
   vertex) some distortion is unavoidable and is NOT a defect -- what matters is
   whether it stays local or propagates away across the mesh.

Report what you actually see, including "nothing wrong" if that is the case. Do
not restate the numbers; they are above."""


def _needs_mesh(session):
    """A refusal when nothing is loaded, or ``None`` when there is."""
    if not session.loaded:
        return {'ok': False,
                'reason': 'no mesh is loaded -- call rhino_pull to take one '
                          'from the Rhino document, or load_mesh to read one '
                          'from a file'}
    return None


def _reading(session):
    """The quality dict and its prose, together."""
    metrics = session.quality()
    told = describe(session.mesh, metrics, thresholds())
    return metrics, told


def _improved(before, after):
    """Which of the three got better, and whether all of them did.

    The same three ``framefield.relax`` gates on. A pass that improves one and
    wrecks another has not improved the mesh, so ``all_improved`` is the number
    that matters and it is reported separately from the individual flags.
    """
    keys = ('min_angle', 'max_angle', 'aspect_max')
    higher_is_better = {'min_angle': True, 'max_angle': False, 'aspect_max': False}
    flags = {}
    for key in keys:
        old, new = (before or {}).get(key), (after or {}).get(key)
        if old is None or new is None:
            flags[key] = None
            continue
        flags[key] = (new > old) if higher_is_better[key] else (new < old)
    decided = [v for v in flags.values() if v is not None]
    unchanged = all(
        (before or {}).get(k) == (after or {}).get(k) for k in keys)
    return flags, (bool(decided) and all(decided)), unchanged


def _outcome(session, action, before, **detail):
    """Measure, record and report a step that has just been applied."""
    after, told = _reading(session)
    # Whatever the numbers say, vertices were moved through: the last picture
    # no longer shows the mesh in hand. Marking here rather than in each tool
    # means a pass added later cannot forget to.
    session.mark_changed()
    flags, all_improved, unchanged = _improved(before, after)
    if unchanged:
        # Nothing moved, so "did it improve" has no useful answer. Reporting
        # False here reads as "this made it worse, undo it" and buys a wasted
        # step restoring positions that are already identical.
        all_improved = None
    entry = session.record(action, before=before, after=after,
                           improved=flags, all_improved=all_improved, **detail)
    payload = {'ok': True, 'action': action, 'step': entry['step'],
               'before': before, 'after': after,
               'improved': flags, 'all_improved': all_improved,
               'unchanged': unchanged,
               'reading': told['reading'], 'verdict': told['verdict'],
               'undo_available': session.can_undo()}
    payload.update(detail)
    if unchanged:
        payload['note'] = ('nothing moved. Trying the same pass again with '
                           'different numbers rarely helps -- read '
                           'guidance://smoothing before spending another step.')
    elif all_improved is False:
        payload['note'] = ('this did NOT improve all three of min angle, max '
                           'angle and aspect. A snapshot was taken before it; '
                           'call undo to take it back.')
    return payload


# ==============================================================================
# reading
# ==============================================================================

@tool(
    'inspect',
    'Measure the mesh in hand and say what is wrong with it, in numbers and in '
    'prose, and OPTIONALLY draw it. Returns min/max corner angle, worst aspect '
    'ratio, share_below (the fraction of corners under the low-angle line), the '
    'count of irregular interior vertices and poles, a handle locating the '
    'worst face, and a reading with a verdict. READ share_below BEFORE choosing '
    'a smoother: it is what separates one bad face from a generally slack mesh, '
    'and therefore what decides between smooth_region and a global pass. Note '
    'that irregular vertices and poles are NOT defects -- their count is fixed '
    'by the shape of the domain, and a change in it means the topology moved. '
    'SET image=true TO SEE THE MESH. The numbers cannot show whether the mesh '
    'left its boundary, ignored a guide, dropped a point feature, or has a '
    'polyedge that kinks -- only a picture can, so leave image off while '
    'iterating and turn it on when you need to judge those. YOU MUST CALL THIS '
    'WITH image=true ONCE MORE BEFORE YOU FINISH: rhino_push and save_mesh '
    'refuse until the current mesh has actually been looked at.',
    properties={
        'low_angle': {
            'type': 'number',
            'description': 'Degrees below which a corner counts toward '
                           'share_below. Default 20.'},
        'image': {
            'type': 'boolean',
            'description': 'Draw the mesh with its input boundaries, guides and '
                           'point features, and return it as an image. Costs '
                           'roughly 850 tokens at the default size. Off by '
                           'default so iterating is cheap; REQUIRED once on the '
                           'final mesh before it can be delivered.'},
        'size': {
            'type': 'integer',
            'description': 'Image width and height in pixels. Default 800, max '
                           '1400. Cost runs with width times height, so raise '
                           'it only to read a detail you could not otherwise.'},
    },
    read_only=True, idempotent=True, title='Inspect the mesh')
def _t_inspect(session, low_angle=None, image=False, size=None):
    refusal = _needs_mesh(session)
    if refusal:
        return refusal
    metrics = session.quality(low_angle=low_angle)
    told = describe(session.mesh, metrics, thresholds())
    payload = {'ok': True, 'quality': metrics, 'reading': told['reading'],
               'verdict': told['verdict'], 'bands': told['bands'],
               'worst_face_at': told['worst_face_at'],
               'inputs': {'walls': len(session.walls),
                          'guides': len(session.guides),
                          'point_features': len(session.points)},
               'singularities': _singularity_handles(session.mesh)}

    if image:
        size = render.DEFAULT_SIZE if size is None else int(size)
        try:
            encoded = render.render_png_base64(session, width=size, height=size)
        except Exception as exc:
            payload['image_error'] = 'could not draw the mesh: {}: {}'.format(
                type(exc).__name__, exc)
            encoded = None
        if encoded:
            # Only a rendered AND delivered picture counts as having looked.
            session.mark_seen()
            payload['image_legend'] = render.LEGEND
            payload['look_at_this'] = VISUAL_CHECKLIST
            payload['_image'] = {'data': encoded, 'mimeType': 'image/png'}
            if not session.walls:
                payload['image_caveat'] = (
                    'No input boundary curves are loaded, so nothing is drawn '
                    'in orange and check 1 cannot be made. Pull from Rhino, or '
                    'pass walls to load_mesh, if the boundary matters.')
            if not session.points:
                payload.setdefault('image_caveat', '')
                payload['image_caveat'] = (payload['image_caveat'] + ' No input '
                    'point features are loaded, so check 2 does not apply.').strip()

    payload['state'] = session.state()
    return payload


def _singularity_handles(mesh, limit=24):
    """Where the irregular vertices are, as handles. Truncated, and says so."""
    from .handle import select as _select
    try:
        keys, _ = _select(mesh, {'kind': 'singularities', 'rings': 0})
    except Exception:
        return {'count': None, 'at': []}
    at = sorted(vertex_handle(mesh.vertex_coordinates(key)) for key in keys)
    out = {'count': len(at), 'at': at[:limit]}
    if len(at) > limit:
        out['truncated'] = 'showing {} of {}'.format(limit, len(at))
    return out


@tool(
    'history',
    'Every step taken in this session, with the quality before and after each '
    'one and any remarks attached to it. Use this to see what has already been '
    'tried before spending another step, and after a context compaction to '
    'recover what you were doing. A step an undo took back is still listed, '
    'with undone=true -- it is what was tried and rejected.',
    read_only=True, idempotent=True, title='Session history')
def _t_history(session):
    steps = []
    for entry in session.history:
        steps.append({
            'step': entry['step'],
            'action': entry['action'],
            'layer': entry.get('layer'),
            'undone': bool(entry.get('undone')),
            'before': entry.get('before'),
            'after': entry.get('after'),
            'all_improved': entry.get('all_improved'),
            'remarks': session.remarks_for(entry['step']),
        })
    return {'ok': True, 'steps': steps, 'state': session.state(),
            'remarks': session.remarks}


@tool(
    'remark',
    'Attach a note to the session: why you chose a setting, what you tried and '
    'undid, what is still wrong. Remarks come back from history and are carried '
    'into the corpus by save_example, so a remark written now is what a later '
    'session reads. Record what the numbers cannot -- the reasoning, the '
    'rejected attempt, and especially anything still wrong that smoothing '
    'cannot fix. See guidance://remarks for what is worth keeping.',
    properties={
        'text': {'type': 'string', 'description': 'The note.'},
        'about': {
            'type': 'string',
            'description': 'Optional vertex handle this is about, e.g. '
                           '"v:8.500,3.250,0.000", from a report.'},
    },
    required=('text',), title='Leave a remark')
def _t_remark(session, text, about=None):
    if not (text or '').strip():
        return {'ok': False, 'reason': 'a remark needs some text'}
    if about is not None and session.loaded:
        from .handle import resolve
        key, how = resolve(session.mesh, about)
        if key is None:
            return {'ok': False, 'reason': 'cannot anchor the remark: ' + how}
    remark = session.add_remark(text.strip(), about=about)
    return {'ok': True, 'remark': remark, 'total': len(session.remarks)}


@tool(
    'snapshot',
    'Remember where every vertex is, so a later undo can come back here. The '
    'two ungated smoothing tools snapshot themselves, so this is only needed to '
    'mark a point you chose -- before a sequence of passes, say. Cheap: a '
    'snapshot is a map of positions, not a copy of the mesh.',
    properties={'label': {'type': 'string',
                          'description': 'What this point is, for the report.'}},
    title='Snapshot positions')
def _t_snapshot(session, label=''):
    refusal = _needs_mesh(session)
    if refusal:
        return refusal
    session.snapshot(label or 'manual snapshot')
    return {'ok': True, 'label': label or 'manual snapshot',
            'undo_depth': session.state()['undo_depth']}


@tool(
    'undo',
    'Put every vertex back where the last snapshot had it, and mark the '
    'dense-mesh steps taken since as undone (coarse-layout steps are '
    "coarse_undo's business and are left alone). Use it when a pass reports "
    'all_improved false. Refuses if the topology has changed since the '
    'snapshot, because a map of positions cannot restore that -- no tool in '
    'this version changes dense topology, so that refusal should not occur.',
    destructive=True, title='Undo to last snapshot')
def _t_undo(session):
    ok, detail = session.undo()
    if not ok:
        return {'ok': False, 'reason': detail}
    metrics, told = _reading(session)
    return {'ok': True, 'restored': detail, 'quality': metrics,
            'reading': told['reading'], 'verdict': told['verdict'],
            'state': session.state()}


# ==============================================================================
# improving
# ==============================================================================

@tool(
    'relax',
    'Global relaxation with an accept/reject gate -- TRY THIS FIRST. Walks a '
    'schedule of settings gentlest-first ((20, 0.3), then (50, 0.5), then '
    '(100, 0.5)) and keeps a result ONLY if min angle, max angle and max aspect '
    'all improve together; a setting that trades one against another is '
    'rejected. Because of that gate it cannot make the mesh worse. '
    'accepted=null means no setting improved all three at once, so the mesh was '
    'left exactly as it was -- that is a PASS, not a failure, and retrying with '
    'different numbers rarely changes it. Leave pin_singularities off: pinning '
    'the irregular vertices is what stops this pass working, because they are '
    'the vertices that must move for the elements around them to open up.',
    properties={
        'seams': {
            'type': 'string', 'enum': ['free', 'slide', 'fixed'],
            'description': "What patch-seam vertices may do. 'free' is the "
                           "default and the only one measured to help. 'slide' "
                           "never hurts and never helps. 'fixed' relaxes patch "
                           "interiors only."},
        'corner_angle': {
            'type': 'number',
            'description': 'Degrees of kink past which a boundary vertex is a '
                           'domain corner and is held. Default 30.'},
        'pin_singularities': {
            'type': 'boolean',
            'description': 'Hold irregular vertices and poles. Off by default; '
                           'turning it on is what stops this pass working.'},
    },
    title='Relax (gated)')
def _t_relax(session, seams='free', corner_angle=30.0, pin_singularities=False):
    refusal = _needs_mesh(session)
    if refusal:
        return refusal
    if seams not in ('free', 'slide', 'fixed'):
        return {'ok': False,
                'reason': "seams must be 'free', 'slide' or 'fixed'; got "
                          '{!r}'.format(seams)}
    before = session.quality()
    session.snapshot('before relax')
    try:
        report = relax_mesh(session.mesh, corner_angle=float(corner_angle),
                            seams=seams,
                            pin_singularities=bool(pin_singularities))
    except Exception as exc:
        return {'ok': False,
                'reason': 'relax failed: {}: {}'.format(type(exc).__name__, exc)}

    accepted = report.get('accepted')
    payload = _outcome(session, 'relax', before,
                       accepted=accepted, seams=seams,
                       pin_singularities=bool(pin_singularities),
                       sliding=report.get('sliding'),
                       pinned=report.get('pinned'))
    if accepted is None:
        payload['note'] = ('accepted is null: no setting in the schedule '
                           'improved all three measures at once, so the mesh '
                           'was left untouched. This is a pass. Reach for a '
                           'local pass, or report that the layout is the limit '
                           '-- do not retry relax with other numbers.')
    else:
        payload['note'] = 'accepted the setting (kmax, damping) = {}'.format(
            tuple(accepted))
    return payload


@tool(
    'smooth_boundary_constrained',
    'The workhorse: interior vertices relax freely, boundary vertices SLIDE '
    'along the walls, and kink vertices stay pinned. Reach for this when the '
    'mesh is generally slack rather than locally broken -- when share_below is '
    'above zero. UNLIKE relax THIS HAS NO GATE: it applies what it is told and '
    'can make the mesh worse, so it snapshots first and reports all_improved; '
    'undo if that is false. Without walls (a mesh loaded from file rather than '
    'pulled from Rhino) it falls back to the mesh boundary, which means the '
    'boundary can only get smoother, never truer -- a faceted boundary stays '
    'faceted.',
    properties={
        'kmax': {'type': 'integer',
                 'description': 'Iterations. Default 100. Cost is linear and '
                                'returns diminish past a few hundred.'},
        'damping': {'type': 'number',
                    'description': 'Between 0 and 1. Default 0.5.'},
        'algorithm': {
            'type': 'string', 'enum': ['centroid', 'area', 'centerofmass'],
            'description': "Default 'centroid'. Use 'area' on a mesh meant to "
                           'be graded: centroid equalises edge lengths, which '
                           'fights the grading.'},
        'corner_angle': {
            'type': 'number',
            'description': 'Degrees of kink past which a boundary vertex is a '
                           'corner and is pinned. Default 30. DEGREES here, '
                           'converted internally.'},
        'fix_corners': {'type': 'boolean',
                        'description': 'Pin the kink vertices. Default true.'},
    },
    title='Smooth, boundary sliding')
def _t_smooth_boundary(session, kmax=100, damping=0.5, algorithm='centroid',
                       corner_angle=30.0, fix_corners=True):
    refusal = _needs_mesh(session)
    if refusal:
        return refusal
    if algorithm not in ('centroid', 'area', 'centerofmass'):
        return {'ok': False,
                'reason': "algorithm must be 'centroid', 'area' or "
                          "'centerofmass'; got {!r}".format(algorithm)}
    before = session.quality()
    session.snapshot('before smooth_boundary_constrained')
    try:
        constraints = boundary_constrained_smoothing(
            session.mesh,
            curves=session.walls or None,
            kmax=int(kmax), damping=float(damping), algorithm=algorithm,
            # The library takes radians here; the tool takes degrees, because
            # handing a model radians is a reliable way to get a wrong call.
            corner_angle=math.radians(float(corner_angle)),
            fix_corners=bool(fix_corners))
    except Exception as exc:
        return {'ok': False,
                'reason': 'smoothing failed: {}: {}'.format(
                    type(exc).__name__, exc)}
    return _outcome(session, 'smooth_boundary_constrained', before,
                    kmax=int(kmax), damping=float(damping),
                    algorithm=algorithm, corner_angle=float(corner_angle),
                    constrained_vertices=len(constraints or {}),
                    walls_used=len(session.walls))


@tool(
    'smooth_region',
    'Smooth ONE PART of the mesh, leaving the rest fixed. Reach for this when '
    'share_below is zero or near it and min_angle is still bad: that '
    'combination means a single bad face, and a global pass would move the '
    'whole mesh to fix it. Area-weighted rather than centroid, because centroid '
    'equalises edge lengths and fights the grading a frame-field mesh is '
    'supposed to have. The damping tapers to zero over blend rings around the '
    'region, so it blends instead of leaving a crease -- do NOT set blend to 0, '
    'because a fully relaxed region against a fixed ring is itself a kink. Has '
    'no gate, so it snapshots first and reports all_improved.',
    properties={
        'region': {
            'type': 'object',
            'description': 'What to smooth. One of: {"kind":"worst_faces",'
                           '"count":3} -- usually the right choice; '
                           '{"kind":"near_point","point":[x,y,z],"radius":r}; '
                           '{"kind":"singularities","rings":2}; '
                           '{"kind":"handles","handles":["v:1.000,2.000,0.000"]}; '
                           '{"kind":"interior"}; {"kind":"all"}.'},
        'kmax': {'type': 'integer', 'description': 'Iterations. Default 50.'},
        'damping': {'type': 'number',
                    'description': 'Between 0 and 1. Default 0.5.'},
        'blend': {'type': 'integer',
                  'description': 'Rings over which the damping falls to zero '
                                 'outside the region. Default 3. Never 0.'},
    },
    required=('region',), title='Smooth a region')
def _t_smooth_region(session, region, kmax=50, damping=0.5, blend=3):
    refusal = _needs_mesh(session)
    if refusal:
        return refusal
    try:
        keys, note = select(session.mesh, region)
    except SelectorError as exc:
        return {'ok': False, 'reason': str(exc)}
    if not keys:
        return {'ok': False,
                'reason': 'that region selected no vertices ({})'.format(note)}
    if int(blend) < 1:
        return {'ok': False,
                'reason': 'blend must be at least 1 -- a fully relaxed region '
                          'against a fixed ring leaves a crease, which is the '
                          'kind of defect this tool exists to remove'}

    before = session.quality()
    session.snapshot('before smooth_region')
    try:
        weights = smoothing_region(
            session.mesh, keys, kmax=int(kmax), damping=float(damping),
            blend=int(blend),
            # None means the region's own boundary vertices slide along the
            # walls instead of being dragged inward.
            constraints=None)
    except Exception as exc:
        return {'ok': False,
                'reason': 'region smoothing failed: {}: {}'.format(
                    type(exc).__name__, exc)}
    return _outcome(session, 'smooth_region', before,
                    region=describe_selector(region), selected=note,
                    core_vertices=len(keys), moved_vertices=len(weights or {}),
                    kmax=int(kmax), damping=float(damping), blend=int(blend))


def _guide_points(curve):
    """Plain points from a ``Polyline`` or an already-bare point list."""
    return getattr(curve, 'points', curve)


@tool(
    'smooth_guides',
    'Attach the mesh to its guide curves and smooth. For each guide curve '
    'pulled from Rhino, chooses the longest run of one polyedge that already '
    'follows it -- not built up vertex by vertex, CHOSEN, so it cannot fold a '
    'face the way a radius-based pick can -- moves it onto the guide, and '
    'runs a constrained smoothing pass with every attached vertex, plus the '
    'boundary, held. Every guide is selected BEFORE any is attached, so the '
    'result does not depend on guide order. A boundary vertex is never moved '
    'onto a guide -- at most it slides along its own wall, since moving it '
    'would take the wall with it. UNLIKE relax THIS HAS NO GATE: it applies '
    'what it is told and can make the mesh worse, so it snapshots first and '
    'reports all_improved. Refuses if no guide curves are loaded.',
    properties={
        'tolerance_factor': {
            'type': 'number',
            'description': 'How far off a guide a chain vertex may sit, as a '
                           'multiple of the mesh mean edge length. Default '
                           '2.0.'},
        'max_angle': {
            'type': 'number',
            'description': "Degrees the polyedge may run off the guide's "
                           'tangent and still be selected -- this is what '
                           'cuts a chain where it veers away near the ends. '
                           'Default 30. Pass 90 to turn this gate off and '
                           'select on distance alone.'},
        'hold': {
            'type': 'string', 'enum': ['fixed', 'sliding'],
            'description': "How an attached INTERIOR vertex is held. "
                           "'fixed' (default) pins it where it lands; "
                           "'sliding' re-projects it every iteration so it "
                           'may travel along the guide but never leave it.'},
        'boundary': {
            'type': 'string', 'enum': ['sliding', 'fixed'],
            'description': "What an UNATTACHED boundary vertex does. "
                           "'sliding' (default) lets it slide along the "
                           "walls; 'fixed' pins it."},
        'kmax': {'type': 'integer',
                 'description': 'Iterations. Default 100.'},
        'damping': {'type': 'number',
                    'description': 'Between 0 and 1. Default 0.5.'},
    },
    title='Smooth to guide curves')
def _t_smooth_guides(session, tolerance_factor=2.0, max_angle=30.0,
                     hold='fixed', boundary='sliding', kmax=100, damping=0.5):
    refusal = _needs_mesh(session)
    if refusal:
        return refusal
    if not session.guides:
        return {'ok': False,
                'reason': 'no guide curves are loaded -- rhino_pull reads '
                          'them from the document; without any there is '
                          'nothing to attach to'}
    if hold not in ('fixed', 'sliding'):
        return {'ok': False,
                'reason': "hold must be 'fixed' or 'sliding'; got "
                          '{!r}'.format(hold)}
    if boundary not in ('sliding', 'fixed'):
        return {'ok': False,
                'reason': "boundary must be 'sliding' or 'fixed'; got "
                          '{!r}'.format(boundary)}

    mesh = session.mesh
    try:
        polyedges = collect_polyedges(mesh)
    except Exception as exc:
        return {'ok': False,
                'reason': 'could not collect polyedges: {}: {}'.format(
                    type(exc).__name__, exc)}
    average = mean_edge_length(mesh)
    tolerance = float(tolerance_factor) * average

    proposals, refused = [], []
    for index, guide in enumerate(session.guides):
        guide_curve = GuideCurve(_guide_points(guide))
        selected, info = guide_chain(mesh, guide_curve, tolerance=tolerance,
                                     max_angle=float(max_angle),
                                     polyedges=polyedges)
        if not selected:
            refused.append({'guide': index, 'reason': info.get('reason', '')})
            continue
        proposals.append((index, guide_curve, selected))

    if not proposals:
        return {'ok': False,
                'reason': 'no guide attached a chain -- widen tolerance_factor '
                          'or max_angle',
                'refused': refused}

    before = session.quality()
    session.snapshot('before smooth_guides')

    attached, moved = {}, set()
    for index, guide_curve, selected in proposals:
        moves, constraints = attach_chain(mesh, selected, guide_curve, hold=hold)
        for vertex, xyz in moves.items():
            mesh.vertex_attributes(vertex, 'xyz', xyz)
        attached.update(constraints)
        moved.update(moves)

    if boundary == 'sliding':
        constraints = automated_boundary_constraints(
            mesh, curves=session.walls or None)
        fixed = None
    else:
        constraints = {}
        fixed = [v for loop in mesh.vertices_on_boundaries() for v in loop
                if v not in attached]
    constraints.update(attached)

    try:
        constrained_smoothing(mesh, kmax=int(kmax), damping=float(damping),
                              constraints=constraints, algorithm='area',
                              fixed=fixed)
    except Exception as exc:
        return {'ok': False,
                'reason': 'smoothing failed: {}: {}'.format(
                    type(exc).__name__, exc)}

    return _outcome(session, 'smooth_guides', before,
                    guides_attached=len(proposals), guides_refused=refused,
                    vertices_moved_onto_guides=len(moved),
                    vertices_constrained=len(constraints),
                    hold=hold, boundary=boundary,
                    kmax=int(kmax), damping=float(damping))


@tool(
    'relax_fdm',
    'Force-density relaxation: fixed vertices held, everything else finds a '
    'minimal-tension shape under uniform force density, with boundary edges '
    'weighted q_factor times heavier than interior ones so the outline holds '
    'its shape. Needs the compas_fd package; a clean refusal comes back if it '
    "is not installed. THIS PASS HAS NO GATE, and there is currently no way "
    "to steer it with custom constraints or loads -- the underlying "
    "function's own 'constraints' argument is not wired up and is not "
    'exposed here, so do not expect anything except the fixed set to hold '
    'its place. Snapshots first and reports all_improved; undo if that is '
    'false.',
    properties={
        'fixed': {
            'type': 'string', 'enum': ['corners', 'boundary', 'manual'],
            'description': "Which vertices are held. 'corners' (default) -- "
                           "just the boundary kinks. 'boundary' -- every "
                           "boundary vertex. 'manual' -- exactly "
                           "fixed_vertices."},
        'fixed_vertices': {
            'type': 'array', 'items': {'type': 'string'},
            'description': "Vertex handles to hold. Only used when "
                           "fixed='manual'."},
        'q_factor': {
            'type': 'number',
            'description': 'How much heavier a boundary edge is weighted '
                           'than an interior one. Default 100.'},
    },
    title='Relax (force density)')
def _t_relax_fdm(session, fixed='corners', fixed_vertices=None, q_factor=100.0):
    refusal = _needs_mesh(session)
    if refusal:
        return refusal
    if fixed not in ('corners', 'boundary', 'manual'):
        return {'ok': False,
                'reason': "fixed must be 'corners', 'boundary' or 'manual'; "
                          'got {!r}'.format(fixed)}
    resolved_fixed = []
    if fixed == 'manual':
        if not fixed_vertices:
            return {'ok': False,
                    'reason': "fixed='manual' needs a non-empty "
                              'fixed_vertices list'}
        missing = []
        for item in fixed_vertices:
            key, how = resolve(session.mesh, item)
            if key is None:
                missing.append(item)
            else:
                resolved_fixed.append(key)
        if missing:
            return {'ok': False,
                    'reason': 'could not place {} of {} handles: {}'.format(
                        len(missing), len(fixed_vertices),
                        ', '.join(missing[:3]))}

    before = session.quality()
    session.snapshot('before relax_fdm')
    try:
        relaxation(session.mesh, fixed=fixed, fixed_vertices=resolved_fixed,
                  q_factor=float(q_factor))
    except ImportError as exc:
        return {'ok': False,
                'reason': 'relax_fdm needs the compas_fd package, which is '
                          'not installed: {}'.format(exc)}
    except Exception as exc:
        return {'ok': False,
                'reason': 'relax_fdm failed: {}: {}'.format(
                    type(exc).__name__, exc)}
    return _outcome(session, 'relax_fdm', before, fixed=fixed,
                    q_factor=float(q_factor))
