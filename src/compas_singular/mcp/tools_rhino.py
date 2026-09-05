"""**Pulling geometry out of a live Rhino document, and pushing a mesh back.**

The server is standalone: everything else in this package works with Rhino
closed, on a mesh read from a file. These three tools are the exception, and all
three go through :mod:`~compas_singular.mcp.bridge.spool` -- a request written to
disk, taken by ``CMD_mcp_link`` on Rhino's main thread, answered on disk. No
socket, no thread, and nothing here imports a Rhino module.

**A pull is free; a push is not.** ``rhino_pull`` and ``rhino_status`` read and
change nothing, so they are marked read-only. ``rhino_push`` writes into a
document a person has open, which is why it is marked destructive and open-world
and why the link keeps the previous mesh on a ``::MCP::Before`` layer.

**Never round-trip a mesh mid-session.** ``bake_mesh`` goes through Rhino's
``Point3f``, so a push followed by a pull comes back single precision -- about
2e-6 on a 20-unit plate, measured, and recorded at the top of
``rhino/coarse_curves.py``. The server holds the authoritative double-precision
mesh and Rhino holds a copy for display. Pull once, improve as many times as it
takes, push at the end.

**A timeout is a refusal with a cause, not a hang.** Rhino may be inside a
command, in which case the link defers rather than touching the document
underneath the person using it. The wait returns and says which of the two it
was -- no link attached, or attached and busy -- because the fix is different.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .bridge import spool
from .bridge import wire
from .describe import describe
from .library import thresholds
from .registry import tool


__all__ = []

#: The layer a dense quad mesh lives on, in the layer scheme every ``CMD_``
#: command already uses.
DEFAULT_LAYER = 'QuadMesh'


def _ask(verb, args=None, timeout=30.0):
    """One request to Rhino, and its answer. Never raises, never hangs."""
    try:
        request_id = spool.post(verb, args or {})
    except OSError as exc:
        return {'ok': False,
                'reason': 'could not write to the spool directory {}: '
                          '{}'.format(spool.default_directory(), exc)}
    return spool.wait(request_id, timeout=float(timeout))


def unseen_refusal(session):
    """A refusal if the current mesh has not been rendered since it changed.

    The quality numbers cannot see a mesh that has left its boundary, ignored a
    guide, or dropped a point feature, and they cannot see a polyedge that
    kinks. Every one of those is obvious in a picture and none of them is in
    ``mesh_quality``. So a result is not delivered until somebody has looked at
    it: one extra call, against shipping a mesh that measures well and is
    visibly wrong.

    Shared by ``rhino_push`` and ``save_mesh`` so the rule, and its wording, is
    one thing rather than two that drift.
    """
    if session.visually_current:
        return None
    return {'ok': False,
            'reason': 'this mesh has changed since it was last drawn, so '
                      'nothing has actually looked at what is about to be '
                      'delivered. Call inspect with image=true first, judge the '
                      'boundaries, the point features, the guides and the '
                      'continuity of the polyedges, and then deliver.',
            'fix': {'tool': 'inspect', 'arguments': {'image': True}}}


@tool(
    'rhino_status',
    'Whether Rhino is listening, which document is open, and how many requests '
    'are queued. Call this FIRST if a pull or push has failed: it separates '
    '"no link attached" -- the person must run CMD_mcp_link in Rhino -- from '
    '"attached but busy", which means Rhino is inside a command and the link is '
    'deferring rather than touching the document underneath them. Reads '
    'nothing from the document and changes nothing.',
    read_only=True, idempotent=True, open_world=True, title='Rhino link status')
def _t_rhino_status(session, timeout=5.0):
    state = spool.link_state()
    queued = len(spool.pending())
    if state is None:
        return {'ok': True, 'attached': False, 'queued': queued,
                'spool': spool.default_directory(),
                'reading': 'No link is attached. Run CMD_mcp_link in Rhino to '
                           'attach one. Requests posted meanwhile will queue '
                           'and be answered when it attaches.'}
    # A live link answers a ping; a stale heartbeat means Rhino went away
    # without detaching, so say that rather than reporting it as attached.
    reply = _ask('ping', timeout=timeout)
    alive = bool(reply.get('ok'))
    return {'ok': True, 'attached': alive, 'queued': queued,
            'document': state.get('document'),
            'seconds_since_heartbeat': round(state.get('age', 0.0), 1),
            'stale': state.get('stale'),
            'spool': spool.default_directory(),
            'reading': ('Attached to {} and answering.'.format(
                state.get('document') or 'an unnamed document')
                if alive else
                'A link is registered but did not answer a ping. Rhino is '
                'either inside a command or has been closed without '
                'detaching.')}


@tool(
    'rhino_pull',
    'Take the mesh and the boundary curves out of the open Rhino document into '
    'this session, replacing whatever was held. Reads the outer and inner '
    'boundaries, the guide curves, the point features, and the mesh on the '
    'named layer. The '
    'boundaries matter: they are the walls that boundary vertices SLIDE along '
    'during smoothing, and without them a faceted boundary can only get '
    'smoother, never truer. Reads the document and changes nothing in it. '
    'Discards the current undo history, since positions from one mesh cannot be '
    'applied to another.',
    properties={
        'layer': {
            'type': 'string',
            'description': 'Layer holding the mesh. Default "QuadMesh" (the '
                           'dense mesh). Use "Mesh" for the coarse layout.'},
        'spacing': {
            'type': 'number',
            'description': 'How finely a CURVED boundary is sampled into a '
                           'polyline. Default 0.125. Sample finer than the '
                           'element size: these points become the wall, so the '
                           'boundary is only as smooth as what is handed in.'},
        'timeout': {'type': 'number',
                    'description': 'Seconds to wait for Rhino. Default 30.'},
    },
    read_only=True, open_world=True, title='Pull from Rhino')
def _t_rhino_pull(session, layer=DEFAULT_LAYER, spacing=0.125, timeout=30.0):
    reply = _ask('pull', {'layer': layer, 'spacing': float(spacing)},
                 timeout=timeout)
    if not reply.get('ok'):
        return {'ok': False, 'reason': reply.get('reason', 'the pull failed'),
                'link': reply.get('link'), 'timed_out': reply.get('timed_out')}

    result = reply.get('result') or {}
    payload = result.get('mesh')
    if not payload:
        return {'ok': False,
                'reason': 'Rhino answered but there is no mesh on layer '
                          '{!r}.'.format(layer)}
    try:
        mesh = wire.mesh_from_wire(payload)
    except wire.WireError as exc:
        return {'ok': False,
                'reason': 'the mesh Rhino sent could not be rebuilt: {}'.format(exc)}

    outer = result.get('outer') or []
    inners = result.get('inners') or []
    walls = ([outer] if outer else []) + list(inners)
    session.adopt(mesh, walls=walls, guides=result.get('guides') or [],
                  points=result.get('points') or [],
                  source={'kind': 'rhino', 'layer': layer,
                          'document': result.get('document')})

    metrics = session.quality()
    told = describe(session.mesh, metrics, thresholds())
    session.record('rhino_pull', before=None, after=metrics, layer=layer,
                   walls=len(walls))
    out = {'ok': True, 'layer': layer,
           'document': result.get('document'),
           'vertices': mesh.number_of_vertices(),
           'faces': mesh.number_of_faces(),
           'walls': len(walls), 'guides': len(session.guides),
           'point_features': len(session.points),
           'quality': metrics, 'reading': told['reading'],
           'verdict': told['verdict'], 'worst_face_at': told['worst_face_at']}
    if not walls:
        out['warning'] = ('no boundary curves came back, so the boundary '
                          'vertices have nothing true to slide along. Smoothing '
                          'can still run, but it will only make the boundary '
                          'smoother, not truer.')
    return out


@tool(
    'rhino_push',
    'Bake the improved mesh back into the open Rhino document, on the given '
    'layer. WRITES TO A DOCUMENT SOMEBODY HAS OPEN: the mesh that was there is '
    'moved to a "::MCP::Before" sub-layer first so the two can be compared, and '
    'the write is one named undo step so it can be taken back in Rhino. It will '
    'not fire while Rhino is inside a command -- it waits. Push only when the '
    'mesh is actually better than it was pulled; if it is not, say so and push '
    'nothing. Do not pull again afterwards to check: Rhino stores mesh vertices '
    'in single precision, so a round trip loses about 2e-6 per coordinate and '
    'the copy here is the accurate one. REFUSES until the mesh being pushed '
    'has been drawn and looked at -- call inspect with image=true first.',
    properties={
        'layer': {'type': 'string',
                  'description': 'Layer to bake onto. Default "QuadMesh".'},
        'timeout': {'type': 'number',
                    'description': 'Seconds to wait for Rhino. Default 30.'},
    },
    destructive=True, open_world=True, title='Push to Rhino')
def _t_rhino_push(session, layer=DEFAULT_LAYER, timeout=30.0):
    if not session.loaded:
        return {'ok': False, 'reason': 'no mesh is loaded, so there is nothing '
                                       'to push'}
    refusal = unseen_refusal(session)
    if refusal:
        return refusal
    payload = wire.mesh_to_wire(session.mesh)
    reply = _ask('push', {'layer': layer, 'mesh': payload}, timeout=timeout)
    if not reply.get('ok'):
        return {'ok': False, 'reason': reply.get('reason', 'the push failed'),
                'link': reply.get('link'), 'timed_out': reply.get('timed_out')}
    result = reply.get('result') or {}
    session.record('rhino_push', before=None, after=session.quality(),
                   layer=layer)
    return {'ok': True, 'layer': layer,
            'document': result.get('document'),
            'vertices': len(payload['vertices']),
            'faces': len(payload['faces']),
            'kept_previous_on': result.get('before_layer'),
            'undo_record': result.get('undo_record'),
            'reading': 'Baked {} faces onto {}. The mesh that was there is on '
                       '{}.'.format(len(payload['faces']), layer,
                                    result.get('before_layer') or 'no backup layer')}
