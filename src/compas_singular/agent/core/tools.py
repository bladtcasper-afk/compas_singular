"""**The tool surface: every stage, as plain functions returning plain dicts.**

One function per thing a caller can do, each taking a
:class:`~compas_singular.agent.core.session.MeshEditSession` first and keyword
arguments after, each returning ``{'ok': bool, ...}``. No SDK, no CAD, no
prompts, no prints. :mod:`compas_singular.agent.runner` wraps these for a model;
a test drives the identical functions with a scripted plan and no key.

**The schemas live here, next to the functions they describe, not in the
runner.** A JSON schema is a dict; keeping it beside the implementation is what
stops the two drifting, and it costs no dependency.

**The descriptions carry the measurements.** That is not documentation
misplaced -- it is the working part. An automated caller handed ``relax`` with
no context will turn it on for a closed guide and destroy the layout (a
2.5-radius hoop in an 8x8 square: 8 singularities and 13 patches at 15.4
degrees, against 0 singularities and 1 patch at 20.0 with ``relax=True``). It
will copy ``guide_band=3.0`` out of ``12_cables``, which is the worst setting
measured. It will read ``target_length`` as the element size when at the field
stage it is the background triangulation spacing. Every one of those is cheap to
prevent in a sentence and expensive to discover in a loop.

**Four tools are marked** ``confirm``. They are the ones that destroy work:
``solve_field`` and ``set_guides`` discard the layout and everything under it,
``commit`` welds and renumbers and may drop drawn curves, and
``revert_to_coarse`` throws away every dense edit. A driver should put those in
front of a person, or at least in front of a budget.

**Nothing here calls** ``quad_mesh``. See ``rebuild.py``.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ...datastructures.mesh.smoothing import automated_boundary_constraints
from ...datastructures.mesh.smoothing import constrained_smoothing
from ...editing.guide_chain import attach_chain
from ...editing.guide_chain import chain_quality
from ...editing.guide_chain import guide_chain
from ...editing.guide_chain import mean_edge_length
from ...framefield import cache as _cache
from .address import strip_address
from .digest import digest as _digest
from .inputs import check_inputs as _check_inputs
from .rebuild import DensifyRefused
from .session import StageError


__all__ = ['ToolSpec', 'TOOLS', 'tool', 'call', 'tool_names', 'schemas']


class ToolSpec(object):
    """One callable tool: the function, what it is for, and its arguments."""

    def __init__(self, name, function, description, schema, stage=None,
                 confirm=False):
        self.name = name
        self.function = function
        self.description = description
        self.schema = schema
        #: Stages this tool is defined at. ``None`` means any.
        self.stage = tuple(stage) if stage else None
        #: Destroys work. A driver should confirm before calling.
        self.confirm = confirm

    def __repr__(self):
        return '<ToolSpec {}>'.format(self.name)


TOOLS = {}


def tool(name, description, properties=None, required=(), stage=None,
         confirm=False):
    """Register a tool. ``properties`` is the JSON-schema body of its arguments."""
    def wrap(function):
        schema = {
            'type': 'object',
            'properties': dict(properties or {}),
            'required': list(required),
            'additionalProperties': False,
        }
        TOOLS[name] = ToolSpec(name, function, description, schema, stage,
                               confirm)
        return function
    return wrap


def tool_names():
    return sorted(TOOLS)


def schemas():
    """``[{name, description, input_schema}, ...]`` -- ready for a tool API."""
    return [{'name': spec.name, 'description': spec.description,
             'input_schema': spec.schema}
            for spec in (TOOLS[n] for n in tool_names())]


def call(session, name, **kwargs):
    """Run a tool by name. Unknown names and wrong stages are refusals, not raises."""
    spec = TOOLS.get(name)
    if spec is None:
        return {'ok': False, 'reason': 'no tool named {!r}'.format(name),
                'available': tool_names()}
    if spec.stage and session.stage not in spec.stage:
        return {'ok': False,
                'reason': '{} is only defined at stage(s) {} -- the session is '
                          'at {!r}'.format(name, ', '.join(spec.stage),
                                           session.stage)}
    try:
        return spec.function(session, **kwargs)
    except TypeError as exc:
        # Python names the IMPLEMENTATION in the message (``_t_move_corner()``),
        # which is not a name the caller has ever seen or could call. Swap in
        # the tool's own name so the correction is obvious.
        detail = str(exc).replace(spec.function.__name__, name)
        return {'ok': False,
                'reason': 'bad arguments for {}: {}. Accepted: {}'.format(
                    name, detail,
                    ', '.join(sorted(spec.schema['properties'])) or 'none')}


def _fail(reason, **extra):
    out = {'ok': False, 'reason': reason}
    out.update(extra)
    return out


def _done(**detail):
    out = {'ok': True, 'reason': ''}
    out.update(detail)
    return out


# ======================================================================
# stage 0 -- the field
# ======================================================================

@tool('check_inputs',
      'Measure how finely the boundaries and guides are sampled, BEFORE solving. '
      'A cross field folds every angle into +/-45 degrees, so a curve that turns '
      'more than 45 degrees between consecutive points is read as a CORNER and '
      'cannot be read any other way. The measured case: an ellipse hole 4.2 by '
      '0.5 sampled at 48 points turns 57.7 degrees at each tip, loses a full '
      'period of boundary winding, and poincare_hopf still reports ok=True at '
      'every resolution because the field is faithfully solving the wrong '
      'problem. The only fix is to resample the curve finer -- no solver '
      'parameter changes it. Run this first; it is free.')
def _t_check_inputs(session):
    inputs = session.inputs()
    if not inputs['complete']:
        return _fail('the session has no outer boundary recorded, so the input '
                     'cannot be checked. Construct MeshEditSession with outer=.')
    result = _check_inputs(outer=inputs['outer'], inners=inputs['inners'],
                           guides=inputs['guides'])
    return _done(**result)


@tool('inspect_field',
      'The solved field: Poincare-Hopf, singularity count and sign, tracing '
      'report, symmetry. Read arm_mismatch first -- when a singularity does not '
      'produce exactly 4-index launch directions the layout comes out with the '
      'wrong number of patches and nothing about the coarse mesh says so. Note '
      'poincare_hopf.ok is NOT a quality verdict: it compares against the '
      'measured boundary winding, so it is True on a field that has faithfully '
      'solved an aliased input. Use check_inputs for that.')
def _t_inspect_field(session):
    data = _digest(session)
    return _done(field=data.get('field', {}), warnings=data.get('warnings', []))


@tool('solve_field',
      'Re-solve the field with different parameters and adopt the result. '
      'DESTROYS the coarse layout, every edit on it, the densities and the dense '
      'mesh -- a solve returns a NEW decomposition rather than mutating this one, '
      'so nothing below it survives. Call discarded_by first.\n'
      'relax: use True on a guided domain with OPEN guides -- it removes the '
      'magnitude ridge that pushes singularities off a hard guide, and is the '
      'only setting where real load paths beat walls alone. Do NOT use it with a '
      'CLOSED guide: a 2.5-radius hoop in an 8x8 square goes from 8 '
      'singularities / 13 patches / 15.4 degrees to 0 singularities / 1 patch / '
      '20.0 degrees, i.e. the layout collapses to the unguided one.\n'
      'guide_band: leave unset. The default (= target_length) measures 18.8 '
      'degrees against 26.7 for the 3.0 used in 12_cables, which is the worst '
      'setting tested.\n'
      'symmetry: leave as auto. On an asymmetric domain it is a no-op.\n'
      'target_length here is the BACKGROUND TRIANGULATION SPACING, not the quad '
      'size. Element size is set later with set_target_length or '
      'set_strip_density.',
      properties={
          'relax': {'type': 'boolean',
                    'description': 'Diffusion + normalisation solver. See above.'},
          'target_length': {'type': 'number',
                            'description': 'Background triangulation spacing.'},
          'guide_weight': {'type': 'number'},
          'guide_band': {'type': 'number',
                         'description': 'Leave unset unless you have a reason.'},
          'mode': {'type': 'string', 'enum': ['perpendicular', 'tangent']},
          'use_cache': {'type': 'boolean',
                        'description': 'Reuse an identical previous solve. '
                                       'Default true; a sweep should keep it on.'},
      },
      stage=('field', 'coarse', 'dense'), confirm=True)
def _t_solve_field(session, use_cache=True, **params):
    inputs = session.inputs()
    if not inputs['complete']:
        return _fail('no outer boundary recorded -- construct MeshEditSession '
                     'with outer= before re-solving.')

    kwargs = dict(inputs['params'])
    kwargs.update({k: v for k, v in params.items() if v is not None})
    cost = session.discarded_by('adopt_decomposition')

    try:
        decomposition = _cache.solve(
            inputs['outer'], inputs['inners'] or None,
            guides=inputs['guides'] or None, build=True,
            cache=None if use_cache else False, **kwargs)
    except Exception as exc:                                # noqa: BLE001
        return _fail('solve failed: {}: {}'.format(type(exc).__name__, exc))

    session.adopt_decomposition(decomposition, solve_params=kwargs)
    data = _digest(session)
    return _done(discarded=cost, params=kwargs, field=data.get('field', {}),
                 warnings=data.get('warnings', []))


@tool('set_guides',
      'Replace the guide curves -- cables, force lines -- and re-solve. Same '
      'destruction as solve_field: everything below the field goes. Each guide '
      'is a list of [x, y, z] points. Check the sampling first: a guide that '
      'turns more than 45 degrees between points is read as a corner.',
      properties={
          'guides': {
              'type': 'array',
              'description': 'One entry per curve; each a list of [x, y, z].',
              'items': {'type': 'array',
                        'items': {'type': 'array',
                                  'items': {'type': 'number'},
                                  'minItems': 2, 'maxItems': 3}},
          },
      },
      required=('guides',), stage=('field', 'coarse', 'dense'), confirm=True)
def _t_set_guides(session, guides):
    inputs = session.inputs()
    if not inputs['complete']:
        return _fail('no outer boundary recorded -- cannot re-solve.')
    curves = [[list(p) for p in curve] for curve in guides]
    check = _check_inputs(guides=curves)
    cost = session.discarded_by('adopt_decomposition')
    try:
        decomposition = _cache.solve(
            inputs['outer'], inputs['inners'] or None, guides=curves or None,
            build=True, **inputs['params'])
    except Exception as exc:                                # noqa: BLE001
        return _fail('solve failed: {}: {}'.format(type(exc).__name__, exc))
    session.adopt_decomposition(decomposition, guides=curves)
    return _done(discarded=cost, guides=len(curves), sampling=check)


# ======================================================================
# stage 1 -- the coarse layout
# ======================================================================

@tool('enter_coarse',
      'Build the coarse layout from the field, or return to the one already in '
      'hand. Idempotent -- an edit in progress is never discarded by this.',
      stage=('field', 'coarse'))
def _t_enter_coarse(session):
    editor = session.enter_coarse()
    return _done(faces=editor.mesh.number_of_faces(),
                 vertices=editor.mesh.number_of_vertices())


@tool('move_corner',
      'Move one corner of the coarse layout. A corner on the layout boundary is '
      'projected back onto the nearest domain wall, so nudging the edge does not '
      'eat the outline; an interior corner goes exactly where it is told. The '
      'edit is NOT finished until commit -- until then the moved edges are '
      'straight chords and the field alignment is not applied.',
      properties={
          'label': {'type': 'string',
                    'description': 'A vertex label, e.g. "v7", from the digest.'},
          'x': {'type': 'number'}, 'y': {'type': 'number'},
          'project': {'type': 'boolean',
                      'description': 'Hold a boundary corner on its wall. '
                                     'Default true; turn off only deliberately.'},
      },
      required=('label', 'x', 'y'), stage=('coarse',))
def _t_move_corner(session, label, x, y, project=True):
    book = session.book()
    vkey = book.vertex(label) if book else None
    if vkey is None:
        return _fail('no vertex labelled {!r} on this layout'.format(label))
    z = session.mesh.vertex_coordinates(vkey)[2]
    ok, _notes = session.coarse_editor.move_vertex(vkey, [x, y, z], project=project)
    if not ok:
        return _fail(session.coarse_editor.last_reason)
    session.record('move_corner', label=label)
    moved = session.mesh.vertex_coordinates(vkey)
    session.refresh_book()
    return _done(label=label, moved_to=moved)


@tool('divide',
      'Cut the layout with a drawn line, arc or polyline: every coarse edge it '
      'crosses is split and every patch it passes through is split in two. The '
      'rule the cut is checked against is that every patch of the RESULT has '
      'four sides. Two things follow: a patch must be entered and left through '
      'OPPOSITE sides (through adjacent ones the split makes a triangle, and '
      'that is refused outright), and the cut must reach a wall at both ends. A '
      'cut that stops in the middle is not wrong, it is short -- set extend=true '
      'and it is continued along the strip to the wall. The drawn shape is kept: '
      'an arc stays an arc in the final mesh. Planned on a copy; a refusal costs '
      'nothing.',
      properties={
          'points': {'type': 'array',
                     'description': 'The curve, as [x, y, z] points.',
                     'items': {'type': 'array', 'items': {'type': 'number'},
                               'minItems': 2, 'maxItems': 3}},
          'extend': {'type': 'boolean',
                     'description': 'Continue a short cut along the strip to the wall.'},
      },
      required=('points',), stage=('coarse',))
def _t_divide(session, points, extend=True):
    pts = [list(p) + [0.0] * (3 - len(p)) for p in points]
    ok, _notes = session.coarse_editor.divide(pts, extend=extend)
    if not ok:
        return _fail(session.coarse_editor.last_reason)
    session.record('divide', points=len(pts), extend=bool(extend))
    carry = session.refresh_book()
    return _done(faces=session.mesh.number_of_faces(),
                 cut=dict(session.coarse_editor.last_cut or {}),
                 relabelled={'kept': len(carry.get('kept', {})),
                             'lost': len(carry.get('lost', []))})


@tool('plan_strip_deletion',
      'What deleting the strip through an edge would cost, WITHOUT doing it. A '
      'deletion is never local: it welds the two sides together so surrounding '
      'corners move, it can take collateral strips with it, and it can collapse '
      'a boundary. No cheap test predicts every failure -- on one measured '
      'layout 10 of 20 strips left a non-manifold result and poles were not the '
      'discriminator -- so this performs the deletion on a copy and reports what '
      'came out. Always call this before remove_strip.',
      properties={'label': {'type': 'string',
                            'description': 'An edge label, e.g. "e12".'}},
      required=('label',), stage=('coarse',))
def _t_plan_strip_deletion(session, label):
    book = session.book()
    edge = book.edge(label) if book else None
    if edge is None:
        return _fail('no edge labelled {!r} on this layout'.format(label))
    plan = session.coarse_editor.plan_strip_deletion(edge)
    return _done(label=label, plan=plan)


@tool('remove_strip',
      'Delete the strip through an edge. A quad layout cannot lose a single '
      'edge -- removing one merges two patches into a hexagon -- so the unit of '
      'deletion is the whole strip, wall to wall. This COLLAPSES a band and welds '
      'its two sides; it does not dissolve the line and merge the patches either '
      'side, and no such operation exists. Call plan_strip_deletion first.',
      properties={'label': {'type': 'string'}},
      required=('label',), stage=('coarse',))
def _t_remove_strip(session, label):
    book = session.book()
    edge = book.edge(label) if book else None
    if edge is None:
        return _fail('no edge labelled {!r} on this layout'.format(label))
    ok, _notes = session.coarse_editor.remove_strip(edge)
    if not ok:
        return _fail(session.coarse_editor.last_reason)
    session.record('remove_strip', label=label)
    carry = session.refresh_book()
    return _done(faces=session.mesh.number_of_faces(),
                 deletion=dict(session.coarse_editor.last_deletion or {}),
                 relabelled={'kept': len(carry.get('kept', {})),
                             'lost': len(carry.get('lost', []))})


@tool('reset_coarse',
      'Undo every uncommitted edit on the layout, back to the last commit. Does '
      'not undo a commit -- nothing does.',
      stage=('coarse',))
def _t_reset_coarse(session):
    session.coarse_editor.reset()
    session.record('reset_coarse')
    session.refresh_book()
    return _done(faces=session.mesh.number_of_faces())


@tool('commit',
      'Hand the edited layout back to the decomposition. THIS is what makes an '
      'edit worth anything: the moved edges are re-matched to the separatrices '
      'they were traced from and those curves are warped onto the new corners, '
      'so the layout keeps its field alignment. Densifying without committing '
      'gives straight chords between the moved corners. It also RENUMBERS every '
      'key -- 2 of 8 vertex keys survived on a measured disc -- so labels are '
      'reissued afterwards and anything held across it must be re-read. A bad '
      'layout is refused rather than silently replaced.',
      stage=('coarse',), confirm=True)
def _t_commit(session):
    layout, notes = session.coarse_editor.commit()
    ok = layout is not None
    if not ok:
        return _fail(notes.get('error', session.coarse_editor.last_reason),
                     notes=notes)
    session.record('commit', **{k: v for k, v in notes.items()
                                if isinstance(v, (int, float, str))})
    carry = session.refresh_book()
    return _done(notes=notes,
                 relabelled={'kept': len(carry.get('kept', {})),
                             'moved': len(carry.get('moved', [])),
                             'lost': len(carry.get('lost', []))})


# ======================================================================
# stages 2 and 3 -- densify, and the densities
# ======================================================================

@tool('rebuild_dense',
      'Densify the layout and move to the dense stage. Per-strip densities set '
      'earlier are re-applied automatically. If the result is not usable this '
      'REFUSES and leaves the layout exactly as it was -- unlike quad_mesh, '
      'which falls back to a triangulation and in doing so replaces the layout '
      'and clears the edited flag, losing the whole edit.',
      properties={
          'target_length': {'type': 'number',
                            'description': 'Base element size. This one IS the '
                                           'quad size.'},
          'density': {'type': 'integer',
                      'description': 'Uniform subdivisions per coarse edge. '
                                     'Overrides target_length.'},
      },
      stage=('field', 'coarse', 'dense'))
def _t_rebuild_dense(session, target_length=None, density=None):
    if session.stage == 'dense':
        session.revert_to_coarse()
    try:
        _dense, report = session.rebuild_dense(target_length=target_length,
                                               density=density)
    except DensifyRefused as exc:
        return _fail(exc.reason, stage=exc.stage, metrics=exc.metrics,
                     kept='the layout is untouched -- nothing was replaced')
    except Exception as exc:                                # noqa: BLE001
        return _fail('{}: {}'.format(type(exc).__name__, exc))
    session.refresh_book()
    data = _digest(session)
    return _done(faces=report['faces'], route=report['route'],
                 unmatched=len(report['density_report']['unmatched']),
                 notes=report['notes'], dense=data.get('dense', {}))


@tool('set_field_aware',
      'Whether patch INTERIORS are integrated from the field (on, the default) '
      'or filled with a Coons blend that never consults it (off). This affects '
      'interiors ONLY -- the coarse edges follow separatrices either way, '
      'through edges_to_curves -- so turning it off is a smaller change than the '
      'name suggests. Takes effect at the next rebuild_dense.',
      properties={'enabled': {'type': 'boolean'}},
      required=('enabled',))
def _t_set_field_aware(session, enabled):
    session.decomposition.field_aware = bool(enabled)
    session.record('set_field_aware', enabled=bool(enabled))
    return _done(field_aware=bool(enabled),
                 note='takes effect at the next rebuild_dense')


@tool('list_strips',
      'Every strip of the coarse layout with its current density, its address '
      'and a representative point. The address is geometric, because strip keys '
      'renumber on every edit and a remembered one silently names a different '
      'strip.',
      stage=('coarse', 'dense'))
def _t_list_strips(session):
    live = (session.coarse_editor.mesh if session.coarse_editor
            else session.decomposition.mesh)
    if live is None:
        return _fail('no coarse layout')

    # On a COPY. ``collect_strips`` does not populate ``strips_density``, so
    # reading a density off a layout that has not been through a base pass
    # raises ``KeyError`` from inside ``get_strip_density`` -- the failure
    # ``edit_coarse``'s docstring warns about, several frames from its cause.
    # Applying the base pass here would also mutate a layout the caller is still
    # editing, so it happens on a copy and only the numbers come back.
    work = live.copy()
    work.collect_strips()
    try:
        work.set_strips_density_target(session.target_length)
    except Exception as exc:                                # noqa: BLE001
        return _fail('could not derive densities at target_length {}: {}: '
                     '{}'.format(session.target_length, type(exc).__name__, exc))

    out = []
    for skey in work.strips():
        address = strip_address(work, skey)
        if address is None:
            continue
        derived = work.get_strip_density(skey)
        override = session.densities.get(address)
        out.append({
            'address': list(address),
            'density': derived if override is None else override,
            'from': 'target_length' if override is None else 'override',
            'edges': len(list(work.strip_edges(skey))),
            'faces': len(work.strip_faces(skey)),
        })
    return _done(strips=out, count=len(out),
                 target_length=session.target_length,
                 note='density is what the next rebuild_dense would use: derived '
                      'from target_length unless an override was set')


@tool('set_strip_density',
      'Set the subdivision count of ONE strip, by its address from list_strips. '
      'Takes effect at the next rebuild_dense, and survives it -- the address '
      'is geometric, so it still names the right strip after a layout edit.',
      properties={
          'address': {'type': 'array', 'items': {'type': 'number'},
                      'minItems': 3, 'maxItems': 3,
                      'description': 'From list_strips.'},
          'density': {'type': 'integer', 'minimum': 1},
      },
      required=('address', 'density'), stage=('coarse', 'dense'))
def _t_set_strip_density(session, address, density):
    if int(density) < 1:
        return _fail('a density below 1 leaves a strip with no faces and the '
                     'densification refuses it')
    key = tuple(round(float(v), 3) + 0.0 for v in address)
    session.densities[key] = int(density)
    session.record('set_strip_density', address=key, density=int(density))
    return _done(address=list(key), density=int(density),
                 pending=len(session.densities),
                 note='takes effect at the next rebuild_dense')


@tool('set_target_length',
      'The base element size used where no per-strip density overrides it. This '
      'is the QUAD SIZE, not the background triangulation spacing that '
      'solve_field takes.',
      properties={'target_length': {'type': 'number', 'exclusiveMinimum': 0}},
      required=('target_length',))
def _t_set_target_length(session, target_length):
    if target_length <= 0:
        return _fail('target_length must be positive')
    session.target_length = float(target_length)
    session.record('set_target_length', target_length=float(target_length))
    return _done(target_length=float(target_length),
                 note='takes effect at the next rebuild_dense')


@tool('set_face_target',
      'Set every strip density so the dense mesh lands near a total face count. '
      'A convenience over set_target_length when the caller is thinking in '
      'elements rather than in size. Clears per-strip overrides.',
      properties={'faces': {'type': 'integer', 'minimum': 1}},
      required=('faces',), stage=('coarse', 'dense'))
def _t_set_face_target(session, faces):
    coarse = (session.coarse_editor.mesh if session.coarse_editor
              else session.decomposition.mesh)
    if coarse is None:
        return _fail('no coarse layout')
    coarse.collect_strips()
    try:
        coarse.set_mesh_density_face_target(int(faces))
    except Exception as exc:                                # noqa: BLE001
        return _fail('{}: {}'.format(type(exc).__name__, exc))
    session.densities = {}
    for skey in coarse.strips():
        address = strip_address(coarse, skey)
        if address is not None:
            session.densities[address] = coarse.get_strip_density(skey)
    session.record('set_face_target', faces=int(faces))
    return _done(requested=int(faces), strips=len(session.densities),
                 note='takes effect at the next rebuild_dense')


# ======================================================================
# stage 4 -- global smoothing
# ======================================================================

@tool('smooth',
      'Relax the dense mesh globally, with the outline and the patch seams '
      'sliding on their own curves and the singularities free. The safest tool '
      'here: it keeps a result ONLY if the minimum angle, the maximum angle and '
      'the aspect ratio all improve together, so a mesh with nothing to gain '
      'comes back untouched and accepted is null. seams="free" is the only '
      'setting measured to help; "slide" never hurts and never helps. Leave '
      'pin_singularities off -- pinning them is what stops this pass working.',
      properties={
          'seams': {'type': 'string', 'enum': ['free', 'slide', 'fixed']},
          'pin_singularities': {'type': 'boolean'},
      },
      stage=('dense',))
def _t_smooth(session, seams='free', pin_singularities=False):
    mesh = session.mesh
    try:
        report = session.decomposition.smooth_quad_mesh(
            mesh, seams=seams, pin_singularities=bool(pin_singularities))
    except Exception as exc:                                # noqa: BLE001
        return _fail('{}: {}'.format(type(exc).__name__, exc))
    session.record('smooth', seams=seams, accepted=str(report.get('accepted')))
    session.refresh_book()
    return _done(accepted=report.get('accepted'), before=report.get('before'),
                 after=report.get('after'), sliding=report.get('sliding'),
                 pinned=report.get('pinned'),
                 note=('no setting improved all three metrics together; the mesh '
                       'is unchanged' if report.get('accepted') is None else ''))


# ======================================================================
# stage 5 -- smoothing to a guide
# ======================================================================

@tool('select_guide_chain',
      'Choose the run of mesh vertices that follows a guide curve, WITHOUT '
      'moving anything. Two gates, doing different jobs: distance asks whether a '
      'vertex is near the guide, angle asks whether its line is still going the '
      'same way -- a polyedge that runs beside the guide then veers off keeps '
      'passing the distance test for a vertex or two after it has stopped '
      'following, and those are the ones whose projection bends the mesh. Always '
      'measure the result before attaching.',
      properties={
          'guide': {'type': 'array',
                    'description': 'The curve, as [x, y, z] points.',
                    'items': {'type': 'array', 'items': {'type': 'number'},
                              'minItems': 2, 'maxItems': 3}},
          'tolerance': {'type': 'number',
                        'description': 'Distance gate. Default is a multiple of '
                                       'the mean edge length.'},
          'max_angle': {'type': 'number',
                        'description': 'Angle gate in degrees. 90 turns it off.'},
          'boundary': {'type': 'string', 'enum': ['anchor', 'exclude', 'free'],
                       'description': '"anchor" is the working setting. "free" '
                                      'lets the outline itself be selected and '
                                      'exists only to measure what that costs.'},
      },
      required=('guide',), stage=('dense',))
def _t_select_guide_chain(session, guide, tolerance=None, max_angle=None,
                          boundary='anchor'):
    pts = [list(p) + [0.0] * (3 - len(p)) for p in guide]
    kwargs = {'boundary': boundary}
    if tolerance is not None:
        kwargs['tolerance'] = tolerance
    if max_angle is not None:
        kwargs['max_angle'] = max_angle
    try:
        chain, info = guide_chain(session.mesh, pts, **kwargs)
    except Exception as exc:                                # noqa: BLE001
        return _fail('{}: {}'.format(type(exc).__name__, exc))
    if not chain:
        return _fail(info.get('reason', 'no chain follows that guide'), info=info)
    quality = chain_quality(session.mesh, chain, pts)
    session._pending_chain = (chain, pts)
    return _done(chain=len(chain), info=info, quality=quality,
                 mean_edge=mean_edge_length(session.mesh),
                 gate='boundary_interior must be 0 -- a chain running along the '
                      'outline takes the outline with it when projected')


@tool('attach_guide_chain',
      'Move the selected chain onto its guide and smooth the mesh with it held '
      'there. A BOUNDARY vertex is never moved onto the guide -- at most it '
      'slides along the outline, because moving it would take the outline of the '
      'building with it. hold="fixed" pins an interior vertex where it lands; '
      '"sliding" re-projects it every iteration so it may travel along the guide '
      'but never leave it. Run select_guide_chain first.',
      properties={
          'hold': {'type': 'string', 'enum': ['fixed', 'sliding']},
          'boundary': {'type': 'string', 'enum': ['sliding', 'fixed'],
                       'description': 'What the rest of the outline may do.'},
          'kmax': {'type': 'integer', 'minimum': 1},
          'damping': {'type': 'number', 'minimum': 0.0, 'maximum': 1.0},
      },
      stage=('dense',))
def _t_attach_guide_chain(session, hold='fixed', boundary='sliding', kmax=100,
                          damping=0.5):
    pending = getattr(session, '_pending_chain', None)
    if not pending:
        return _fail('nothing selected -- call select_guide_chain first')
    chain, guide = pending
    mesh = session.mesh

    moves, attached = attach_chain(mesh, chain, guide, hold=hold)
    for vertex, xyz in moves.items():
        mesh.vertex_attributes(vertex, 'xyz', xyz)

    if boundary == 'sliding':
        constraints = automated_boundary_constraints(mesh)
        fixed = None
    else:
        constraints = {}
        fixed = [v for loop in mesh.vertices_on_boundaries() for v in loop]
    # An attached BOUNDARY vertex is already held on its own outline by
    # attach_chain -- the same constraint the sliding treatment would give it --
    # so this never takes a vertex off a wall.
    constraints.update(attached)
    if fixed:
        fixed = [v for v in fixed if v not in attached]

    try:
        constrained_smoothing(mesh, kmax=int(kmax), damping=float(damping),
                              constraints=constraints, algorithm='area',
                              fixed=fixed)
    except Exception as exc:                                # noqa: BLE001
        return _fail('{}: {}'.format(type(exc).__name__, exc))

    quality = chain_quality(mesh, chain, guide)
    session._pending_chain = None
    session.record('attach_guide_chain', hold=hold, moved=len(moves))
    session.refresh_book()
    return _done(attached=len(attached), moved=len(moves),
                 slid=len(chain) - len(moves), quality=quality)


# ======================================================================
# stage 6 -- the dense mesh
# ======================================================================

@tool('plan_dense_strip_deletion',
      'What removing the line through an edge of the DENSE mesh would cost, '
      'without doing it. Same three consequences as on the layout -- welded '
      'sides, collateral strips, collapsed boundaries -- performed on a copy.',
      properties={'label': {'type': 'string'}},
      required=('label',), stage=('dense',))
def _t_plan_dense(session, label):
    book = session.book()
    edge = book.edge(label) if book else None
    if edge is None:
        return _fail('no edge labelled {!r} on this mesh'.format(label))
    return _done(label=label,
                 plan=session.dense_editor.plan_strip_deletion(edge))


@tool('remove_line',
      'Delete the strip through an edge of the dense mesh, welding the two sides '
      'together. Refused on a mesh with any non-quad face: collect_strip walks '
      'face_opposite_edge and a pole fan has none. Call '
      'plan_dense_strip_deletion first. Note dense edits cannot be replayed onto '
      'a regenerated layout -- rebuild_dense discards them.',
      properties={'label': {'type': 'string'},
                  'preserve_boundaries': {
                      'type': 'boolean',
                      'description': 'Pre-split the strips that would otherwise '
                                     'let a boundary collapse.'}},
      required=('label',), stage=('dense',))
def _t_remove_line(session, label, preserve_boundaries=False):
    book = session.book()
    edge = book.edge(label) if book else None
    if edge is None:
        return _fail('no edge labelled {!r} on this mesh'.format(label))
    ok = session.dense_editor.remove_line(
        edge, preserve_boundaries=bool(preserve_boundaries))
    if not ok:
        return _fail(session.dense_editor.last_reason)
    session.record('remove_line', label=label)
    carry = session.refresh_book()
    return _done(deletion=dict(session.dense_editor.last_deletion or {}),
                 faces=session.mesh.number_of_faces(),
                 relabelled={'lost': len(carry.get('lost', []))})


@tool('add_line',
      'Grow a strip along the polyedge through an edge of the dense mesh: one '
      'new row of elements, wall to wall or all the way round. A strip cannot '
      'run part-way -- anything less leaves a five-sided face -- so a polyedge '
      'that neither closes nor ends on the boundary at both ends is refused. The '
      'new vertices are created on top of the ones they replace, so the mesh is '
      'relaxed to open the strip.',
      properties={'label': {'type': 'string'},
                  'relax': {'type': 'boolean'}},
      required=('label',), stage=('dense',))
def _t_add_line(session, label, relax=True):
    book = session.book()
    edge = book.edge(label) if book else None
    if edge is None:
        return _fail('no edge labelled {!r} on this mesh'.format(label))
    ok = session.dense_editor.add_line(edge, relax=bool(relax))
    if not ok:
        return _fail(session.dense_editor.last_reason)
    session.record('add_line', label=label)
    carry = session.refresh_book()
    return _done(addition=dict(session.dense_editor.last_addition or {}),
                 faces=session.mesh.number_of_faces(),
                 relabelled={'lost': len(carry.get('lost', []))})


@tool('move_vertex',
      'Move one vertex of the dense mesh. A boundary vertex is projected back '
      'onto the domain wall. Topology is unchanged, so this is the one operation '
      'that still works on a mesh with poles.',
      properties={'label': {'type': 'string'},
                  'x': {'type': 'number'}, 'y': {'type': 'number'},
                  'project': {'type': 'boolean'}},
      required=('label', 'x', 'y'), stage=('dense',))
def _t_move_vertex(session, label, x, y, project=True):
    book = session.book()
    vkey = book.vertex(label) if book else None
    if vkey is None:
        return _fail('no vertex labelled {!r} on this mesh'.format(label))
    z = session.mesh.vertex_coordinates(vkey)[2]
    ok = session.dense_editor.move_vertex(vkey, [x, y, z], project=project)
    if not ok:
        return _fail(session.dense_editor.last_reason)
    session.record('move_vertex', label=label)
    session.refresh_book()
    return _done(label=label, moved_to=session.mesh.vertex_coordinates(vkey))


# ======================================================================
# session control
# ======================================================================

@tool('inspect',
      'The full state of the session: the stage, the field, the layout or the '
      'dense mesh with its element quality, the recent repair notes and the '
      'warnings. Free, and safe to call at any point.',
      properties={'max_labels': {'type': 'integer', 'minimum': 0}})
def _t_inspect(session, max_labels=None):
    kwargs = {} if max_labels is None else {'max_labels': max_labels}
    return _done(**_digest(session, **kwargs))


@tool('discarded_by',
      'What a backward move would destroy, before it is made. The two backward '
      'moves are revert_to_coarse (drops the dense mesh and every dense edit) '
      'and adopt_decomposition (a re-solve, which drops the layout, its edits, '
      'the densities and the dense mesh).',
      properties={'transition': {'type': 'string',
                                 'enum': ['revert_to_coarse',
                                          'adopt_decomposition']}},
      required=('transition',))
def _t_discarded_by(session, transition):
    return _done(**session.discarded_by(transition))


@tool('revert_to_coarse',
      'Go back from the dense mesh to the layout. DESTROYS every dense edit; a '
      'hand-edited dense mesh has no layout that reproduces it. Call '
      'discarded_by first.',
      stage=('dense',), confirm=True)
def _t_revert_to_coarse(session):
    # ``revert_to_coarse`` RAISES StageError when there is no layout to go back
    # to -- an adopted dense mesh has none. ``call`` only converts TypeError, so
    # without this the exception escapes into the loop and ends the run over
    # what is merely a refusal.
    try:
        return _done(discarded=session.revert_to_coarse())
    except StageError as exc:
        return _fail(str(exc))


@tool('snapshot',
      'Push the whole session onto the undo stack so undo can come back here. '
      'Not free -- it copies the field -- so take one before something risky '
      'rather than every step.',
      properties={'label': {'type': 'string'}})
def _t_snapshot(session, label=''):
    depth = session.snapshot(label)
    return _done(depth=depth, label=label)


@tool('undo',
      'Restore the last snapshot. Only reaches back to a snapshot that was '
      'explicitly taken.')
def _t_undo(session):
    ok, label = session.undo()
    if not ok:
        return _fail(label)
    session.refresh_book()
    return _done(restored=label, stage=session.stage)
