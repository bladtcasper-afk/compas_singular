"""**What the session looks like right now, in numbers and in words.**

The text half of the observation an automated caller works from; ``render.py``
is the other half. Two rules shape it.

**The useful facts are different at every stage**, so this is not one summary
with fields that happen to be empty. At the field stage the question is whether
the field is a faithful reading of the input -- Poincare-Hopf, the singularity
count, and above all ``arm_mismatch``. At the coarse stage it is the patch
count and whether anything is not a quad. At the dense stage it is element
quality, and the route that produced it.

**Two numbers here look like verdicts and are not:**

* ``poincare_hopf.ok`` is ``True`` on a field that has faithfully solved the
  WRONG problem -- an aliased boundary arc gives ``ok=True``, residual 0, at
  every resolution and under both solvers, because it compares against the
  MEASURED boundary winding rather than against ``4*chi``. It says the unwrapping
  is sound, nothing more. ``inputs.check_inputs`` is what catches the other
  thing, and it has to be run on the input, before the solve.
* ``route`` is structural, never quality. A mesh on route ``'field'`` can carry
  a 180-degree angle; ``12_cables``' ring cable did. So ``hard_floor`` is
  reported next to it, always, and the two are separate lines because they
  answer separate questions.

**And one that is easy to over-read:** ``densify.guarded`` counts patches that
rejected the field-integrated interior and kept the Coons one. On an UNGUIDED
domain that is the ordinary outcome and means nothing -- there was no alignment
to buy. Only where a guide was supposed to reach that patch is it a complaint.
The digest says which case it is rather than leaving the caller to guess.

**Labels are only listed when there are few enough to be useful.** A dense mesh
has thousands of vertices; dumping their addresses is a large payload that no
caller can act on. Past ``max_labels`` the digest says how many there are and
stops, and a caller that needs a specific one asks for it by position.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ...framefield.quality import hard_floor


__all__ = ['digest', 'render_text']


#: Above this many, labels are counted rather than listed.
MAX_LABELS = 40

#: How many of the most recent repair notes to carry. They accumulate over a
#: session and the old ones describe a layout that no longer exists.
NOTE_TAIL = 6


def _field_section(session):
    d = session.decomposition
    field = d.field
    out = {}

    try:
        report = field.report()
    except Exception as exc:                                # noqa: BLE001
        return {'error': '{}: {}'.format(type(exc).__name__, exc)}

    ph = {k: report[k] for k in
          ('interior_index_sum', 'boundary_winding', 'singular_faces', 'ok')
          if k in report}
    out['poincare_hopf'] = ph
    out['poincare_hopf_means'] = (
        'the angle unwrapping is consistent with the measured boundary winding. '
        'It does NOT mean the boundary was sampled finely enough -- run '
        'check_inputs for that.')
    out['solver'] = {
        'relaxed': report.get('relaxed'),
        'iterations': report.get('iterations'),
        'residual': report.get('residual'),
        'min_magnitude': report.get('min_magnitude'),
    }

    singular = field.singularities()
    out['singularities'] = {
        'count': len(singular),
        'positive': sum(1 for _f, k in singular if k > 0),
        'negative': sum(1 for _f, k in singular if k < 0),
    }

    trace = d.trace_report or {}
    mismatch = trace.get('arm_mismatch') or []
    out['arm_mismatch'] = len(mismatch)
    if mismatch:
        out['arm_mismatch_detail'] = [
            'singularity {}: {} launch directions, expected {}'.format(*m)
            for m in mismatch[:6]]
        out['arm_mismatch_means'] = (
            'the layout will come out with the WRONG NUMBER OF PATCHES and '
            'nothing about the coarse mesh says so. This is the field-stage '
            'warning that matters.')
    out['trace_reasons'] = dict(trace.get('reasons') or {})
    out['separatrices'] = len(d.polylines or [])
    out['snap_report'] = dict(d.snap_report or {})
    out['symmetry'] = repr(d.symmetry) if d.symmetry is not None else None
    return out


def _coarse_section(session, max_labels):
    mesh = session.mesh
    out = {}
    sides = [len(mesh.face_vertices(f)) for f in mesh.faces()]
    out['faces'] = len(sides)
    out['quads'] = sum(1 for n in sides if n == 4)
    out['non_quads'] = sorted(set(n for n in sides if n != 4))
    out['poles'] = len(mesh.attributes.get('face_pole') or {})
    out['vertices'] = mesh.number_of_vertices()
    out['edges'] = mesh.number_of_edges()
    # Read off the EDITOR, not the decomposition: a commit writes into the
    # layout it was given and no longer marks the decomposition as edited.
    editor = session.coarse_editor
    out['edited'] = bool(editor is not None and editor.committed)
    out['user_curves'] = len(editor.curves) if editor is not None else 0

    book = session.book()
    if book is not None:
        counts = book.counts()
        out['label_counts'] = counts
        total = sum(counts.values())
        if total <= max_labels:
            out['labels'] = {
                'faces': {k: list(v) for k, v in sorted(book.faces.items())},
                'vertices': {k: list(v) for k, v in sorted(book.vertices.items())},
            }
        else:
            out['labels_note'] = (
                '{} labels -- too many to list. Ask for the ones near a point, '
                'or work from the picture.'.format(total))
    return out


def _dense_section(session, low_angle=None):
    mesh = session.mesh
    d = session.decomposition
    out = {}

    try:
        metrics = d.quality(mesh=mesh, low_angle=low_angle)
    except Exception as exc:                                # noqa: BLE001
        return {'error': '{}: {}'.format(type(exc).__name__, exc)}

    out['quality'] = {k: metrics[k] for k in (
        'faces', 'poles', 'min_angle', 'max_angle', 'aspect_max',
        'share_below', 'irregular_interior', 'worst_face') if k in metrics}
    out['route'] = metrics.get('route')
    out['coverage'] = metrics.get('coverage')
    out['coarse_faces'] = metrics.get('coarse_faces')

    ok, why = hard_floor(metrics)
    out['hard_floor'] = {'ok': ok, 'reason': why}
    if not ok:
        out['hard_floor_means'] = (
            'this element is DEGENERATE, not merely poor. The mesh is still '
            'returned -- quad_mesh always returns a mesh -- but nothing '
            'downstream should assume otherwise.')

    stats = d.densify_stats or {}
    guarded = stats.get('guarded')
    if guarded is not None:
        out['densify'] = {'guarded': guarded, 'patches': stats.get('patches')}
        out['densify_means'] = (
            '{} patch(es) kept the Coons interior instead of the '
            'field-integrated one. {}'.format(
                guarded,
                'A guide was supposed to reach those -- this is a complaint.'
                if session.guides else
                'There are no guides, so this is the ordinary outcome and means '
                'nothing.'))

    sides = [len(mesh.face_vertices(f)) for f in mesh.faces()]
    out['non_quads'] = sorted(set(n for n in sides if n != 4))
    editor = session.dense_editor
    if editor is not None:
        out['strips_available'] = editor.strips_available()
        if not out['strips_available']:
            out['strips_available_means'] = (
                'add_line and remove_line are refused: collect_strip walks '
                'face_opposite_edge and a pole fan has none. move_vertex still '
                'works.')
    return out


def digest(session, max_labels=MAX_LABELS, notes=NOTE_TAIL, low_angle=None):
    """**The whole state, as a dict.** Stage-dependent -- see the module docstring.

    Parameters
    ----------
    session : MeshEditSession
    max_labels : int, optional
        Above this many, labels are counted rather than listed.
    notes : int, optional
        How many of the most recent ``repair_notes`` to carry.
    low_angle : float, optional
        Threshold for ``share_below``, passed to ``quality``.

    Returns
    -------
    dict
        ``state`` (the session's own summary), ``field`` / ``coarse`` / ``dense``
        as the stage allows, ``notes``, ``warnings`` and ``carry``.
    """
    d = session.decomposition
    out = {'state': session.state()}

    out['field'] = _field_section(session)
    if session.stage in ('coarse', 'dense') and session.mesh is not None:
        if session.stage == 'coarse':
            out['coarse'] = _coarse_section(session, max_labels)
        else:
            out['dense'] = _dense_section(session, low_angle=low_angle)
            if session.coarse_editor is not None:
                out['coarse'] = {
                    'faces': session.coarse_editor.mesh.number_of_faces(),
                    'edited': bool(d._edited),
                }

    out['notes'] = list(d.repair_notes[-notes:]) if d.repair_notes else []
    try:
        out['warnings'] = list(d.warnings())
    except Exception as exc:                                # noqa: BLE001
        out['warnings'] = ['warnings() failed: {}: {}'.format(
            type(exc).__name__, exc)]

    carry = session.last_carry
    if carry:
        out['carry'] = {'kept': len(carry.get('kept', {})),
                        'moved': len(carry.get('moved', [])),
                        'lost': len(carry.get('lost', []))}
    return out


# ======================================================================
# words
# ======================================================================

def _fmt(value):
    if isinstance(value, float):
        return '{:.4g}'.format(value)
    return str(value)


def render_text(data, indent=0):
    """A digest as plain text. What actually goes into a prompt.

    Nested dicts are indented; lists of scalars are joined. Deliberately not
    JSON: the ``_means`` keys are sentences meant to be read, and quoting them
    into JSON wastes tokens and reads worse.
    """
    lines = []
    pad = '  ' * indent
    for key, value in data.items():
        # Not every dict here is keyed by string: ``edit_notes['sides']`` maps
        # a side COUNT to how many patches had it, so the keys are ints. This
        # renders whatever it is given rather than assuming.
        key = key if isinstance(key, str) else str(key)
        if key.endswith('_means'):
            lines.append('{}  ({})'.format(pad, value))
            continue
        if isinstance(value, dict):
            if not value:
                continue
            lines.append('{}{}:'.format(pad, key))
            lines.append(render_text(value, indent + 1))
        elif isinstance(value, (list, tuple)):
            if not value:
                continue
            if all(not isinstance(v, (dict, list, tuple)) for v in value):
                lines.append('{}{}: {}'.format(
                    pad, key, ', '.join(_fmt(v) for v in value)))
            else:
                lines.append('{}{}:'.format(pad, key))
                for item in value:
                    lines.append('{}  - {}'.format(pad, _fmt(item)))
        else:
            lines.append('{}{}: {}'.format(pad, key, _fmt(value)))
    return '\n'.join(line for line in lines if line.strip())
