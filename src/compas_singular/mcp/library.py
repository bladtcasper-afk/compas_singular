"""**The knowledge the server serves: guidance, workflows, and worked examples.**

``agent`` keeps the equivalent in ``runner/prompt.py``, as Python string
constants assembled at import time. That makes it byte-stable, which is worth
something for prompt caching, and completely unchangeable without editing the
package, which is worth rather less. Here it is Markdown and JSON on disk, read
per request, exposed over MCP's ``resources`` and ``prompts``. Retuning what
counts as a good mesh is an edit to a text file.

Three kinds of thing live here:

``guidance/*.md``
    Prose the client can read: what good quality looks like, when to reach for
    which smoother, how to write a remark. Also read by
    :mod:`~compas_singular.mcp.describe`, so the thresholds a reading uses and
    the thresholds a model is told about cannot drift apart.
``prompts/*.md``
    Whole workflows, offered through ``prompts/list``.
``sessions/*.json``
    Finished sessions with a verdict on them. ``agent`` writes a transcript and
    never reads one back; this is the half that was missing.

**Nothing here fails loudly.** A missing directory, an unreadable file or a
malformed example is a gap in the library, not a broken server -- the tools must
still work with an empty one, because on a fresh checkout that is what there is.
Every loader returns a default and logs nothing to stdout.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import json
import os

from .describe import DEFAULT_THRESHOLDS


__all__ = [
    'LIBRARY_ENVVAR',
    'library_root',
    'thresholds',
    'list_resources',
    'read_resource',
    'list_prompts',
    'get_prompt',
    'save_example',
    'load_examples',
    'fingerprint',
    'recall',
]


#: Points the library somewhere else -- a working copy, or a project's own.
LIBRARY_ENVVAR = 'COMPAS_SINGULAR_MCP_LIBRARY'

GUIDANCE_SCHEME = 'guidance://'
EXAMPLE_SCHEME = 'example://'
SESSION_SCHEME = 'session://'


def library_root():
    """The library directory: the environment override, or the shipped one."""
    override = (os.environ.get(LIBRARY_ENVVAR) or '').strip()
    if override:
        return override
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), 'library')


def _folder(name):
    return os.path.join(library_root(), name)


def _read(path):
    try:
        with open(path, 'r', encoding='utf-8') as stream:
            return stream.read()
    except (IOError, OSError, ValueError, UnicodeDecodeError):
        return None


def _names_in(folder, suffix):
    try:
        return sorted(n[:-len(suffix)] for n in os.listdir(folder)
                      if n.endswith(suffix))
    except OSError:
        return []


def _summarise(text, limit=160):
    """The first real sentence of a Markdown file, for a resource listing."""
    for line in (text or '').splitlines():
        line = line.strip()
        if not line or line.startswith('#') or line.startswith('>'):
            continue
        line = line.replace('**', '').replace('`', '')
        return line[:limit]
    return ''


def thresholds():
    """Quality bands, from ``thresholds.json``, over the built-in defaults.

    A partial file is merged per metric, so overriding ``min_angle`` alone does
    not silently drop the bands for everything else.
    """
    merged = dict((name, dict(limits))
                  for name, limits in DEFAULT_THRESHOLDS.items())
    raw = _read(os.path.join(library_root(), 'thresholds.json'))
    if not raw:
        return merged
    try:
        loaded = json.loads(raw)
    except ValueError:
        return merged
    if not isinstance(loaded, dict):
        return merged
    for name, limits in loaded.items():
        if isinstance(limits, dict):
            merged.setdefault(name, {}).update(limits)
    return merged


def _unsafe(name):
    """Whether a resource or prompt name would escape the library folder.

    The name comes off the wire, so it is not trusted. Built from ``os.sep``
    rather than a literal separator so it is right on both platforms.
    """
    if not name or name.startswith('.'):
        return True
    return any(sep and sep in name for sep in ('/', os.sep, os.altsep))


# ==============================================================================
# resources
# ==============================================================================

def list_resources():
    """Everything readable, as MCP resource descriptors.

    ``session://current`` is NOT listed here -- it is live state, so the server
    adds it. The library only knows about files.
    """
    out = []
    guidance = _folder('guidance')
    for name in _names_in(guidance, '.md'):
        text = _read(os.path.join(guidance, name + '.md')) or ''
        out.append({'uri': GUIDANCE_SCHEME + name,
                    'name': name.replace('_', ' '),
                    'description': _summarise(text),
                    'mimeType': 'text/markdown'})
    for record in load_examples():
        out.append({'uri': EXAMPLE_SCHEME + record['name'],
                    'name': 'example: ' + record['name'],
                    'description': _example_line(record),
                    'mimeType': 'text/markdown'})
    return out


def read_resource(uri):
    """The body behind a ``guidance://`` or ``example://`` uri. ``None`` if absent."""
    if uri.startswith(GUIDANCE_SCHEME):
        name = uri[len(GUIDANCE_SCHEME):]
        if _unsafe(name):
            return None                       # never escape the library folder
        return _read(os.path.join(_folder('guidance'), name + '.md'))
    if uri.startswith(EXAMPLE_SCHEME):
        name = uri[len(EXAMPLE_SCHEME):]
        for record in load_examples():
            if record['name'] == name:
                return render_example(record)
        return None
    return None


# ==============================================================================
# prompts
# ==============================================================================

def list_prompts():
    folder = _folder('prompts')
    out = []
    for name in _names_in(folder, '.md'):
        text = _read(os.path.join(folder, name + '.md')) or ''
        out.append({'name': name,
                    'description': _summarise(text),
                    'arguments': []})
    return out


def get_prompt(name, arguments=None):
    """One workflow, as a single user message.

    Returns ``None`` when there is no such prompt, which the protocol turns into
    an invalid-params error naming it.
    """
    if _unsafe(name):
        return None
    text = _read(os.path.join(_folder('prompts'), name + '.md'))
    if text is None:
        return None
    return {'description': _summarise(text),
            'messages': [{'role': 'user',
                          'content': {'type': 'text', 'text': text}}]}


# ==============================================================================
# worked examples
# ==============================================================================

def fingerprint(metrics, walls=0, faces=None):
    """The shape of a PROBLEM, for matching one session against another.

    Deliberately coarse, and deliberately not text. The corpus will be small and
    what repeats in it is geometry: how many boundary loops, roughly how big,
    and which defect dominated. Matching on the wording of an instruction would
    rank a differently-shaped mesh above an identically-shaped one.
    """
    metrics = metrics or {}
    faces = faces if faces is not None else metrics.get('faces') or 0
    aspect = metrics.get('aspect_max')
    dominant = 'none'
    limits = thresholds()
    from .describe import band
    worst_band, dominant = 'good', 'none'
    for metric in ('min_angle', 'max_angle', 'aspect_max'):
        if metric not in limits:
            continue
        rating = band(metrics.get(metric), limits[metric])
        order = {'good': 0, 'usable': 1, 'poor': 2, 'unusable': 3, 'unknown': 0}
        if order.get(rating, 0) > order.get(worst_band, 0):
            worst_band, dominant = rating, metric
    return {'loops': int(walls),
            'faces_band': _faces_band(faces),
            'dominant': dominant,
            'severity': worst_band,
            'aspect_max': aspect}


def _faces_band(faces):
    faces = faces or 0
    for edge, name in ((50, 'tiny'), (300, 'small'), (2000, 'medium')):
        if faces < edge:
            return name
    return 'large'


def _score(a, b):
    """How alike two fingerprints are. Higher is closer."""
    score = 0
    if a.get('loops') == b.get('loops'):
        score += 3
    if a.get('faces_band') == b.get('faces_band'):
        score += 2
    if a.get('dominant') == b.get('dominant'):
        score += 4
    if a.get('severity') == b.get('severity'):
        score += 1
    return score


def load_examples():
    """Every stored session. Malformed files are skipped, not raised on."""
    folder = _folder('sessions')
    out = []
    for name in _names_in(folder, '.json'):
        raw = _read(os.path.join(folder, name + '.json'))
        if not raw:
            continue
        try:
            record = json.loads(raw)
        except ValueError:
            continue
        if not isinstance(record, dict):
            continue
        record.setdefault('name', name)
        out.append(record)
    return out


def _example_line(record):
    """One line describing a stored example, for a listing."""
    before = record.get('before') or {}
    after = record.get('after') or {}
    return '{} -- {} steps, min angle {} -> {}, verdict {}'.format(
        record.get('instruction') or record.get('title')
        or 'no instruction recorded',
        len(record.get('steps') or []),
        _round(before.get('min_angle')), _round(after.get('min_angle')),
        record.get('verdict', '?'))


def _round(value):
    try:
        return '{:.1f}'.format(float(value))
    except (TypeError, ValueError):
        return '?'


def render_example(record):
    """A stored session as Markdown, which is what a model reads best."""
    lines = ['# {}'.format(record.get('title') or record.get('name')
                            or 'example'), '']
    verdict = record.get('verdict', '?')
    lines.append('**Verdict: {}**'.format(verdict))
    if record.get('instruction'):
        lines.append('')
        lines.append('Asked: {}'.format(record['instruction']))
    before, after = record.get('before') or {}, record.get('after') or {}
    lines += ['', '| | min angle | max angle | aspect |', '|---|---|---|---|',
              '| before | {} | {} | {} |'.format(
                  _round(before.get('min_angle')), _round(before.get('max_angle')),
                  _round(before.get('aspect_max'))),
              '| after | {} | {} | {} |'.format(
                  _round(after.get('min_angle')), _round(after.get('max_angle')),
                  _round(after.get('aspect_max')))]
    steps = record.get('steps') or []
    if steps:
        lines += ['', '## What was done', '']
        for i, step in enumerate(steps):
            detail = step.get('arguments') or {}
            shown = ', '.join('{}={}'.format(k, v) for k, v in sorted(detail.items()))
            lines.append('{}. `{}({})`{}'.format(
                i + 1, step.get('action', '?'), shown,
                ' -- ' + step['note'] if step.get('note') else ''))
    remarks = record.get('remarks') or []
    if remarks:
        lines += ['', '## Remarks', '']
        lines += ['- {}'.format(r if isinstance(r, str) else r.get('text', ''))
                  for r in remarks]
    if record.get('lesson'):
        lines += ['', '## Lesson', '', record['lesson']]
    return '\n'.join(lines)


def save_example(record, directory=None):
    """Write a finished session into the corpus. Returns the path.

    The name is slugged and collisions get a numeric suffix, so saving twice
    under the same name keeps both rather than overwriting the earlier one --
    the older run is evidence too.
    """
    folder = directory or _folder('sessions')
    if not os.path.isdir(folder):
        os.makedirs(folder)
    stem = _slug(record.get('name') or 'session')
    name, index = stem, 2
    while os.path.exists(os.path.join(folder, name + '.json')):
        name = '{}-{}'.format(stem, index)
        index += 1
    record = dict(record)
    record['name'] = name
    path = os.path.join(folder, name + '.json')
    with open(path, 'w', encoding='utf-8') as stream:
        json.dump(record, stream, indent=2, sort_keys=True, default=str)
    return path


def _slug(text):
    out = []
    for character in str(text).strip().lower():
        if character.isalnum():
            out.append(character)
        elif character in ' -_' and out and out[-1] != '-':
            out.append('-')
    return (''.join(out).strip('-') or 'session')[:60]


def recall(target, limit=3, verdicts=('good', 'acceptable', 'bad')):
    """The stored examples closest to ``target``, best match first.

    Parameters
    ----------
    target : dict
        A fingerprint, from :func:`fingerprint`.
    limit : int, optional
    verdicts : sequence, optional
        Which verdicts to include. A bad example is instructive -- it is what
        NOT to repeat -- so it is included by default.

    Returns
    -------
    list[dict]
        ``{name, score, verdict, summary, uri}``, so the caller can read the
        full text through ``resources/read`` if a summary is not enough.
    """
    out = []
    for record in load_examples():
        if record.get('verdict') not in verdicts:
            continue
        stored = record.get('fingerprint') or {}
        out.append({'name': record['name'],
                    'score': _score(target or {}, stored),
                    'verdict': record.get('verdict', '?'),
                    'fingerprint': stored,
                    'summary': _example_line(record),
                    'uri': EXAMPLE_SCHEME + record['name']})
    out.sort(key=lambda item: (-item['score'], item['name']))
    return out[:max(1, int(limit))]
