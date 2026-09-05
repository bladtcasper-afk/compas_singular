"""**The loop: an instruction in, an edited mesh and a transcript out.**

Provider-free by construction -- it imports :mod:`backend` for the protocol and
nothing else, and every message it builds is in the neutral format that module
defines. Point it at a different :class:`~.backend.Backend` and nothing here
changes.

**Three policies live here, and each exists because the obvious default is
wrong.**

*Destructive tools are refused unless somebody said otherwise.* Four of the 29
destroy work that cannot be recovered -- ``solve_field`` and ``set_guides``
discard the layout and everything under it, ``commit`` welds and renumbers, and
``revert_to_coarse`` throws away every dense edit. ``confirm=False``, the
default, turns those into refusals carrying an explanation the model can act on;
it does not hide them. Approving them is a decision for whoever starts the run.

*Images are attached at transitions, not every step.* A render after each of
forty steps is mostly forty pictures of the same mesh, and images dominate the
context. The default sends one when the stage changes or when the model asks to
look.

*A failed tool is a message, not an exception.* Every refusal comes back as a
``tool_result`` marked as an error, and the loop continues. The reasons in this
library are specific -- "entered and left through adjacent sides", "10 of 20
strips leave a broken result" -- and are worth more to the model than a stack
trace is to anyone.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import base64

from ..core import tools as _tools
from ..core.digest import digest as _digest
from ..core.digest import render_text
from ..core.render import render_png
from .prompt import system_prompt


__all__ = ['SessionResult', 'run_session']


#: How many steps before the loop stops on its own.
DEFAULT_MAX_STEPS = 40


class SessionResult(object):
    """What a run produced.

    Attributes
    ----------
    mesh : Mesh or None
        The mesh at whatever stage the session ended in.
    stage : str
    transcript : list[dict]
        Every tool call: ``step``, ``name``, ``arguments``, ``approved``,
        ``ok``, ``reason``.
    text : str
        What the model said last.
    digest : dict
        The final state.
    images : list[bytes]
        Every PNG rendered during the run, in order.
    stop_reason : str
        ``'end_turn'`` -- the model finished. ``'max_steps'`` -- the loop did.
        Anything else came from the backend.
    usage : dict
        Accumulated token counts, when the backend reports them.
    steps : int
    """

    def __init__(self, mesh=None, stage='', transcript=None, text='',
                 digest=None, images=None, stop_reason='', usage=None, steps=0):
        self.mesh = mesh
        self.stage = stage
        self.transcript = list(transcript or [])
        self.text = text
        self.digest = digest or {}
        self.images = list(images or [])
        self.stop_reason = stop_reason
        self.usage = dict(usage or {})
        self.steps = steps

    def calls(self, ok=None):
        """Transcript entries, optionally filtered by outcome."""
        if ok is None:
            return list(self.transcript)
        return [e for e in self.transcript if bool(e['ok']) is bool(ok)]

    def summary(self):
        """One line. What ran, what stuck, where it ended."""
        good = len(self.calls(ok=True))
        return ('{} step(s), {} of {} tool call(s) succeeded, ended at stage '
                '{!r} ({})'.format(self.steps, good, len(self.transcript),
                                   self.stage, self.stop_reason))

    def __repr__(self):
        return '<SessionResult {}>'.format(self.summary())


# ======================================================================


def _approve(confirm, name, arguments):
    """``(allowed, reason)`` for one call. Only destructive tools are gated."""
    spec = _tools.TOOLS.get(name)
    if spec is None or not spec.confirm:
        return True, ''
    if confirm is True:
        return True, ''
    if callable(confirm):
        try:
            return bool(confirm(name, dict(arguments))), 'declined by the caller'
        except Exception as exc:                            # noqa: BLE001
            return False, 'the confirmation callback raised {}: {}'.format(
                type(exc).__name__, exc)
    return False, (
        '{} destroys work that cannot be recovered, and this run was started '
        'without permission to do that. Nothing was changed. Either achieve the '
        'goal without it, or stop and report that it is needed and why.'.format(
            name))


def _result_text(name, result):
    """A tool result as the text the model reads."""
    body = {k: v for k, v in result.items()
            if k not in ('ok', 'reason') and v not in (None, '', [], {})}
    head = 'OK' if result.get('ok') else 'REFUSED'
    reason = result.get('reason') or ''
    lines = ['{}: {}{}'.format(name, head, ' -- ' + reason if reason else '')]
    rendered = render_text(body)
    if rendered:
        lines.append(rendered)
    return '\n'.join(lines)


def _image_block(session):
    """A PNG of the current stage as a neutral image block, or ``None``."""
    try:
        png = render_png(session, width=760, height=760)
    except Exception:                                       # noqa: BLE001
        # A picture is a convenience. Never let one stop the run.
        return None, None
    return ({'type': 'image', 'media_type': 'image/png',
             'data': base64.b64encode(png).decode('ascii')}, png)


def run_session(session, instruction, backend=None, confirm=False,
                max_steps=DEFAULT_MAX_STEPS, images='transitions',
                on_event=None, extra_prompt=None):
    """**Drive the session from an instruction.**

    Parameters
    ----------
    session : MeshEditSession
        Built over a solved decomposition. Mutated in place -- the caller keeps
        it and reads the mesh off it afterwards.
    instruction : str
        What to do, in words.
    backend : Backend, optional
        Defaults to :class:`~.backend.AnthropicBackend`, which needs the SDK and
        credentials. Pass a :class:`~.backend.ScriptedBackend` to drive the loop
        without either.
    confirm : bool or callable, optional
        ``False`` (default) refuses the four destructive tools with an
        explanation. ``True`` allows them. A callable ``(name, arguments) ->
        bool`` is asked each time -- this is where a Rhino prompt goes.
    max_steps : int, optional
    images : {'transitions', 'always', 'never'}, optional
        When to attach a render. ``'transitions'`` sends one when the stage
        changes or the model calls ``inspect``.
    on_event : callable, optional
        Called with ``{'type', 'message', ...}`` as the run proceeds. ``print``
        headless; a Rhino prompt in a command.
    extra_prompt : str, optional
        Appended to the system prompt. Keep it stable within a run.

    Returns
    -------
    SessionResult
    """
    if backend is None:
        from .backend import AnthropicBackend
        backend = AnthropicBackend()

    def emit(kind, message, **detail):
        if on_event is not None:
            payload = {'type': kind, 'message': message}
            payload.update(detail)
            on_event(payload)

    system = system_prompt(extra_prompt)
    schemas = _tools.schemas()

    opening = [{'type': 'text', 'text': '{}\n\nThe session is at stage {!r}. '
                                       'Current state:\n\n{}'.format(
                                           instruction, session.stage,
                                           render_text(_digest(session)))}]
    rendered = []
    if images != 'never':
        block, png = _image_block(session)
        if block is not None:
            opening.append(block)
            rendered.append(png)

    messages = [{'role': 'user', 'content': opening}]

    transcript = []
    usage = {}
    text = ''
    stop_reason = 'max_steps'
    steps = 0

    while steps < max_steps:
        steps += 1
        reply = backend.send(system, schemas, messages)

        for key, value in (reply.usage or {}).items():
            usage[key] = usage.get(key, 0) + value
        if reply.text:
            text = reply.text
            emit('text', reply.text)

        if not reply.tool_calls:
            stop_reason = reply.stop_reason or 'end_turn'
            break

        assistant = []
        if reply.text:
            assistant.append({'type': 'text', 'text': reply.text})
        for call in reply.tool_calls:
            assistant.append({'type': 'tool_call', 'id': call['id'],
                              'name': call['name'],
                              'arguments': call['arguments']})
        messages.append({'role': 'assistant', 'content': assistant})

        stage_before = session.stage
        looked = False
        results = []

        for call in reply.tool_calls:
            name = call['name']
            arguments = call['arguments']
            allowed, why = _approve(confirm, name, arguments)

            if not allowed:
                result = {'ok': False, 'reason': why}
                emit('refused', '{} not approved'.format(name), tool=name)
            else:
                result = _tools.call(session, name, **arguments)
                emit('tool', '{}: {}'.format(
                    name, 'ok' if result.get('ok') else result.get('reason', '')),
                    tool=name, ok=bool(result.get('ok')))

            transcript.append({'step': steps, 'name': name,
                               'arguments': dict(arguments),
                               'approved': allowed,
                               'ok': bool(result.get('ok')),
                               'reason': result.get('reason', '')})
            results.append({'type': 'tool_result', 'id': call['id'],
                            'content': _result_text(name, result),
                            'is_error': not result.get('ok')})
            if name == 'inspect':
                looked = True

        # EVERY result in ONE user message. Splitting them across messages
        # silently trains the model out of calling tools in parallel.
        content = list(results)

        want_image = (images == 'always'
                      or (images == 'transitions'
                          and (session.stage != stage_before or looked)))
        if want_image and session.mesh is not None:
            block, png = _image_block(session)
            if block is not None:
                content.append(block)
                rendered.append(png)

        messages.append({'role': 'user', 'content': content})

    emit('done', 'finished: {}'.format(stop_reason), stop_reason=stop_reason)

    return SessionResult(
        mesh=session.mesh, stage=session.stage, transcript=transcript,
        text=text, digest=_digest(session), images=rendered,
        stop_reason=stop_reason, usage=usage, steps=steps)
