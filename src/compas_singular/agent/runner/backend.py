"""**The seam between the loop and whatever model is driving it.**

:class:`Backend` is three things: a message format, one method, and a reply
object. Everything above it -- the loop, the gating, the image policy, the
transcript -- is written against those and nothing else, so swapping the model
means writing one class and changing nothing else.

**Why a hand-written loop and not the SDK's tool runner.**
``client.beta.messages.tool_runner`` is less code and is the right default for
most agents. It is the wrong default here for two reasons. It is
Anthropic-SDK-specific by construction, so adopting it would put the LOOP on the
provider side of this seam -- exactly the thing this module exists to prevent.
And it derives tool schemas from function signatures, while ``core.tools``
hand-writes them because their descriptions carry measured numbers that a
signature cannot express (which ``relax`` setting destroys a closed guide's
layout, which ``guide_band`` is the worst tested). Deriving them would either
lose that or duplicate it.

**The message format is neutral, not Anthropic's.** A message is
``{'role': 'user'|'assistant', 'content': [block, ...]}`` and a block is one of

* ``{'type': 'text', 'text': str}``
* ``{'type': 'tool_call', 'id': str, 'name': str, 'arguments': dict}``
* ``{'type': 'tool_result', 'id': str, 'content': str, 'is_error': bool}``
* ``{'type': 'image', 'media_type': 'image/png', 'data': <base64 str>}``

Translating those into a provider's own blocks is the backend's whole job. It is
mechanical, and it is the only part of the system that is.

**What does NOT port**, and is therefore configured on the backend rather than
passed down from the loop: adaptive thinking, ``output_config.effort``, prompt
caching and context editing. A backend for a provider without them simply does
without; the loop never asks.

Two backends ship here. :class:`AnthropicBackend` is the real one.
:class:`ScriptedBackend` replays a fixed plan and is what the offline test uses
-- it is not a mock bolted on afterwards, it is the second implementation, and
having two is what keeps the seam honest.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


__all__ = ['Reply', 'Backend', 'AnthropicBackend', 'ScriptedBackend',
           'DEFAULT_MODEL']


#: Opus, deliberately. The layout rules here are geometric and mostly about
#: what NOT to do -- a cut that leaves a pentagon, a strip that cannot go, a
#: density that will not densify -- and a cheaper model spends the difference on
#: refused calls. Override per run if a job is simple enough to want it.
DEFAULT_MODEL = 'claude-opus-5'


class Reply(object):
    """One turn back from a model, in the loop's own terms.

    Attributes
    ----------
    text : str
        Whatever the model said in prose this turn.
    tool_calls : list[dict]
        ``{'id', 'name', 'arguments'}`` each. Empty when the model is done.
    stop_reason : str
        ``'tool_use'``, ``'end_turn'``, ``'max_tokens'``, ``'refusal'``, ...
    usage : dict
        Whatever the provider reports. Free-form; only ever displayed.
    raw : object
        The provider's own response, for a caller that needs to dig.
    """

    def __init__(self, text='', tool_calls=None, stop_reason='end_turn',
                 usage=None, raw=None):
        self.text = text or ''
        self.tool_calls = list(tool_calls or [])
        self.stop_reason = stop_reason
        self.usage = dict(usage or {})
        self.raw = raw

    def __repr__(self):
        return '<Reply {} tool call(s), stop={}>'.format(
            len(self.tool_calls), self.stop_reason)


class Backend(object):
    """What the loop needs from a model. Implement :meth:`send` and nothing else."""

    #: Shown in transcripts and errors.
    name = 'backend'

    def send(self, system, tools, messages):
        """One turn.

        Parameters
        ----------
        system : str
            The system prompt. Stable for the whole run -- a backend that caches
            should cache on it.
        tools : list[dict]
            ``{'name', 'description', 'input_schema'}`` from
            ``core.tools.schemas()``. Plain JSON Schema; also stable.
        messages : list[dict]
            The conversation so far, in the neutral format above.

        Returns
        -------
        Reply
        """
        raise NotImplementedError


# ======================================================================
# Anthropic
# ======================================================================

class AnthropicBackend(Backend):
    """The Claude backend.

    Parameters
    ----------
    model : str, optional
        Defaults to :data:`DEFAULT_MODEL`.
    effort : str, optional
        ``'low'`` to ``'max'``. ``'xhigh'`` by default: the model has to reason
        about a geometric rule set before acting, and the refused calls a lower
        effort produces cost more than the thinking saved.
    max_tokens : int, optional
    thinking : bool, optional
        Adaptive thinking. On by default.
    cache : bool, optional
        Put a cache breakpoint on the system prompt and the tool list. Both are
        byte-stable for a whole run and together are a few thousand tokens, so
        this is close to free and worth it from the second turn onwards.
    context_editing : bool, optional
        Let the server clear old tool results once the conversation grows. On by
        default because this loop attaches images, which dominate the context.
    client : optional
        A pre-built ``anthropic.Anthropic``. Otherwise one is constructed with
        no arguments, which picks up ``ANTHROPIC_API_KEY``, ``ANTHROPIC_AUTH_TOKEN``
        or an ``ant auth login`` profile in that order.

    Notes
    -----
    **No** ``strict: true``. It requires every property of a tool's schema to
    appear in ``required``, and 12 of the 29 tools here have genuinely optional
    arguments. Validation is not lost: ``core.tools.call`` catches ``TypeError``
    and returns a refusal naming the bad argument, which the model can read and
    correct -- better feedback than a wire-level rejection.
    """

    name = 'anthropic'

    #: Beta for ``clear_tool_uses_20250919``.
    CONTEXT_BETA = 'context-management-2025-06-27'

    def __init__(self, model=None, effort='xhigh', max_tokens=64000,
                 thinking=True, cache=True, context_editing=True, client=None,
                 api_key=None, key_file=None):
        if client is None:
            from .keys import resolve_api_key
            # ``required=False`` matters here: with nothing passed the client is
            # built zero-arg, which is what lets an `ant auth login` profile
            # work. Raising for a missing variable would break that route.
            api_key = resolve_api_key(
                key=api_key, key_file=key_file,
                env=('ANTHROPIC_API_KEY', 'ANTHROPIC_AUTH_TOKEN'),
                required=False)
            try:
                import anthropic
            except ImportError:
                raise ImportError(
                    'AnthropicBackend needs the Anthropic SDK, which is not '
                    'installed in this interpreter. Install it with '
                    '"pip install anthropic". Everything else in '
                    'compas_singular.agent works without it -- use '
                    'ScriptedBackend to drive the loop with no model at all.')
            client = (anthropic.Anthropic(api_key=api_key) if api_key
                      else anthropic.Anthropic())
        self.client = client
        self.model = model or DEFAULT_MODEL
        self.effort = effort
        self.max_tokens = max_tokens
        self.thinking = thinking
        self.cache = cache
        self.context_editing = context_editing

    # ------------------------------------------------------------------
    # translation
    # ------------------------------------------------------------------

    @staticmethod
    def _block(block):
        """One neutral block as an Anthropic block."""
        kind = block['type']
        if kind == 'text':
            return {'type': 'text', 'text': block['text']}
        if kind == 'tool_call':
            return {'type': 'tool_use', 'id': block['id'],
                    'name': block['name'], 'input': block.get('arguments') or {}}
        if kind == 'tool_result':
            out = {'type': 'tool_result', 'tool_use_id': block['id'],
                   'content': block.get('content', '')}
            if block.get('is_error'):
                out['is_error'] = True
            return out
        if kind == 'image':
            return {'type': 'image',
                    'source': {'type': 'base64',
                               'media_type': block.get('media_type', 'image/png'),
                               'data': block['data']}}
        raise ValueError('unknown block type {!r}'.format(kind))

    def _messages(self, messages):
        return [{'role': m['role'],
                 'content': [self._block(b) for b in m['content']]}
                for m in messages]

    def _tools(self, tools):
        out = [dict(t) for t in tools]
        if out and self.cache:
            # One breakpoint at the end of the tool list covers tools+system,
            # which render before messages and never change during a run.
            out[-1] = dict(out[-1], cache_control={'type': 'ephemeral'})
        return out

    # ------------------------------------------------------------------

    def send(self, system, tools, messages):
        system_blocks = [{'type': 'text', 'text': system}]
        if self.cache:
            system_blocks[0]['cache_control'] = {'type': 'ephemeral'}

        kwargs = {
            'model': self.model,
            'max_tokens': self.max_tokens,
            'system': system_blocks,
            'tools': self._tools(tools),
            'messages': self._messages(messages),
            'output_config': {'effort': self.effort},
        }
        if self.thinking:
            kwargs['thinking'] = {'type': 'adaptive'}

        betas = []
        if self.context_editing:
            betas.append(self.CONTEXT_BETA)
            kwargs['context_management'] = {
                'edits': [{'type': 'clear_tool_uses_20250919'}]}
        if betas:
            kwargs['betas'] = betas

        # Streamed because max_tokens is large and a non-streamed request that
        # size can outlive the SDK's HTTP timeout. get_final_message gives back
        # the whole thing, so nothing downstream has to know.
        endpoint = self.client.beta.messages if betas else self.client.messages
        with endpoint.stream(**kwargs) as stream:
            response = stream.get_final_message()

        text = []
        calls = []
        for block in response.content:
            if block.type == 'text':
                text.append(block.text)
            elif block.type == 'tool_use':
                calls.append({'id': block.id, 'name': block.name,
                              'arguments': dict(block.input or {})})

        usage = {}
        if getattr(response, 'usage', None) is not None:
            for field in ('input_tokens', 'output_tokens',
                          'cache_read_input_tokens',
                          'cache_creation_input_tokens'):
                value = getattr(response.usage, field, None)
                if value is not None:
                    usage[field] = value

        stop = response.stop_reason
        if stop == 'refusal':
            # Populated only on a refusal; guard before reading it.
            details = getattr(response, 'stop_details', None)
            text.append('[the model declined: {}]'.format(
                getattr(details, 'explanation', None) or 'no explanation given'))

        return Reply(text='\n'.join(t for t in text if t), tool_calls=calls,
                     stop_reason=stop, usage=usage, raw=response)


# ======================================================================
# scripted
# ======================================================================

class ScriptedBackend(Backend):
    """Replays a fixed plan. No SDK, no key, no tokens.

    The second implementation of :class:`Backend`, and the reason the loop can
    be tested at all: every branch of the gating, the image policy, the error
    handling and the transcript is reachable without a model.

    Parameters
    ----------
    script : list
        One entry per turn. A **list of dicts** ``{'name', 'arguments'}`` is a
        turn that calls those tools (all of them at once, which is also how a
        real parallel call arrives). A **string** is a turn that says that and
        stops. Anything left when the script runs out behaves as a final turn.
    final : str, optional
        What to say when the script is exhausted.

    Attributes
    ----------
    seen : list
        The ``messages`` list as it was on each :meth:`send`, so a test can
        assert what the loop actually sent -- that tool results came back in ONE
        message, that an image was attached where the policy says.
    """

    name = 'scripted'

    def __init__(self, script=None, final='done'):
        self.script = list(script or [])
        self.final = final
        self.seen = []
        self.turn = 0

    def send(self, system, tools, messages):
        self.seen.append([{'role': m['role'],
                           'content': [dict(b) for b in m['content']]}
                          for m in messages])
        self.turn += 1

        if not self.script:
            return Reply(text=self.final, stop_reason='end_turn')

        step = self.script.pop(0)
        if isinstance(step, str):
            return Reply(text=step, stop_reason='end_turn')

        calls = []
        for i, call in enumerate(step):
            calls.append({'id': 'call_{}_{}'.format(self.turn, i),
                          'name': call['name'],
                          'arguments': dict(call.get('arguments') or {})})
        return Reply(text=step[0].get('say', '') if step else '',
                     tool_calls=calls, stop_reason='tool_use')
