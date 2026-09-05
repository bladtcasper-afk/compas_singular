"""**Google Gemini as the model behind the loop.**

The second real backend, and the proof that the seam in :mod:`backend` is worth
having: nothing in ``loop.py``, ``core/`` or the CLI changes to run the whole
29-tool pipeline on a different provider.

**It renders images, so the model can see its own work.** That is the reason a
text-only route was not enough here: the loop's job is partly to judge what it
just did, and the digest -- face counts, angles, warnings -- says whether a mesh
is *valid*, not whether it looks right. So this backend carries the PNGs
``render.py`` produces, and a vision-capable model is a requirement rather than
a nicety.

**Written against ``client.models.generate_content``, not ``client.interactions``.**
Both exist in ``google-genai``; ``interactions`` is an untyped passthrough in the
installed version (its ``create`` takes a raw ``body``), while
``models.generate_content`` is fully typed, is what ``types.Content``,
``types.Tool`` and ``types.FunctionDeclaration`` are built for, and takes the
whole conversation on every call -- which is what this loop does anyway.

**Three things are genuinely different from the Anthropic backend, and each one
is handled here rather than leaking upwards.**

*The model's own turns must be resent exactly as they came back.* A ``Part``
carries a ``thought_signature``, and a turn rebuilt from its text and its tool
calls would drop it. So this backend keeps its OWN native history and appends
``response.candidates[0].content`` verbatim, rather than translating the loop's
reconstruction of the assistant turn. That makes it **stateful for the length of
a run**: use one backend per :func:`~..loop.run_session` call.

*A function response needs the function's NAME, and the loop's tool_result block
only carries the id.* Anthropic correlates on ``tool_use_id`` alone. So the ids
handed out with each call are remembered here and looked up when the results
come back.

*Schemas go through ``parameters_json_schema``, not ``parameters``.* The latter
is the restricted OpenAPI ``Schema`` type; the former takes raw JSON Schema,
which is what ``core.tools`` already emits. Only ``additionalProperties`` is
stripped -- the one keyword Google documents as unsupported, and one Gemini
would not enforce either way.

A key comes from ``GEMINI_API_KEY`` or ``GOOGLE_API_KEY``, or is passed in.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import base64
import os
import warnings

from .backend import Backend
from .backend import Reply
from .keys import resolve_api_key


__all__ = ['GeminiBackend', 'DEFAULT_GEMINI_MODEL', 'repair_ssl_env']


def repair_ssl_env():
    """Work around a ``SSL_CERT_FILE`` that points at a file which is not there.

    ``google-genai`` honours ``SSL_CERT_FILE`` explicitly -- it says so in a
    comment, because ``httpx`` otherwise ignores it -- and hands the value
    straight to :func:`ssl.create_default_context`. A missing path therefore
    fails as a bare ``FileNotFoundError`` from inside ``ssl.py``, naming no file
    and mentioning neither the variable nor the SDK.

    **Which is exactly what conda does on Git Bash.** Its
    ``etc/conda/activate.d/openssl_activate.sh`` exports
    ``$CONDA_PREFIX/ssl/cacert.pem`` -- the UNIX layout -- while on Windows the
    bundle is at ``$CONDA_PREFIX/Library/ssl/cacert.pem``, which is what the
    ``.bat`` and ``.ps1`` scripts beside it correctly use. So activating an
    environment in MINGW64 breaks every library that respects the variable,
    ``google-genai`` included, and leaves no clue that it did.

    This repairs the variable for THIS PROCESS only, and only when it is already
    broken: a value that resolves to a real file is never touched, so a
    deliberate corporate trust store still wins. The Anthropic backend needs
    none of this, which is why it lives here rather than in ``backend.py``.

    Returns
    -------
    str or None
        What the variable was changed to, or ``None`` if nothing needed doing.
    """
    current = os.environ.get('SSL_CERT_FILE')
    if not current or os.path.exists(current):
        return None

    replacement = None
    prefix = os.environ.get('CONDA_PREFIX')
    if prefix:
        candidate = os.path.join(prefix, 'Library', 'ssl', 'cacert.pem')
        if os.path.exists(candidate):
            replacement = candidate
    if replacement is None:
        try:
            import certifi
            if os.path.exists(certifi.where()):
                replacement = certifi.where()
        except ImportError:
            pass

    if replacement is None:
        # Nothing better to offer. Drop it and let the SDK fall back to its own
        # default, which is certifi -- still better than a guaranteed crash.
        del os.environ['SSL_CERT_FILE']
        warnings.warn(
            'SSL_CERT_FILE pointed at {!r}, which does not exist, and no '
            'replacement bundle was found. It has been unset for this process '
            'so the SDK can fall back to its own default.'.format(current))
        return None

    os.environ['SSL_CERT_FILE'] = replacement
    warnings.warn(
        'SSL_CERT_FILE pointed at {!r}, which does not exist -- conda sets that '
        'path on Git Bash even though the Windows bundle is elsewhere. Using '
        '{!r} for this process instead. To fix it for good, put '
        '"export SSL_CERT_FILE=$CONDA_PREFIX/Library/ssl/cacert.pem" in your '
        '~/.bashrc.'.format(current, replacement))
    return replacement


#: A Flash model: fast, vision-capable, and the tier Google AI Studio's free
#: quota normally covers. Which models your key can reach is worth confirming in
#: AI Studio -- the free allowance is per-model and changes.
DEFAULT_GEMINI_MODEL = 'gemini-3.7-flash'

#: Tried in order when the chosen model is unavailable. Two things shape it,
#: both measured against a real free-tier key on 2026-08-31:
#:
#: * **quota is per MODEL**, so switching genuinely helps. With 3.7, 3.6, 3.5
#:   and ``flash-latest`` all answering 429, the ``-lite`` models answered
#:   normally -- they have their own bucket. That is why the list ends with
#:   lite variants rather than more full-size Flash models.
#: * **a model can be listed and still not callable.** ``gemini-2.5-flash``
#:   comes back from ``models.list()`` and then answers 404 "no longer available
#:   to new users" on the first request. It was in this list and is now out.
#:
#: Set ``fallback_models=[]`` to fail instead of switching.
GEMINI_FALLBACK_MODELS = ('gemini-3.6-flash', 'gemini-3.5-flash',
                          'gemini-3.5-flash-lite', 'gemini-3.1-flash-lite')

#: HTTP statuses worth retrying: rate limit, and the three server-side ones.
#: 503 is what "this model is currently experiencing high demand" arrives as.
RETRY_STATUS = (429, 500, 502, 503, 504)

#: Attempts per request, including the first. The SDK does NOT retry by default
#: -- every field of ``HttpRetryOptions`` is None unless set -- so without this a
#: single transient 503 ends a run and loses everything the session has done.
RETRY_ATTEMPTS = 5


def _clean_schema(schema):
    """A tool schema Gemini will accept.

    ``parameters_json_schema`` takes raw JSON Schema, so almost nothing has to
    change. ``additionalProperties`` is dropped because Google documents the
    accepted subset as ``type``/``properties``/``description``/``enum``/
    ``required``/``items`` and it buys nothing here -- Gemini does not enforce
    it, and ``core.tools.call`` already rejects an unexpected argument by name.
    """
    if not isinstance(schema, dict):
        return schema
    out = {}
    for key, value in schema.items():
        if key == 'additionalProperties':
            continue
        if isinstance(value, dict):
            out[key] = _clean_schema(value)
        elif isinstance(value, list):
            out[key] = [_clean_schema(v) for v in value]
        else:
            out[key] = value
    return out


class GeminiBackend(Backend):
    """Gemini, through ``google-genai``.

    Parameters
    ----------
    model : str, optional
        Defaults to :data:`DEFAULT_GEMINI_MODEL`.
    api_key : str, optional
        Otherwise the client reads ``GEMINI_API_KEY`` or ``GOOGLE_API_KEY``.
    max_output_tokens : int, optional
    client : optional
        A pre-built ``google.genai.Client``, mostly for tests.

    Notes
    -----
    **One backend per run.** It holds the native conversation, because the
    model's turns have to be resent unmodified -- see the module docstring.
    :meth:`reset` clears it if a caller really wants to reuse the object.
    """

    name = 'gemini'

    def __init__(self, model=None, api_key=None, key_file=None,
                 max_output_tokens=8192, client=None,
                 fallback_models=GEMINI_FALLBACK_MODELS,
                 retry_attempts=RETRY_ATTEMPTS):
        self.fallback_models = list(fallback_models or [])
        self.retry_attempts = retry_attempts
        #: Models this run gave up on, in order. Reported so a caller can say
        #: which model actually produced the result.
        self.switched_from = []
        if client is None:
            # Before the client, not after: the SSL context is built inside the
            # constructor, so a broken SSL_CERT_FILE fails there. See
            # :func:`repair_ssl_env`.
            repair_ssl_env()
            # ``required=False``: with neither a file nor an environment
            # variable, let the SDK do its own resolution rather than raising
            # here on a caller who has configured it some other way.
            api_key = resolve_api_key(
                key=api_key, key_file=key_file,
                env=('GEMINI_API_KEY', 'GOOGLE_API_KEY'), required=False)
            try:
                from google import genai
            except ImportError:
                raise ImportError(
                    'GeminiBackend needs the Google GenAI SDK, which is not '
                    'installed in this interpreter. Install it with '
                    '"pip install google-genai". Everything else in '
                    'compas_singular.agent works without it -- use '
                    'ScriptedBackend to drive the loop with no model at all.')
            # Retry at the TRANSPORT, where the SDK can honour Retry-After and
            # back off properly, rather than re-sending a whole conversation
            # from up here. Every HttpRetryOptions field defaults to None, i.e.
            # no retries at all, so this has to be spelled out.
            from google.genai import types as _types
            http_options = _types.HttpOptions(
                retry_options=_types.HttpRetryOptions(
                    attempts=retry_attempts,
                    initial_delay=1.0,
                    max_delay=30.0,
                    exp_base=2.0,
                    jitter=1.0,
                    http_status_codes=list(RETRY_STATUS)))
            client = (genai.Client(api_key=api_key, http_options=http_options)
                      if api_key else genai.Client(http_options=http_options))
        self.client = client
        self.model = model or DEFAULT_GEMINI_MODEL
        self.max_output_tokens = max_output_tokens
        self.reset()

    def reset(self):
        """Forget the conversation. Call between runs if the object is reused."""
        self._native = []          # list[types.Content], Gemini's own history
        self._consumed = 0         # loop messages already folded in
        self._call_names = {}      # call id -> tool name, for the responses
        #: Models that answered 429 or 503 this run. Never asked again -- see
        #: :meth:`_generate`.
        self._exhausted = set()

    # ------------------------------------------------------------------
    # translation
    # ------------------------------------------------------------------

    def _part(self, block, types):
        """One neutral block as a Gemini ``Part``, or ``None`` to skip it."""
        kind = block['type']
        if kind == 'text':
            return types.Part(text=block['text'])
        if kind == 'image':
            # The neutral format carries base64; from_bytes wants the bytes.
            return types.Part.from_bytes(
                data=base64.b64decode(block['data']),
                mime_type=block.get('media_type', 'image/png'))
        if kind == 'tool_result':
            call_id = block['id']
            return types.Part(function_response=types.FunctionResponse(
                id=call_id,
                # Gemini keys a response by the function's NAME as well as its
                # id, and the loop's block does not carry one -- see the module
                # docstring.
                name=self._call_names.get(call_id, 'unknown'),
                response={'ok': not block.get('is_error'),
                          'result': block.get('content', '')}))
        if kind == 'tool_call':
            # Only reached if a caller hand-builds an assistant turn; the normal
            # path stores the model's own Content instead.
            return types.Part.from_function_call(
                name=block['name'], args=block.get('arguments') or {})
        return None

    def _fold(self, messages, types):
        """Add any loop messages not yet in the native history."""
        for message in messages[self._consumed:]:
            if message['role'] == 'assistant':
                # Already held verbatim from the response it came from. Never
                # rebuilt: a reconstructed turn loses its thought_signature.
                continue
            parts = [p for p in (self._part(b, types) for b in message['content'])
                     if p is not None]
            if parts:
                self._native.append(types.Content(role='user', parts=parts))
        self._consumed = len(messages)

    def _tools(self, tools, types):
        if not tools:
            return None
        declarations = [
            types.FunctionDeclaration(
                name=t['name'], description=t['description'],
                parameters_json_schema=_clean_schema(t['input_schema']))
            for t in tools]
        return [types.Tool(function_declarations=declarations)]

    # ------------------------------------------------------------------

    @staticmethod
    def _is_overloaded(exc):
        """Is this the model being saturated rather than something being wrong?

        503 UNAVAILABLE ("this model is currently experiencing high demand") and
        429 RESOURCE_EXHAUSTED are the two that another model would not have.
        """
        blob = '{}: {}'.format(type(exc).__name__, exc)
        return any(marker in blob for marker in (
            'UNAVAILABLE', '503', 'RESOURCE_EXHAUSTED', '429',
            'high demand', 'overloaded'))

    @staticmethod
    def _is_missing_model(exc):
        """Is this model simply not callable by this key?

        A model can appear in ``models.list()`` and still answer 404 on the
        first request -- ``gemini-2.5-flash`` returns "no longer available to
        new users". Listing is not permission, so the only way to find out is to
        ask, and the right response is to move to the next model rather than end
        the run.
        """
        blob = '{}: {}'.format(type(exc).__name__, exc)
        return any(marker in blob for marker in (
            'NOT_FOUND', '404', 'no longer available', 'is not found',
            'not supported'))

    @classmethod
    def _should_switch(cls, exc):
        """Would ANOTHER model plausibly succeed where this one failed?

        Saturation, quota and a model this key cannot call are all specific to
        the model. A bad schema or a bad key is not: those fail identically
        everywhere, so switching would spend three more calls learning nothing.
        """
        return cls._is_overloaded(exc) or cls._is_missing_model(exc)

    @staticmethod
    def _is_quota(exc):
        """Quota exhausted, as opposed to the model merely being busy.

        429 RESOURCE_EXHAUSTED means the KEY is out of allowance, and on the
        free tier that is a per-model, per-minute or per-day budget. Waiting a
        few seconds does not fix it, so it is worth saying differently from a
        503 spike -- and worth never asking that model again this run.
        """
        blob = '{}: {}'.format(type(exc).__name__, exc)
        return 'RESOURCE_EXHAUSTED' in blob or '429' in blob or 'quota' in blob

    def _generate(self, types, config):
        """One request, switching models if this one is saturated or out of quota.

        The transport already retried with backoff by the time this sees an
        error, so reaching here means the model is not merely spiking.

        **A model that failed once is never asked again this run.** Without that
        memory the candidate list is rebuilt from ``[self.model] + fallbacks``
        every call, so after switching to the third model a failure there sends
        it back to the second -- which is already exhausted. Measured on a real
        free-tier key: 3.7 -> 3.6 -> 3.5 -> back to 3.6 -> 2.5, two of those
        calls spent re-confirming a known dead end and each one costing quota.
        """
        candidates = [m for m in [self.model] + list(self.fallback_models)
                      if m not in self._exhausted]
        # Preserve order, drop duplicates.
        seen, ordered = set(), []
        for model in candidates:
            if model not in seen:
                seen.add(model)
                ordered.append(model)

        if not ordered:
            raise RuntimeError(
                'every model has been exhausted this run ({}). On the free tier '
                'this is usually a per-minute or per-day quota rather than load '
                '-- wait, or use a key with more allowance.'.format(
                    ', '.join(sorted(self._exhausted))))

        last = None
        for index, model in enumerate(ordered):
            try:
                response = self.client.models.generate_content(
                    model=model, contents=list(self._native), config=config)
            except Exception as exc:                        # noqa: BLE001
                last = exc
                if not self._should_switch(exc):
                    raise
                self._exhausted.add(model)
                self.switched_from.append(model)
                if index == len(ordered) - 1:
                    raise
                warnings.warn(
                    '{} is {} ({}). Switching to {} for the rest of this run, '
                    'and not asking {} again.'.format(
                        model,
                        'out of quota' if self._is_quota(exc)
                        else ('not available to this key'
                              if self._is_missing_model(exc) else 'unavailable'),
                        str(exc)[:110], ordered[index + 1], model))
                continue
            if model != self.model:
                # Stay switched: going back would only hit the same wall, and a
                # conversation that changes model every turn is harder to reason
                # about than one that changed once.
                self.model = model
            return response
        raise last

    def send(self, system, tools, messages):
        from google.genai import types

        self._fold(messages, types)

        config = types.GenerateContentConfig(
            system_instruction=system,
            tools=self._tools(tools, types),
            max_output_tokens=self.max_output_tokens)

        response = self._generate(types, config)

        candidates = getattr(response, 'candidates', None) or []
        if not candidates:
            # A safety block or an empty completion. Say so rather than looping
            # on an empty turn.
            feedback = getattr(response, 'prompt_feedback', None)
            return Reply(
                text='[the model returned no candidates: {}]'.format(
                    feedback or 'no reason given'),
                stop_reason='refusal', raw=response)

        content = candidates[0].content
        if content is not None:
            # Verbatim -- signatures and all.
            self._native.append(content)

        text = []
        calls = []
        for i, part in enumerate(getattr(content, 'parts', None) or []):
            if getattr(part, 'text', None) and not getattr(part, 'thought', False):
                text.append(part.text)
            call = getattr(part, 'function_call', None)
            if call is not None:
                # Gemini does not always supply an id; the loop needs one to
                # match the result back, so mint a stable one when it is absent.
                call_id = call.id or 'gem_{}_{}'.format(len(self._native), i)
                self._call_names[call_id] = call.name
                calls.append({'id': call_id, 'name': call.name,
                              'arguments': dict(call.args or {})})

        usage = {}
        meta = getattr(response, 'usage_metadata', None)
        if meta is not None:
            for src, dst in (('prompt_token_count', 'input_tokens'),
                             ('candidates_token_count', 'output_tokens'),
                             ('thoughts_token_count', 'thinking_tokens'),
                             ('cached_content_token_count',
                              'cache_read_input_tokens')):
                value = getattr(meta, src, None)
                if value is not None:
                    usage[dst] = value

        finish = getattr(candidates[0], 'finish_reason', None)
        stop = 'tool_use' if calls else 'end_turn'
        if not calls and finish is not None and 'MAX_TOKENS' in str(finish):
            stop = 'max_tokens'

        return Reply(text='\n'.join(t for t in text if t), tool_calls=calls,
                     stop_reason=stop, usage=usage, raw=response)
