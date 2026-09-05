"""Drive the whole runner with ScriptedBackend -- no SDK, no key, no tokens.

Every branch of the loop is reachable without a model: the confirm gate both
ways, the image policy, a bad tool name, bad arguments, max_steps, and the
transcript. If this passes, the only thing the live backend can get wrong is the
wire format.

Run:
    C:/Users/Casper/anaconda3/envs/singular312/python.exe test_runner_offline.py
"""
import os
import sys
from math import cos, pi, sin

REPO = r"C:\Users\Casper\libraries\carbcomn\compas_singular\compas_singular"
for path in (os.path.join(REPO, "src"),):
    if path not in sys.path:
        sys.path.insert(0, path)

from compas_singular.framefield.decomposition import FieldDecomposition   # noqa: E402
from compas_singular.agent.core import MeshEditSession                   # noqa: E402
from compas_singular.agent.runner import ScriptedBackend, run_session    # noqa: E402
from compas_singular.agent.runner.prompt import system_prompt            # noqa: E402


FAILURES = []


def check(name, condition, detail=""):
    print("  {} {}{}".format("ok  " if condition else "FAIL", name,
                             (" -- " + str(detail)) if detail else ""))
    if not condition:
        FAILURES.append(name)


def square(size=10.0):
    return [[0.0, 0.0, 0.0], [size, 0.0, 0.0], [size, size, 0.0],
            [0.0, size, 0.0], [0.0, 0.0, 0.0]]


def circle(cx, cy, r, n=32):
    pts = [[cx + r * cos(2 * pi * i / n), cy + r * sin(2 * pi * i / n), 0.0]
           for i in range(n)]
    return pts + pts[:1]


OUTER = square()
INNERS = [circle(5.0, 5.0, 2.0)]


def fresh():
    d = FieldDecomposition.from_boundary(OUTER, INNERS, target_length=1.0)
    return MeshEditSession(d, outer=OUTER, inners=INNERS, target_length=2.0)


def calls(*names, **kwargs):
    """One scripted turn calling these tools with no arguments."""
    return [{"name": n, "arguments": kwargs.get(n, {})} for n in names]


def blocks(messages, kind):
    return [b for m in messages for b in m["content"] if b["type"] == kind]


# ======================================================================

print("\nthe prompt")

prompt = system_prompt()
check("builds", len(prompt) > 1500, "{} chars".format(len(prompt)))
for phrase in ("FOUR RULES", "DELETION IS A STRIP", "THREE STATES",
               "check_inputs", "poincare_hopf"):
    check("  carries {!r}".format(phrase), phrase in prompt)
check("takes an extra section",
      "SITE RULE" in system_prompt(extra="SITE RULE: keep the north wall."))
check("is stable across calls", system_prompt() == prompt)


# ======================================================================

print("\na plain run: field -> coarse -> dense")

s = fresh()
backend = ScriptedBackend(script=[
    calls("check_inputs"),
    calls("enter_coarse"),
    calls("rebuild_dense", rebuild_dense={"target_length": 2.0}),
    calls("smooth"),
    "Densified at 2.0 and smoothed.",
])
r = run_session(s, "make a quad mesh", backend=backend, confirm=True)

check("ran to the model's own end", r.stop_reason == "end_turn", r.stop_reason)
check("kept the closing text", "Densified" in r.text, r.text[:40])
check("every scripted call ran", len(r.transcript) == 4, len(r.transcript))
check("all succeeded", len(r.calls(ok=True)) == 4,
      [(e["name"], e["reason"]) for e in r.calls(ok=False)])
check("ended at dense", r.stage == "dense", r.stage)
check("mesh came back", r.mesh is not None and r.mesh.number_of_faces() > 0,
      r.mesh.number_of_faces() if r.mesh else None)
check("digest came back", "state" in r.digest, sorted(r.digest)[:4])
check("summary reads", "4 of 4" in r.summary(), r.summary())

order = [e["name"] for e in r.transcript]
check("transcript is in order",
      order == ["check_inputs", "enter_coarse", "rebuild_dense", "smooth"],
      order)
check("transcript keeps arguments",
      r.transcript[2]["arguments"] == {"target_length": 2.0},
      r.transcript[2]["arguments"])


# ======================================================================

print("\nthe confirm gate")

s = fresh()
backend = ScriptedBackend(script=[calls("enter_coarse"), calls("commit"), "done"])
r = run_session(s, "commit it", backend=backend, confirm=False)
commit = [e for e in r.transcript if e["name"] == "commit"][0]
check("a destructive tool is refused by default", commit["ok"] is False)
check("  and marked unapproved", commit["approved"] is False)
check("  with a reason the model can act on",
      "destroys work" in commit["reason"] and "Nothing was changed" in commit["reason"],
      commit["reason"][:70])
results = blocks(backend.seen[-1], "tool_result")
check("  and it reached the model as an error",
      any(b["is_error"] for b in results), [b["is_error"] for b in results])

s = fresh()
backend = ScriptedBackend(script=[calls("enter_coarse"), calls("commit"), "done"])
r = run_session(s, "commit it", backend=backend, confirm=True)
commit = [e for e in r.transcript if e["name"] == "commit"][0]
check("confirm=True lets it through", commit["ok"] is True, commit["reason"])

asked = []


def ask(name, arguments):
    asked.append(name)
    return name != "commit"


s = fresh()
backend = ScriptedBackend(script=[calls("enter_coarse"), calls("commit"), "done"])
r = run_session(s, "commit it", backend=backend, confirm=ask)
check("a callable is consulted", asked == ["commit"], asked)
check("  and its refusal is honoured",
      [e for e in r.transcript if e["name"] == "commit"][0]["ok"] is False)

s = fresh()
backend = ScriptedBackend(script=[calls("enter_coarse"), "done"])
r = run_session(s, "look", backend=backend, confirm=False)
check("a harmless tool is never gated",
      r.transcript[0]["approved"] is True and r.transcript[0]["ok"] is True)


def explode(name, arguments):
    raise RuntimeError("callback is broken")


s = fresh()
backend = ScriptedBackend(script=[calls("enter_coarse"), calls("commit"), "done"])
r = run_session(s, "commit it", backend=backend, confirm=explode)
commit = [e for e in r.transcript if e["name"] == "commit"][0]
check("a raising callback refuses rather than crashing the run",
      commit["ok"] is False and "RuntimeError" in commit["reason"],
      commit["reason"][:60])


# ======================================================================

print("\nbad calls keep the run alive")

# Three DIFFERENT failure modes, each reached deliberately: an unknown name, a
# wrong stage, and -- only once the stage is right, or the stage gate would
# shadow it -- genuinely bad arguments.
s = fresh()
backend = ScriptedBackend(script=[
    [{"name": "no_such_tool", "arguments": {}}],
    calls("smooth"),                                   # dense-only, at field
    calls("enter_coarse"),
    [{"name": "move_corner", "arguments": {"nonsense": 1}}],   # now at coarse
    calls("list_strips"),
    "recovered",
])
r = run_session(s, "try things", backend=backend, confirm=True)

check("the run survived every bad call", r.stop_reason == "end_turn",
      r.stop_reason)
check("an unknown tool is a refusal",
      r.transcript[0]["ok"] is False and "no tool named" in r.transcript[0]["reason"],
      r.transcript[0]["reason"][:50])
check("a wrong-stage call is a refusal",
      r.transcript[1]["ok"] is False
      and "only defined at stage" in r.transcript[1]["reason"],
      r.transcript[1]["reason"][:60])
check("bad arguments are a refusal -- at the RIGHT stage, so the stage gate "
      "cannot shadow it",
      r.transcript[3]["ok"] is False
      and "bad arguments" in r.transcript[3]["reason"],
      r.transcript[3]["reason"][:70])
check("  naming the tool and the bad argument",
      "move_corner" in r.transcript[3]["reason"]
      and "nonsense" in r.transcript[3]["reason"],
      r.transcript[3]["reason"][:70])
check("and a good call after them all still ran",
      r.transcript[4]["ok"] is True, r.transcript[4]["reason"])


# ======================================================================

print("\nparallel calls come back in ONE message")

s = fresh()
backend = ScriptedBackend(script=[
    [{"name": "enter_coarse", "arguments": {}},
     {"name": "check_inputs", "arguments": {}}],
    "done",
])
r = run_session(s, "both", backend=backend, confirm=True)
check("both ran", len(r.transcript) == 2, len(r.transcript))
last = backend.seen[-1]
user_msgs = [m for m in last if m["role"] == "user"]
results = [b for b in user_msgs[-1]["content"] if b["type"] == "tool_result"]
check("both results are in the final user message", len(results) == 2,
      len(results))


# ======================================================================

print("\nthe image policy")

# The SAME script under all three policies, so the counts can be compared
# rather than each read on its own. One stage change (enter_coarse), one look
# (inspect), and three calls that are neither.
SCRIPT = [calls("enter_coarse"), calls("inspect"), calls("list_strips"),
          calls("list_strips"), calls("list_strips"), "done"]


def images_for(mode):
    return run_session(fresh(), "look around",
                       backend=ScriptedBackend(script=list(SCRIPT)),
                       images=mode, confirm=True).images


transitions = images_for("transitions")
always = images_for("always")
never = images_for("never")

check("images are real PNGs",
      all(png[:8] == b"\x89PNG\r\n\x1a\n" for png in transitions))
check("'transitions' sends one opening + stage change + look",
      len(transitions) == 3, len(transitions))
check("'always' sends one per step as well as the opening",
      len(always) == 6, len(always))
check("  so the policy demonstrably suppresses images",
      len(transitions) < len(always),
      "{} vs {}".format(len(transitions), len(always)))
check("'never' sends none", never == [], len(never))

backend = ScriptedBackend(script=list(SCRIPT))
run_session(fresh(), "look", backend=backend, images="never", confirm=True)
check("  and none reached the model",
      not [b for m in backend.seen[-1] for b in m["content"]
           if b["type"] == "image"])


# ======================================================================

print("\nmax_steps")

s = fresh()
backend = ScriptedBackend(script=[calls("inspect")] * 20)
r = run_session(s, "keep looking", backend=backend, max_steps=3, confirm=True)
check("stopped at the limit", r.steps == 3, r.steps)
check("  and said so", r.stop_reason == "max_steps", r.stop_reason)
check("  having run three calls", len(r.transcript) == 3, len(r.transcript))
check("  and still returns a usable result", r.mesh is None or r.stage == "field")


# ======================================================================

print("\nevents")

seen = []
s = fresh()
backend = ScriptedBackend(script=[calls("enter_coarse"), calls("commit"),
                                  "finished"])
r = run_session(s, "go", backend=backend, confirm=False,
                on_event=lambda e: seen.append(e))
kinds = [e["type"] for e in seen]
check("tool events fire", "tool" in kinds, kinds)
check("refusals fire their own event", "refused" in kinds, kinds)
check("the closing text fires", "text" in kinds, kinds)
check("done fires last", kinds[-1] == "done", kinds[-3:])
check("  carrying the stop reason", seen[-1]["stop_reason"] == "end_turn",
      seen[-1])


# ======================================================================

print("\nwhat the backend actually received")

s = fresh()
backend = ScriptedBackend(script=[calls("enter_coarse"), "done"])
run_session(s, "an instruction that should appear verbatim", backend=backend,
            confirm=True)

first = backend.seen[0]
check("the first message is from the user", first[0]["role"] == "user")
opening = first[0]["content"][0]
check("  and carries the instruction",
      "an instruction that should appear verbatim" in opening["text"])
check("  and the opening digest", "stage" in opening["text"])

second = backend.seen[1]
assistant = [m for m in second if m["role"] == "assistant"]
check("the assistant turn was echoed back", len(assistant) == 1, len(assistant))
tool_calls = [b for b in assistant[0]["content"] if b["type"] == "tool_call"]
check("  with the tool call", len(tool_calls) == 1 and
      tool_calls[0]["name"] == "enter_coarse", tool_calls)
check("blocks are the neutral format, not a provider's",
      all(b["type"] in ("text", "tool_call", "tool_result", "image")
          for m in second for b in m["content"]),
      sorted(set(b["type"] for m in second for b in m["content"])))


# ======================================================================

print("\nthe SDK is not required")

import compas_singular.agent.runner as runner_pkg                        # noqa: E402
check("importing the runner does not import anthropic",
      "anthropic" not in sys.modules)
check("  and AnthropicBackend is exported anyway",
      hasattr(runner_pkg, "AnthropicBackend"))
try:
    runner_pkg.AnthropicBackend()
    made = True
    message = ""
except ImportError as exc:
    made = False
    message = str(exc)
except Exception as exc:                                    # noqa: BLE001
    made = False
    message = "{}: {}".format(type(exc).__name__, exc)
if not made and "pip install anthropic" in message:
    check("constructing it without the SDK says how to fix it", True,
          message[:60])
else:
    check("constructing it without the SDK says how to fix it",
          made, "unexpected: {}".format(message[:80]))


# ======================================================================

# ======================================================================
# The one part of the LIVE path that can be checked without credentials: what
# AnthropicBackend actually puts on the wire. A fake client captures the kwargs.
# ======================================================================

print("\nthe Anthropic request shape")


class _Block(object):
    def __init__(self, **kw):
        self.__dict__.update(kw)


class _Stream(object):
    def __init__(self, message):
        self.message = message

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def get_final_message(self):
        return self.message


class _Messages(object):
    def __init__(self, sink):
        self.sink = sink

    def stream(self, **kwargs):
        self.sink.append(kwargs)
        return _Stream(_Block(
            content=[_Block(type="text", text="looking"),
                     _Block(type="tool_use", id="t1", name="inspect", input={})],
            stop_reason="tool_use",
            usage=_Block(input_tokens=10, output_tokens=5,
                         cache_read_input_tokens=0,
                         cache_creation_input_tokens=100)))


class _FakeClient(object):
    def __init__(self):
        self.sent = []
        self.beta = _Block(messages=_Messages(self.sent))
        self.messages = _Messages(self.sent)


from compas_singular.agent.runner import AnthropicBackend                # noqa: E402
from compas_singular.agent.core import tools as _core_tools              # noqa: E402

fake = _FakeClient()
b = AnthropicBackend(client=fake)
reply = b.send("SYSTEM", _core_tools.schemas(),
               [{"role": "user", "content": [
                   {"type": "text", "text": "hello"},
                   {"type": "image", "media_type": "image/png", "data": "AAAA"}]}])

check("no key was needed to build the request", True)
sent = fake.sent[-1]
check("model is Opus by default", sent["model"] == "claude-opus-5", sent["model"])
check("adaptive thinking, not budget_tokens",
      sent["thinking"] == {"type": "adaptive"} and "budget_tokens" not in str(sent["thinking"]),
      sent["thinking"])
check("effort is inside output_config",
      sent["output_config"] == {"effort": "xhigh"}, sent.get("output_config"))
check("context editing is on, with its beta",
      "context-management-2025-06-27" in sent["betas"]
      and sent["context_management"]["edits"][0]["type"] == "clear_tool_uses_20250919",
      sent.get("context_management"))
check("the system prompt is cached",
      sent["system"][0]["cache_control"] == {"type": "ephemeral"},
      sent["system"][0].get("cache_control"))
check("the tool list carries the breakpoint on its LAST entry",
      "cache_control" in sent["tools"][-1]
      and not any("cache_control" in t for t in sent["tools"][:-1]),
      [i for i, t in enumerate(sent["tools"]) if "cache_control" in t])
check("no tool is marked strict",
      not any(t.get("strict") for t in sent["tools"]))
check("all 29 tools were sent", len(sent["tools"]) == len(_core_tools.tool_names()),
      len(sent["tools"]))

blocks_sent = sent["messages"][0]["content"]
check("a text block became an Anthropic text block",
      blocks_sent[0] == {"type": "text", "text": "hello"}, blocks_sent[0])
check("an image block became a base64 source block",
      blocks_sent[1]["type"] == "image"
      and blocks_sent[1]["source"]["type"] == "base64"
      and blocks_sent[1]["source"]["media_type"] == "image/png",
      blocks_sent[1])

check("the reply was parsed back into the neutral shape",
      reply.tool_calls == [{"id": "t1", "name": "inspect", "arguments": {}}],
      reply.tool_calls)
check("  keeping the text", reply.text == "looking", reply.text)
check("  and the usage", reply.usage.get("cache_creation_input_tokens") == 100,
      reply.usage)

nc = AnthropicBackend(client=_FakeClient(), cache=False, context_editing=False,
                      thinking=False, model="claude-sonnet-5")
nc.send("S", _core_tools.schemas()[:2], [{"role": "user",
                                          "content": [{"type": "text", "text": "x"}]}])
sent = nc.client.sent[-1]
check("cache can be turned off",
      "cache_control" not in sent["system"][0]
      and not any("cache_control" in t for t in sent["tools"]))
check("context editing can be turned off",
      "betas" not in sent and "context_management" not in sent, sorted(sent))
check("thinking can be turned off", "thinking" not in sent)
check("the model is a parameter", sent["model"] == "claude-sonnet-5", sent["model"])

# A tool_result round trip, since that is the block the loop sends most.
rt = AnthropicBackend(client=_FakeClient())
rt.send("S", [], [{"role": "user", "content": [
    {"type": "tool_result", "id": "t1", "content": "REFUSED", "is_error": True}]}])
block = rt.client.sent[-1]["messages"][0]["content"][0]
check("a tool_result carries tool_use_id and is_error",
      block["type"] == "tool_result" and block["tool_use_id"] == "t1"
      and block["is_error"] is True, block)


# ======================================================================
# Gemini, the same way: a fake client captures what goes on the wire. The types
# are REAL google.genai types, so this also checks the constructions are legal.
# ======================================================================

print("\nthe Gemini request shape")

try:
    from google.genai import types as gtypes
    HAVE_GENAI = True
except ImportError:
    HAVE_GENAI = False
    print("  (google-genai not installed -- skipped)")

if HAVE_GENAI:
    from compas_singular.agent.runner import GeminiBackend                # noqa: E402
    from compas_singular.agent.runner.backend_gemini import _clean_schema  # noqa: E402

    class _GModels(object):
        def __init__(self, sink, reply_parts):
            self.sink = sink
            self.reply_parts = reply_parts
            self.returned = []

        def generate_content(self, model, contents, config):
            self.sink.append({"model": model, "contents": list(contents),
                              "config": config})
            content = gtypes.Content(role="model", parts=self.reply_parts)
            # Kept so the test can assert the backend stored THIS object, not a
            # reconstruction -- that is what preserves thought_signature.
            self.returned.append(content)
            cand = _Block(content=content, finish_reason="STOP")
            return _Block(candidates=[cand], prompt_feedback=None,
                          usage_metadata=_Block(prompt_token_count=11,
                                                candidates_token_count=4,
                                                thoughts_token_count=2,
                                                cached_content_token_count=0))

    class _GClient(object):
        def __init__(self, reply_parts):
            self.sent = []
            self.models = _GModels(self.sent, reply_parts)

    reply_parts = [
        gtypes.Part(text="having a look"),
        gtypes.Part(function_call=gtypes.FunctionCall(
            id=None, name="inspect", args={})),
    ]
    gclient = _GClient(reply_parts)
    gb = GeminiBackend(client=gclient)

    reply = gb.send("SYSTEM", _core_tools.schemas(),
                    [{"role": "user", "content": [
                        {"type": "text", "text": "hello"},
                        {"type": "image", "media_type": "image/png",
                         "data": "iVBORw0KGgo="}]}])

    sent = gclient.sent[-1]
    check("model defaults to a Flash model", "flash" in sent["model"],
          sent["model"])
    check("the system prompt goes in system_instruction",
          sent["config"].system_instruction == "SYSTEM")
    check("all 29 tools became function declarations",
          len(sent["config"].tools[0].function_declarations)
          == len(_core_tools.tool_names()),
          len(sent["config"].tools[0].function_declarations))
    decl = sent["config"].tools[0].function_declarations[0]
    check("schemas go through parameters_json_schema",
          decl.parameters_json_schema is not None and decl.parameters is None,
          type(decl.parameters_json_schema).__name__)
    check("  with additionalProperties stripped",
          not any("additionalProperties" in str(d.parameters_json_schema)
                  for d in sent["config"].tools[0].function_declarations))
    check("  but the properties kept",
          "properties" in decl.parameters_json_schema,
          sorted(decl.parameters_json_schema))

    parts = sent["contents"][0].parts
    check("the user turn has role 'user'", sent["contents"][0].role == "user")
    check("a text block became a text part", parts[0].text == "hello")
    check("an image block became inline_data with its mime type",
          parts[1].inline_data is not None
          and parts[1].inline_data.mime_type == "image/png",
          parts[1].inline_data.mime_type if parts[1].inline_data else None)
    check("  and was base64-DECODED to bytes",
          isinstance(parts[1].inline_data.data, bytes)
          and parts[1].inline_data.data.startswith(b"\x89PNG"),
          repr(parts[1].inline_data.data[:6]))

    check("the reply parsed into the neutral shape",
          len(reply.tool_calls) == 1 and reply.tool_calls[0]["name"] == "inspect",
          reply.tool_calls)
    check("  minting an id when Gemini supplies none",
          bool(reply.tool_calls[0]["id"]), reply.tool_calls[0]["id"])
    check("  keeping the text", reply.text == "having a look", reply.text)
    check("  and mapping usage to the common names",
          reply.usage.get("input_tokens") == 11
          and reply.usage.get("thinking_tokens") == 2, reply.usage)
    check("  stop_reason is tool_use when it called something",
          reply.stop_reason == "tool_use", reply.stop_reason)

    # Turn two: the model's own Content must come back verbatim, and the
    # tool_result must find the function NAME from the remembered id.
    call_id = reply.tool_calls[0]["id"]
    gb.send("SYSTEM", _core_tools.schemas(), [
        {"role": "user", "content": [{"type": "text", "text": "hello"}]},
        {"role": "assistant", "content": [
            {"type": "text", "text": "having a look"},
            {"type": "tool_call", "id": call_id, "name": "inspect",
             "arguments": {}}]},
        {"role": "user", "content": [
            {"type": "tool_result", "id": call_id, "content": "OK",
             "is_error": False}]}])

    contents = gclient.sent[-1]["contents"]
    roles = [c.role for c in contents]
    check("the model turn is in history as role 'model'", "model" in roles, roles)
    model_turn = contents[roles.index("model")]
    check("  and is the SAME object the API returned, not a rebuild",
          model_turn is gclient.models.returned[0],
          "identity holds" if model_turn is gclient.models.returned[0]
          else "a rebuilt turn would drop thought_signature")
    response_part = contents[-1].parts[0]
    check("the tool result became a function_response",
          response_part.function_response is not None)
    check("  carrying the function NAME looked up from the id",
          response_part.function_response.name == "inspect",
          response_part.function_response.name)
    check("  and the id", response_part.function_response.id == call_id)
    check("  with ok/result in the response dict",
          response_part.function_response.response == {"ok": True, "result": "OK"},
          response_part.function_response.response)

    check("an assistant turn is never rebuilt into history",
          roles.count("model") == 1, roles)

    check("_clean_schema leaves a plain schema alone",
          _clean_schema({"type": "object",
                         "properties": {"x": {"type": "number"}}})
          == {"type": "object", "properties": {"x": {"type": "number"}}})
    check("  and strips additionalProperties at any depth",
          "additionalProperties" not in str(_clean_schema(
              {"type": "object", "additionalProperties": False,
               "properties": {"a": {"type": "object",
                                    "additionalProperties": False}}})))

    check("no key was needed for any of this", True)

    # Saturation. A free-tier 503 on the newest model is routine, and losing a
    # whole edited session to one would be absurd.
    print("\n  when the model is saturated")

    class _Saturating(object):
        def __init__(self, fail_for, exc):
            self.fail_for, self.exc, self.tried = fail_for, exc, []

        def generate_content(self, model, contents, config):
            self.tried.append(model)
            if model in self.fail_for:
                raise self.exc
            content = gtypes.Content(
                role="model", parts=[gtypes.Part(text="ok from " + model)])
            return _Block(candidates=[_Block(content=content,
                                             finish_reason="STOP")],
                          prompt_feedback=None, usage_metadata=None)

    class _SatClient(object):
        def __init__(self, models):
            self.models = models

    _msg = [{"role": "user", "content": [{"type": "text", "text": "hi"}]}]

    import warnings as _w                                                # noqa: E402
    sat = _Saturating({"gemini-3.7-flash"},
                      RuntimeError("503 UNAVAILABLE. high demand"))
    b = GeminiBackend(client=_SatClient(sat))
    with _w.catch_warnings(record=True) as caught:
        _w.simplefilter("always")
        r = b.send("S", [], _msg)
    check("a 503 switches to the next model rather than failing",
          r.text == "ok from gemini-3.6-flash", r.text)
    check("  trying them in order", sat.tried
          == ["gemini-3.7-flash", "gemini-3.6-flash"], sat.tried)
    check("  and staying switched", b.model == "gemini-3.6-flash", b.model)
    check("  recording what it gave up on",
          b.switched_from == ["gemini-3.7-flash"], b.switched_from)
    check("  loudly", len(caught) == 1 and "Switching to" in str(caught[0].message))

    bad = _Saturating({"gemini-3.7-flash"},
                      RuntimeError("400 INVALID_ARGUMENT bad schema"))
    b2 = GeminiBackend(client=_SatClient(bad))
    try:
        b2.send("S", [], _msg)
        raised = False
    except RuntimeError:
        raised = True
    check("a non-overload error raises instead of switching", raised)
    check("  without wasting calls on other models",
          bad.tried == ["gemini-3.7-flash"], bad.tried)

    b_all = GeminiBackend(client=_SatClient(_Saturating(set(), None)))
    every = set([b_all.model] + list(b_all.fallback_models))
    allsat = _Saturating(every, RuntimeError("503 UNAVAILABLE"))
    b3 = GeminiBackend(client=_SatClient(allsat))
    with _w.catch_warnings(record=True):
        _w.simplefilter("ignore")
        try:
            b3.send("S", [], _msg)
            raised = False
        except RuntimeError:
            raised = True
    check("everything saturated still raises", raised)
    check("  having tried every fallback exactly once",
          sorted(allsat.tried) == sorted(every), allsat.tried)

    b4 = GeminiBackend(client=_SatClient(_Saturating(set(), None)),
                       fallback_models=[])
    check("fallbacks can be turned off", b4.fallback_models == [])

    # A model that failed must never be asked again this run. Without that
    # memory the candidate list is rebuilt every call and cycles back to a
    # known-dead model -- measured on a real free-tier key as
    # 3.7 -> 3.6 -> 3.5 -> back to 3.6 -> 2.5, each wasted call costing quota.
    quota = _Saturating({"gemini-3.7-flash", "gemini-3.6-flash"},
                        RuntimeError("429 RESOURCE_EXHAUSTED quota"))
    b5 = GeminiBackend(client=_SatClient(quota))
    with _w.catch_warnings(record=True) as caught:
        _w.simplefilter("always")
        b5.send("S", [], _msg)
    check("quota failures are named as quota, not load",
          any("out of quota" in str(c.message) for c in caught),
          [str(c.message)[:50] for c in caught][:1])
    check("  and it settles on the first live model",
          b5.model == "gemini-3.5-flash", b5.model)

    # A model can be LISTED and still answer 404 "no longer available to new
    # users" -- gemini-2.5-flash does exactly that. That is model-specific, so
    # it must switch like an overload does, not end the run.
    gone = _Saturating({"gemini-3.7-flash"}, RuntimeError(
        "404 NOT_FOUND. This model models/gemini-3.7-flash is no longer "
        "available to new users."))
    b404 = GeminiBackend(client=_SatClient(gone))
    with _w.catch_warnings(record=True) as caught:
        _w.simplefilter("always")
        r404 = b404.send("S", [], _msg)
    check("a 404 for a missing model switches rather than failing",
          r404.text == "ok from gemini-3.6-flash", r404.text)
    check("  and says the model is not available to this key",
          any("not available to this key" in str(c.message) for c in caught),
          [str(c.message)[:60] for c in caught][:1])

    check("gemini-2.5-flash is no longer a fallback -- it 404s on a real key",
          "gemini-2.5-flash" not in b404.fallback_models,
          b404.fallback_models)

    quota.fail_for.add("gemini-3.5-flash")
    quota.tried = []
    with _w.catch_warnings(record=True):
        _w.simplefilter("ignore")
        b5.send("S", [], _msg)
    check("a second failure does NOT retry already-dead models",
          "gemini-3.7-flash" not in quota.tried
          and "gemini-3.6-flash" not in quota.tried, quota.tried)
    check("  going straight to the next live one",
          quota.tried == ["gemini-3.5-flash", "gemini-3.5-flash-lite"],
          quota.tried)
    check("  and remembering every one that died",
          b5._exhausted == {"gemini-3.7-flash", "gemini-3.6-flash",
                            "gemini-3.5-flash"}, sorted(b5._exhausted))

    # The SSL_CERT_FILE trap. conda's Git Bash activation script exports the
    # UNIX path ($CONDA_PREFIX/ssl/cacert.pem) on Windows, where the bundle is
    # under Library/ssl -- and google-genai hands the variable straight to
    # ssl.create_default_context, so it dies with a bare FileNotFoundError.
    import os as _os                                                     # noqa: E402
    import warnings as _warnings                                         # noqa: E402
    from compas_singular.agent.runner.backend_gemini import repair_ssl_env  # noqa: E402

    _held = _os.environ.get("SSL_CERT_FILE")
    _held_prefix = _os.environ.get("CONDA_PREFIX")
    try:
        broken = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)),
                               "no-such-bundle-should-not-exist.pem")
        _os.environ["SSL_CERT_FILE"] = broken
        with _warnings.catch_warnings(record=True) as caught:
            _warnings.simplefilter("always")
            fixed = repair_ssl_env()
        check("a missing SSL_CERT_FILE is repaired",
              fixed is None or _os.path.exists(fixed), fixed)
        check("  and the replacement is a real file",
              _os.path.exists(_os.environ.get("SSL_CERT_FILE", fixed or "")),
              _os.environ.get("SSL_CERT_FILE"))
        check("  loudly, not silently", len(caught) == 1, len(caught))
        check("  naming the variable and the fix",
              "SSL_CERT_FILE" in str(caught[0].message)
              and "bashrc" in str(caught[0].message))

        good = _os.environ["SSL_CERT_FILE"]
        with _warnings.catch_warnings(record=True) as caught:
            _warnings.simplefilter("always")
            again = repair_ssl_env()
        check("a VALID SSL_CERT_FILE is never touched",
              again is None and _os.environ["SSL_CERT_FILE"] == good)
        check("  and says nothing", len(caught) == 0, len(caught))

        del _os.environ["SSL_CERT_FILE"]
        check("an unset SSL_CERT_FILE is left unset",
              repair_ssl_env() is None and "SSL_CERT_FILE" not in _os.environ)
    finally:
        _os.environ.pop("SSL_CERT_FILE", None)
        if _held is not None:
            _os.environ["SSL_CERT_FILE"] = _held
        _os.environ.pop("CONDA_PREFIX", None)
        if _held_prefix is not None:
            _os.environ["CONDA_PREFIX"] = _held_prefix


# ======================================================================
# Keys from a file. Rhino reads its environment once at launch, so a key
# exported in a terminal never reaches a running Rhino -- hence a path.
# ======================================================================

print("\nkeys from a file")

import os as _o                                                          # noqa: E402
import tempfile as _tf                                                   # noqa: E402
from compas_singular.agent.runner import KeyNotFound                     # noqa: E402
from compas_singular.agent.runner import read_key_file, resolve_api_key  # noqa: E402

_dir = _tf.mkdtemp(prefix="agent_keys_")


def keyfile(name, text):
    path = _o.path.join(_dir, name)
    with open(path, "w") as handle:
        handle.write(text)
    return path


check("reads a bare key",
      read_key_file(keyfile("bare.txt", "AIza-secret\n")) == "AIza-secret")
check("reads NAME=value",
      read_key_file(keyfile("named.txt", "GEMINI_API_KEY=AIza-named\n"))
      == "AIza-named")
check("skips comments and blank lines",
      read_key_file(keyfile("commented.txt",
                            "# created 2026-08-29\n\n  \nAIza-after\n"))
      == "AIza-after")
check("strips quotes and whitespace",
      read_key_file(keyfile("quoted.txt", '  "AIza-quoted"  \n'))
      == "AIza-quoted")
check("takes only the first real line",
      read_key_file(keyfile("two.txt", "AIza-first\nAIza-second\n"))
      == "AIza-first")

for name, text, why in (("empty.txt", "", "an empty file"),
                        ("only-comments.txt", "# nothing here\n\n",
                         "a file with only comments")):
    try:
        read_key_file(keyfile(name, text))
        raised, message = False, ""
    except KeyNotFound as exc:
        raised, message = True, str(exc)
    check("{} raises KeyNotFound".format(why), raised, message[:60])
    check("  naming the file", name in message, message[:70])

missing = _o.path.join(_dir, "not-here.txt")
try:
    read_key_file(missing)
    raised, message = False, ""
except KeyNotFound as exc:
    raised, message = True, str(exc)
check("a missing file raises KeyNotFound", raised)
check("  naming the path", "not-here.txt" in message, message[:70])

good = keyfile("good.txt", "AIza-from-file\n")
_o.environ["_AGENT_TEST_KEY"] = "AIza-from-env"
try:
    check("explicit key wins over everything",
          resolve_api_key(key="AIza-explicit", key_file=good,
                          env=("_AGENT_TEST_KEY",)) == "AIza-explicit")
    check("a file wins over the environment",
          resolve_api_key(key_file=good, env=("_AGENT_TEST_KEY",))
          == "AIza-from-file")
    check("the environment is the fallback",
          resolve_api_key(env=("_AGENT_TEST_KEY",)) == "AIza-from-env")
    check("later env names are tried in order",
          resolve_api_key(env=("_AGENT_NOPE", "_AGENT_TEST_KEY"))
          == "AIza-from-env")

    # The one that matters: a path that was GIVEN but is wrong must raise, not
    # quietly serve a stale environment key. Otherwise you edit the file, see no
    # change, and have no idea why.
    try:
        resolve_api_key(key_file=missing, env=("_AGENT_TEST_KEY",))
        fell_through = True
    except KeyNotFound:
        fell_through = False
    check("a bad key_file raises rather than falling back to the environment",
          not fell_through)

    check("nothing found and required=False returns None",
          resolve_api_key(env=("_AGENT_NOPE",), required=False) is None)
    try:
        resolve_api_key(env=("_AGENT_NOPE",))
        raised = False
    except KeyNotFound as exc:
        raised, message = True, str(exc)
    check("nothing found and required=True raises", raised)
    check("  naming the variable it looked for", "_AGENT_NOPE" in message,
          message[:70])
finally:
    _o.environ.pop("_AGENT_TEST_KEY", None)

# Both backends accept key_file, using their fake clients so nothing dials out.
gb = GeminiBackend(client=_GClient([gtypes.Part(text="hi")]),
                   key_file=good) if HAVE_GENAI else None
check("GeminiBackend accepts key_file", gb is not None or not HAVE_GENAI)

ab = AnthropicBackend(client=_FakeClient(), key_file=good)
check("AnthropicBackend accepts key_file", ab is not None)

import shutil as _sh                                                     # noqa: E402
_sh.rmtree(_dir, ignore_errors=True)


print("\n{} check(s) failed{}".format(
    len(FAILURES), (": " + ", ".join(FAILURES)) if FAILURES else ""))
sys.exit(1 if FAILURES else 0)
