"""The wire: framing, dispatch, notifications, errors, and the stdout guard.

Verification step B. Knows nothing about meshes -- it drives the protocol with a
fake handler, so a failure here is a protocol failure and nothing else.
"""
import io
import json
import sys

from _harness import check, report

from compas_singular.mcp import protocol


class Fake(object):
    def list_tools(self):
        return [{"name": "echo", "description": "echo",
                 "inputSchema": {"type": "object"}}]

    def call_tool(self, name, arguments):
        if name == "boom":
            raise RuntimeError("handler exploded")
        if name == "refuse":
            return {"ok": False, "reason": "not today"}
        return {"ok": True, "got": arguments}

    def list_resources(self):
        return [{"uri": "guidance://quality", "name": "quality"}]

    def read_resource(self, uri):
        return "body of " + uri if uri.startswith("guidance://") else None


class Noisy(Fake):
    def call_tool(self, name, arguments):
        print("THIS PRINT WOULD CORRUPT THE STREAM")
        return {"ok": True}


d = protocol.Dispatcher(Fake(), instructions="hello")

print("\n=== 1. initialize negotiates and advertises honestly ===")
out = d.handle({"jsonrpc": "2.0", "id": 1, "method": "initialize",
                "params": {"protocolVersion": "2024-11-05"}})["result"]
check("echoes a version we speak", out["protocolVersion"] == "2024-11-05")
check("carries serverInfo", "name" in out["serverInfo"])
check("carries instructions", out["instructions"] == "hello")
check("advertises only what the handler has",
      set(out["capabilities"]) == {"tools", "resources"}, out["capabilities"])
fallback = d.handle({"jsonrpc": "2.0", "id": 1, "method": "initialize",
                     "params": {"protocolVersion": "1999-01-01"}})["result"]
check("falls back on a version we do not speak",
      fallback["protocolVersion"] == protocol.PREFERRED_VERSION)

print("\n=== 2. a notification is never answered ===")
check("no reply", d.handle({"jsonrpc": "2.0",
                            "method": "notifications/initialized"}) is None)
check("but the state is recorded", d.initialized is True)
check("an unknown notification is also silent",
      d.handle({"jsonrpc": "2.0", "method": "notifications/whatever"}) is None)

print("\n=== 3. a refusing tool is a RESULT, not an error ===")
ok = d.handle({"jsonrpc": "2.0", "id": 2, "method": "tools/call",
               "params": {"name": "echo", "arguments": {"a": 1}}})["result"]
check("a working call is not an error", ok["isError"] is False)
check("the payload round-trips",
      json.loads(ok["content"][0]["text"])["got"] == {"a": 1})
ref = d.handle({"jsonrpc": "2.0", "id": 3, "method": "tools/call",
                "params": {"name": "refuse"}})
check("a refusal sets isError", ref["result"]["isError"] is True)
check("but is still a result, so the model reads the reason",
      "error" not in ref and "not today" in ref["result"]["content"][0]["text"])

print("\n=== 4. a crashing handler does not kill the session ===")
saved, sys.stderr = sys.stderr, io.StringIO()
boom = d.handle({"jsonrpc": "2.0", "id": 4, "method": "tools/call",
                 "params": {"name": "boom"}})
trace = sys.stderr.getvalue()
sys.stderr = saved
check("it becomes an internal error", boom["error"]["code"] == protocol.INTERNAL_ERROR)
check("naming the exception", "exploded" in boom["error"]["message"])
check("with the traceback on stderr", "RuntimeError" in trace)

print("\n=== 5. protocol faults are JSON-RPC errors ===")
check("unknown method",
      d.handle({"jsonrpc": "2.0", "id": 5, "method": "nope"})["error"]["code"]
      == protocol.METHOD_NOT_FOUND)
check("an unadvertised capability is method-not-found, not a crash",
      d.handle({"jsonrpc": "2.0", "id": 6, "method": "prompts/list"})["error"]["code"]
      == protocol.METHOD_NOT_FOUND)
check("tools/call with no name",
      d.handle({"jsonrpc": "2.0", "id": 7, "method": "tools/call",
                "params": {}})["error"]["code"] == protocol.INVALID_PARAMS)
check("arguments that are not an object",
      d.handle({"jsonrpc": "2.0", "id": 8, "method": "tools/call",
                "params": {"name": "echo", "arguments": 5}})["error"]["code"]
      == protocol.INVALID_PARAMS)
check("a missing resource",
      d.handle({"jsonrpc": "2.0", "id": 9, "method": "resources/read",
                "params": {"uri": "nope://x"}})["error"]["code"]
      == protocol.INVALID_PARAMS)
check("a message that is not an object",
      d.handle("banana")["error"]["code"] == protocol.INVALID_REQUEST)

print("\n=== 6. the stdio loop: framing, blanks, junk, batches ===")
lines = [
    json.dumps({"jsonrpc": "2.0", "id": 1, "method": "initialize", "params": {}}),
    json.dumps({"jsonrpc": "2.0", "method": "notifications/initialized"}),
    "",
    "this is not json",
    json.dumps({"jsonrpc": "2.0", "id": 2, "method": "tools/list"}),
    json.dumps([{"jsonrpc": "2.0", "method": "notifications/initialized"}]),
    json.dumps({"jsonrpc": "2.0", "id": 3, "method": "ping"}),
]
sink = io.StringIO()
protocol.serve(protocol.Dispatcher(Fake()),
               stdin=io.StringIO("\n".join(lines) + "\n"), stdout=sink,
               guard=False)
replies = [json.loads(line) for line in sink.getvalue().splitlines() if line.strip()]
check("one line per message, all valid JSON", len(replies) == 4, len(replies))
check("messages are newline-delimited, not Content-Length framed",
      "Content-Length" not in sink.getvalue())
check("the notification produced no reply", [r.get("id") for r in replies] == [1, None, 2, 3],
      [r.get("id") for r in replies])
check("junk is a parse error and the loop survives",
      replies[1]["error"]["code"] == protocol.PARSE_ERROR)
check("an all-notification batch produces nothing", replies[3]["id"] == 3)
check("ping answers empty", replies[3]["result"] == {})

print("\n=== 7. the stdout guard: a stray print cannot corrupt the stream ===")
sink = io.StringIO()
saved, sys.stderr = sys.stderr, io.StringIO()
protocol.serve(protocol.Dispatcher(Noisy()),
               stdin=io.StringIO(json.dumps(
                   {"jsonrpc": "2.0", "id": 1, "method": "tools/call",
                    "params": {"name": "chatty"}}) + "\n"),
               stdout=sink)
captured = sys.stderr.getvalue()
sys.stderr = saved
check("nothing but protocol reached stdout", "CORRUPT" not in sink.getvalue())
check("every stdout line still parses",
      all(json.loads(l) for l in sink.getvalue().splitlines() if l.strip()))
check("the print was redirected to stderr", "CORRUPT" in captured)
check("sys.stdout is restored afterwards", sys.stdout is not sys.stderr)

sys.exit(report("test_protocol"))
