# MCP server tests

Offline. No Rhino, no model, no network, no key.

```
C:/Users/Casper/.rhinocode/py39-rh8/python.exe run_all.py
```

Use **Rhino 8's own CPython**. It is the interpreter the server runs on and the
one the link runs on inside Rhino, so a name that resolves here resolves there.
`PYTHONPATH` must include the site-env:

```
C:/Users/Casper/.rhinocode/py39-rh8/site-envs/brg-csd-rjJylh8p
```

Each file runs in its own process. That is deliberate: `test_rhino_side` installs
fake `Rhino` and `rhinoscriptsyntax` modules, and compas decides whether it is
running inside Rhino by looking for exactly those. Sharing an interpreter would
let the stub leak into the tests that must run as if Rhino were absent.

| File | What it covers |
|---|---|
| `test_spool.py` | The file protocol: atomic writes, the rename-as-claim, ordering, timeouts that refuse instead of hanging, crash litter. |
| `test_wire.py` | The mesh encoding: exact round trip, poles as points, plain JSON with no compas types, malformed payloads refused where the fault is. |
| `test_protocol.py` | JSON-RPC: negotiation, notifications never answered, a refusing tool as a result rather than an error, and the stdout guard. |
| `test_library.py` | Guidance, prompts, threshold merging, path-escape refusal, and fingerprint ranking. |
| `test_tools_offline.py` | Every tool through the dispatcher, on a real mesh. |
| `test_rhino_side.py` | `CMD_mcp_link` with Rhino stubbed -- above all the four gates that keep Rhino usable. |
| `test_roundtrip.py` | Pull, improve, push, with a thread playing the link. Covers `tools_rhino`'s success path. |

## What these cannot check

Two things need a real Rhino, and are steps E and E2 of the plan:

- **Whether `RhinoApp.Idle` fires often enough.** It can go quiet when Rhino has
  nothing to redraw. Attach the link, post 20 requests, and count how many drain
  untouched. If it drops them, swap the pump for a blocking `serve()` loop --
  `handle()` does not change.
- **Whether the gates hold in practice.** Draw a polyline, run `CMD_densities`
  end to end, and Ctrl+Z your own edit while requests are in flight. A `push`
  arriving mid-command must defer, your selection must survive, and the bake must
  appear in the undo stack as one step named "MCP push".
