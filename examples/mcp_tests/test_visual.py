"""The renderer, the image content block, and the look-before-you-deliver gate.

The gate is the part worth testing hardest. A picture the model never looks at
is worth nothing, and a gate that can be walked past is not a gate.
"""
import base64
import json
import os
import struct
import sys
import tempfile
import zlib

from _harness import check, data_path, library_dir, report, spool_dir

spool_dir()
library_dir()

from compas_singular.datastructures import CoarsePseudoQuadMesh
from compas_singular.datastructures.mesh.smoothing import mesh_boundary_polylines
from compas_singular.mcp import render
from compas_singular.mcp.server import build
from compas_singular.mcp.session import MeshSession

coarse = CoarsePseudoQuadMesh.from_json(
    data_path("coarse_quad_mesh_british_museum_poles.json"))
coarse.collect_strips()
coarse.set_strips_density(3)
coarse.densification()
MESH = coarse.get_quad_mesh()
WALLS = [[list(p) for p in pl.points] for pl in mesh_boundary_polylines(MESH)]
POLES = [list(MESH.vertex_coordinates(v)) for v in MESH.poles()]


def png_size(data):
    """Width and height straight out of the IHDR chunk."""
    assert data[:8] == b"\x89PNG\r\n\x1a\n"
    return struct.unpack(">II", data[16:24])


print("\n=== 1. it produces a real PNG ===")
session = MeshSession().adopt(MESH, walls=WALLS, points=POLES)
data = render.render_png(session, 400, 400)
check("bytes come back", bool(data), len(data or b""))
check("with a PNG signature", data[:8] == b"\x89PNG\r\n\x1a\n")
check("at the size asked for", png_size(data) == (400, 400), png_size(data))
check("and the IDAT decompresses", bool(zlib.decompress(
    data[data.index(b"IDAT") + 4:-12])))
check("it ends with IEND", data[-8:-4] == b"IEND")
check("base64 encodes for the wire",
      base64.b64decode(render.render_png_base64(session, 200, 200)) is not None)

print("\n=== 2. sizes are clamped, not obeyed blindly ===")
check("a huge request is capped", png_size(render.render_png(session, 9000, 9000))
      == (render.MAX_SIZE, render.MAX_SIZE))
check("a tiny one is floored", png_size(render.render_png(session, 1, 1)) == (160, 160))

print("\n=== 3. an empty session draws nothing rather than a white square ===")
check("nothing to draw returns None", render.render_png(MeshSession()) is None)
check("and so does the base64 form", render.render_png_base64(MeshSession()) is None)

print("\n=== 4. THE SIGNAL: orange appears only when the mesh leaves its wall ===")


def orange_pixels(sess):
    """How much of the wall colour is visible. The whole point of the drawing."""
    data = render.render_png(sess, 300, 300)
    raw = zlib.decompress(data[data.index(b"IDAT") + 4:-12])
    width, height = png_size(data)
    stride = width * 3 + 1
    wall = render.COLORS["wall"]
    seen = 0
    for y in range(height):
        row = raw[y * stride + 1:(y + 1) * stride]
        for x in range(width):
            r, g, b = row[x * 3], row[x * 3 + 1], row[x * 3 + 2]
            if (abs(r - wall[0]) < 40 and abs(g - wall[1]) < 40
                    and abs(b - wall[2]) < 40):
                seen += 1
    return seen


on_wall = orange_pixels(MeshSession().adopt(MESH, walls=WALLS, points=POLES))
check("a mesh sitting on its boundary hides the wall almost entirely",
      on_wall < 200, "{} orange pixels".format(on_wall))

drifted = CoarsePseudoQuadMesh.from_json(
    data_path("coarse_quad_mesh_british_museum_poles.json"))
drifted.collect_strips()
drifted.set_strips_density(3)
drifted.densification()
drifted = drifted.get_quad_mesh()
cx = sum(drifted.vertex_coordinates(v)[0] for v in drifted.vertices()) / drifted.number_of_vertices()
cy = sum(drifted.vertex_coordinates(v)[1] for v in drifted.vertices()) / drifted.number_of_vertices()
for v in drifted.vertices_on_boundary():
    x, y, z = drifted.vertex_coordinates(v)
    drifted.vertex_attributes(v, "xyz", [x + (cx - x) * 0.10, y + (cy - y) * 0.10, z])
off_wall = orange_pixels(MeshSession().adopt(drifted, walls=WALLS, points=POLES))
check("a mesh pulled off its boundary exposes it", off_wall > on_wall * 5,
      "{} -> {} orange pixels".format(on_wall, off_wall))

print("\n=== 5. inspect: the image is opt-in, and arrives as an image block ===")
dispatcher = build()
handler = dispatcher.handler
live = handler.session
live.adopt(MESH, walls=WALLS, points=POLES, source={"kind": "file", "name": "x"})


def call_raw(_tool, **arguments):
    return dispatcher.handle({"jsonrpc": "2.0", "id": 1, "method": "tools/call",
                              "params": {"name": _tool, "arguments": arguments}})


def call(_tool, **arguments):
    return json.loads(call_raw(_tool, **arguments)["result"]["content"][0]["text"])


plain = call_raw("inspect")
check("without image=true there is only text",
      [c["type"] for c in plain["result"]["content"]] == ["text"])
check("and no base64 leaks into the payload",
      "_image" not in plain["result"]["content"][0]["text"])

shown = call_raw("inspect", image=True)
kinds = [c["type"] for c in shown["result"]["content"]]
check("with image=true an image block is appended", kinds == ["text", "image"], kinds)
image = shown["result"]["content"][1]
check("declared as a PNG", image["mimeType"] == "image/png")
check("carrying real PNG bytes",
      base64.b64decode(image["data"])[:8] == b"\x89PNG\r\n\x1a\n")
body = json.loads(shown["result"]["content"][0]["text"])
check("the base64 is NOT also in the text", "_image" not in body)
check("a legend comes with it", "orange" in body["image_legend"])
check("and the four checks are spelled out", all(
    word in body["look_at_this"] for word in
    ("BOUNDARIES", "POINT FEATURES", "GUIDES", "POLYEDGE")))
check("the inputs it can be judged against are reported",
      body["inputs"] == {"walls": 2, "guides": 0, "point_features": 4},
      body["inputs"])
check("size is honoured",
      png_size(base64.b64decode(
          call_raw("inspect", image=True, size=300)["result"]["content"][1]["data"]))
      == (300, 300))

print("\n=== 6. THE GATE: nothing is delivered unlooked-at ===")
tmp = tempfile.mkdtemp()
out = os.path.join(tmp, "out.json")

fresh = build()
fresh.handler.session.adopt(MESH, walls=WALLS, points=POLES)


def on(_tool, **kw):
    reply = fresh.handle({"jsonrpc": "2.0", "id": 1, "method": "tools/call",
                          "params": {"name": _tool, "arguments": kw}})
    return json.loads(reply["result"]["content"][0]["text"])


check("a freshly loaded mesh has not been looked at",
      fresh.handler.session.visually_current is False)
refused = on("save_mesh", path=out)
check("save_mesh refuses", refused.get("ok") is False)
check("naming the fix", refused.get("fix", {}).get("arguments") == {"image": True},
      refused.get("fix"))
check("and saying why in words", "image=true" in refused["reason"])
check("nothing was written", not os.path.exists(out))
check("rhino_push refuses too",
      on("rhino_push", timeout=0.2).get("reason", "").startswith("this mesh has changed"))

check("inspect WITHOUT the image does not satisfy the gate",
      on("inspect") and fresh.handler.session.visually_current is False)
on("inspect", image=True)
check("inspect WITH the image does", fresh.handler.session.visually_current is True)
check("save_mesh now works", on("save_mesh", path=out).get("ok") is True)
check("and the file is there", os.path.isfile(out))

print("\n=== 7. any change to the mesh invalidates the last look ===")
os.remove(out)
on("relax")
check("a smoothing pass makes the picture stale",
      fresh.handler.session.visually_current is False)
check("so delivery is refused again", on("save_mesh", path=out).get("ok") is False)
on("inspect", image=True)
on("undo")
check("an undo moves every vertex, so it is stale too",
      fresh.handler.session.visually_current is False)
on("inspect", image=True)
check("looking again clears it", on("save_mesh", path=out).get("ok") is True)
check("state reports it for anyone who asks",
      on("inspect")["state"]["visually_current"] is True)

sys.exit(report("test_visual"))
