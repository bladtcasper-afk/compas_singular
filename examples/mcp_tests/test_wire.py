"""The mesh encoding that crosses to Rhino. Part of verification step A.

The wire is plain JSON with no compas types in it, and that is load-bearing
rather than tidy: it is what lets the link survive another command purging
``compas_singular`` out of ``sys.modules``.
"""
import json
import sys

from _harness import check, data_path, report

from compas_singular.datastructures import CoarsePseudoQuadMesh, CoarseQuadMesh
from compas_singular.mcp.bridge import wire

print("\n=== 1. a quad mesh round-trips exactly ===")
mesh = CoarseQuadMesh.from_json(
    data_path("coarse_quad_mesh_british_museum.json"))
payload = wire.mesh_to_wire(mesh)
check("vertices are a list of triples",
      all(len(v) == 3 for v in payload["vertices"]))
check("faces index INTO that list, not vertex keys",
      all(all(isinstance(i, int) and 0 <= i < len(payload["vertices"]) for i in f)
          for f in payload["faces"]))
back = wire.mesh_from_wire(payload)
check("vertex count survives", back.number_of_vertices() == mesh.number_of_vertices())
check("face count survives", back.number_of_faces() == mesh.number_of_faces())
check("coordinates are not rounded",
      wire.mesh_to_wire(back)["vertices"] == payload["vertices"])
check("it comes back a QuadMesh", type(back).__name__ == "QuadMesh", type(back).__name__)

print("\n=== 2. poles travel as POINTS, because keys renumber ===")
poled = CoarsePseudoQuadMesh.from_json(
    data_path("coarse_quad_mesh_british_museum_poles.json"))
payload = wire.mesh_to_wire(poled)
check("poles are on the wire", len(payload["poles"]) == len(list(poled.poles())),
      len(payload["poles"]))
check("as coordinates, not keys",
      all(len(p) == 3 and isinstance(p[0], float) for p in payload["poles"]))
back = wire.mesh_from_wire(payload)
check("a payload with poles rebuilds as a pseudo-quad mesh",
      type(back).__name__ == "CoarsePseudoQuadMesh", type(back).__name__)
check("and keeps them all", len(list(back.poles())) == len(list(poled.poles())))

print("\n=== 3. the payload is plain JSON: no compas type crosses ===")
text = json.dumps(payload)
check("it serialises with the stdlib alone", isinstance(text, str))
check("and comes back identical", json.loads(text) == payload)
check("with no type markers in it", "dtype" not in text and "compas" not in text)

print("\n=== 4. a malformed payload fails where the fault is ===")
for bad, want in (
        ({"vertices": [[0, 0, 0]], "faces": [[0, 1, 2, 3]]}, "refers to vertex"),
        ({"vertices": [], "faces": []}, "no vertices"),
        ({"faces": []}, "needs"),
        ({"vertices": [[0, 0, 0]], "faces": [[0, 1]]}, "at least 3"),
        ("nope", "not a dict")):
    try:
        wire.mesh_from_wire(bad)
        check("refuses {}".format(want), False, "it did not raise")
    except wire.WireError as exc:
        check("refuses: {}".format(want), want in str(exc), str(exc)[:56])

print("\n=== 5. curves coerce from whatever the library hands around ===")
from compas.geometry import Point, Polyline
check("a Polyline", wire.points_of(Polyline([[0, 0, 0], [1, 1, 0]]))
      == [[0.0, 0.0, 0.0], [1.0, 1.0, 0.0]])
check("a list of Points", wire.points_of([Point(0, 0, 0), Point(1, 1, 0)])
      == [[0.0, 0.0, 0.0], [1.0, 1.0, 0.0]])
check("bare triples", wire.points_of([[0, 0, 0], [1, 1, 0]])
      == [[0.0, 0.0, 0.0], [1.0, 1.0, 0.0]])
check("several at once", len(wire.curves_to_wire([[[0, 0, 0], [1, 0, 0]]])) == 1)
check("none at all", wire.curves_to_wire(None) == [])

sys.exit(report("test_wire"))
