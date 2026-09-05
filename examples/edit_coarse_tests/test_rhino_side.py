"""Exercise the RHINO half without opening Rhino, by stubbing rhinoscriptsyntax.

Run under Rhino 8's own CPython -- same interpreter as in-app, so an import or
a name that resolves here resolves there. This catches the class of bug that
killed CoarseLayoutObject before: a name bound inside a try/except ImportError
block that only fails at use.
"""
import os
import sys
import types

REPO = r"C:\Users\Casper\libraries\carbcomn\compas_singular\compas_singular"
COMMANDS = r"C:\Users\Casper\libraries\carbcomn\compas_singular\compas_singular\rhino_plugin\commands"
for path in (os.path.join(REPO, "src"), os.path.join(REPO, "examples", "New approach"),
             COMMANDS):
    if path not in sys.path:
        sys.path.insert(0, path)

FAILURES = []


def check(name, condition, detail=""):
    print("  {} {}{}".format("ok  " if condition else "FAIL", name,
                             (" -- " + str(detail)) if detail else ""))
    if not condition:
        FAILURES.append(name)


# ----------------------------------------------------------------------
print("\n=== 1. the editing core imports with NO Rhino present ===")
try:
    from compas_singular.editing import CoarseLayoutEditor
    ok = True
    err = ""
except Exception as e:                                   # noqa: BLE001
    ok, err = False, "{}: {}".format(type(e).__name__, e)
check("compas_singular.editing imports", ok, err)
check("CoarseLayoutEditor is a class", ok and isinstance(CoarseLayoutEditor, type))

print("\n=== 2. mesh_ui imports without Rhino and says so on use ===")
from compas_singular.rhino import mesh_ui
check("mesh_ui imports outside Rhino", True)
check("mesh_ui knows Rhino is absent", mesh_ui.RHINO is False, mesh_ui.RHINO)
try:
    mesh_ui.ensure_layer("X")
    raised = ""
except RuntimeError as e:
    raised = str(e)
check("it raises a message that names the alternative",
      "compas_singular.editing" in raised, raised[:70])

# ----------------------------------------------------------------------
print("\n=== 4. now WITH a stubbed Rhino, every mesh_ui name resolves ===")


class _Filter(object):
    point = 1
    curve = 4


class _StubRs(object):
    """Enough rhinoscriptsyntax to walk the code paths, and no more."""

    filter = _Filter()

    def __init__(self):
        self.layers = set()
        self.objects = {}
        self.messages = []
        self._next = 0

    # layers
    def IsLayer(self, name):
        return name in self.layers

    def AddLayer(self, name=None, parent=None, color=None):
        full = "{}::{}".format(parent, name) if parent else name
        self.layers.add(full)
        return full

    def LayerLocked(self, name, value=None):
        return False

    # objects
    def _add(self, kind):
        self._next += 1
        guid = "{}-{}".format(kind, self._next)
        self.objects[guid] = kind
        return guid

    def AddPoint(self, point):
        return self._add("point")

    def AddLine(self, a, b):
        return self._add("line")

    def AddPolyline(self, points):
        return self._add("polyline")

    def IsObject(self, guid):
        return guid in self.objects

    def IsObjectLocked(self, guid):
        return False

    def UnlockObject(self, guid):
        pass

    def LockObject(self, guid):
        pass

    def ObjectLayer(self, guid, layer=None):
        return layer

    def ObjectColor(self, guid, color=None):
        return color

    def DeleteObjects(self, guids):
        for guid in guids:
            self.objects.pop(guid, None)

    def EnableRedraw(self, value):
        pass

    def MessageBox(self, message, style, title):
        self.messages.append((title, message))
        return 1

    def GetString(self, message, default=None, options=None):
        return default


class _StubSc(object):
    class _Views(object):
        def Redraw(self):
            pass

    class _Doc(object):
        ModelAbsoluteTolerance = 0.001
        Views = None

    def __init__(self):
        self.doc = self._Doc()
        self.doc.Views = self._Views()


class _Point3d(object):
    """Stands in for Rhino.Geometry.Point3d, which is a .NET type.

    It has to be stubbed too: mesh_ui binds it inside the same
    ``try: import ... except ImportError`` block as ``rs``, so replacing only
    ``rs`` leaves it undefined -- which is the exact failure mode this file
    exists to catch, and it caught it.
    """

    def __init__(self, x, y, z=0.0):
        self.X, self.Y, self.Z = x, y, z

    def DistanceTo(self, other):
        return ((self.X - other.X) ** 2 + (self.Y - other.Y) ** 2) ** 0.5


stub_rs = _StubRs()
mesh_ui.rs = stub_rs
mesh_ui.sc = _StubSc()
mesh_ui.Point3d = _Point3d
mesh_ui.RHINO = True

from math import cos, pi, sin                                            # noqa: E402
from compas_singular.framefield.decomposition import FieldDecomposition  # noqa: E402

square = [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 10.0, 0.0], [0.0, 10.0, 0.0],
          [0.0, 0.0, 0.0]]
circle = [[5.0 + 1.5 * cos(2 * pi * i / 24), 5.0 + 1.5 * sin(2 * pi * i / 24), 0.0]
          for i in range(24)]
circle = circle + circle[:1]

d = FieldDecomposition.from_boundary(square, [circle], target_length=1.0)
d.decomposition_mesh()
editor = CoarseLayoutEditor(d)

scene = mesh_ui.PickableMesh("TopologyProblem::Skeleton::TempEdit")
scene.draw(editor.mesh, edge_shape=editor.edge_shape, special=editor.mesh.poles())
check("draw made one object per corner",
      len(scene.guid_vertices) == editor.mesh.number_of_vertices(),
      "{} vs {}".format(len(scene.guid_vertices), editor.mesh.number_of_vertices()))
check("draw made one object per edge",
      len(scene.guid_edges) == editor.mesh.number_of_edges(),
      "{} vs {}".format(len(scene.guid_edges), editor.mesh.number_of_edges()))
check("the layer path was created",
      "TopologyProblem::Skeleton::TempEdit" in stub_rs.layers, stub_rs.layers)

state = scene.unlock()
scene.relock(state)
check("unlock/relock round trip", True)

drawn = len(stub_rs.objects)
scene.clear()
check("clear removes everything it drew", len(stub_rs.objects) == 0, drawn)

mesh_ui.refuse("Add line", "nope")
check("refuse goes to a dialog, not a print", stub_rs.messages == [("Add line", "nope")],
      stub_rs.messages)

answer = mesh_ui.ask("next", ["Yes", "No"], "Yes")
check("ask lowercases the answer", answer == "yes", answer)

print("\n=== 5. edge_shape drives the curved-edge draw ===")
# An edge carrying a shape must come back as a POLYLINE, not a line -- drawing
# an added arc as its chord makes a cut that worked look like one that snapped
# straight. Registering the shape directly, because a cut across a straight
# square produces none and the check would then pass without exercising it.
u, v = list(editor.mesh.edges())[0]
pa = editor.mesh.vertex_coordinates(u)
pb = editor.mesh.vertex_coordinates(v)
bulge = [(pa[0] + pb[0]) / 2.0, (pa[1] + pb[1]) / 2.0 + 0.4, 0.0]
editor.curves[editor._curve_key(pa, pb)] = [list(pa), bulge, list(pb)]

check("edge_shape finds the registered shape", editor.edge_shape(u, v) is not None)
check("edge_shape finds it the other way round too",
      editor.edge_shape(v, u) is not None)
shaped = [(a, b) for a, b in editor.mesh.edges()
          if editor.edge_shape(a, b) is not None]
check("exactly one edge is shaped", len(shaped) == 1, len(shaped))

scene.draw(editor.mesh, edge_shape=editor.edge_shape)
kinds = sorted(stub_rs.objects[g] for g in scene.guid_edges)
check("the shaped edge is drawn as a polyline, the rest as lines",
      kinds.count("polyline") == 1 and set(kinds) == {"line", "polyline"},
      dict((k, kinds.count(k)) for k in set(kinds)))
scene.clear()

print("\n" + "=" * 60)
if FAILURES:
    print("FAILED ({}): {}".format(len(FAILURES), ", ".join(FAILURES)))
    sys.exit(1)
print("ALL CHECKS PASSED")
