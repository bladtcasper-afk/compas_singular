#! python3
# r: compas
# r: google-genai
"""**Ask an AI to do the meshing work, from inside Rhino.**

    reads   Outer / Inner / Guides                   the domain, as every step does
            Mesh + Poles                             the coarse layout, if started there
            TopologyProblem::QuadMesh                the dense mesh, if started there
    writes  TopologyProblem::QuadMesh                the result
            Mesh + Poles + Skeleton::Polylines       the layout, if one was committed
            TopologyProblem::QuadMesh::AI::Before    the mesh it started from
    prints  every tool call as it happens

**ASKS WHERE TO START.** ``Field`` generates a layout from the field, ignoring
whatever is drawn -- what to use on a fresh document. ``Coarse`` picks up the
layout on ``Mesh``, so the model continues from the patches you have already
edited. ``Dense`` picks up the mesh on ``QuadMesh`` and edits it directly; note
a dense mesh has no layout that reproduces it, so the model cannot go back to
the coarse stage from there.

This is a DRIVER, not an implementation. The session, the tools, the state
machine and the loop all live in ``compas_singular.agent`` and are tested
without Rhino in ``examples/agent_tests``. This file supplies the document, the
prompts and the bake, exactly as ``CMD_edit_coarse_mesh`` supplies picks for
``CoarseEditor``. Nothing here decides anything about a mesh.

**The seam is the decomposition object, not geometry.** ``FieldDecomposition``
is Rhino-free, so it crosses into the agent layer unchanged and no Rhino type
ever reaches it. That is why the same run is reproducible from a terminal with
``python -m compas_singular.agent.cli`` -- same session, same tools, different
front end.

**FOUR TOOLS DESTROY WORK, AND YOU ARE ASKED ABOUT EACH ONE.**
``solve_field`` and ``set_guides`` discard the layout, its edits, the densities
and the dense mesh, because a re-solve returns a NEW decomposition and nothing
below it survives. ``commit`` welds and renumbers. ``revert_to_coarse`` throws
away every dense edit. The command runs with a confirmation callback, so each
one becomes a dialog; answering No is not a failure -- the refusal goes back to
the model with the reason and it carries on without that tool.

**Nothing is written to the document until the run finishes.** The mesh the run
started from is kept on ``::AI::Before`` so the result can be compared with it,
and a run that produces nothing usable leaves ``QuadMesh`` untouched.

**TWO THINGS TO SET, BOTH AT THE TOP OF THIS FILE.** ``BACKEND`` picks the
model provider and ``API_KEY_FILE`` says where the key is. Everything else is
already configured.

**The key comes from a FILE, not the environment, and that is not a
preference.** Rhino reads its environment once, when it starts, so a variable
exported in a terminal never reaches a running Rhino and changing a key would
mean restarting the application. A path in a text file can be edited between
runs. Keep that file OUTSIDE both repositories -- ``compas_singular`` pushes to
a public fork, and a key committed once is a key that has to be revoked.

**The SDK has to be in RHINO's interpreter**, which is a different install from
the conda env the tests use. The ``# r:`` lines at the top of this file are what
put it there: Rhino's script editor resolves them into the site-env on first
run. ``# r: google-genai`` covers the Gemini backend; add ``# r: anthropic``
beside it to use ``BACKEND = "anthropic"``. If the SDK is missing the command
says so and stops without touching anything.
"""

# Temporary import of compas_singular development library. MUST come before any
# compas_singular import: it fixes sys.path and purges a stale copy from
# sys.modules. ``compas_singular.agent`` is purged with everything else because
# the purge matches the whole package -- which is exactly why the agent layer
# lives under ``src/compas_singular`` rather than beside it.
from CMD_start import import_compas_singular
import_compas_singular()

import rhinoscriptsyntax as rs

from compas_singular.datastructures import QuadMesh
from compas_singular.rhino.helpers.helpers import bake_mesh
from compas_singular.rhino.helpers.helpers import bake_points
from compas_singular.rhino.helpers.helpers import bake_polylines
from compas_singular.rhino.helpers.helpers import clear_layer
from compas_singular.rhino.helpers.helpers import read_boundaries
from compas_singular.rhino.helpers.helpers import read_coarse
from compas_singular.rhino.helpers.helpers import read_mesh
from CMD_start import ensure_paths
from CMD_start import get_decomposition
from CMD_start import get_settings


# ======================================================================
# SET THESE TWO
# ======================================================================

#: WHICH MODEL PROVIDER. "gemini" or "anthropic". Whichever you pick, its SDK
#: must be on a ``# r:`` line at the top of this file so Rhino installs it.
BACKEND = "gemini"

#: WHERE YOUR API KEY IS -- a plain text file OUTSIDE this repository, holding
#: the key on one line. ``KEY=value`` and ``#`` comment lines are both fine.
#: Set to None to use the environment instead (GEMINI_API_KEY / GOOGLE_API_KEY,
#: or ANTHROPIC_API_KEY), but note Rhino only sees variables that existed when
#: it started.
API_KEY_FILE = r"C:\Users\Casper\.secrets\gemini_api_key.txt"

#: Model name, or None for the backend's default -- Gemini Flash, or Claude
#: Opus. On the Gemini free tier the NEWEST model is the busiest and is the
#: first to answer 503 "high demand", so if runs keep stalling set this to an
#: older one: "gemini-3.6-flash", "gemini-3.5-flash" or "gemini-2.5-flash" are
#: all reachable and all capable of this work. Leaving it None is fine -- the
#: backend retries with backoff and then switches down the list on its own.
MODEL = None

# ======================================================================

ROOT = "TopologyProblem"
QUADMESH_LAYER = ROOT + "::QuadMesh"
AI_LAYER = QUADMESH_LAYER + "::AI"
BEFORE_LAYER = AI_LAYER + "::Before"

#: The COARSE layout, as CMD_edit_coarse_mesh writes it and CMD_quad_mesh reads
#: it back. Baked again here when the run changed the layout, or the next
#: command would regenerate from a stale one and lose the edit.
MESH_LAYER = "Mesh"
POLES_LAYER = "Poles"
POLYLINE_LAYER = "Skeleton::Polylines"

#: Steps before the loop stops itself. A step is one turn, which may carry
#: several tool calls.
MAX_STEPS = 40


def ensure_layer(path, color=None):
    """Create a ``::`` layer path, parents first, and return the full path.

    ``rs.AddLayer`` does not create intermediate parents, and ``bake_mesh``'s own
    layer creation hard-codes ``parent="TopologyProblem"``, which is wrong for a
    layer two levels down.

    **``name`` is the LEAF, never the full path.** Passing the ``::`` path as the
    name AND a parent makes Rhino create a layer literally called
    "TopologyProblem::QuadMesh" nested inside "TopologyProblem" -- duplicated
    layers at every level, and the real path never comes into existence. The
    colour goes through ``AddLayer`` for the same reason: a later
    ``rs.LayerColor(path, ...)`` raises ``ValueError: ... does not exist in
    LayerTable`` on a path that was never made.

    This is the same implementation as ``CMD_edit_quad_mesh.ensure_layer``.
    """
    parts = path.split("::")
    for i in range(len(parts)):
        name = "::".join(parts[:i + 1])
        if not rs.IsLayer(name):
            rs.AddLayer(name=parts[i],
                        parent="::".join(parts[:i]) if i else None,
                        color=color if i == len(parts) - 1 else None)
    return path


def confirm(name, arguments):
    """One dialog per destructive tool. ``True`` allows it."""
    what = {
        "solve_field": "Re-solve the field.\n\nThis DISCARDS the coarse layout, "
                       "every edit on it, the densities and the dense mesh -- a "
                       "re-solve returns a new decomposition and nothing below "
                       "it survives.",
        "set_guides": "Replace the guide curves and re-solve.\n\nSame as "
                      "re-solving: everything below the field is discarded.",
        "commit": "Commit the edited layout.\n\nThis is what makes the edit keep "
                  "its field alignment. It also renumbers every key, so labels "
                  "are reissued afterwards. It cannot be undone.",
        "revert_to_coarse": "Go back to the coarse layout.\n\nThis DESTROYS "
                            "every edit made to the dense mesh. A hand-edited "
                            "dense mesh has no layout that reproduces it.",
    }.get(name, "{} destroys work that cannot be recovered.".format(name))

    detail = ""
    if arguments:
        detail = "\n\narguments:\n" + "\n".join(
            "  {} = {}".format(k, v) for k, v in sorted(arguments.items()))

    print("  asking about {}".format(name))
    # 4 = Yes/No, 32 = question icon. 6 is Yes.
    answer = rs.MessageBox(what + detail + "\n\nAllow it?", 4 | 32,
                           "AI edit: {}".format(name))
    allowed = answer == 6
    print("    {}".format("allowed" if allowed else "refused by you"))
    return allowed


def on_event(event):
    """Progress, on the command line and in the prompt."""
    kind = event["type"]
    if kind == "tool":
        print("  {} {}".format("." if event.get("ok") else "x", event["message"]))
        rs.Prompt(event["message"][:80])
    elif kind == "refused":
        print("  ! {}".format(event["message"]))
    elif kind == "text":
        print("\n{}\n".format(event["message"]))
    elif kind == "done":
        rs.Prompt("")


# ----------------------------------------------------------------------
# run
# ----------------------------------------------------------------------

def main():
    # Late imports, so a missing SDK is a sentence rather than a traceback, and
    # so this module can be imported for its helpers without paying for them.
    # ``ensure_paths`` first: Rhino re-initialises sys.path between runs while
    # keeping sys.modules warm, so a late import on a second run fails with a
    # bare ModuleNotFoundError that looks nothing like a path problem.
    ensure_paths()
    from compas_singular.agent.core import MeshEditSession
    from compas_singular.agent.runner import run_session

    from compas_singular.agent.runner import KeyNotFound

    try:
        if BACKEND == "gemini":
            from compas_singular.agent.runner import GeminiBackend
            backend = GeminiBackend(model=MODEL, key_file=API_KEY_FILE)
        elif BACKEND == "anthropic":
            from compas_singular.agent.runner import AnthropicBackend
            backend = AnthropicBackend(model=MODEL, key_file=API_KEY_FILE)
        else:
            rs.MessageBox(
                "BACKEND is set to {!r}, which is not a backend. Edit the top "
                "of CMD_ai_edit.py and set it to \"gemini\" or "
                "\"anthropic\".\n\nNothing was changed.".format(BACKEND),
                0 | 16, "AI edit")
            return
    except KeyNotFound as exc:
        # A wrong PATH and a wrong KEY are different problems with different
        # fixes, so they get different dialogs.
        rs.MessageBox(
            "The API key file could not be read.\n\n{}\n\nAPI_KEY_FILE is set "
            "at the top of CMD_ai_edit.py and currently points at:\n  {}\n\n"
            "Put the key on one line in that file, or change the path. "
            "Nothing was changed.".format(exc, API_KEY_FILE),
            0 | 16, "AI edit")
        return
    except ImportError as exc:
        package = "google-genai" if BACKEND == "gemini" else "anthropic"
        rs.MessageBox(
            "{}\n\nThis command needs the {} SDK inside RHINO's own Python, "
            "which is a separate install from the conda environment the tests "
            "use. Add a line\n\n    # r: {}\n\nto the top of this file and run "
            "it again -- Rhino resolves those into its site-env.\n\nNothing "
            "was changed.".format(exc, package, package), 0 | 16, "AI edit")
        return
    except Exception as exc:
        rs.MessageBox("Could not start the model client.\n\n{}: {}\n\nNothing "
                      "was changed.".format(type(exc).__name__, exc),
                      0 | 16, "AI edit")
        return

    instruction = rs.StringBox(
        "What should the AI do to the mesh?",
        "Densify the layout and even out the elements around the holes.",
        "AI edit")
    if not instruction:
        print("Cancelled -- nothing was changed.")
        return

    # WHERE TO START. The three states the session has, offered as the three
    # things actually present in a document at this point in the workflow.
    start = (rs.GetString(
        "Start from what?", "Coarse",
        ["Field", "Coarse", "Dense"]) or "").lower()
    if not start:
        print("Cancelled -- nothing was changed.")
        return

    settings = get_settings()

    print("reading the document...")
    outer, inners, guides, _poles = read_boundaries(
        spacing=settings["triangulation_spacing"])

    print("getting the field (cached unless the document changed)...")
    decomposition = get_decomposition()

    session = MeshEditSession(
        decomposition, outer=outer, inners=inners, guides=guides,
        solve_params={"target_length": settings["triangulation_spacing"]},
        target_length=settings["target_length"])

    # ``Field`` leaves the session where it is: the layout is generated from the
    # field when the model asks for it, which is what happens with no document
    # mesh at all. The other two adopt what is already drawn.
    if start.startswith("c"):
        try:
            coarse, poles = read_coarse()
        except Exception as exc:
            rs.MessageBox(
                "No coarse layout to start from.\n\n{}\n\nRun CMD_coarse_mesh "
                "and CMD_edit_coarse_mesh first, or start from Field.\n\n"
                "Nothing was changed.".format(exc), 0 | 16, "AI edit")
            return
        # The same two lines CMD_quad_mesh uses to adopt a document layout.
        # ``strict=False`` because a layout a person has been editing may not
        # densify yet, and refusing to start on it would be unhelpful -- the
        # model can be asked to repair it.
        decomposition.edit_coarse(coarse, poles=poles or None, strict=False)
        session.enter_coarse()
        print("starting from the coarse layout on '{}': {} patch(es)".format(
            MESH_LAYER, session.mesh.number_of_faces()))

    elif start.startswith("d"):
        try:
            dense = read_mesh(QUADMESH_LAYER)
        except Exception as exc:
            rs.MessageBox(
                "No dense mesh to start from.\n\n{}\n\nRun CMD_quad_mesh "
                "first, or start from Coarse or Field.\n\nNothing was "
                "changed.".format(exc), 0 | 16, "AI edit")
            return
        # ``mesh_to_compas`` gives a plain Mesh; the strip and polyedge grammar
        # lives on QuadMesh, so it is rebuilt. Rhino stores a triangle as a quad
        # with a repeated corner, so the duplicate is dropped -- otherwise the
        # face is degenerate and every valence through it is wrong.
        vertices, faces = dense.to_vertices_and_faces()
        faces = [[v for i, v in enumerate(f) if v != f[i - 1]] for f in faces]
        session.adopt_dense(QuadMesh.from_vertices_and_faces(vertices, faces))
        print("starting from the dense mesh on '{}': {} faces".format(
            QUADMESH_LAYER, session.mesh.number_of_faces()))
        print("  note: a dense mesh has no layout that reproduces it, so the "
              "model cannot go back to the coarse stage from here.")

    # This session's starting point, so the result can be compared with it.
    # Refreshed every run rather than written once: a backup kept from before a
    # regeneration would be a mesh of a different layout, and stale is worse
    # than absent.
    if session.decomposition.dense is not None:
        ensure_layer(BEFORE_LAYER, (170, 170, 170))
        bake_mesh(session.decomposition.dense, BEFORE_LAYER)
        rs.LayerVisible(BEFORE_LAYER, False)
        print("starting mesh kept on '{}'".format(BEFORE_LAYER))

    # ``effort`` is an Anthropic knob and GeminiBackend has none, so it is read
    # defensively rather than assumed -- this line used to raise AttributeError
    # before a Gemini run had started.
    effort = getattr(backend, "effort", None)
    print("\nmodel: {} ({}){}. Up to {} steps.".format(
        backend.model, backend.name,
        ", effort {}".format(effort) if effort else "", MAX_STEPS))
    print("You will be asked before any tool that destroys work.\n")

    try:
        result = run_session(session, instruction, backend=backend,
                             confirm=confirm, max_steps=MAX_STEPS,
                             images="transitions", on_event=on_event)
    except Exception as exc:
        name = type(exc).__name__
        blob = "{}: {}".format(name, exc)
        if any(marker in blob for marker in (
                "Authentication", "authentication_error", "API key", "api_key",
                "PERMISSION_DENIED", "UNAUTHENTICATED", "401", "403")):
            where = ("the file API_KEY_FILE points at:\n  {}".format(API_KEY_FILE)
                     if API_KEY_FILE else
                     "the environment Rhino was started from -- note Rhino only "
                     "sees variables that existed when it launched")
            rs.MessageBox(
                "The API rejected the credentials.\n\nThe key was read from "
                "{}\n\nCheck the key itself is valid and has quota. Nothing "
                "was changed.\n\n{}".format(where, blob[:200]), 0 | 16,
                "AI edit")
        else:
            rs.MessageBox("The run failed.\n\n{}\n\nNothing was "
                          "changed.".format(blob), 0 | 16, "AI edit")
        return

    print("\n" + result.summary())
    # If the backend had to move off a saturated model, say which one answered.
    switched = getattr(backend, "switched_from", None)
    if switched:
        print("  {} was unavailable; the run finished on {}".format(
            ", ".join(switched), backend.model))
    for entry in result.calls(ok=False):
        print("  refused: {} -- {}".format(entry["name"], entry["reason"][:110]))

    # WHAT TO DELIVER is not the same question as WHAT STAGE IT ENDED IN. A run
    # that densifies well and then steps back to the layout on its last turn --
    # measured, after 37 tool calls of density work -- used to bake nothing at
    # all. ``last_dense`` survives a revert precisely so that result is not lost.
    deliver = result.mesh if result.stage == "dense" else None
    if deliver is None and getattr(session, "last_dense", None) is not None:
        deliver = session.last_dense
        print("the run ended at stage '{}', so the last dense mesh it produced "
              "is what gets baked".format(result.stage))

    if deliver is None:
        rs.MessageBox(
            "The run finished at stage '{}' and never produced a dense mesh, so "
            "nothing was baked.\n\n{}\n\nWhat it said:\n{}".format(
                result.stage, result.summary(), (result.text or "")[:400]),
            0 | 48, "AI edit")
        return

    quality = (result.digest.get("dense") or {}).get("quality") or {}
    floor = (result.digest.get("dense") or {}).get("hard_floor") or {}
    if quality:
        print("quality: min angle {:.1f}, max angle {:.1f}, aspect {:.2f}".format(
            quality.get("min_angle", 0.0), quality.get("max_angle", 0.0),
            quality.get("aspect_max", 0.0)))
    if not floor.get("ok", True):
        print("QUALITY FAILURE: {}".format(floor.get("reason")))
    for warning in result.digest.get("warnings") or []:
        print("warning: {}".format(warning))

    # ensure_layer FIRST, at every bake site. The helpers create a missing layer
    # themselves and both ways of doing it are wrong for a nested path:
    # ``bake_mesh`` hard-codes ``parent="TopologyProblem"``, and
    # ``bake_points``/``bake_polylines`` pass the whole "::" path as a layer
    # NAME. Neither shows up until the layer happens not to exist yet.
    ensure_layer(QUADMESH_LAYER)
    clear_layer(QUADMESH_LAYER, clean_sublayers=False)
    guid = bake_mesh(deliver, QUADMESH_LAYER)

    # VERIFY, rather than trust the call returned. bake_mesh raises only when
    # AddMesh fails outright; it says nothing about whether the object ended up
    # somewhere you can see. A layer left switched off -- by an earlier command,
    # or by a hidden sublayer's parent -- makes a perfectly good bake invisible,
    # which is indistinguishable from no mesh at all.
    rs.LayerVisible(QUADMESH_LAYER, True)
    if rs.IsLayerLocked(QUADMESH_LAYER):
        rs.LayerLocked(QUADMESH_LAYER, False)
    on_layer = rs.ObjectsByLayer(QUADMESH_LAYER) or []
    print("baked {} faces to '{}' -- {} object(s) now on that layer".format(
        deliver.number_of_faces(), QUADMESH_LAYER, len(on_layer)))
    if not on_layer:
        print("  WARNING: the bake returned {} but the layer is empty. Check "
              "whether another layer of the same name exists.".format(guid))
    rs.SelectObject(guid)
    rs.Command("_Zoom _Selected", echo=False)
    rs.UnselectAllObjects()

    # THE LAYOUT TOO, when the run committed one. Without this the document is
    # left inconsistent: 'Mesh' still holds the layout from before the edit, and
    # the next CMD_quad_mesh reads THAT and regenerates -- silently throwing the
    # model's layout work away while the dense mesh on screen still shows it.
    # The three layers written here are exactly the three CMD_quad_mesh reads.
    committed = session.coarse_editor.mesh if session.coarse_editor else None
    if committed is not None and decomposition.edit_notes:
        ensure_layer(MESH_LAYER)
        ensure_layer(POLES_LAYER)
        ensure_layer(POLYLINE_LAYER)
        bake_mesh(committed, MESH_LAYER)
        pole_keys = committed.poles() if hasattr(committed, "poles") else []
        bake_points([committed.vertex_coordinates(v) for v in pole_keys],
                    POLES_LAYER)
        # ``decomposition.polylines``, NOT decomposition_polylines(): that
        # method reassigns self.polylines and drops the user curves ``commit``
        # published, so a line the model drew would densify as a chord.
        ribs = decomposition.polylines or decomposition.decomposition_polylines()
        _guids, skipped = bake_polylines(ribs, POLYLINE_LAYER)
        print("layout committed: baked {} patch(es) to '{}', {} pole(s), "
              "{} separatrix/ices".format(committed.number_of_faces(),
                                          MESH_LAYER, len(pole_keys),
                                          len(_guids)))
        if skipped:
            # Not fatal, and it must not be: this runs AFTER the meshes are
            # written, so raising would leave the document half updated.
            print("  {} rib(s) skipped -- Rhino refused them".format(skipped))
    elif committed is not None:
        print("layout unchanged -- '{}' left as it was".format(MESH_LAYER))

    summary = "{}\n\n{} faces".format(result.summary(),
                                      deliver.number_of_faces())
    if quality:
        summary += ", min angle {:.1f} deg, aspect {:.2f}".format(
            quality.get("min_angle", 0.0), quality.get("aspect_max", 0.0))
    if not floor.get("ok", True):
        summary += "\n\nQUALITY FAILURE: {}".format(floor.get("reason"))
    if result.text:
        summary += "\n\nWhat it said:\n{}".format(result.text[:600])
    rs.MessageBox(summary, 0 | 64, "AI edit")


# ``main()`` is called only under this guard because every other CMD_ file that
# imports from this one would otherwise run it on import -- the same reason
# CMD_start guards ``reset_project``.
if __name__ == "__main__":
    main()
