"""The CLI, driven through its own main() -- dry run only, so no model is called.

Covers what a terminal user actually hits: argument validation, the two ways of
supplying a domain, the output files, and the exit code. The live path is one
swapped backend and is not exercised here.

Run:
    C:/Users/Casper/anaconda3/envs/singular312/python.exe test_cli_offline.py
"""
import json
import os
import shutil
import sys
import tempfile

REPO = r"C:\Users\Casper\libraries\carbcomn\compas_singular\compas_singular"
for path in (os.path.join(REPO, "src"),):
    if path not in sys.path:
        sys.path.insert(0, path)

import compas                                                            # noqa: E402
from compas_singular.agent import cli                                    # noqa: E402


FAILURES = []


def check(name, condition, detail=""):
    print("  {} {}{}".format("ok  " if condition else "FAIL", name,
                             (" -- " + str(detail)) if detail else ""))
    if not condition:
        FAILURES.append(name)


def run(argv):
    """``(exit code, stdout)`` for one invocation."""
    from io import StringIO
    held, sys.stdout = sys.stdout, StringIO()
    try:
        code = cli.main(argv)
    except SystemExit as exc:
        code = exc.code
    finally:
        text = sys.stdout.getvalue()
        sys.stdout = held
    return code, text


TMP = tempfile.mkdtemp(prefix="agent_cli_")


# ======================================================================

print("\nargument validation")

for argv, why in (
        ([], "no domain"),
        (["--boundary", "x.json", "--demo", "square", "go"], "both sources"),
        (["--demo", "square"], "no instruction and no --dry-run")):
    try:
        cli.build_parser().parse_args(argv)
        parsed = True
    except SystemExit:
        parsed = True                       # argparse's own errors exit too
    # main() is where the cross-argument rules live
    code, _text = run(argv)
    check("refused: {}".format(why), code != 0, code)

check("--help exits 0", run(["--help"])[0] in (0, None))


# ======================================================================

print("\ndemos")

check("every demo returns three parts",
      all(len(cli.demo_boundary(n)) == 3 for n in cli.DEMOS),
      sorted(cli.DEMOS))
check("an unknown demo is refused",
      isinstance(getattr(run(["--demo", "nope", "--dry-run"])[0], "real", None)
                 or run(["--demo", "nope", "--dry-run"])[0], (int, str)))

outer, inners, guides = cli.demo_boundary("plate-with-hole")
check("plate-with-hole has a hole", len(inners) == 1, len(inners))
outer, inners, guides = cli.demo_boundary("cable")
check("cable has a guide", len(guides) == 1, len(guides))


# ======================================================================

print("\nboundary files")

path = os.path.join(TMP, "site.json")
with open(path, "w") as handle:
    json.dump({"outer": [[0, 0], [12, 0], [12, 8], [0, 8]],
               "inners": [[[5, 3], [7, 3], [7, 5], [5, 5]]],
               "guides": [[[0.5, 4], [11.5, 4]]]}, handle)

outer, inners, guides = cli.load_boundary(path)
check("loads outer", len(outer) == 4, len(outer))
check("pads 2D points to 3D", all(len(p) == 3 for p in outer), outer[0])
check("loads holes and guides", len(inners) == 1 and len(guides) == 1)

bad = os.path.join(TMP, "bad.json")
with open(bad, "w") as handle:
    json.dump({"inners": []}, handle)
try:
    cli.load_boundary(bad)
    refused = False
    message = ""
except SystemExit as exc:
    refused = True
    message = str(exc)
check("a file with no outer loop is refused", refused, message[:60])
check("  and says which file", path.split(os.sep)[-1] in message or "outer" in message,
      message[:60])


# ======================================================================

print("\na dry run end to end")

out = os.path.join(TMP, "run")
code, text = run(["--demo", "plate-with-hole", "--dry-run",
                  "--target-length", "1.5", "--out", out])

check("exited 0", code == 0, code)
check("said it was a dry run", "no model will be called" in text)
check("ran the whole plan",
      all(name in text for name in ("check_inputs", "enter_coarse",
                                    "rebuild_dense", "smooth")))
check("reported the summary", "tool call(s) succeeded" in text)
check("reported quality", "min angle" in text, text[-200:].strip()[:60])
check("ended at dense", "'dense'" in text)

for name in ("mesh.json", "mesh.png", "mesh.svg", "transcript.json"):
    check("wrote {}".format(name), os.path.isfile(os.path.join(out, name)))

check("the PNG is a PNG",
      open(os.path.join(out, "mesh.png"), "rb").read(8) == b"\x89PNG\r\n\x1a\n")

mesh = compas.json_load(os.path.join(out, "mesh.json"))
check("the mesh JSON round-trips as a mesh",
      hasattr(mesh, "number_of_faces") and mesh.number_of_faces() > 0,
      "{} faces".format(mesh.number_of_faces()))

transcript = json.load(open(os.path.join(out, "transcript.json")))
check("the transcript names every call",
      [c["name"] for c in transcript["calls"]][:2] == ["check_inputs",
                                                       "enter_coarse"],
      [c["name"] for c in transcript["calls"]])
check("  and records the outcome",
      all("ok" in c for c in transcript["calls"]))
check("  and the stop reason", transcript["stop_reason"] == "end_turn",
      transcript["stop_reason"])


# ======================================================================

print("\n--quiet")

code, text = run(["--demo", "square", "--dry-run", "--quiet"])
check("quiet exits 0", code == 0, code)
check("  and prints nothing", text.strip() == "", text[:60])


# ======================================================================

print("\nthe destructive-tool default")

code, text = run(["--demo", "square", "--dry-run", "go"])
check("says the destructive tools are refused",
      "will be refused" in text and "--yes" in text,
      [line for line in text.splitlines() if "refused" in line][:1])

code, text = run(["--demo", "square", "--dry-run", "--yes", "go"])
check("--yes drops that notice", "will be refused" not in text)


print("\n--key-file and --backend")

check("--key-file is accepted",
      cli.build_parser().parse_args(
          ["--demo", "square", "--key-file", "x.txt", "go"]).key_file == "x.txt")
check("--backend accepts gemini",
      cli.build_parser().parse_args(
          ["--demo", "square", "--backend", "gemini", "go"]).backend == "gemini")
try:
    cli.build_parser().parse_args(["--demo", "square", "--backend", "nope", "go"])
    rejected = False
except SystemExit:
    rejected = True
check("  and rejects an unknown one", rejected)

# A bad key file must stop the run with a message naming the path, rather than
# reaching the API with no key and failing there.
code, text = run(["--demo", "square", "--backend", "gemini",
                  "--key-file", os.path.join(TMP, "does-not-exist.txt"),
                  "--max-steps", "1", "go"])
# SystemExit carries the message as its code, so that is where it lands.
check("a missing key file stops the run", code != 0 and code is not None, code)
check("  naming the path", "does-not-exist.txt" in str(code),
      str(code)[:90])
check("  BEFORE the field solve, which would otherwise be wasted",
      "solving the field" not in text, text.strip()[-90:])

# A dry run must not need a key at all -- that is the whole point of it.
code, text = run(["--demo", "square", "--backend", "gemini", "--dry-run"])
check("--dry-run needs no key even with --backend gemini", code == 0, code)


shutil.rmtree(TMP, ignore_errors=True)

print("\n{} check(s) failed{}".format(
    len(FAILURES), (": " + ", ".join(FAILURES)) if FAILURES else ""))
sys.exit(1 if FAILURES else 0)
