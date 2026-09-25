"""Run every MCP test in one go. No Rhino, no model, no network.

    python run_all.py

Use Rhino 8's own CPython -- the same interpreter the server runs on and the
same one the link runs on inside Rhino::

    C:/Users/Casper/.rhinocode/py39-rh8/python.exe run_all.py

Each test runs in its OWN process, deliberately: ``test_rhino_side`` installs
fake ``Rhino`` and ``rhinoscriptsyntax`` modules into ``sys.modules``, and
compas decides whether it is inside Rhino by looking for exactly those. Sharing
an interpreter would let that stub leak into the tests that must run as if Rhino
were absent.
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = ("test_spool.py", "test_wire.py", "test_protocol.py",
         "test_library.py", "test_tools_offline.py", "test_rhino_side.py",
         "test_visual.py", "test_roundtrip.py")


def main():
    failed = []
    for name in TESTS:
        print("\n" + "#" * 62)
        print("# {}".format(name))
        print("#" * 62)
        result = subprocess.run([sys.executable, os.path.join(HERE, name)],
                                cwd=HERE)
        if result.returncode:
            failed.append(name)
    print("\n" + "#" * 62)
    if failed:
        print("# {} of {} FAILED: {}".format(len(failed), len(TESTS),
                                             ", ".join(failed)))
        return 1
    print("# all {} test files passed".format(len(TESTS)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
