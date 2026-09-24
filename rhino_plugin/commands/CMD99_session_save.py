#! python3

# r: compas
# r: pydantic

"""**Save the whole project to one JSON file**: settings, domain, layout, field
and dense mesh.

The file is a ``SingularSession``, so a plain Python script -- no Rhino -- opens
it with ``SingularSession.load(path)``, edits it, and ``dump``s it again;
``CMD_session_open`` brings it back. The document keeps its own copy whatever
happens here: this is an export, not the save of the ``.3dm``.
"""
import os

import rhinoscriptsyntax as rs

from compas_singular.rhino.session import RhinoSession


def main():
    session = RhinoSession.current()
    doc = session.doc
    folder = os.path.dirname(doc.Path) if doc.Path else None
    name = (os.path.splitext(doc.Name)[0] if doc.Name else "project") + ".json"
    path = rs.SaveFileName("Save the project", "JSON (*.json)|*.json||", folder, name)
    if not path:
        print("Cancelled -- nothing saved.")
        return
    session.dump(path)
    items = [item for item in session.ITEMS if getattr(session, item) is not None]
    print("project saved to {}: settings, {}".format(path, ", ".join(items) or "no items yet"))


if __name__ == "__main__":
    main()
