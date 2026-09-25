#! python3

# r: compas
# r: pydantic

"""Reset the project: empty the ``TopologyProblem`` layers and recreate them.

The session is reset with them: every item goes, the settings stay. Layer names
live in ``compas_singular.rhino.project``; no command imports another. Edit the
settings with ``CMD_settings``.
"""
import rhinoscriptsyntax as rs

from compas_singular.rhino.project import LAYER_DATA
from compas_singular.rhino.project import ROOT
from compas_singular.rhino.session import RhinoSession


def reset_project():
    """Empty the project layers and recreate them. This command's actual job.

    Under the ``__name__`` guard below. Other commands used to import this file,
    which ran its body -- so before the guard existed, the first command of a Rhino
    session DELETED the user's whole ``TopologyProblem`` layer as a side effect of
    an import. Nothing imports this file any more; the guard stays anyway.
    """
    # Before the purge below, which deletes the session's anchor with every other
    # object under ROOT.
    settings = RhinoSession.current().settings

    #Setup of project folder
    project_folder = ROOT
    if rs.IsLayer(project_folder):
        if rs.IsLayer("Default"):
            rs.CurrentLayer("Default")
        else:
            print("Will not be able to clear problem folder. Set another mayer as active.")

        sublayers = rs.LayerChildren(project_folder)

        for layer in sublayers:
            rs.PurgeLayer(layer)
        objects = rs.ObjectsByLayer(project_folder)
        rs.DeleteObjects(objects)
    else:
        rs.AddLayer(name=project_folder)

    #Setup of subfolders
    for name, color in LAYER_DATA.values():
        rs.AddLayer(name=name, color=color)

    # An empty session with the settings kept.
    session = RhinoSession.current()
    session.clear(*session.ITEMS)
    session.settings = settings
    session.record("Start")
    print(settings.model_dump())

    # Say so. A reset that silently stops happening -- if a future Rhino ran
    # this file under a name other than "__main__" -- would look like a command
    # that did nothing, and the absence of this line is how you would spot it.
    print("CMD_start: project '{}' reset.".format(project_folder))
    print("note: this deletes every input, layout and mesh under '{}' -- nothing "
          "from an earlier session survives.".format(project_folder))
    print("next: CMD_boundary_selection to define a new domain.")
    return project_folder


if __name__ == "__main__":
    reset_project()

#Some general informational links
#Filter numbers: https://developer.rhino3d.com/api/rhinoscript/selection_methods/filterobjects.htm
