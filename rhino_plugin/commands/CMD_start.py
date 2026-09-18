#! python3

# r: compas

"""Reset the project: empty the ``TopologyProblem`` layers and recreate them.

Settings, layer names and side-car paths live in ``compas_singular.rhino.project``;
no command imports another. Edit the settings with ``CMD_settings``.
"""

# ----------------------------------------------------------------------
# Development bootstrap -- delete once compas_singular is installed into Rhino's
# Python. MUST run before any compas_singular import.
#
# Rhino keeps its interpreter alive between runs: it re-initialises ``sys.path``
# but keeps ``sys.modules``. So the source goes on the path every run, and every
# ``compas_singular`` module is dropped so this run imports fresh ones.
# ``framefield`` must be dropped WITH everything else (2026-08-28): kept, it holds
# the previous run's ``Mesh`` classes, ``isinstance`` is silently false, and
# pickling fails with "not the same object as ...Mesh".
import sys
SINGULAR_SRC = r"C:\Users\Casper\libraries\carbcomn\compas_singular\compas_singular\src"
if SINGULAR_SRC not in sys.path:
    sys.path.insert(0, SINGULAR_SRC)
if not getattr(sys, "compas_singular_keep_modules", False):  # set by headless tests
    for _mod in list(sys.modules):
        if _mod == "compas_singular" or _mod.startswith("compas_singular."):
            del sys.modules[_mod]
# ----------------------------------------------------------------------

import rhinoscriptsyntax as rs

from compas_singular.rhino.project import LAYER_DATA
from compas_singular.rhino.project import ROOT
from compas_singular.rhino.project import get_settings
from compas_singular.rhino.project import set_settings


def reset_project():
    """Empty the project layers and recreate them. This command's actual job.

    Under the ``__name__`` guard below. Other commands used to import this file,
    which ran its body -- so before the guard existed, the first command of a Rhino
    session DELETED the user's whole ``TopologyProblem`` layer as a side effect of
    an import. Nothing imports this file any more; the guard stays anyway.
    """
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

    # Write the settings back, so a document that has none gets the defaults
    # stored on it. Settings the document already has are kept.
    settings = set_settings(get_settings())
    print(settings)

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
