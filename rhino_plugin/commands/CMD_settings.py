#! python3
# r: compas
"""Show and edit this document's settings.

    reads   document user text "settings"      over the defaults
    writes  document user text "settings"      on every change

The settings travel with the ``.3dm``. Their defaults and meaning are in
``compas_singular.rhino.project.DEFAULT_SETTINGS``; this command is only the menu.
Several other commands still change single settings in passing
(``CMD_boundary_selection``: spacing, guide alignment, relax; ``CMD_densities``:
the density target; ``CMD_quad_mesh``: field_aware) -- all through the same store,
so what is shown here is always what they will use.
"""

# Development bootstrap -- delete once compas_singular is installed into Rhino's
# Python. MUST run before any compas_singular import: Rhino resets sys.path between
# runs but keeps sys.modules, so put the source on the path and drop a stale copy
# (see CMD_start for why every module, framefield included, has to go).
import sys
SINGULAR_SRC = r"C:\Users\Casper\libraries\carbcomn\compas_singular\compas_singular\src"
if SINGULAR_SRC not in sys.path:
    sys.path.insert(0, SINGULAR_SRC)
if not getattr(sys, "compas_singular_keep_modules", False):  # set by headless tests
    for _mod in list(sys.modules):
        if _mod == "compas_singular" or _mod.startswith("compas_singular."):
            del sys.modules[_mod]

import rhinoscriptsyntax as rs

from compas_singular.rhino.project import DEFAULT_SETTINGS
from compas_singular.rhino.project import get_settings
from compas_singular.rhino.project import set_settings


def show(settings):
    print("settings:")
    for key in sorted(settings):
        marker = "" if settings[key] == DEFAULT_SETTINGS.get(key) else "   (changed)"
        print("    {:<22} {!r}{}".format(key, settings[key], marker))


def choose(message, current, strings):
    """One of ``strings``, lower-cased, or None on Esc."""
    answer = rs.GetString(message=message, defaultString=current, strings=strings)
    return answer.lower() if answer else None


def main():
    settings = get_settings()
    show(settings)

    while True:
        section = rs.GetString(
            message="Edit settings",
            defaultString="Done",
            strings=["Background_Spacing", "Guide_Allignment", "Density_Mode",
                     "Target_Length", "Target_Density", "Field_Aware", "Symmetry",
                     "Relax", "Defaults", "Done"])
        # rs.GetString returns the option as typed or clicked; compare lower-case.
        section = (section or "Done").lower()
        changed = True

        if section == "background_spacing":
            value = rs.GetReal("Background triangulation spacing (NOT the quad size)",
                               settings["triangulation_spacing"], 1e-3)
            if value:
                settings["triangulation_spacing"] = value
        elif section == "guide_allignment":
            answer = choose("Elements run ALONG (tangent) or ACROSS (perpendicular) the guides?",
                            settings["guide_allignment"], ["tangent", "perpendicular"])
            if answer:
                settings["guide_allignment"] = answer
        elif section == "density_mode":
            answer = choose("Strips with no picked density: by target LENGTH or fixed DENSITY?",
                            settings["density_mode"], ["length", "density"])
            if answer:
                settings["density_mode"] = answer
        elif section == "target_length":
            value = rs.GetReal("Target quad edge length", settings["target_length"], 1e-3)
            if value:
                settings["target_length"] = value
        elif section == "target_density":
            value = rs.GetInteger("Elements across every strip", int(settings["target_density"]), 1)
            if value:
                settings["target_density"] = value
        elif section == "field_aware":
            answer = choose("Integrate patch interiors from the field?",
                            "On" if settings["field_aware"] else "Off", ["On", "Off"])
            if answer:
                settings["field_aware"] = answer == "on"
        elif section == "symmetry":
            # JSON null is off -- see resolve_symmetry.
            current = "Off" if settings.get("symmetry") is None else "Auto"
            answer = choose("Detect and use symmetry?", current, ["Auto", "Off"])
            if answer:
                settings["symmetry"] = "auto" if answer == "auto" else None
        elif section == "relax":
            # Stored as "auto", True or False -- the values CMD_boundary_selection writes.
            current = settings.get("relax", "auto")
            default = "Auto" if str(current).lower() == "auto" else ("On" if current else "Off")
            answer = choose("Relax the field solve? Auto = on when guides exist",
                            default, ["Auto", "On", "Off"])
            if answer:
                settings["relax"] = "auto" if answer == "auto" else (answer == "on")
        elif section == "defaults":
            if rs.GetString("Reset every setting to its default?", "No", ["Yes", "No"]) == "Yes":
                settings = dict(DEFAULT_SETTINGS)
            else:
                changed = False
        else:
            break

        if changed:
            set_settings(settings)
            show(settings)


if __name__ == "__main__":
    main()
