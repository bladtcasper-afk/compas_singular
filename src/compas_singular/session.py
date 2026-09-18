import rhinoscriptsyntax as rs  # type: ignore

from compas.scene import Scene
from compas_rhino.scene import RhinoSceneObject
from compas_rv.settings import RVSettings
from compas_session.session import Session

class SingularSession(Session):
    