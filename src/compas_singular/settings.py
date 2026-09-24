"""The project settings: one typed model, stored in the session.

Flat on purpose: eight values do not need groups. The field names are the keys
the Rhino commands have always read, so ``settings.model_dump()`` is the dict
they expect. The one rename, ``guide_allignment`` -> ``guide_alignment``, still
reads the old spelling, so a document or file written before it loads.

A pydantic model, like ``compas_session.Settings`` and RhinoVAULT's
``RVSettings``: typed, with defaults, and editable field by field in a form.
"""
from typing import Optional
from typing import Union

from pydantic import AliasChoices
from pydantic import BaseModel
from pydantic import Field


__all__ = ['Settings']


class Settings(BaseModel):

    #: Background triangulation spacing. NOT the quad size -- the field is
    #: solved on this. ``None`` takes thesis eq. 4.1, 0.04 times the domain's
    #: bounding-box diagonal, in both routes; a number overrides it. Finer means
    #: a better field and a slower solve; 0.3 on a 20 x 14 m plate is 2 438
    #: vertices and 2.1 s. Rhino-side lengths: ``project.resolve_spacing``.
    triangulation_spacing: Optional[float] = None

    #: 'tangent' makes elements run ALONG the guides, 'perpendicular' across.
    guide_alignment: str = Field(
        'tangent', validation_alias=AliasChoices('guide_alignment', 'guide_allignment'))

    #: Target quad edge length for strips with no explicit density.
    target_length: float = 0.5

    #: Elements across every strip with no explicit density, when
    #: ``density_mode`` is ``"density"``.
    target_density: int = 5

    #: Which global rule step 5 last applied: ``"length"`` or ``"density"``.
    #: Only consulted when a layout has no densities saved on it.
    density_mode: str = 'length'

    #: Whether patch INTERIORS are integrated from the field. Real on the field
    #: route; measured to COST quality on a skeleton layout, whose patch edges
    #: do not follow the field -- aspect 1.99 -> 3.20 on a square with a cable.
    field_aware: bool = True

    #: FRAME-FIELD route only: the symmetry group the field is solved under.
    #: ``"auto"`` detects it; ``None`` (JSON null) disables it. The skeleton
    #: route ignores it.
    field_symmetry: Optional[str] = 'auto'

    #: ``"auto"`` relaxes when there are guides; or a bool.
    relax: Union[bool, str] = 'auto'
