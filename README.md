# compas_singular

Welcome to **compas_singular**, a Python library tackles topology finding of patterns, particularly singularities in structured quad meshes.
Based on the [PhD research of Robin Oval](https://blockresearchgroup.github.io/compas_singular/latest/07_publications.html), this library implements data structures for quad mesh modelling, several algorithms for topological exploration and interface with Rhino3D/Grasshopper3D.

**compas_singular** is based on the **COMPAS** framework is an open-source, Python-based framework for computational research and collaboration in architecture, engineering, digital fabrication and construction.

## Getting Started

### Installation

From the root of this repository (the folder containing `pyproject.toml`):

```bash
# runtime only
pip install -e .

# runtime + the 3D viewer used by the examples
pip install -e .[viewer]
```

The `-e` (editable) install puts a `.pth` file pointing at `src/` on your
`sys.path`, so `from compas_singular.datastructures import CoarseQuadMesh`
works from any directory and any script, and your edits to `src/` take effect
immediately without reinstalling.

Prefer conda? A ready-made environment (Python 3.12 + compas 2.15.1 + viewer +
an editable install of this package) is described in `environment.yml`:

```bash
conda env create -f environment.yml
conda activate compas-singular
```

Verify with:

```bash
python -c "import compas_singular; print(compas_singular.__file__)"
```

If that prints `None` instead of a path ending in
`src/compas_singular/__init__.py`, you are not running an installed copy —
Python found the bare repository folder and treated it as an empty namespace
package. Run one of the install commands above.

Dependencies live in `requirements.txt` (core) and `requirements-viewer.txt`
(viewer extra); `pyproject.toml` reads both, so those files are the single
source of truth. The Rhino/Grasshopper modules under
`src/compas_singular/rhino/` additionally need `compas_rhino` and are only
importable inside Rhino.

Upstream docs:
[Installation](https://blockresearchgroup.github.io/compas_singular/latest/01_getting_started.html#installation-1)
[Updating](https://blockresearchgroup.github.io/compas_singular/latest/01_getting_started.html#updates-1)
[Rhino](https://blockresearchgroup.github.io/compas_singular/latest/01_getting_started.html#rhino-1)

## First Steps

See the [overview](https://blockresearchgroup.github.io/compas_singular/latest/02_overview.html) of the algorithms and tools currently proposed and try out the [examples](https://blockresearchgroup.github.io/compas_singular/latest/03_examples.html).

## Questions and feedback

## Issue tracker

## Contributing

## Changelog

## License

**compas_singular** is [released under the MIT license](https://blockresearchgroup.github.io/compas_singular/latest/05_license.html).
