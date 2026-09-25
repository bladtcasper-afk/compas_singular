# compas_singular

Welcome to **compas_singular**, a Python library tackles topology finding of patterns, particularly singularities in structured quad meshes.
Based on the [PhD research of Robin Oval](https://bladtcasper-afk.github.io/compas_singular/latest/publications.html), this library implements data structures for quad mesh modelling, several algorithms for topological exploration and interface with Rhino3D/Grasshopper3D.

**compas_singular** is based on the **COMPAS** framework is an open-source, Python-based framework for computational research and collaboration in architecture, engineering, digital fabrication and construction.

## Getting Started

### Installation

Python 3.9 or later. The package is not on PyPI yet; install it from GitHub:

```bash
pip install "compas_singular @ git+https://github.com/bladtcasper-afk/compas_singular.git@dev"
```

Optional extras: `[viewer]` for the 3D viewer used by the examples, `[fd]` for
constrained smoothing (compas_fd), and `[rhino]` for the Rhino 8 scene objects
(compas_rui).

To work on the code, install it editable from the root of this repository (the
folder containing `pyproject.toml`):

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

Dependencies live in `requirements.txt` (core) and `requirements-viewer.txt`,
`requirements-fd.txt` and `requirements-rhino.txt` (the extras);
`pyproject.toml` reads them, so those files are the single source of truth.
The modules under `src/compas_singular/rhino/` additionally need
`compas_rhino` and are only importable inside Rhino.

### Quick start

A rectangle with a round hole, meshed by both routes:

```python
import math

from compas_singular.algorithms import SkeletonDecomposition
from compas_singular.framefield import FieldDecomposition

outer = [[-6, -4, 0], [6, -4, 0], [6, 4, 0], [-6, 4, 0]]
hole = [[1.6 * math.cos(2 * math.pi * i / 64), 1.6 * math.sin(2 * math.pi * i / 64), 0] for i in range(64)]

# Skeleton route: boundary -> coarse quad layout -> densities -> dense quad mesh.
coarse = SkeletonDecomposition.from_boundary(outer, inner_boundaries=[hole], target_length=0.5).coarse_mesh()
coarse.collect_strips()
coarse.set_strips_density_target(0.5)
mesh = coarse.densify()
print(mesh.number_of_faces(), 'quads')

# Frame-field route: the same domain, the layout traced from a cross field.
decomposition = FieldDecomposition.from_boundary(outer, inner_boundaries=[hole], target_length=0.5)
mesh = decomposition.quad_mesh(target_length=0.5)
print(mesh.number_of_faces(), 'quads on route', decomposition.route())
```

[`examples/GitHub/`](https://github.com/bladtcasper-afk/compas_singular/tree/dev/examples/GitHub) has thirteen worked examples, from this workflow to layout
editing, patterns, guides, symmetry and dual blocks.

### Rhino 8

The Rhino commands are in `rhino_plugin/commands/`, numbered in workflow order
(`CMD00_start` to `CMD09_edit_quad_mesh`). They run on Rhino 8's CPython 3.9
and need `compas_singular` installed in Rhino's Python. Until the package is
published, install it with Rhino's interpreter (Rhino closed):

```bash
%USERPROFILE%\.rhinocode\py39-rh8\python.exe -m pip install "compas_singular[rhino,fd] @ git+https://github.com/bladtcasper-afk/compas_singular.git@dev"
```

[`markdowns/RHINO_PLUGIN.md`](https://github.com/bladtcasper-afk/compas_singular/blob/dev/markdowns/RHINO_PLUGIN.md) describes the workflow and the layers it creates.

The [Rhino page](https://bladtcasper-afk.github.io/compas_singular/latest/tutorial/rhino.html) of the documentation lists the commands.

## First Steps

The [documentation](https://bladtcasper-afk.github.io/compas_singular/latest/) has an [overview](https://bladtcasper-afk.github.io/compas_singular/latest/tutorial/overview.html) of the workflow, the [examples](https://bladtcasper-afk.github.io/compas_singular/latest/examples.html) with pictures, and the [API reference](https://bladtcasper-afk.github.io/compas_singular/latest/api.html).

## Questions and feedback

Open an issue on the [issue tracker](https://github.com/bladtcasper-afk/compas_singular/issues).

## Issue tracker

Bugs and feature requests: <https://github.com/bladtcasper-afk/compas_singular/issues>.
Please include a script that reproduces the problem.

## Contributing

See [CONTRIBUTING.md](https://github.com/bladtcasper-afk/compas_singular/blob/dev/CONTRIBUTING.md).

## Changelog

See [CHANGELOG.md](https://github.com/bladtcasper-afk/compas_singular/blob/dev/CHANGELOG.md). What this fork changed relative to the
upstream [BRG-research/compas_singular](https://github.com/BRG-research/compas_singular)
is in [markdowns/CHANGES_VS_UPSTREAM.md](https://github.com/bladtcasper-afk/compas_singular/blob/dev/markdowns/CHANGES_VS_UPSTREAM.md).

## License

**compas_singular** is [released under the MIT license](LICENSE).
