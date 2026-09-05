"""Draw what the automatic selection picks: the mesh, the guide, and the chain.

The table in :mod:`eval_guide_chain` says the chains are clean long polyedge lines. This
is where you look at them and agree or not. Guides are drawn in blue, the chain as a thick
red line over the mesh, and the meshes are laid out in a row.

    python examples/guide_chain_tests/view_guide_chain.py
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import sys

from compas.colors import Color
from compas.geometry import Polyline
from compas.geometry import Translation

from compas_singular.editing import GuideCurve
from compas_singular.editing import chain_quality
from compas_singular.editing import collect_polyedges
from compas_singular.editing import guide_chain

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from build_meshes import cases  # noqa: E402


def main():
    from compas_viewer import Viewer

    viewer = Viewer()
    offset = 0.0

    for name, mesh, guides in cases():
        polyedges = collect_polyedges(mesh)
        width = mesh.aabb().xsize
        shift = Translation.from_vector([offset, 0, 0])

        drawn = mesh.copy()
        drawn.transform(shift)
        group = viewer.scene.add_group(name=name)
        group.add(drawn, facecolor=Color.from_hex('#E8E8E8'), show_lines=True)

        for guide_name, points in guides:
            guide = GuideCurve(points)
            chain, info = guide_chain(mesh, guide, polyedges=polyedges)

            curve = Polyline([list(point) for point in guide.points])
            curve.transform(shift)
            group.add(curve, linecolor=Color.from_hex('#4C72B0'), linewidth=2,
                      name='{} guide'.format(guide_name))

            if not chain:
                print('{:<12} {:<9} {}'.format(name, guide_name, info['reason']))
                continue

            line = Polyline([mesh.vertex_coordinates(vertex) for vertex in chain])
            line.transform(shift)
            group.add(line, linecolor=Color.from_hex('#C44E52'), linewidth=6,
                      name='{} chain'.format(guide_name))

            quality = chain_quality(mesh, chain, guide)
            print('{:<12} {:<9} {:>3} vertices, {:>4.0f}% of the guide, alignment {:.2f}'.format(
                name, guide_name, quality['n'], 100 * quality['coverage'],
                quality['alignment']))

        offset += width + 3.0

    total = offset - 3.0
    viewer.renderer.camera.target = [0.5 * total, 5.0, 0.0]
    viewer.renderer.camera.position = [0.5 * total, 4.0, 1.1 * total]
    viewer.show()


if __name__ == '__main__':
    main()
