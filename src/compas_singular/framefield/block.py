"""**The singularity as a BLOCK** -- the local operator, plus the global dual.

The local operator MOVED. It now lives in ``compas_singular.blocks``, because
it stopped being a frame-field experiment the moment a picked point could drive
it: see :func:`compas_singular.blocks.block_points`. Everything it exported is
re-exported here unchanged, so ``from compas_singular.framefield import block
as B`` and every ``B.block_at`` / ``B.blockable`` / ``B.head_strips`` call
still resolves.

``quad_dual`` stays here. It is the GLOBAL alternative -- every primal vertex
becomes a face at once, nothing propagates, but every block in the plate moves
by half a block and on a curved wall the dual boundary chords the outline and
loses area. The two routes answer the same question at opposite scales and are
worth reading side by side; ``23_singularity_block.py`` does exactly that.

Their restrictions are complementary, which is the useful thing to remember:

    local  refuses when a POLE sits in the joint's own 1-ring, or when a repair
           strip has no wall to run to. A pole ELSEWHERE in the mesh is fine.

    global refuses nothing but must not be pointed at a mesh with poles at all
           -- a pole is already a face, so dualising trades it the wrong way
           and turns each one into a valence-3 joint.
"""
from compas.datastructures import Mesh

from compas.tolerance import TOL

from compas_singular.blocks import *            # noqa: F401 F403
from compas_singular import blocks as _blocks


__all__ = list(_blocks.__all__) + ['quad_dual']


# ----------------------------------------------------------------------------
# the global alternative, for comparison
# ----------------------------------------------------------------------------

def quad_dual(mesh):
    """The GLOBAL dual of an all-quad mesh, closed on the boundary with quads.

    Every primal vertex becomes a face, so every singularity becomes a block at
    once and nothing has to propagate -- but every block in the plate moves by
    half a block, and on a CURVED wall the dual boundary chords the outline and
    loses area. Lifted from ``FloorExample vault.py``; see that file for why the
    corner case (valence 2) is the only one that keeps its own vertex.

    Do not use on a mesh with poles: a pole is already a face, so dualising
    trades it the wrong way and turns it into a valence-3 joint. Use
    :func:`compas_singular.blocks.block_points` there instead -- it only needs
    the joint's own 1-ring to be quads, so poles elsewhere do not stop it.
    """
    index, points = {}, []

    def add(xyz):
        gk = TOL.geometric_key(xyz)
        if gk not in index:
            index[gk] = len(points)
            points.append([float(c) for c in xyz])
        return index[gk]

    faces = []
    for v in mesh.vertices():
        fkeys = [f for f in mesh.vertex_faces(v, ordered=True) if f is not None]
        if not mesh.is_vertex_on_boundary(v):
            faces.append([add(mesh.face_centroid(f)) for f in fkeys])
            continue
        nbrs = mesh.vertex_neighbors(v, ordered=True)
        ring = [add(mesh.edge_midpoint((v, nbrs[0])))]
        ring += [add(mesh.face_centroid(f)) for f in fkeys]
        ring.append(add(mesh.edge_midpoint((v, nbrs[-1]))))
        if len(nbrs) == 2:
            ring.append(add(mesh.vertex_coordinates(v)))
        faces.append(ring)

    dual = Mesh.from_vertices_and_faces(points, faces)
    dual.unify_cycles()
    return dual
