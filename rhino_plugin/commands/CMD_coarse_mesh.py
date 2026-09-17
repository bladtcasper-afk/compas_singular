#! python3

# r: compas

#Temporary import of compas_singular development library
from CMD_start import import_compas_singular
import_compas_singular()


import rhinoscriptsyntax as rs
import scriptcontext as sc
import Rhino

import compas_rhino as cr
from compas_rhino.conversions import point_to_rhino

from compas_singular.algorithms import SkeletonDecomposition
from compas_singular.framefield.decomposition import FieldDecomposition
from compas_singular.rhino.helpers.helpers import clear_layer, read_boundary_loops, bake_edge_curves, bake_polylines, read_boundaries
from compas_singular.rhino.coarse_curves import coarse_edges_to_curves, snap_corners_to_walls
from CMD_start import get_settings, resolve_relax, resolve_symmetry
from CMD_start import cache_path, COARSE_CACHE, FIELD_CACHE

#The walls the layout's boundary edges densify ALONG, as a multiple of the
#background spacing. Finer than the background, because these points are the
#boundary from here on: the final mesh follows the real curve only as closely as
#these follow it.
#
#(A CURVED input is sampled at the background spacing itself -- ``read_boundaries``
#takes it directly. An arc or a NURBS curve is usable at all only because of
#that; ``curve_to_compas_polyline`` raised ConversionError and the command died
#here.)
WALL_SAMPLING_FACTOR = 0.25


#Coarse mesh generation

#New approach: frame-field decomposition -> all-quad mesh
#See examples/New approach/README.md -- "the two lines that matter":
#  d = FieldDecomposition.from_boundary(outline, guides=cables)
#  mesh = d.quad_mesh(target_length=0.5)

def coarse_from_field(settings, outer, inners, guides):
    """Frame-field decomposition. The mirror of ``coarse_from_skeleton`` above.

    Solved fresh every time this step runs, and the only step that solves: steps
    4 and 6 read the field back from the ``field.json`` side-car written below.
    The user is choosing to (re)build the coarse mesh, so the extra computation
    when nothing changed is an accepted cost for a call site that reads the same,
    directly, as ``SkeletonDecomposition.from_boundary`` below.
    """
    decomposition = FieldDecomposition.from_boundary(
        outer, inner_boundaries=inners, guides=guides,
        mode=settings["guide_allignment"],
        target_length=settings["triangulation_spacing"],
        relax=resolve_relax(settings, guides),
        symmetry=resolve_symmetry(settings))
    print(decomposition.field.report())
    # Print the group: a smaller group than the drawing suggests means the
    # layout is ALLOWED to be less symmetric, and nothing downstream complains.
    group = decomposition.symmetry
    if group is None or group.trivial:
        print("symmetry: none detected -- the layout is under no obligation to have any")
    else:
        print("symmetry: order {} about ({:.3f}, {:.3f}) -- {}".format(
            len(group), group.centre[0], group.centre[1], ", ".join(group.names())))

    skeleton = decomposition.decomposition_polylines()
    coarse_mesh = decomposition.coarse_mesh()

    #Which route built this layout, so step 6 does not have to guess. It rides
    #in ``attributes`` and so travels with the side-car JSON.
    coarse_mesh.attributes["route"] = "field"

    return coarse_mesh, skeleton, decomposition.get_field()

#From Skeleton
def coarse_from_skeleton(settings, outer, inners, line_features, point_features):
    """Medial-axis decomposition. The mirror of ``coarse_from_field`` above.

    ``from_boundary`` owns the whole front half: it resamples the walls so no
    segment is longer than ``spacing`` -- the medial axis is read off a
    Delaunay triangulation of the boundary POINTS, so a four-point square gives
    a two-triangle mesh and a caricature layout -- triangulates them, and
    remembers the domain it was given. ``coarse_mesh`` then takes its poles from
    that domain rather than being handed them a second time.

    ``spacing`` is the BACKGROUND spacing, not the quad size; the quad size
    is ``settings["spacing"]`` and is applied in step 5.

    **Guides are passed as polyline features**, so here they CUT the domain and
    the layout is built around them: each guide becomes a chain of mesh edges,
    a guide ending on a wall a three-valent boundary vertex, a free end a
    singularity. The field route instead steers the layout onto a guide without
    cutting along it. See ``HOW_IT_WORKS.md`` section 5 for what still does not
    work.
    """
    decomposition = SkeletonDecomposition.from_boundary(
        outer,
        inner_boundaries=inners,
        polyline_features=line_features,
        point_features=point_features)

    coarse_mesh = decomposition.coarse_mesh()
    # ``decomposition.polylines`` rather than a second ``decomposition_polylines()``
    # call: that method reassigns the list, and on the field route it would
    # discard any user-drawn curves published into it. ``coarse_mesh`` has
    # already filled it.
    skeleton = decomposition.polylines

    coarse_mesh.attributes["route"] = "skeleton"

    #No field on this route. Step 6 can still solve one and densify with it --
    #a layout and a field only have to agree on the DOMAIN, not on an origin --
    #but nothing here has one to hand over.
    return coarse_mesh, skeleton, None


def main():
    settings = get_settings()
    #Finding an checking of boundaries
    outer, inners, guides, point_features = read_boundaries(
        spacing=settings["triangulation_spacing"])
    wall_sampling = settings["triangulation_spacing"] * WALL_SAMPLING_FACTOR

    coarse_mesh, skeleton, field = None, None, None
    while True:
        mode = rs.GetString(message="Choose coarse mesh generation", defaultString="Skeleton", strings=["FrameField", "Skeleton", "Exit"])
        mode = (mode or "Exit").lower()

        if mode == "skeleton":
            coarse_mesh, skeleton, field = coarse_from_skeleton(
                settings, outer, inners, guides, point_features)
            break
        elif mode == "framefield":
            coarse_mesh, skeleton, field = coarse_from_field(settings, outer, inners, guides)
            break
        elif mode == "exit":
            break
        else:
            print("Unrecognised mode: {}".format(mode))

    
    #Bake results
    if coarse_mesh and skeleton:
        rs.EnableRedraw(False)
        clear_layer("Skeleton", clean_sublayers = True) #Clean Layer before adding objects
        
        #BEFORE the bake, or the baked mesh and the edge curves disagree about where
        #the corners are. A boundary corner of the layout IS a point of the domain
        #boundary, but nothing upstream puts it there -- the background
        #triangulation places it, and on a curved wall that means on a CHORD of the
        #wall. An arc is pinned to its two corners, so those corners would be the
        #only points of the dense boundary not on the wall. Measured on a plate with
        #three r=0.7 holes near the wall: worst corner 0.106 off its own hole, 15.1%
        #of the radius. Snapping costs nothing -- same patches, same face count,
        #quality unchanged to a few tenths of a degree.
        #
        #NOT ``decomposition.edges_to_curves()``, which does the same two steps for a
        #caller that has no document. Here the walls are RE-READ from the Rhino
        #curves at WALL_SAMPLING, four times finer than the background the
        #decomposition was built on -- and on a curved wall that is the difference
        #between the true curve and a resampling of the background's own chords.
        #The decomposition cannot beat that; it only has what it was handed.
        outer_loop, inner_loops = read_boundary_loops(wall_sampling)
        snapped, worst = snap_corners_to_walls(
            coarse_mesh, loops=[outer_loop] + inner_loops)
        if snapped:
            print("snapped {} boundary corner(s) onto their wall, worst {:.4f}".format(
                snapped, worst))

        rs.AddLayer(name="Skeleton", parent="TopologyProblem")

        layer = rs.AddLayer(name="Poles", parent="Skeleton")
        vkeys = coarse_mesh.poles()
        pts = [coarse_mesh.vertex_coordinates(v) for v in vkeys]
        points = [point_to_rhino(pt) for pt in pts]
        guids = rs.AddPoints(points)
        for guid in guids:
            rs.ObjectLayer(guid, layer=layer)

        #``rs.AddPolyline`` RAISES -- 'Unable to add polyline to document' -- on a
        #rib Rhino will not take, and it judges that at the DOCUMENT tolerance, not
        #at 1e-9: a separatrix that doubles back on itself by half a tolerance is
        #geometrically fine and still unbakeable. That used to kill the command
        #here, which cost the MESH as well, because the mesh is baked below this.
        #``bake_polylines`` cleans each rib and guards the add.
        layer = rs.AddLayer(name="Polylines", parent="Skeleton")
        _guids, skipped = bake_polylines(skeleton, layer)
        print("separatrices: {} baked".format(len(_guids)))
        if skipped:
            print("  {} rib(s) skipped -- Rhino refused them (coincident points at "
                  "the document tolerance, or fewer than two)".format(skipped))

        layer = rs.AddLayer(name="Mesh", parent="Skeleton")
        # face_vertices gives vertex KEYS, the vertex list is positional, and the
        # two only agree while the keys are 0..n-1. Repair and weld leave gaps, and
        # the baked faces would then reference the wrong corners -- silently.
        index = {v: i for i, v in enumerate(coarse_mesh.vertices())}
        vertices = [coarse_mesh.vertex_attributes(v, "xyz") for v in coarse_mesh.vertices()]
        faces = [[index[v] for v in coarse_mesh.face_vertices(f)] for f in coarse_mesh.faces()]
        guid = rs.AddMesh(vertices, faces)
        rs.ObjectLayer(guid, layer=layer)

        #The layout WITH its curvature. A Rhino mesh has straight edges, so on a
        #curved domain the mesh above is drawn cutting the corner off its own
        #boundary -- which reads as the workflow having lost the curve when only
        #the display has. These polylines are what CMD_quad_mesh densifies along,
        #rebuilt there from the same document data rather than read back from here.
        curves, tally = coarse_edges_to_curves(
            coarse_mesh, loops=[outer_loop] + inner_loops, polylines=skeleton)
        # AddLayer returns the FULL '::' path, which is what every rs call below
        # needs -- a nested layer's short name is not resolvable on its own.
        edge_layer = rs.AddLayer(name="EdgeCurves", parent="Skeleton", color=(0, 120, 200))
        _guids, skipped = bake_edge_curves(curves.values(), edge_layer)
        print("coarse edges: {}".format(tally))
        if skipped:
            print("  {} edge curve(s) Rhino refused -- those edges densify as "
                  "chords".format(skipped))
        if tally["chord"]:
            print("  {} edge(s) densify as a straight chord -- interior edges with "
                  "no traced branch, which is expected, or a boundary corner that "
                  "is not on a wall, which is not".format(tally["chord"]))
        if tally.get("wall_missed"):
            print("  {} boundary edge(s) did NOT get their wall arc -- on a hole "
                  "this is what makes one circle come out round and the next a "
                  "polygon".format(tally["wall_missed"]))

        #The side-cars: everything the bake above could not carry.
        #
        #``rs.AddMesh`` stores vertices and faces and nothing else, so the
        #layout's own knowledge -- strips, poles, the route that built it --
        #dies at the bake and every later step re-derives it. ``save_to_json``
        #carries the whole ``attributes`` dict, and the field cannot be baked at
        #all. Written LAST, so a failure in the bake above never leaves a
        #side-car describing a layout the document does not have.
        #
        #``curves`` -- already computed above for the EdgeCurves bake -- is
        #stored on the layout too, so a later step (CMD_quad_mesh, or a reload
        #of this side-car) can just ask ``coarse.edges_to_curves()`` instead of
        #re-deriving it from the document every time.
        coarse_mesh.set_edges_to_curves(curves)
        coarse_mesh.set_global_face_pattern("ortho")
        coarse_mesh.save_to_json(cache_path(COARSE_CACHE))
        print("layout cached: {}".format(cache_path(COARSE_CACHE, create=False)))
        if field is not None:
            field.save_to_json(cache_path(FIELD_CACHE))
            print("field cached:  {}".format(cache_path(FIELD_CACHE, create=False)))
        print("note: this replaces any previous layout -- densities (CMD_densities) and "
              "patterns (CMD_dense_pattern) set on it do not carry over.")
        print("next: CMD_edit_coarse_mesh to hand-edit the layout, then CMD_densities and "
              "CMD_quad_mesh.")

        rs.EnableRedraw(True)


if __name__ == "__main__":
    main()
