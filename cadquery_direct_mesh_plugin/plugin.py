from OCP.TopLoc import TopLoc_Location
from OCP.BRep import BRep_Tool
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP import GCPnts, BRepAdaptor
import cadquery as cq


# vertices: list[tuple[float, float, float]] | list["cadquery.occ_impl.geom.Vector"],
# triangles_by_solid_by_face: list[list[tuple[int, int, int]]],
# material_tags: list[str],
# h5m_filename: str = "dagmc.h5m",
# implicit_complement_material_tag: str | None = None,


def to_mesh(
    self,
    imprint=True,
    tolerance=0.1,
    angular_tolerance=0.1,
    include_brep_edges=False,
    include_brep_vertices=False,
):
    """
    Converts an assembly to a custom mesh format defined by the CadQuery team.
    :param imprint: Whether or not the assembly should be impr
    :type bar: int
    """

    # To keep track of the vertices and triangles in the mesh
    vertices = []
    solids = []
    solid_face_triangle = []
    imprinted_assembly = None
    imprinted_solids_with_orginal_ids = None
    solid_colors = []
    solid_locs = []
    solid_brep_edge_segments = []
    solid_brep_vertices = []

    # Imprinted assemblies end up being compounds, whereas you have to step through each of the
    # parts in an assembly and extract the solids.
    if imprint:
        # Imprint the assembly and process it as a compound
        (
            imprinted_assembly,
            imprinted_solids_with_orginal_ids,
        ) = cq.occ_impl.assembly.imprint(self)
        solids.append(imprinted_assembly)

        # Keep track of the colors and location of each of the solids
        solid_colors.append((0.5, 0.5, 0.5, 1.0))
        solid_locs.append(cq.Location())
    else:
        # Step through every child in the assembly and save their solids
        for child in self.children:
            # Make sure we end up with a base shape
            obj = child.obj
            if type(child.obj).__name__ == "Workplane":
                solids.append(obj.val())
            else:
                solids.append(obj)

            # Use the color set for the assembly component, or use a default color
            if child.color:
                solid_colors.append(child.color.toTuple())
            else:
                solid_colors.append((0.5, 0.5, 0.5, 1.0))

            # Keep track of the location of each of the solids
            solid_locs.append(child.loc)

    # Step through all of the collected solids and their respective faces to get the vertices
    solid_idx = 0
    for solid in solids:
        # Reset this each time so that we get the correct number of faces per solid
        face_triangles = []

        for face in solid.Faces():
            # Location information of the face to place the vertices and edges correctly
            loc = TopLoc_Location()

            # Perform the tessellation
            BRepMesh_IncrementalMesh(face.wrapped, tolerance, True, angular_tolerance)
            face_mesh = BRep_Tool.Triangulation_s(face.wrapped, loc)

            # If this is not an imprinted assembly, override the location of the triangulation
            if not imprint:
                loc = solid_locs[solid_idx].wrapped

            # Save the transformation so that we can place vertices in the correct locations later
            Trsf = loc.Transformation()

            # Step through the triangles of the face
            cur_triangles = []
            for i in range(1, face_mesh.NbTriangles() + 1):
                triangle_vertex_indices = []

                # Get the current triangle and its index vertices
                cur_tri = face_mesh.Triangle(i)
                idx_1, idx_2, idx_3 = cur_tri.Get()

                # Save the vertices
                tri_verts = []
                tri_verts.append(face_mesh.Node(idx_1))
                tri_verts.append(face_mesh.Node(idx_2))
                tri_verts.append(face_mesh.Node(idx_3))

                for vert in tri_verts:
                    # Apply the assembly location transformation to each vertex
                    v_trsf = vert.Transformed(Trsf)
                    temp_tris = (v_trsf.X(), v_trsf.Y(), v_trsf.Z())

                    # Handle duplicate vertices
                    if temp_tris in vertices:
                        triangle_vertex_indices.append(vertices.index(temp_tris))
                    else:
                        # Append the vertices for this face and triangle
                        vertices.append(temp_tris)

                        # The vertex we just added is the one we need to track
                        triangle_vertex_indices.append(len(vertices) - 1)

                cur_triangles.append(triangle_vertex_indices)

            # Save this triangle for the current face
            face_triangles.append(cur_triangles)

        solid_face_triangle.append(face_triangles)

        # If the caller wants to track edges, include them
        if include_brep_edges:
            # If this is not an imprinted assembly, override the location of the edges
            loc = TopLoc_Location()
            if not imprint:
                loc = solid_locs[solid_idx].wrapped

            # Save the transformation so that we can place vertices in the correct locations later
            Trsf = loc.Transformation()

            # Add CadQuery-reported edges
            current_segments = []
            for edge in solid.edges():
                # We need to handle different kinds of edges differently
                gt = edge.geomType()

                # Line edges are just point to point
                if gt == "LINE":
                    start = edge.startPoint().toPnt()
                    end = edge.endPoint().toPnt()

                    # Apply the assembly location transformation to each vertex
                    start_trsf = start.Transformed(Trsf)
                    located_start = (start_trsf.X(), start_trsf.Y(), start_trsf.Z())
                    end_trsf = end.Transformed(Trsf)
                    located_end = (end_trsf.X(), end_trsf.Y(), end_trsf.Z())

                    # Save the start and end points for the edge
                    current_segments.append([located_start, located_end])
                # If dealing with some sort of arc, discretize it into individual lines
                elif (
                    gt == "CIRCLE"
                    or gt == "ARC"
                    or gt == "SPLINE"
                    or gt == "BSPLINE"
                    or gt == "ELLIPSE"
                ):
                    # Discretize the curve
                    disc = GCPnts.GCPnts_TangentialDeflection(
                        BRepAdaptor.BRepAdaptor_Curve(edge.wrapped),
                        tolerance,
                        angular_tolerance,
                    )

                    # Add each of the discretized sections to the edge list
                    if disc.NbPoints() > 1:
                        for i in range(2, disc.NbPoints() + 1):
                            p_0 = disc.Value(i - 1)
                            p_1 = disc.Value(i)

                            # Apply the assembly location transformation to each vertex
                            p_0_trsf = p_0.Transformed(Trsf)
                            located_p_0 = (p_0_trsf.X(), p_0_trsf.Y(), p_0_trsf.Z())
                            p_1_trsf = p_1.Transformed(Trsf)
                            located_p_1 = (p_1_trsf.X(), p_1_trsf.Y(), p_1_trsf.Z())

                            # Save the start and end points for the edge
                            current_segments.append([located_p_0, located_p_1])

            solid_brep_edge_segments.append(current_segments)

        # Add CadQuery-reported vertices, if requested
        if include_brep_vertices:
            # If this is not an imprinted assembly, override the location of the edges
            loc = TopLoc_Location()
            if not imprint:
                loc = solid_locs[solid_idx].wrapped

            # Save the transformation so that we can place vertices in the correct locations later
            Trsf = loc.Transformation()

            current_vertices = []
            for vertex in solid.vertices():
                p = BRep_Tool.Pnt_s(vertex.wrapped)

                # Apply the assembly location transformation to each vertex
                p_trsf = p.Transformed(Trsf)
                located_p = (p_trsf.X(), p_trsf.Y(), p_trsf.Z())
                current_vertices.append(located_p)

            solid_brep_vertices.append(current_vertices)

        solid_idx += 1

    return {
        "vertices": vertices,
        "solid_face_triangle_vertex_map": solid_face_triangle,
        "solid_colors": solid_colors,
        "solid_brep_edge_segments": solid_brep_edge_segments,
        "solid_brep_vertices": solid_brep_vertices,
    }


# Patch the function(s) into the Workplane class
cq.Assembly.toMesh = to_mesh
cq.Assembly.to_mesh = to_mesh
