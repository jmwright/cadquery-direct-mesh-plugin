from OCP.TopLoc import TopLoc_Location
from OCP.BRep import BRep_Tool
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP import GCPnts, BRepAdaptor
from OCP.TopAbs import TopAbs_REVERSED, TopAbs_IN
from OCP.gp import gp_Pnt, gp_Vec
import cadquery as cq


def to_mesh(
    self,
    imprint=True,
    tolerance=0.1,
    angular_tolerance=0.1,
    scale_factor=1.0,
    include_brep_edges=False,
    include_brep_vertices=False,
    parallel: bool = True,
):
    """
    Converts an assembly to a custom mesh format defined by the CadQuery team.

    :param imprint: Whether or not the assembly should be imprinted
    :param tolerance: Tessellation tolerance for mesh generation
    :param angular_tolerance: Angular tolerance for tessellation
    :param include_brep_edges: Whether to include BRep edge segments
    :param include_brep_vertices: Whether to include BRep vertices
    :param parallel: If True, OCCT will use parallel processing to mesh the assembly. Default is True.
    """

    # To keep track of the vertices and triangles in the mesh
    vertices = []
    vertex_map = {}
    solids = []
    solid_face_triangle = {}
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

        # Extract the solids from the imprinted assembly because we should not mesh the compound
        for solid in imprinted_assembly.Solids():
            solids.append(solid)

        # Keep track of the colors and location of each of the solids
        solid_colors.append((0.5, 0.5, 0.5, 1.0))
        solid_locs.append(cq.Location())
    else:
        # Step through every child in the assembly and save their solids
        for child in self.children:
            # Make sure we end up with a base shape
            obj = child.obj
            if isinstance(obj, cq.Workplane):
                val = obj.val()
                if isinstance(val, cq.Solid) or isinstance(val, cq.Compound):
                    solids.append(val)
            elif isinstance(obj, cq.Solid):
                solids.append(obj)
            elif isinstance(obj, cq.Compound):
                solids.append(obj)
            else:
                continue

            # Use the color set for the assembly component, or use a default color
            if child.color:
                solid_colors.append(child.color.toTuple())
            else:
                solid_colors.append((0.5, 0.5, 0.5, 1.0))

            # Keep track of the location of each of the solids
            solid_locs.append(child.loc)

    # Solid and face IDs need to be unique unless they are a shared face
    solid_idx = 1  # We start at 1 to mimic gmsh
    face_idx = 1  # We start at id of 1 to mimic gmsh

    # Step through all of the collected solids and their respective faces to get the vertices
    for solid in solids:
        # Reset this each time so that we get the correct number of faces per solid
        face_triangles = {}

        # Perform the tessellation
        BRepMesh_IncrementalMesh(
            solid.wrapped, tolerance, False, angular_tolerance, parallel
        )

        # Walk through all the faces
        for face in solid.Faces():
            # Figure out if the face has a reversed orientation so we can handle the triangles accordingly
            is_reversed = False
            if face.wrapped.Orientation() == TopAbs_REVERSED:
                is_reversed = True

            # Location information of the face to place the vertices and edges correctly
            loc = TopLoc_Location()

            # Retrieve the triangulation
            face_mesh = BRep_Tool.Triangulation_s(face.wrapped, loc)

            # If this is not an imprinted assembly, override the location of the triangulation
            if not imprint:
                loc = solid_locs[solid_idx - 1].wrapped

            # Save the transformation so that we can place vertices in the correct locations later
            Trsf = loc.Transformation()

            # Pre-process all vertices from the face mesh for better performance
            face_vertices = {}  # Map from face mesh node index to global vertex index
            for node_idx in range(1, face_mesh.NbNodes() + 1):
                node = face_mesh.Node(node_idx)
                v_trsf = node.Transformed(Trsf)
                vertex_coords = (
                    v_trsf.X() * scale_factor,
                    v_trsf.Y() * scale_factor,
                    v_trsf.Z() * scale_factor,
                )

                # Use dictionary for O(1) lookup instead of O(n) list operations
                if vertex_coords in vertex_map:
                    face_vertices[node_idx] = vertex_map[vertex_coords]
                else:
                    global_vertex_idx = len(vertices)
                    vertices.append(vertex_coords)
                    vertex_map[vertex_coords] = global_vertex_idx
                    face_vertices[node_idx] = global_vertex_idx

            # Step through the triangles of the face
            cur_triangles = []
            for i in range(1, face_mesh.NbTriangles() + 1):
                # Get the current triangle and its index vertices
                cur_tri = face_mesh.Triangle(i)
                idx_1, idx_2, idx_3 = cur_tri.Get()

                # Look up pre-processed vertex indices - O(1) operation
                if is_reversed:
                    triangle_vertex_indices = [
                        face_vertices[idx_1],
                        face_vertices[idx_3],
                        face_vertices[idx_2],
                    ]
                else:
                    triangle_vertex_indices = [
                        face_vertices[idx_1],
                        face_vertices[idx_2],
                        face_vertices[idx_3],
                    ]

                cur_triangles.append(triangle_vertex_indices)

            # Save this triangle for the current face
            face_triangles[face_idx] = cur_triangles

            # Move to the next face
            face_idx += 1

        solid_face_triangle[solid_idx] = face_triangles

        # If the caller wants to track edges, include them
        if include_brep_edges:
            # If this is not an imprinted assembly, override the location of the edges
            loc = TopLoc_Location()
            if not imprint:
                loc = solid_locs[solid_idx - 1].wrapped

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
                elif gt in ("CIRCLE", "ARC", "SPLINE", "BSPLINE", "ELLIPSE"):
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
                loc = solid_locs[solid_idx - 1].wrapped

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

        # Move to the next solid
        solid_idx += 1

    return {
        "vertices": vertices,
        "solid_face_triangle_vertex_map": solid_face_triangle,
        "solid_colors": solid_colors,
        "solid_brep_edge_segments": solid_brep_edge_segments,
        "solid_brep_vertices": solid_brep_vertices,
        "imprinted_assembly": imprinted_assembly,
        "imprinted_solids_with_orginal_ids": imprinted_solids_with_orginal_ids,
    }


# Patch the function(s) into the Workplane class
cq.Assembly.toMesh = to_mesh
