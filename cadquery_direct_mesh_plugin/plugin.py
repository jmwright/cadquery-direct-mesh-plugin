from OCP.TopLoc import TopLoc_Location
from OCP.BRep import BRep_Tool
from OCP.BRepMesh import BRepMesh_IncrementalMesh
import cadquery as cq


# vertices: list[tuple[float, float, float]] | list["cadquery.occ_impl.geom.Vector"],
# triangles_by_solid_by_face: list[list[tuple[int, int, int]]],
# material_tags: list[str],
# h5m_filename: str = "dagmc.h5m",
# implicit_complement_material_tag: str | None = None,


def to_mesh(
    self, imprint=True, reverse_winding=False, tolerance=0.1, angular_tolerance=0.1
):
    """
    Converts an assembly to a custom mesh format defined by the CadQuery team.
    :param imprint: Whether or not the assembly should be impr
    :type bar: int
    """

    # To keep track of the vertices and triangles in the mesh
    vertices = []
    face_triangles = []
    solids = []
    solid_face_triangle = []
    imprinted_assembly = None
    imprinted_solids_with_orginal_ids = None

    # Imprinted assemblies end up being compounds, whereas you have to step through each of the
    # parts in an assembly and extract the solids.
    if imprint:
        # Imprint the assembly and process it as a compound
        (
            imprinted_assembly,
            imprinted_solids_with_orginal_ids,
        ) = cq.occ_impl.assembly.imprint(self)
        solids.append(imprinted_assembly)
    else:
        # Step through every child in the assembly and save their solids
        for child in self.children:
            # Make sure we end up with a base shape
            obj = child.obj
            if type(child.obj).__name__ == "Workplane":
                solids.append(obj.val())
            else:
                solids.append(obj)

    # Step through all of the collected solids and their respective faces to get the vertices
    for solid in solids:
        for face in solid.faces():
            # Location information of the face to place the vertices and edges correctly
            loc = TopLoc_Location()

            # Perform the tessellation
            BRepMesh_IncrementalMesh(face.wrapped, tolerance, True, angular_tolerance)
            face_mesh = BRep_Tool.Triangulation_s(face.wrapped, loc)

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

    return {"vertices": vertices, "solid_face_triangle_vertex_map": solid_face_triangle}


# Patch the function(s) into the Workplane class
cq.Assembly.toMesh = to_mesh
cq.Assembly.to_mesh = to_mesh
