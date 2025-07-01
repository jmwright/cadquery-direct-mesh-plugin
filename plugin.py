from OCP.TopLoc import TopLoc_Location
from OCP.BRep import BRep_Tool
from OCP.BRepMesh import BRepMesh_IncrementalMesh
import cadquery as cq


# vertices: list[tuple[float, float, float]] | list["cadquery.occ_impl.geom.Vector"],
# triangles_by_solid_by_face: list[list[tuple[int, int, int]]],
# material_tags: list[str],
# h5m_filename: str = "dagmc.h5m",
# implicit_complement_material_tag: str | None = None,


def to_mesh(self, imprint=False, reverse_winding=False, tolerance=0.1, angular_tolerance=0.1):
    """
    Converts an assembly to a custom mesh format defined by the CadQuery team.
    :param imprint: Whether or not the assembly should be impr
    :type bar: int
    """

    # To keep track of the vertices and triangles in the mesh
    vertices = []
    triangles = []
    solids = []
    triangles_by_solid_by_face = []  # : list[list[tuple[int, int, int]]]

    # Imprinted assemblies end up being compounds, whereas you have to step through each of the
    # parts in an assembly and extract the solids.
    if imprint:
        pass
        # Save the compound in the solids list here
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

            # Step through all the triangle vertices
            temp_tris = None
            for i in range(1, face_mesh.NbNodes() + 1):
                v = face_mesh.Node(i)
                v_trsf = v.Transformed(Trsf)
                temp_tris = (v_trsf.X(), v_trsf.Y(), v_trsf.Z())

                # Append the vertices for this face and triangle
                vertices.append(temp_tris)

    print(vertices)

# Patch the function(s) into the Workplane class
cq.Assembly.toMesh = to_mesh
