import pytest
import cadquery as cq
import cadquery_direct_mesh_plugin


@pytest.fixture
def basic_assy():
    """
    Creates a single object assembly which can be shared across tests.
    """

    # Create an assembly with a single cube in it
    assy = cq.Assembly()
    cuboid = cq.Workplane().box(10, 10, 10)
    assy.add(cuboid)

    return assy


@pytest.fixture
def basic_multipart_assy():
    """
    Creates an assembly with multiple objects which can be shared across tests.
    """

    # Create a simple test assembly
    assy = cq.Assembly()
    assy.add(
        cq.Workplane("XY").box(1, 1, 1, centered=False),
        name="box1",
        color=cq.Color(0.0, 0.5, 0.0),
        loc=(cq.Location(0, 0, 0)),
    )
    assy.add(
        cq.Workplane("XY").box(2, 2, 2),
        name="box2",
        color=cq.Color(0.0, 0.0, 0.5),
        loc=cq.Location(3, 3, 3),
    )

    return assy


@pytest.fixture
def cylinder_assy():
    """
    Creates an assembly with a single cylinder in it, which can be shared across tests.
    """

    assy = cq.Assembly()
    cyl = cq.Workplane().cylinder(5.0, 10.0)
    assy.add(cyl)

    return assy


@pytest.fixture
def more_faces_assy():
    """
    Creates a slightly more complex assembly with some more interesting features.
    """

    # Create two cubes with rounded edges
    cube_1 = cq.Workplane().box(10, 10, 10).fillet(3.0)
    cube_2 = cq.Workplane().box(5, 5, 5).fillet(1.0)

    # Put the assembly together
    assy = cq.Assembly()
    assy.add(cube_1)
    assy.add(cube_2, loc=cq.Location(0, 0, 7.5))

    return assy


def test_basic_assembly(basic_assy):
    """
    Tests to make sure a basic assembly will work.
    """

    # Call the main conversion method we want to test
    mesh = basic_assy.toMesh()

    # Make sure we have the correct number of vertices
    assert len(mesh["vertices"]) == 8

    # Make sure we have the correct number of faces
    assert len(mesh["solid_face_triangle_vertex_map"][1]) == 6


def test_basic_multipart_assembly(basic_multipart_assy):
    """
    Tests to make sure basic multi-part assemblies work correctly.
    """

    # Mesh the assemby
    mesh = basic_multipart_assy.toMesh(imprint=False)

    # Make sure we have the correct number of vertices
    assert len(mesh["vertices"]) == 16

    # Make sure that we have the correct number of solids
    assert len(mesh["solid_face_triangle_vertex_map"]) == 2

    # Make sure that each of the solids has the correct number of faces
    assert len(mesh["solid_face_triangle_vertex_map"][1]) == 6


def test_more_faces(more_faces_assy):
    """
    Tests to make sure a slightly more challenging model can be meshed.
    """

    # Convert the model to a mesh
    mesh = more_faces_assy.toMesh(imprint=False)

    # Make sure the mesh has the correct number of vertices
    assert len(mesh["vertices"]) > 9400

    # Make sure that we have the correct number of solids
    assert len(mesh["solid_face_triangle_vertex_map"]) == 2

    # Make sure that each of the solids has the correct number of faces
    assert len(mesh["solid_face_triangle_vertex_map"][1]) == 26

    # Reset and do an imprinted mesh
    mesh = more_faces_assy.toMesh(imprint=True)

    # Make sure the mesh has the correct number of vertices
    assert len(mesh["vertices"]) > 9400 and len(mesh["vertices"]) < 9500

    # Make sure that we have the correct number of solids
    assert len(mesh["solid_face_triangle_vertex_map"]) == 2

    # Make sure that each of the solids has the correct number of faces
    assert len(mesh["solid_face_triangle_vertex_map"][1]) == 27


def test_edge_handling(basic_assy, cylinder_assy):
    """
    Tests to make sure edges can be extracted from an assembly.
    """

    # Call the main conversion method we want to test
    mesh = basic_assy.toMesh(imprint=False, include_brep_edges=True)

    # Make sure we have the correct number of edges
    assert len(mesh["solid_brep_edge_segments"][0]) == 12

    # Convert the cylinder assembly to a mesh
    mesh = cylinder_assy.toMesh(imprint=False, include_brep_edges=True)

    # Make sure we have the correct number of edges
    assert len(mesh["solid_brep_edge_segments"][0]) == 127


def test_vertex_handling(basic_assy, cylinder_assy):
    """
    Tests to make sure vertices can be extracted from an assembly.
    """

    # Call the main conversion method we want to test
    mesh = basic_assy.toMesh(imprint=False, include_brep_vertices=True)

    # Make sure we have the correct number of vertices
    assert len(mesh["solid_brep_vertices"][0]) == 8

    # Convert the cylinder assembly to a mesh
    mesh = cylinder_assy.toMesh(imprint=False, include_brep_vertices=True)

    # Make sure we have the correct number of edges
    assert len(mesh["solid_brep_vertices"][0]) == 2


def test_assembly_material_meshing():
    """
    Makes sure that assembly materials make it into the mesh data structure.
    """

    # Build a basic assembly with two cubes of different materials
    cube_1 = cq.Workplane().box(10, 10, 10)
    cube_2 = cq.Workplane().box(5, 5, 5)
    assy = cq.Assembly()
    assy.add(
        cube_1, name="cube_1", color=cq.Color(0.722, 0.451, 0.2, 1.0), material="copper"
    )
    assy.add(cube_2, name="cube_2", material="steel", loc=cq.Location(0, 0, 5))

    # Add two other objects to increase the test coverage
    assy.add(cq.Workplane().box(5, 5, 5).val(), loc=cq.Location(0, 0, -5))
    assy.add(cq.Workplane().rect(5, 5).val())

    # Mesh the assembly without imprinting
    mesh = assy.toMesh(imprint=False)
    imprinted_mesh = assy.toMesh(imprint=True)

    # Make sure that each mode of meshing has the material in the correct place
    assert mesh["solid_materials"][0] == "copper"
    assert mesh["solid_materials"][1] == "steel"
    assert imprinted_mesh["solid_materials"][0] == "copper"
    assert imprinted_mesh["solid_materials"][1] == "steel"


def test_overlapping_spheres_shared_faces():
    """
    Tests that two overlapping spheres produce shared faces with non-zero
    triangles on both solids after imprinting.

    After imprint, BRepMesh only tessellates a shared face for whichever
    solid is meshed first.  The face-hash reuse logic should detect the
    shared face and copy the triangles to the second solid so that both
    volumes are watertight.
    """
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeSphere
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
    from OCP.gp import gp_Pnt

    d, r = 3.0, 5.0

    sphere1_shape = BRepPrimAPI_MakeSphere(gp_Pnt(-d, 0, 0), r).Shape()
    sphere2_shape = BRepPrimAPI_MakeSphere(gp_Pnt(d, 0, 0), r).Shape()
    crescent_shape = BRepAlgoAPI_Cut(sphere2_shape, sphere1_shape).Shape()

    assy = cq.Assembly()
    assy.add(cq.Workplane().add(cq.Solid(sphere1_shape)))
    assy.add(cq.Workplane().add(cq.Solid(crescent_shape)))

    mesh = assy.toMesh(imprint=True, tolerance=0.1, angular_tolerance=0.1)
    tris = mesh["solid_face_triangle_vertex_map"]

    # Both solids must be present
    assert len(tris) == 2

    # Every face on every solid must have at least one triangle
    for solid_id, face_map in tris.items():
        for face_id, triangles in face_map.items():
            assert len(triangles) > 0, (
                f"solid {solid_id}, face {face_id} has 0 triangles"
            )

    # At least one face_id must appear in both solids (the shared face)
    faces_solid1 = set(tris[1].keys())
    faces_solid2 = set(tris[2].keys())
    shared = faces_solid1 & faces_solid2
    assert len(shared) > 0, "No shared faces detected between the two solids"
