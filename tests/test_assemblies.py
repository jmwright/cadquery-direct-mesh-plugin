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
