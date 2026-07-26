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


def test_imprinted_shared_face_shares_id():
    """
    Makes sure an imprinted face that is shared by two solids is given the same
    face id in both of them.

    Imprinting makes the interface between two touching solids the same
    underlying face. Consumers such as DAGMC exporters build one surface per
    face id, so if each solid reports its own id for that interface they end up
    with two coincident surfaces instead of one surface shared by two solids,
    and the resulting model is not watertight there.
    """

    # Two cuboids sharing a 10 x 10 face
    assy = cq.Assembly()
    assy.add(cq.Workplane().box(10, 10, 10))
    assy.add(cq.Workplane().transformed(offset=(0, 7, 0)).box(10, 4, 10))

    mesh = assy.toMesh(imprint=True)
    face_map = mesh["solid_face_triangle_vertex_map"]

    # Both solids should still report six faces each
    assert len(face_map) == 2
    assert len(face_map[1]) == 6
    assert len(face_map[2]) == 6

    # Collect the solids each face id appears in
    solids_by_face_id = {}
    for solid_id, faces in face_map.items():
        for face_id in faces:
            solids_by_face_id.setdefault(face_id, []).append(solid_id)

    shared = [f for f, solids in solids_by_face_id.items() if len(solids) > 1]

    # Exactly one face is shared, so the two solids use 11 ids rather than 12
    assert len(shared) == 1
    assert len(solids_by_face_id) == 11
    assert sorted(solids_by_face_id[shared[0]]) == [1, 2]
