import cadquery as cq
import cadquery_direct_mesh_plugin
from model_benchmark_zoo import Cuboid


def test_basic_assembly():
    """
    Tests to make sure a basic assembly will work.
    """

    # Create a simple test assembly
    cuboid = Cuboid(width=10)
    assy = cuboid.cadquery_assembly()

    # Call the main conversion method we want to test
    mesh = assy.toMesh()

    # Make sure we have the correct number of vertices
    assert len(mesh["vertices"]) == 8

    # Make sure we have the correct number of faces
    assert len(mesh["solid_face_triangle_vertex_map"][0]) == 6


def test_basic_multipart_assembly():
    """
    Tests to make sure basic multi-part assemblies work correctly.
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

    # Mesh the assemby
    mesh = assy.toMesh(imprint=False)

    # Make sure we have the correct number of vertices
    assert len(mesh["vertices"]) == 16

    # Make sure that we have the correct number of solids
    assert len(mesh["solid_face_triangle_vertex_map"]) == 2

    # Make sure that each of the solids has the correct number of faces
    assert len(mesh["solid_face_triangle_vertex_map"][0]) == 6
