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
