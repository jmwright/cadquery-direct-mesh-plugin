import cadquery as cq
import cadquery_direct_mesh_plugin


def test_basic_assembly():
    """
    Tests to make sure a basic assembly will work.
    """

    # Create a simple test assembly
    assy = cq.Assembly()
    assy.add(cq.Workplane().box(10, 10, 10), name="box_1")

    mesh = assy.toMesh()

    # Make sure we have the correct number of vertices
    assert len(mesh["vertices"]) == 24
