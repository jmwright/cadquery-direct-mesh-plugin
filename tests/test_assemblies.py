import cadquery as cq
import cadquery_direct_mesh_plugin
from model_benchmark_zoo import (
    Cuboid,
    TwoTouchingCuboids,
    Cylinder,
    Circulartorus,
    Ellipticaltorus,
    NestedCylinder,
    Sphere,
    NestedSphere,
    Oktavian,
    SimpleTokamak,
)

# from cad_to_dagmc.core import vertices_to_h5m


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
    assert len(mesh["solid_face_triangle_vertex_map"][1]) == 6


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
    assert len(mesh["solid_face_triangle_vertex_map"][1]) == 6


def test_edge_handling():
    """
    Tests to make sure edges can be extracted from an assembly.
    """

    # Create a simple test assembly
    cuboid = Cuboid(width=10)
    assy = cuboid.cadquery_assembly()

    # Call the main conversion method we want to test
    mesh = assy.toMesh(imprint=False, include_brep_edges=True)

    # Make sure we have the correct number of edges
    assert len(mesh["solid_brep_edge_segments"][0]) == 12

    # Test with a cylinder that has non-linear edges
    cyl = cq.Workplane().cylinder(5.0, 10.0)
    assy = cq.Assembly()
    assy.add(cyl)

    # Convert the cylinder assembly to a mesh
    mesh = assy.toMesh(imprint=False, include_brep_edges=True)

    # Make sure we have the correct number of edges
    assert len(mesh["solid_brep_edge_segments"][0]) == 127


def test_vertex_handling():
    """
    Tests to make sure vertices can be extracted from an assembly.
    """

    # Create a simple test assembly
    cuboid = Cuboid(width=10)
    assy = cuboid.cadquery_assembly()

    # Call the main conversion method we want to test
    mesh = assy.toMesh(imprint=False, include_brep_vertices=True)

    # Make sure we have the correct number of vertices
    assert len(mesh["solid_brep_vertices"][0]) == 8

    # # Test with a cylinder that has non-linear edges
    cyl = cq.Workplane().cylinder(5.0, 10.0)
    assy = cq.Assembly()
    assy.add(cyl)

    # Convert the cylinder assembly to a mesh
    mesh = assy.toMesh(imprint=False, include_brep_vertices=True)

    # Make sure we have the correct number of edges
    assert len(mesh["solid_brep_vertices"][0]) == 2


def test_model_zoo_shapes():
    """
    Walks through all of the model zoo shapes and makes sure they mesh properly.
    """

    #
    # Two touching cuboids
    #
    two = TwoTouchingCuboids(width1=10, width2=4)
    assy = two.cadquery_assembly()
    mesh = assy.toMesh(imprint=True)
    assert len(mesh["vertices"]) == 16  # Number of vertices
    assert len(mesh["solid_face_triangle_vertex_map"]) == 2  # Number of solids

    #
    # Cylinder
    #
    cyl = Cylinder(radius=1, height=100)
    assy = cyl.cadquery_assembly()
    mesh = assy.toMesh(imprint=True)
    assert len(mesh["vertices"]) == 378  # Number of vertices
    assert len(mesh["solid_face_triangle_vertex_map"]) == 1  # Number of solids

    #
    # Nested cylinder
    #
    cyl = NestedCylinder()
    assy = cyl.cadquery_assembly()
    mesh = assy.toMesh(imprint=True)
    assert len(mesh["vertices"]) == 756  # Number of vertices
    assert len(mesh["solid_face_triangle_vertex_map"]) == 2  # Number of solids

    #
    # Sphere
    #
    sphere = Sphere(radius=10)
    assy = sphere.cadquery_assembly()
    mesh = assy.toMesh(imprint=True)
    assert len(mesh["vertices"]) == 4002  # Number of vertices
    assert len(mesh["solid_face_triangle_vertex_map"]) == 1  # Number of solids

    #
    # Nested sphere
    #
    sphere = NestedSphere(radius1=10, radius2=1)
    assy = sphere.cadquery_assembly()
    mesh = assy.toMesh(imprint=True)
    assert len(mesh["vertices"]) == 8004  # Number of vertices
    assert len(mesh["solid_face_triangle_vertex_map"]) == 2  # Number of solids

    #
    # Circular torus
    #
    major_radius = 10
    minor_radius = 4
    torus = Circulartorus(major_radius=major_radius, minor_radius=minor_radius)
    assy = torus.cadquery_assembly()
    mesh = assy.toMesh(imprint=True)
    assert len(mesh["vertices"]) == 15876  # Number of vertices
    assert len(mesh["solid_face_triangle_vertex_map"]) == 1  # Number of solids

    #
    # Elliptical torus
    #
    major_radius = 10
    minor_radius1 = 4
    minor_radius2 = 2
    torus = Ellipticaltorus(
        major_radius=major_radius,
        minor_radius1=minor_radius1,
        minor_radius2=minor_radius2,
    )
    assy = torus.cadquery_assembly()
    mesh = assy.toMesh(imprint=True)
    assert len(mesh["vertices"]) == 7057  # Number of vertices
    assert len(mesh["solid_face_triangle_vertex_map"]) == 1  # Number of solids

    #
    # Oktavian
    #
    oct = Oktavian()
    assy = oct.cadquery_assembly()
    mesh = assy.toMesh(imprint=True)
    assert len(mesh["vertices"]) == 14326  # Number of vertices
    assert len(mesh["solid_face_triangle_vertex_map"]) == 2  # Number of solids

    #
    # Simple Tokamak
    #
    tok = SimpleTokamak(
        radius=500, blanket_thicknesses=100, center_column_thicknesses=50
    )
    assy = tok.cadquery_assembly()
    mesh = assy.toMesh(imprint=True)
    assert len(mesh["vertices"]) == 8517  # Number of vertices
    assert len(mesh["solid_face_triangle_vertex_map"]) == 2  # Number of solids
