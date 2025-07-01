import cadquery as cq
import plugin

# Create a simple test assembly
assy = cq.Assembly()
assy.add(cq.Workplane().box(10, 10, 10), name="box_1")

mesh = assy.toMesh()
