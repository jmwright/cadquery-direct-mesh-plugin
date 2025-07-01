# cadquery-direct-mesh-plugin

This plugin is meant to allow users to go directly from a CadQuery assembly to a mesh without relying on Gmsh.

## Goals

* Should not require gmsh or an external meshing library.
* As a baseline, it should work with [this method](https://github.com/fusion-energy/cad_to_dagmc/blob/5db93e0e34ab75da54e76b750efec2f718d08888/src/cad_to_dagmc/core.py#L66) from `cad_to_dagmc`.
* Prevent duplication of vertices between shared edges between faces.
* Allow imprinting of assemblies as an option.
* Needs to respect surface sense, both via vertex winding and possibly via a moab tag.
* Should be able to handle all of the test cases outlined in the [model benchmark zoo](https://github.com/fusion-energy/model_benchmark_zoo).
* If duplicate faces are found, needs to handle the flipping of the surface sense.

## Stretch Goals

* Ability to minimize the number of facets generated.
* Generate volume meshes.
* Support [pip-installable moab](https://github.com/shimwell/wheels).

## Related projects

* [CAD_to_OpenMC](https://github.com/united-neux/CAD_to_OpenMC)
* [Fork of CAD_to_OpenMC](https://github.com/united-neux/CAD_to_OpenMC)
* [assembly-mesh-plugin](https://github.com/CadQuery/assembly-mesh-plugin)
