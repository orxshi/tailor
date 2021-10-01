# Tailor

Tailor is a spatially load balancing flow solver which can operate on three-dimensional moving overset meshes.

## Why use overset mesh technique

The left-most figure below is a helicopter model consisting of six components: A fuselage, a hub and four blades. It is difficult to generate a single structured mesh around all the components. It is possible to generate a single unstructured mesh around all the components, however, it is usually desired to have structured mesh especially in the boundary layer of components. Even if the single structured mesh is generated with great difficulties, in moving body problems, the single structured would need to be re-generated. It is possible to avoid mesh re-generation by allowing meshes to deform. However, careful mesh deformation techniques should be applied in order to avoid excessive mesh deformation which causes lower solution accuracy.

![](https://github.com/orxshi/tailor/blob/main/images/helicopter.png)
![](https://github.com/orxshi/tailor/blob/main/images/hub_mesh_outline.png)
![](https://github.com/orxshi/tailor/blob/main/images/blade_mesh_outline.png)

## How overset mesh technique works

In overset mesh technique, each component is meshed independently. In the case of helicopter, a spherical mesh can be generated for the hub as shown in center figure above. Note that, only outline of the mesh shown for clarity. Similarly, a cylindrical mesh can be generated for a blade as shown in right-most figure above. For the fuselage, a spherical mesh can be used which would also act as a background mesh containing all other meshes. The component meshes, per se, are not useful since data cannot transmit across. The assembler determines the overlapping cells and forms donor-receptor pairs to transmit data. In addition, the assembler determines the invalid cells which overlap invalid regions of space such as solids/holes and exclude those from the mesh-system. In general, the assembler categorizes cells so that the flow solver operates different tasks based on the cell type. For example, there is a category called field/computational cell in which flow solver solves the Euler equations.

## Overset mesh technique versus sliding mesh technique

In sliding mesh technique meshes cannot overlap but only slide. It is possible to simulate the helicopter model shown above with sliding mesh technique. However, it is impossible to add a tail rotor to the model in sliding mesh technique but there is no such limitation in overset mesh technique.

## Partitioning and load balancing

<img src="https://github.com/orxshi/tailor/blob/main/images/partitioning.png" width="400" align="left" />

The figure on the left shows the flowchart for partitioning of mesh-system. In parallel computing environment, data, which is the mesh-system in this case, need to be decomposed into partitions on which processors work on. Tailor expects n partitions of each component mesh to be read by n processors. At this point, processors contain partitions/mesh-blocks that are not overlapping. Since data is transmitted across donor-receptor pairs, the pairs are brought together with geometric/spatial partitioning to avoid repetitive inter-processor communication. The geometric and graph partitioning continues back-and-forth until balanced distribution of partitions to processors is possible.

<br clear="left"/>

## Flow solver
Tailor solves three-dimensional Euler equations. At each interface a Riemann problem is approximately solved with either explicit/implicit Roe solver or explicit HLLC solver. I am working on implicit HLLC solver.

## Documentation

I am preparing a readthedocs documentation to explain the theory and the code. Code documentation will take some time though since a lot of code clean up is needed.

## Test cases

I will provide parallel performance and validation results on pitching transonic airfoil, isentropic vortex, ONERA M6 and helicopter test cases.

## Required libraries
* Few [Boost](https://www.boost.org/) libraries such as
    * Boost MPI for parallelization
    * Boost Serialization to save & restore data.
    * Boost Program_options for reading configuration files.
* [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) for load balancing.
* [Gmsh](https://gmsh.info/) for mesh generation in msh format.
* [amgcl](https://github.com/ddemidov/amgcl) for solution of linear system of equations if implicit formulation is used.

## Tested on
* GNU C++ compiler 8.3.0 and 9.2.0.
