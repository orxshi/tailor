# Tailor

Tailor is a spatially load balancing flow solver which can operate on three-dimensional moving overset meshes.

# Features
* Can solve moving mesh problems.
* Uses overset mesh technique so that each component has its own mesh.
* Has lower inter-processor communication thanks to geometric/spatial partitioning.
* Makes load balance thanks to hybrid graph/geometric partitioning.
* Can solve the three-dimensional Euler equations.

# Installation

Clone the repository.

`git clone https://github.com/orxshi/tailor.git`

Make a build directory and go to the build directory.

`mkdir build && cd "$_"`

Execute CMake with locations of Boost and METIS if they are not in system directory. The following arguments are just examples.

`cmake .. -DBOOST_ROOT=~/boost -DMETIS_DIR=~/metis`

Build the code.

`make`

The output is the dynamic library `libtailor.so` which can be linked to a driver. Sample driver files are available in test cases. Driver executables can be generated by building driver files with CMakeLists.txt.

# Why use overset mesh technique

The left-most figure below is a helicopter model consisting of six components: A fuselage, a hub and four blades. It is difficult to generate a single structured mesh around all the components. It is possible to generate a single unstructured mesh around all the components, however, it is usually desired to have structured mesh especially in the boundary layer of components. Even if the single structured mesh is generated with great difficulties, in moving body problems, the single structured would need to be re-generated. It is possible to avoid mesh re-generation by allowing meshes to deform. However, careful mesh deformation techniques should be applied in order to avoid excessive mesh deformation which causes lower solution accuracy.

![](https://github.com/orxshi/tailor/blob/main/images/helicopter.png)
![](https://github.com/orxshi/tailor/blob/main/images/hub_mesh_outline.png)
![](https://github.com/orxshi/tailor/blob/main/images/blade_mesh_outline.png)

# How overset mesh technique works

In overset mesh technique, each component is meshed independently. In the case of helicopter, a spherical mesh can be generated for the hub as shown in center figure above. Note that, only outline of the mesh shown for clarity. Similarly, a cylindrical mesh can be generated for a blade as shown in right-most figure above. For the fuselage, a spherical mesh can be used which would also act as a background mesh containing all other meshes. The component meshes, per se, are not useful since data cannot transmit across. The assembler determines the overlapping cells and forms donor-receptor pairs to transmit data. In addition, the assembler determines the invalid cells which overlap invalid regions of space such as solids/holes and exclude those from the mesh-system. In general, the assembler categorizes cells so that the flow solver operates different tasks based on the cell type. For example, there is a category called field/computational cell in which flow solver solves the Euler equations. Categorization involves geometric search operations with Alternating Digital Trees and stencil walking.

# Overset mesh technique versus sliding mesh technique

In sliding mesh technique meshes cannot overlap but only slide. It is possible to simulate the helicopter model shown above with sliding mesh technique. However, it is impossible to add a tail rotor to the model in sliding mesh technique but there is no such limitation in overset mesh technique.

# Partitioning and load balancing

<img src="https://github.com/orxshi/tailor/blob/main/images/partitioning.png" width="400" align="left" />

The figure on the left shows the flowchart for partitioning of mesh-system. In parallel computing environment, data, which is the mesh-system in this case, need to be decomposed into partitions on which processors work on. Tailor expects n partitions of each component mesh to be read by n processors. At this point, processors contain partitions/mesh-blocks that are not overlapping. Since data is transmitted across donor-receptor pairs, the pairs are brought together with geometric/spatial partitioning to avoid repetitive inter-processor communication. Cells are registered to an octree based on their spatial locations. The octree is refined until balanced distribution of partitions to processors is possible. At each refinement step, the octree is converted to a graph which in turn passed on to METIS to find the optimum load distribution.

<br clear="left"/>

# Flow solver
Tailor solves three-dimensional Euler equations with finite-volume formulation. At each cell-face a Riemann problem is approximately solved with either explicit/implicit Roe solver or explicit HLLC solver. To admit discontinuities, integral form of the Euler equations are discretized instead of differential form. Mathematical formulation is based on the book of [Toro](https://www.springer.com/gp/book/9783540252023) where flux is calculated in axis-aligned frame of reference (x-split form). Conservative variables and flux are rotated back-on-forth for consistent flux calculation.

Up to second order spatial and temporal order of accuracy is possible. Gradients are computed by solving a least-squares problem and limited with the Venkatakrishnan multiplier. Time derivative can be discretized with forward and backward time levels and 4-stage Runge-Kutta.

It is also possible to use dual-stepping to convert an unsteady flow problem to a steady one. This method is useful to lower time step size or apply acceleration techniques for steady-state solution.

It is known that HLLC is possitivity-preserving however, it is quite time consuming to derive HLLC Jacobian for the case of implicit Euler equations. Currently, the discretized equations can be solved explicitly if interface flux is computed with HLLC. Only Roe's Riemann solver can be used if implicit solution is desired. I am working on implicit HLLC solver by adapting from [Batten et. al.](https://www.sciencedirect.com/science/article/pii/S0021999197957934).

Tailor can operate on both structured and unstructured meshes but does not make use of the structured indexing in the case of structured mesh.

# Moving mesh

Overset mesh technique shines in moving mesh problems. Component meshes are free to *rigidly* move. The flow solver can operate on moving meshes thanks to Arbitrary Lagrangian Eulerian (ALE) formulation. Mesh velocity is added to Eulerian velocity to accomodate for mesh motion.

# Documentation

I am preparing a [readthedocs](https://tailor.readthedocs.io/en/latest/index.html) documentation to explain the theory and the code. Code documentation will take some time though since a lot of code-clean-up is needed.
# Test cases

I will provide validation results on the following test cases.
* transonic flow past NACA 0012 airfoil (2D, single static mesh)
* flow around oscillating NACA 0012 airfoil (2D, overset moving mesh)
* isentropic vortex (2D, single mesh, to test spatial and temporal accuracy)
* isentropic vortex (2D, overset mesh, to test how well vortex passes through different meshes)
* flow past ONERA M6 wing (3D, single mesh)
* Rotor-fuselage interaction in ROBIN helicopter configuration (3D, overset moving mesh) 

I will also provide parallel performance results on helicopter configuration.

# Dependencies
* Linux operating system (I have tested on Ubuntu and CentOS and haven't tried other operating systems).
* C++ compiler (I tested only on GNU C++ compiler 8.3.0 and 9.2.0).
* Few [Boost](https://www.boost.org/) libraries such as
    * Boost MPI for parallelization
    * Boost Serialization to save & restore data.
    * Boost Program_options for reading configuration files.
* [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) for load balancing.
* [Gmsh](https://gmsh.info/) for mesh generation in msh format.
* [amgcl](https://github.com/ddemidov/amgcl) for solution of linear system of equations if implicit formulation is used.
* [CMake](https://cmake.org/) for build file generation.

# Warning

Tailor is in active development. You will not be happy if you dig into source files or try to execute test cases.
