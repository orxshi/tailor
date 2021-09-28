# tailor

tailor is a flow solver which can operate on overset meshes.

## Why use overset grid technique

Below is a helicopter model consisting of six components: A fuselage, a hub and four blades. It is difficult to generate a single structured grid around all the components. It is possible to generate a single unstructured grid around all the components, however, it is usually desired to have structured mesh especially in the boundary layer of components. Even if the single structured mesh is generated with great difficulties, in moving body problems, the single structured would need to be re-generated. It is possible to avoid mesh re-generation by allowing meshes to deform. However, careful mesh deformation techniques should be applied in order to avoid excessive mesh deformation which causes lower solution accuracy.

## How overset mesh technique works


In overset mesh technique, overlapping mesh-cells need to be in the processors either temporarily or permanently. In this code, the overlapping mesh-cells are stored in the same partitions. Graph partitioners such as METIS, partitions each component mesh

## Features
* Three-dimensional meshes
* 
3D Overset grid assembler and the Euler equations solver.
Parallellization with OpenMPI.
Spatial load balancing.

# Required libraries
* Boost MPI for parallelization
* Boost Serialization to save & restore data.
* Boost Program_options for reading configuration files.
* METIS for load balancing.

# Tested on
* GNU C++ compiler 
