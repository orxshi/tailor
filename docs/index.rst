Tailor documentation
====================

Tailor is a spatially load balancing flow solver which can operate on three-dimensional moving overset meshes.

Why use overset grid technique
------------------------------

The left-most figure below is a helicopter model consisting of six components: A fuselage, a hub and four blades. It is difficult to generate a single structured grid around all the components. It is possible to generate a single unstructured grid around all the components, however, it is usually desired to have structured mesh especially in the boundary layer of components. Even if the single structured mesh is generated with great difficulties, in moving body problems, the single structured would need to be re-generated. It is possible to avoid mesh re-generation by allowing meshes to deform. However, careful mesh deformation techniques should be applied in order to avoid excessive mesh deformation which causes lower solution accuracy.

.. image:: ../images/helicopter.png
  :width: 200

.. image:: ../images/hub_mesh_outline.png
  :width: 200

.. image:: ../images/blade_mesh_outline.png
  :width: 200

How overset grid technique works
--------------------------------

In overset grid technique, a mesh is generated independently for each component. In the case of helicopter, a spherical mesh can be generated for the hub as shown in center figure above. Note that, only outline of the mesh shown for clarity. Similarly, a cylindrical mesh can be generated for a blade as shown in right-most figure above. For the fuselage, a spherical mesh can be used which would also act as a background mesh containing all other meshes.

Overset grid technique versus sliding mesh technique
----------------------------------------------------

In sliding mesh technique meshes cannot overlap but only slide. It is possible to simulate the helicopter model shown above with sliding mesh technique. However, it is impossible to add a tail rotor to the model in sliding mesh technique but there is no such limitation in overset mesh technique.

Dependencies
------------

* Few `Boost`_ libraries such as
    * Boost MPI for parallelization
    * Boost Serialization to save & restore data.
    * Boost Program_options for reading configuration files.
* `METIS`_ for load balancing.
* `Gmsh`_ for mesh generation in msh format.
* `amgcl`_ for solution of linear system of equations if implicit formulation is used.

.. _Boost : https://www.boost.org/
.. _METIS : http://glaros.dtc.umn.edu/gkhome/metis/metis/overview
.. _Gmsh : https://gmsh.info/
.. _amgcl : https://github.com/ddemidov/amgcl

Contents
--------

.. toctree::
   :maxdepth: 2

   test_cases
