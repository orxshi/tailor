tailor
======

tailor is a spatially load balancing flow solver which can operate on three-dimensional moving overset meshes.

Why use overset grid technique
------------------------------

The left-most figure below is a helicopter model consisting of six components: A fuselage, a hub and four blades. It is difficult to generate a single structured grid around all the components. It is possible to generate a single unstructured grid around all the components, however, it is usually desired to have structured mesh especially in the boundary layer of components. Even if the single structured mesh is generated with great difficulties, in moving body problems, the single structured would need to be re-generated. It is possible to avoid mesh re-generation by allowing meshes to deform. However, careful mesh deformation techniques should be applied in order to avoid excessive mesh deformation which causes lower solution accuracy.

.. image:: ../images/helicopter.png
  :width: 400
  :alt: Alternative text
