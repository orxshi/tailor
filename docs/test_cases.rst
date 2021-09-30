Test cases
==========

.. _steady-transonic-airfoil:
Steady transonic airfoil
------------------------

A single unstructured mesh is used to solve the Euler equations on at transonic air speed. The unstructured mesh is of O-shape as shown in the figure below. Cell sizes are set with respect to distance from the airfoil surface. The unstructured mesh, in this case, avoids alongated cells which prevents residuals to converge below 1e-6. Therefore, an unstructured mesh is chosen instead of a structured mesh. There x cells in the unstructured mesh.

Boundary conditions on the airfoil and in outer boundary are slip-wall and Riemann, respectively. Riemann boundary condition is based on Riemann invariant equations. Initially flow is set to freestream values everywhere in the domain.

Flow properties are shown in the following table.

.. list-table:: Flow properties
   :header-rows: 0

   * - Mach
     - 0.8
   * - Angle of attack
     - 1.25 :math:`^\circ`
   * - Freestream pressure
     - 0.7
   * - Freestream density
     - 1
   * - Ratio of specific heats
     - 1.4

.. list-table:: Solver parameters
   :header-rows: 0

   * - Rieman solver
     - Roe
   * - Formulation
     - implicit
   * - Spatial order of accuracy
     - 2
   * - Temporal order of accuracy
     - 1
   * - CFL
     - 10

The code works even when CFL is greater than 10 however, residuals do not converge below 1e-2 in that case.

Figure below shows convergence history.
Figure below shows pressure coefficients at the airfoil surface.

:ref:`steady-transonic-airfoil`


