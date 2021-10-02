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
     - :math:`1.25^\circ`
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

Figures below show the mesh for the airfoil.

.. image:: ../images/transonic-airfoil-mesh-far.png
  :width: 200

.. image:: ../images/transonic-airfoil-mesh-mid.png
  :width: 200

.. image:: ../images/transonic-airfoil-mesh-close.png
  :width: 200

Figure below shows convergence history.

.. image:: ../images/transonic-airfoil-converge.png
  :width: 200

Figure below shows pressure coefficients at the airfoil surface.

.. image:: ../images/transonic-airfoil-pc.png
  :width: 200

Reference 1: Manzano, Luis, Jason Lassaline, and David Zingg. "A Newton-Krylov algorithm for the Euler equations using unstructured grids." 41st Aerospace Sciences Meeting and Exhibit. 2003.

Reference 2: https://su2code.github.io/tutorials/Inviscid_2D_Unconstrained_NACA0012/


See :download:`pdf <https://github.com/orxshi/tailor/blob/main/test/airfoil_static_single_mesh/msh/NACA0012_O.geo>`

Settings are read from the settings.ini: ::

   [tailor]
   mesh_folder = msh/32/unstr/NACA0012_O
   save_folder = sv-unsteady
   save_interval = 4
   max_time_step = 1
   save = false
   load = false
   load_folder = sv-unsteady3418
   profiler = false
   solver = true
   compute-pres-coef = true

   [general]
   pseudo3D = true

   [assembler]
   can-rebalance = false
   make-load-balance = true
   load-estim = 2
   merge-bins = true
   donor-search-algo = 1
   print-map = true
   print-repart-info = true
   print-imbalance = true
   print-pre-vtk = true

   [solver]
   increase_cfl = false
   cfl_multiplier = 2.0
   cfl_ratio = 10
   can-rebalance = true
   steady = true
   use-local-time-step = true
   dt = 1e-3
   tol = 1e-12
   sorder = 2
   torder = 1
   cfl = 10
   dual-ts = false
   riemann-solver = 0
   #temporal_discretization = runge_kutta_4
   temporal_discretization = backward_euler
   #temporal_discretization = forward_euler
   load-estim = 2
   max-time-step = 40000
   make-load-balance = true
   print-repart-info = false
   print-imbalance = false
   repart-ratio = 1000
   print-vtk-only-last-step = true

   [linear-solver]
   #max-iteration = 1000
   #max-restart = 10
   #abs-error = 1e-15
   #rel-error = 1e-15
   print-error = false

   [loadcalc]
   arearep = 0

   [loadmap]
   refine-tol = 10
   print-dev = false
   printlm = false

   [partition]
   print-cell-dist = false
   print-bin-dist = false
   print-mesh-system-size = false


   

Oscillating airfoil
-------------------
