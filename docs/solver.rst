Flow solver
===========

The Euler equations model fluid flow problems when inertial or convective force dominates the flow and the effect of viscosity through solid-fluid and fluid-fluid is negligible. Euler equations have the advantage of admitting discontinuous solutions such as shock wave. In the case of Navier-Stokes equations, shock capturing is cumbersome as shock wave has a definite thickness, therefore, the validity of continuum assumptions becomes questionable. In absense of viscous forces, solution of governing equations simplifies without sacrificing important features of flow. In the case of rotor-fuselage interaction, the Euler equations are satisfactory for predicition pressure estimations. However, the Euler equations have drawbacks in capturing rotor wake structure.

Differential form of the time dependent and three dimensional Euler equations is

.. math::

\mathbf{U}_t + \mathbf{F}_x + \mathbf{G}_y + \mathbf{H}_z = 0
