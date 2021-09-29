Flow solver
===========

The Euler equations model fluid flow problems when inertial or convective force dominates the flow and the effect of viscosity through solid-fluid and fluid-fluid is negligible. Differential form of the time dependent and three dimensional Euler equations is

.. math::

   \mathbf{U}_t + \mathbf{F}_x + \mathbf{G}_y + \mathbf{H}_z = 0,

where, as shown below, :math:`\mathbf{U}` is vector of conserved variables and :math:`\mathbf{F}`, :math:`\mathbf{G}` and :math:`\mathbf{H}` are vectors of fluxes in x-, y- and z- directions, respectively.

.. math::

    \mathbf{U}
    =
    \left[
        \begin{matrix}
            \rho \\
            \rho u \\
            \rho v \\
            \rho w \\
            E
        \end{matrix}
    \right],
    \quad
    \mathbf{F}
    =
    \left[
        \begin{matrix}
            \rho u \\
            \rho u^2 + p \\
            \rho uv \\
            \rho uw \\
            u(E+p)
        \end{matrix}
    \right],
    \quad
    \mathbf{G}
    =
    \left[
        \begin{matrix}
            \rho v \\
            \rho uv \\
            \rho v^2 + p \\
            \rho vw \\
            v(E+p)
        \end{matrix}
    \right],
    \quad
    \mathbf{H}
    =
    \left[
        \begin{matrix}
            \rho w \\
            \rho uw \\
            \rho vw \\
            \rho w^2 + p \\
            w(E+p)
        \end{matrix}
    \right].

