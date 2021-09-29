Flow solver
===========

The Euler equations model fluid flow problems when inertial or convective force dominates the flow and the effect of viscosity through solid-fluid and fluid-fluid is negligible. Differential form of the time dependent and three dimensional Euler equations is

.. math::

   \mathbf{U}_t + \mathbf{F}_x + \mathbf{G}_y + \mathbf{H}_z = 0,

where, as shown below, :math:`\mathbf{U}` is vector of conserved variables and :math:`\mathbf{F}`, :math:`\mathbf{G}` and :math:`\mathbf{H}` are vectors of fluxes in x-, y- and z- directions, respectively. In the code, 

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

In order to admit discontinuous solutions such as shocks and contact discontinuity, integral form of Euler equations is used:

.. math::

    {\frac{\partial}{\partial t}} \int_V  \mathbf{U}\, \text{d}V
    +
    \int_A \mathcal{H} \cdot \mathbf{n} \, \text{d}A
    =
    0

where, :math:`V` is a control volume, :math:`A` is area of one of the faces of :math:`V`, :math:`n` is the outward unit vector normal to the area :math:`A` and :math:`\mathcal{H}` is tensor of fluxes :math:`\mathcal{H} = (\mathbf{F}, \mathbf{G}, \mathbf{H})`. The term :math:`\mathcal{H} \cdot \mathbf{n} \, \text{d}A` represents flux component normal to the face :math:`A`, and is found by rotating fluxes :math:`\mathbf{F}`, :math:`\mathbf{G}` and :math:`\mathbf{H}` to face-normal direction with the rotation matrix as shown below.

.. math::

    \mathbf{T}
    =
    \left[
        \begin{matrix}
            1 & 0 & 0 & 0 & 0 \\
            0 & \cos\theta^{(y)} \cos\theta^{(z)} & \cos\theta^{(y)} \sin\theta^{(z)} & \sin\theta^{(y)} & 0 \\
            0 & -\sin\theta^{(z)} & \cos\theta^{(z)} & 0 & 0 \\
            0 & -\sin\theta^{(y)} \cos\theta^{(z)} & -\sin\theta^{(y)} \sin\theta^{(z)} & \cos\theta^{(y)} & 0 \\
            0 & 0 & 0 & 0 & 1
        \end{matrix}
    \right]

where, :math:`\theta^{(y)}` and :math:`\theta^{(z)}` are found from components of unit normal vector :math:`[n_x, n_y, n_z]^T` or explicitly

.. math::

    \begin{aligned}
        n_x &= \cos\theta^{(y)} \cos\theta^{(z)}\\
        n_y &= \cos\theta^{(y)} \sin\theta^{(z)}\\
        n_z &= \sin\theta^{(y)}
    \end{aligned}

such that

.. math::

    \begin{aligned}
        \theta^{(z)} &= \text{atan2}(n_y, n_x)\\
        \theta^{(y)} &= \sin^{-1}(n_z)
    \end{aligned}

:math:`\text{atan2}` is another version of :math:`\tan^{-1}` and by taking two arguments instead of one, it returns a unique value unlike :math:`\tan^{-1}`.
