#include "state.h"

namespace Tailor
{
    State::State(const Vector5& conservative_var, double gamma, const Matrix5& rotation_matrix, const Vector3& face_velocity)
    {
        auto primitive_var = cons_to_prim(conservative_var, gamma);

        primitive_var(1) -= face_velocity(0);
        primitive_var(2) -= face_velocity(1);
        primitive_var(3) -= face_velocity(2);

        primitive_var = rotation_matrix * primitive_var;

        cons = prim_to_cons(primitive_var, gamma);

        rho = primitive_var(0);
        u   = primitive_var(1);
        v   = primitive_var(2);
        w   = primitive_var(3);
        p   = primitive_var(4);
        a   = std::sqrt(gamma * p / rho);
        e   = spec_inter_energy(rho, p, gamma);
        k   = spec_kine_energy(u, v, w);
        E   = total_energy(rho, k, e);
        H   = total_enthalpy(rho, p, E);

        flux = calc_flux(rho, p, u, v, w, H, 0.);
    }

    void State::print() const
    {
        std::cout << "rho: " << rho << std::endl;
        std::cout << "u: " << u << std::endl;
        std::cout << "v: " << v << std::endl;
        std::cout << "w: " << w << std::endl;
        std::cout << "p: " << p << std::endl;
        std::cout << "e: " << e << std::endl;
        std::cout << "k: " << k << std::endl;
        std::cout << "E: " << E << std::endl;
        std::cout << "H: " << H << std::endl;
        std::cout << "a: " << a << std::endl;
        std::cout << "flux(0): " << flux(0) << std::endl;
        std::cout << "flux(1): " << flux(1) << std::endl;
        std::cout << "flux(2): " << flux(2) << std::endl;
        std::cout << "flux(3): " << flux(3) << std::endl;
        std::cout << "flux(4): " << flux(4) << std::endl;
    }
}
