#include "state.h"

namespace Tailor
{
    //State::State(const vararray& conss, double gamma, const varmat& T, double vfn)
    //State::State(const vararray& conss, double gamma, const varmat& T, const vec3<double>& mv, double vfn)
    State::State(const vararray& conss, double gamma, const varmat& T, const vec3<double>& mv, double vfn, const vec3<double>& vf, bool verbose, const vec3<double>& normal, int rank)
    //State::State(const vararray& conss, double gamma, const varmat& T, const MeshCell& mc)
    {
        //NormTang nt(normal);

        //nx = nt.n(0);
        //ny = nt.n(1);
        //nz = nt.n(2);

        //lx = nt.l(0);
        //ly = nt.l(1);
        //lz = nt.l(2);

        //mx = nt.m(0);
        //my = nt.m(1);
        //mz = nt.m(2);

        cons = T * conss;
        //cons = conss;
        assert(cons[0] == conss[0]);
        assert(cons[4] == conss[4]);

        //if (mc.btype() != BouType::interog)
        {
            //auto mv = mc.vgn();

            cons[1] -= mv(0); 
            cons[2] -= mv(1); 
            cons[3] -= mv(2); 
        }

        if (cons[0] < 0.)
        {
        std::cout << "cons[0]: " << cons[0] << std::endl;
        std::cout << "conss[0]: " << conss[0] << std::endl;
        std::cout << "T: " << T << std::endl;
        }
        auto prim = cons_to_prim(cons, gamma);

        //Matrix<TAILOR_N_DIM, TAILOR_N_DIM> RT;
        //RT(0,0) = T(1,1);
        //RT(0,1) = T(1,2);
        //RT(0,2) = T(1,3);

        //RT(1,0) = T(2,1);
        //RT(1,1) = T(2,2);
        //RT(1,2) = T(2,3);

        //RT(2,0) = T(3,1);
        //RT(2,1) = T(3,2);
        //RT(2,2) = T(3,3);

        //Vector<TAILOR_N_DIM> yaxis;
        //Vector<TAILOR_N_DIM> zaxis;
        //yaxis[0] = 0.;
        //yaxis[1] = 1.;
        //yaxis[2] = 0.;

        //zaxis[0] = 0.;
        //zaxis[1] = 0.;
        //zaxis[2] = 1.;

        vararray vff;
        vff[0] = 0.; 
        vff[1] = vf(0); 
        vff[2] = vf(1); 
        vff[3] = vf(2); 
        vff[4] = 0.; 

        auto rotvf = T * vff;
        auto vfnn = rotvf[1];
        auto vfl = rotvf[2];
        auto vfm = rotvf[3];

        prim[1] -= vfnn;
        prim[2] -= vfl;
        prim[3] -= vfm;

        //if ((rotvf[1] - vfn) > TAILOR_ZERO)
        //{
        //    std::cout << "--------------" << std::endl;
        //    std::cout << "rank: " << rank << std::endl;
        //    std::cout << "vfn: " << vfn << std::endl;
        //    std::cout << "rotvf: " << rotvf[1] << std::endl;
        //    std::cout << "normal(0) : " << normal(0) << std::endl;
        //    std::cout << "normal(1) : " << normal(1) << std::endl;
        //    std::cout << "normal(2) : " << normal(2) << std::endl;
        //    std::cout << "vf(0) : " << vf(0) << std::endl;
        //    std::cout << "vf(1) : " << vf(1) << std::endl;
        //    std::cout << "vf(2) : " << vf(2) << std::endl;
        //    std::cout << "normalize(vf)(0) : " << normalize(vf)(0) << std::endl;
        //    std::cout << "normalize(vf)(1) : " << normalize(vf)(1) << std::endl;
        //    std::cout << "normalize(vf)(2) : " << normalize(vf)(2) << std::endl;
        //    std::cout << "thetaz: " << atan2(normal(1), normal(0)) << std::endl;
        //    std::cout << "thetay: " << asin(normal(2)) << std::endl;
        //    std::cout << "thetay: " << atan2(normal(2) * std::cos(atan2(normal(1), normal(0))), std::abs(normal(0))) << std::endl;
        //    std::cout << "============" << std::endl;

        //}
        //assert((rotvf[1] - vfn) < TAILOR_ZERO);

        //auto tmpy = RT * yaxis;
        //auto tmpz = RT * zaxis;

        //vec3<double> ya(tmpy[0], tmpy[1], tmpy[2]);
        //vec3<double> za(tmpz[0], tmpz[1], tmpz[2]);

        //auto vfl = dotp(vf, ya);
        //auto vfm = dotp(vf, za);

        //if (verbose)
        //{
            //std::cout << "ya: " << ya(0) << " " << ya(1) << " " << ya(2) << std::endl;
            //std::cout << "za: " << za(0) << " " << za(1) << " " << za(2) << std::endl;
            //std::cout << "vf: " << vf(0) << " " << vf(1) << " " << vf(2) << std::endl;
        //}

        rho = prim[0];
        u = prim[1];
        v = prim[2];
        w = prim[3];
        p = prim[4];
        a = std::sqrt(gamma * p / rho);
        e = spec_inter_energy(rho, p, gamma);
        k = spec_kine_energy(u, v, w);
        E = total_energy(rho, k, e);
        H = total_enthalpy(rho, p, E);

        cons = prim_to_cons(prim, gamma);

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
        std::cout << "flux(0): " << flux[0] << std::endl;
        std::cout << "flux(1): " << flux[1] << std::endl;
        std::cout << "flux(2): " << flux[2] << std::endl;
        std::cout << "flux(3): " << flux[3] << std::endl;
        std::cout << "flux(4): " << flux[4] << std::endl;
    }
}
