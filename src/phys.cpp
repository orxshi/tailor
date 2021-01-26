#include "phys.h"

namespace Tailor
{
    double spec_inter_energy(double rho, double p, double gamma)
    {
        return (p / ((gamma - 1.) * rho));
    }

    double spec_kine_energy(double u, double v, double w)
    {
        return (0.5 * (std::pow(u, 2.) + std::pow(v, 2.) + std::pow(w, 2.)));
    }

    double total_energy(double rho, double k, double e)
    {
        return (rho * (k + e));
    }

    double total_enthalpy(double rho, double p, double E)
    {
        return ((E + p) / rho);
    }

    double speed_of_sound(double rho, double p, double gamma)
    {
        return std::sqrt(gamma * p / rho);
    }

    vararray prim_to_cons(const vararray& prim, double gamma)
    {
        double k = spec_kine_energy(prim[1], prim[2], prim[3]);
        double e = spec_inter_energy(prim[0], prim[4], gamma); 

        vararray cons;

        cons[0] = prim[0];
        cons[1] = prim[0] * prim[1];
        cons[2] = prim[0] * prim[2];
        cons[3] = prim[0] * prim[3];
        cons[4] = total_energy(prim[0], k, e);

        return cons;
    }
    
    vararray cons_to_prim(const vararray& cons, double gamma)
    {
        vararray prim;

        prim[0] = cons[0];
        prim[1] = cons[1] / prim[0];
        prim[2] = cons[2] / prim[0];
        prim[3] = cons[3] / prim[0];

        double k = spec_kine_energy(prim[1], prim[2], prim[3]);
        double e = cons[4] / prim[0] - k;
        if (e <= 0. || std::isnan(e))
        {
            std::cout << "cons[0]: " << cons[0] << std::endl;
            std::cout << "cons[1]: " << cons[1] << std::endl;
            std::cout << "cons[2]: " << cons[2] << std::endl;
            std::cout << "cons[3]: " << cons[3] << std::endl;
            std::cout << "cons[4]: " << cons[4] << std::endl;

            std::cout << "prim[0]: " << prim[0] << std::endl;
            std::cout << "prim[1]: " << prim[1] << std::endl;
            std::cout << "prim[2]: " << prim[2] << std::endl;
            std::cout << "prim[3]: " << prim[3] << std::endl;

            std::cout << "k: " << k << std::endl;
            std::cout << "E: " << cons[4] / prim[0] << std::endl;

        }
        assert(e > 0. && !std::isnan(e));

        prim[4] = prim[0] * (gamma - 1.) * e;
        
        return prim;
    }

    vararray calc_flux(double rho, double p, double u, double v, double w, double H, double vfn)
    {
        vararray flx;

        double vnet = u - vfn;

        flx[0] = rho * vnet;
        flx[1] = rho * vnet * vnet + p;
        flx[2] = rho * vnet * v;
        flx[3] = rho * vnet * w;
        flx[4] = rho * vnet * H;

        return flx;
    }

    NormTang::NormTang(vec3<double> n): n(n)
    {
        double nx = n(0);
        double ny = n(1);
        double nz = n(2);

        assert(!std::isnan(nx));
        assert(!std::isnan(ny));
        assert(!std::isnan(nz));

        double tempx = ny*ny + nz*nz;
        double tempy = nz*nz + nx*nx;
        double tempz = nx*nx + ny*ny;

        double lx, ly, lz;
        double mx, my, mz;

        if (tempx >= tempy && tempx >= tempz)
        {
            lx = 0.;
            ly = -nz;
            lz =  ny;
        }
        else if (tempy >= tempx && tempy >= tempz)
        {
            lx = -nz;
            ly =  0.;
            lz =  nx;
        }
        else if (tempz >= tempx && tempz >= tempy )
        {
            lx = -ny;
            ly =  nx;
            lz =  0.;
        }
        else
        {
            assert(false);
        }

        double temp = std::sqrt(lx*lx + ly*ly + lz*lz);
        if (temp == 0.)
        {
            std::cout << "nx: " << nx << std::endl;
            std::cout << "ny: " << ny << std::endl;
            std::cout << "nz: " << nz << std::endl;
        }
        assert(temp != 0.);
        lx = lx/temp;
        ly = ly/temp;
        lz = lz/temp;

        mx = ny*lz - nz*ly;
        my = nz*lx - nx*lz;
        mz = nx*ly - ny*lx;

        double abs_n_cross_l = std::sqrt(mx*mx + my*my + mz*mz);
        mx = mx / abs_n_cross_l;
        my = my / abs_n_cross_l;
        mz = mz / abs_n_cross_l;

        l.set_x(lx);
        l.set_y(ly);
        l.set_z(lz);

        m.set_x(mx);
        m.set_y(my);
        m.set_z(mz);

        //l.set_x(ny);
        //l.set_y(-nx);
        //l.set_z(0.);

        //m.set_x(-nz);
        //m.set_y(0.);
        //m.set_z(nx);
    }
}
