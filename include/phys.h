#ifndef POES_PHYS_H
#define POES_PHYS_H

#include "array.h"
#include "vec3.h"

namespace Tailor
{
    double spec_inter_energy(double rho, double p, double gamma);
    double spec_kine_energy(double u, double v, double w);
    double total_energy(double rho, double k, double e);
    double total_enthalpy(double rho, double p, double E);
    double speed_of_sound(double rho, double p, double gamma);
    vararray prim_to_cons(const vararray& prim, double gamma);
    vararray cons_to_prim(const vararray& cons, double gamma);
    //vararray calc_flux(double rho, double p, double u, double v, double w, double H, double qn, double vb, const vec3<double>& n);
    vararray calc_flux(double rho, double p, double u, double v, double w, double H, double vfn);
    
    struct NormTang
    {
        NormTang(vec3<double> n);

        vec3<double> n;
        vec3<double> l;
        vec3<double> m;
    };
}

#endif
