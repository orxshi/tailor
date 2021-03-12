#ifndef POES_PHYS_H
#define POES_PHYS_H

#include "matrix.h"

namespace Tailor
{
    double spec_inter_energy(double rho, double p, double gamma);
    double spec_kine_energy(double u, double v, double w);
    double total_energy(double rho, double k, double e);
    double total_enthalpy(double rho, double p, double E);
    double speed_of_sound(double rho, double p, double gamma);
    Vector5 prim_to_cons(const Vector5& prim, double gamma);
    Vector5 cons_to_prim(const Vector5& cons, double gamma);
    void test_phys(const Vector5& cons, double gamma);
    //Vector5 calc_flux(double rho, double p, double u, double v, double w, double H, double qn, double vb, const Vector3& n);
    Vector5 calc_flux(double rho, double p, double u, double v, double w, double H, double vfn);
    
    struct NormTang
    {
        NormTang(Vector3 n);

        Vector3 n;
        Vector3 l;
        Vector3 m;
    };
}

#endif
