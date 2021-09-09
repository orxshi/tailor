#ifndef AERODYNCOEF_H
#define	AERODYNCOEF_H

#include "matrix.h"

namespace Tailor
{
    struct AeroCoefPara
    {
        double rho_ref;
        double p_ref;
        double u_ref;
        double area_ref;
        double moment_length;
        Vector3 moment_center;

        double dpp() const
        {
            return rho_ref * u_ref * u_ref;
        }
        double dpf() const
        {
            return dpp() * area_ref;
        }
        double dpm() const
        {
            return dpf() * moment_length;
        }
    };

    struct AeroCoef
    {
        std::vector<std::tuple<Vector3, double>> p;
        Vector3 f;
        Vector3 m;
        double thrust;

        void get_force_coef(const Vector3& F, const AeroCoefPara& para)
        {
            for (int i=0; i<TAILOR_N_DIM; ++i)
            {
                f(i) = F(i) / para.dpf();
            }
            thrust = std::sqrt(F(0) * F(0) + F(1) * F(1) + F(2) * F(2)) / para.dpf();
        }
        void get_moment_coef(const Vector3& M, const AeroCoefPara& para)
        {
            for (int i=0; i<TAILOR_N_DIM; ++i)
            {
                m(i) = M(i) / para.dpm();
            }
        }
        void get_pressure_coef(const std::vector<std::tuple<Vector3, double>>& P)
        {
            p = P;
        }
    };
}

#endif
