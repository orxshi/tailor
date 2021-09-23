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
            double d = 0.5 * rho_ref * u_ref * u_ref;
            assert(d != 0.);
            assert(!isnan(d));

            return 0.5 * rho_ref * u_ref * u_ref;
        }
        double dpf() const
        {
            double d = dpp() * area_ref;
            assert(d != 0.);
            assert(!isnan(d));

            return dpp() * area_ref;
        }
        double dpm() const
        {
            double d = dpf() * moment_length;
            assert(d != 0.);
            assert(!isnan(d));

            return dpf() * moment_length;
        }
    };

    struct AeroCoef
    {
        std::vector<std::tuple<Vector3, double>> p; // cnt, pres
        Vector3 F;
        Vector3 M;
        double thrust;

        void compute_coef(const std::vector<std::tuple<Vector3, Vector3, double, double>>& P, const AeroCoefPara& para)
        {
            F = Vector3(0., 0., 0.);
            M = Vector3(0., 0., 0.);
            p.reserve(P.size());

            double p_ref = para.p_ref;
            double dpp = para.dpp();
            double dpf = para.dpf();
            double dpm = para.dpm();
            double moment_length = para.moment_length; 
            Vector3 moment_center = para.moment_center;

            for (int i=0; i<P.size(); ++i)
            {
                auto [cnt, normal, abs_area, pres] = P[i];

                p.push_back(std::make_tuple(cnt, (pres - p_ref) / dpp));
                Vector3 f = normal * pres * abs_area;
                assert(!isnan(f(0)));
                assert(!isnan(f(1)));
                assert(!isnan(f(2)));
                F = F + f;
                M = M + cross(moment_center - moment_length, f);
            }

            F = F / dpf;
            M = M / dpm;

            thrust = std::sqrt(F(0) * F(0) + F(1) * F(1) + F(2) * F(2)) / dpf;
        }
    };
}

#endif
