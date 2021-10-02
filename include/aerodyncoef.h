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

        AeroCoefPara():
            rho_ref(-1.),
            p_ref(-1.),
            u_ref(-1.),
            area_ref(-1.),
            moment_length(-1.),
            moment_center(-1., -1., -1.)
        {
        }
    };

    struct AeroCoef
    {
        std::vector<std::tuple<Vector3, double>> p; // cnt, pres
        Vector3 F;
        Vector3 M;
        double thrust;

        void check(const AeroCoefPara& para, bool compute_pres_coef, bool compute_force_coef, bool compute_moment_coef) const
        {
            if (compute_pres_coef || compute_force_coef || compute_moment_coef)
            {
                assert(para.rho_ref != -1.);
                assert(para.p_ref != -1.);
                assert(para.u_ref != -1.);
            }

            if (compute_force_coef || compute_moment_coef)
            {
                assert(para.area_ref != -1.);
            }

            if (compute_moment_coef)
            {
                assert(para.moment_length != -1.);

                if (para.moment_center(0) == -1.)
                {
                    if (para.moment_center(1) == -1.)
                    {
                        if (para.moment_center(2) == -1.)
                        {
                            assert(false);
                        }
                    }

                }
            }
        }

        void compute_coef(const std::vector<std::tuple<Vector3, Vector3, double, double>>& P, const AeroCoefPara& para, bool compute_pres_coef, bool compute_force_coef, bool compute_moment_coef)
        {
            F = Vector3(0., 0., 0.);
            M = Vector3(0., 0., 0.);
            p.reserve(P.size());

            check(para, compute_pres_coef, compute_force_coef, compute_moment_coef);

            double p_ref = para.p_ref;
            double moment_length = para.moment_length; 
            Vector3 moment_center = para.moment_center;

            if (compute_pres_coef)
            {
                double dpp = para.dpp();

                for (int i=0; i<P.size(); ++i)
                {
                    auto [cnt, normal, abs_area, pres] = P[i];

                    p.push_back(std::make_tuple(cnt, (pres - p_ref) / dpp));
                }
            }

            if (compute_force_coef)
            {
                for (int i=0; i<P.size(); ++i)
                {
                    auto [cnt, normal, abs_area, pres] = P[i];

                    Vector3 f = normal * pres * abs_area;
                    F = F + f;

                    if (compute_moment_coef) 
                    {
                        M = M + cross(cnt - moment_center, f);
                    }
                }

                double dpf = para.dpf();
                F = F / dpf;
                thrust = std::sqrt(F(0) * F(0) + F(1) * F(1) + F(2) * F(2)) / dpf;

                if (compute_moment_coef) 
                {
                    M = M / para.dpm();
                }
            }
        }
    };
}

#endif
