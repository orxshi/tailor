#include "limiter.h"

namespace Tailor
{
    Limiter::Limiter(LimiterType lt): type_(lt)
    {
    }


    LimiterType Limiter::type() const
    {
        return type_;
    }

    std::array<double, NVAR> Limiter::limit(const Mesh& mesh, const MeshCell& mc, const std::array<Vector3, NVAR>& grad)
    {
        assert(type_ != LimiterType::undefined);

        if (type_ == LimiterType::none)
        {
            std::array<double, NVAR> phi;

            phi[0] = 1.;
            phi[1] = 1.;
            phi[2] = 1.;
            phi[3] = 1.;
            phi[4] = 1.;

            return phi;
        }
        else if (type_ == LimiterType::venkatakrishnan)
        {
            return venkatakrishnan(mesh, mc, grad);
        }
        else if (type_ == LimiterType::barth_jespersen)
        {
            return barth_jespersen(mesh, mc, grad);
        }
        else
        {
            assert(false);
        }
    }

    std::array<double, NVAR> Limiter::venkatakrishnan(const Mesh& mesh, const MeshCell& mc, const std::array<Vector3, NVAR>& grad)
    {
        // See page 3 of http://tetra.mech.ubc.ca/ANSLab/publications/michalak2008.pdf

        auto foo = [&] (double x)
        {
            return (std::pow(x, 2.) + 2. * x) / (std::pow(x, 2.) + x + 2.);
        };

        std::array<double, NVAR> phi;
        std::array<double, NVAR> max_dif;
        std::array<double, NVAR> min_dif;

        for (int k=0; k<NVAR; ++k)
        {
            max_dif[k] = TAILOR_BIG_NEG_NUM;
            min_dif[k] = TAILOR_BIG_POS_NUM;
            phi[k] = TAILOR_BIG_POS_NUM;
        }

        for (const MeshFace& face: mc.face())
        {
            const MeshCell* nei = opposing_nei(mesh, face, mc.tag());
            assert(nei != nullptr);

            for (int k=0; k<NVAR; ++k)
            {
                double dif = nei->prim(k) - mc.prim(k);

                max_dif[k] = std::max(max_dif[k], dif);
                min_dif[k] = std::min(min_dif[k], dif);
            }
        }

        for (const MeshFace& face: mc.face())
        {
            const MeshCell* nei = opposing_nei(mesh, face, mc.tag());
            assert(nei != nullptr);

            auto dis = nei->poly().centroid() - mc.poly().centroid();

            if (face.is_boundary())
            {
                dis *= 2.;
            }

            for (int k=0; k<NVAR; ++k)
            {
                double urecon = mc.prim(k) + dot(grad[k], dis);  // unconstrained reconstructed value.

                if (urecon > 0.)
                {
                    phi[k] = std::min(phi[k], foo(max_dif[k] / urecon));
                }
                else if (urecon < 0.)
                {
                    phi[k] = std::min(phi[k], foo(min_dif[k] / urecon));
                }
                else
                {
                    phi[k] = std::min(phi[k], 1.);
                }
            }
        }

        return phi;
    }

    std::array<double, NVAR> Limiter::barth_jespersen(const Mesh& mesh, const MeshCell& mc, const std::array<Vector3, NVAR>& grad)
    {
        // See page 3 of http://tetra.mech.ubc.ca/ANSLab/publications/michalak2008.pdf

        std::array<double, NVAR> phi;
        std::array<double, NVAR> max_dif;
        std::array<double, NVAR> min_dif;

        for (int k=0; k<NVAR; ++k)
        {
            max_dif[k] = TAILOR_BIG_NEG_NUM;
            min_dif[k] = TAILOR_BIG_POS_NUM;
            phi[k] = TAILOR_BIG_POS_NUM;
        }

        for (const MeshFace& face: mc.face())
        {
            const MeshCell* nei = opposing_nei(mesh, face, mc.tag());
            assert(nei != nullptr);

            for (int k=0; k<NVAR; ++k)
            {
                double dif = nei->prim(k) - mc.prim(k);

                max_dif[k] = std::max(max_dif[k], dif);
                min_dif[k] = std::min(min_dif[k], dif);
            }
        }

        for (const MeshFace& face: mc.face())
        {
            const MeshCell* nei = opposing_nei(mesh, face, mc.tag());
            assert(nei != nullptr);

            auto dis = nei->poly().centroid() - mc.poly().centroid();

            if (face.is_boundary())
            {
                dis *= 2.;
            }

            for (int k=0; k<NVAR; ++k)
            {
                double urecon = mc.prim(k) + dot(grad[k], dis);  // unconstrained reconstructed value.

                if (urecon > 0.)
                {
                    phi[k] = std::min(phi[k], std::min(1., max_dif[k] / urecon));
                }
                else if (urecon < 0.)
                {
                    phi[k] = std::min(phi[k], std::min(1., min_dif[k] / urecon));
                }
                else
                {
                    phi[k] = std::min(phi[k], 1.);
                }
            }
        }

        return phi;
    }
}
