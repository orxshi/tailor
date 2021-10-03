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

        if (type_ == LimiterType::venka)
        {
            assert(false);
            return venka(mesh, mc, grad);
        }
        else if (type_ == LimiterType::barth_jespersen)
        {
            return barth_jespersen(mesh, mc);
        }
        else
        {
            assert(false);
        }
    }

    std::array<double, NVAR> Limiter::venka(const Mesh& mesh, const MeshCell& mc, const std::array<Vector3, NVAR>& grad)
    {
        auto phi = [&] (double x, double eps)
        {
            return (std::pow(x, 2.) + 2. * x + eps) / (std::pow(x, 2.) + x + 2. + eps);
        };

        std::array<double, NVAR> ksif;
        std::array<double, NVAR> ksiv;
        double K = 1e6;

        double eps = std::pow(K * mc.poly().volume(), 3.);

        for (int k=0; k<NVAR; ++k)
        {
            ksiv[k] = TAILOR_BIG_POS_NUM;
        }

        for (const MeshFace& face: mc.face())
        {
            const MeshCell* nei = opposing_nei(mesh, face, mc.tag());
            assert(nei != nullptr);

            auto dis = face.face().centroid() - mc.poly().centroid();

            for (int k=0; k<NVAR; ++k)
            {
                double tmp = dot(grad[k], dis);

                if (tmp > 0.)
                {
                    ksif[k] = phi((std::max(nei->prim(k), mc.prim(k) ) - mc.prim(k)) / tmp, eps);
                }
                else if (tmp < 0.)
                {
                    ksif[k] = phi((std::min(nei->prim(k), mc.prim(k) ) - mc.prim(k)) / tmp, eps);
                }
                else
                {
                    ksif[k] = 1.;
                }

                ksiv[k] = std::min(ksiv[k], ksif[k]);
            }
        }

        return ksiv;
    }

    std::array<double, NVAR> Limiter::barth_jespersen(const Mesh& mesh, const MeshCell& mc)
    {
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

            for (int k=0; k<NVAR; ++k)
            {
                double dif = nei->prim(k) - mc.prim(k);

                if (dif > 0.)
                {
                    phi[k] = std::min(phi[k], std::min(1., max_dif[k] / dif));
                }
                else if (dif < 0.)
                {
                    phi[k] = std::min(phi[k], std::min(1., min_dif[k] / dif));
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
