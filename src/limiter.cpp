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

    void Limiter::limit(const Mesh& mesh, const MeshCell& mc, const std::array<vec3<double>, NVAR>& grad)
    {
        assert(type_ != LimiterType::undefined);

        if (type_ == LimiterType::venka)
        {
            data_ = venka(mesh, mc, grad);
        }
        else
        {
            assert(false);
        }
    }

    double Limiter::operator()(int var) const
    {
        return data_[var];
    }

    std::array<double, NVAR> Limiter::venka(const Mesh& mesh, const MeshCell& mc, const std::array<vec3<double>, NVAR>& grad)
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

            auto dis = face.face().centroid() - mc.poly().centroid();

            for (int k=0; k<NVAR; ++k)
            {
                double tmp = dotp(grad[k], dis);

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
}

