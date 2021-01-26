#ifndef POES_LIMITER_H
#define POES_LIMITER_H

#include "mesh.h"
#include "gradient.h"

enum class LimiterType
{
    undefined = -1,
    minmod = 0,
    valbada = 1,
    bj = 2,
    venka = 3,
};

namespace Tailor
{
    class Limiter
    {
        public:

            Limiter(LimiterType);

            LimiterType type() const;
            void limit(const Mesh& mesh, const MeshCell& mc, const std::array<vec3<double>, NVAR>& grad);
            double operator()(int var) const;
            void limit_prim_cons(Mesh& mesh);

        private:

            LimiterType type_;
            std::array<double, NVAR> data_;

            std::array<double, NVAR> venka(const Mesh& mesh, const MeshCell& mc, const std::array<vec3<double>, NVAR>& grad);

    };
}

#endif
