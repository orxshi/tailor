#ifndef TAILOR_GRADIENT_H
#define TAILOR_GRADIENT_H

#include "mesh.h"

namespace Tailor
{
    class Gradient
    {
        public:

            static std::array<Vector3, NVAR> ls_grad(const Mesh& mesh, const MeshCell& mc);
            static void calc_ls_coef(Mesh& mesh);
    };
}

#endif
