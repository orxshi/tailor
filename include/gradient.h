#ifndef POES_GRADIENT_H
#define POES_GRADIENT_H

#include "mesh.h"
#include "meshcell.h"

namespace Tailor
{
    class Gradient
    {
        public:

            std::array<Vector3, NVAR> ls_grad(const Mesh& mesh, const MeshCell& mc);
            void calc_ls_coef(Mesh& mesh);

        private:

            std::array<Vector3, NVAR> data_;
    };
}

#endif
