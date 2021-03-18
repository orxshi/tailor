#ifndef POES_STATE_H
#define POES_STATE_H

#include "phys.h"
#include "meshface.h"

namespace Tailor
{
    struct State
    {
        State(const Vector5& conservative_var, double gamma, const Matrix5& rotation_matrix, const Vector3& face_velocity);

        void print() const;

        double rho, u, v, w, p, e, E, k, H, a;
        Vector5 flux;
        Vector5 cons;
    };
}

#endif
