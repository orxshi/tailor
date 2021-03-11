#ifndef POES_STATE_H
#define POES_STATE_H

#include "phys.h"
#include "meshface.h"

namespace Tailor
{
    struct State
    {
        //State(const Vector5& conss, double gamma, const Matrix5& T, double vfn);
        //State(const Vector5& conss, double gamma, const Matrix5& T, const MeshCell& mc);
        //State(const Vector5& conss, double gamma, const Matrix5& T, const Vector3& mv, double vfn);
        State(const Vector5& conss, double gamma, const Matrix5& T, const Vector3& mv, double vfn, const Vector3& vf, bool verbose, const Vector3& normal, int rank);

        void print() const;

        double rho, u, v, w, p, e, E, k, H, a;
        //double nx, ny, nz;
        //double lx, ly, lz;
        //double mx, my, mz;
        //double qn, ql, qm;
        Vector5 flux;
        Vector5 cons;
    };
}

#endif
