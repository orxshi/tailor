#ifndef POES_STATE_H
#define POES_STATE_H

#include "phys.h"
#include "meshface.h"

namespace Tailor
{
    struct State
    {
        //State(const vararray& conss, double gamma, const varmat& T, double vfn);
        //State(const vararray& conss, double gamma, const varmat& T, const MeshCell& mc);
        //State(const vararray& conss, double gamma, const varmat& T, const vec3<double>& mv, double vfn);
        State(const vararray& conss, double gamma, const varmat& T, const vec3<double>& mv, double vfn, const vec3<double>& vf, bool verbose, const vec3<double>& normal, int rank);

        void print() const;

        double rho, u, v, w, p, e, E, k, H, a;
        //double nx, ny, nz;
        //double lx, ly, lz;
        //double mx, my, mz;
        //double qn, ql, qm;
        vararray flux;
        vararray cons;
    };
}

#endif
