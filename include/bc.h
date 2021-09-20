#ifndef POES_BC_H
#define POES_BC_H

#include "mesh.h"
#include "profiler.h"

namespace Tailor
{
    struct BoundaryCondition
    {
        Freestream fs_;
        bool first_;

        void interog(Mesh& mesh);
        void init_sod(Mesh& mesh);
        void set_bc(Mesh& mesh, Profiler*);
        void empty(Mesh& mesh, MeshCell& mc);
        void farfield(Mesh& mesh);
        void farfield3(Mesh& mesh);
        void farfield4(Mesh& mesh);
        void farfield5(Mesh& mesh);
        void slipwall(Mesh& mesh, MeshCell& mc);
        void symmetry(Mesh& mesh, MeshCell& mc);
        void init_farfield(Mesh& mesh, const Vector5& prim);
        //void set_dirichlet(Mesh& mesh, const Vector5& prim);
        void set_dirichlet(Mesh& mesh);
        void update_fs(const Freestream& fs);
        Vector5 read_dirichlet(const Tag& meshtag);

        BoundaryCondition();

        template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            ar & fs_;
            ar & first_;
        }
    };

    Matrix5 make_rot_matrix(const Vector3& n, bool verbose=false);
}


#endif
