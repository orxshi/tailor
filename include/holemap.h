#ifndef HOLE_MAP_H
#define HOLE_MAP_H

#include "mesh.h"
#include "regular_mesh.h"

namespace Tailor
{
    class HoleMap
    {
        public:
            //HoleMap(const MPI_Comm& comm, const Mesh* wall, const Mesh* dirichlet, const Mesh* farfield, const Mesh* interog, bool pseudo3D);
            HoleMap(const MPI_Comm& comm, const Mesh* mesh, bool pseudo3D, const std::vector<AABB>& hole_aabb);
            bool is_inside_holebin(const vec3<double>& querypoint) const;
            const Tag& tag() const;
            //const AABB& wallaabb() const;
            //bool holeless() const;

        private:
            //void make(const Mesh* wall, const Mesh* dirichlet, const Mesh* farfield, const Mesh* interog, bool pseudo3D);
            void make(const Mesh* mesh, bool pseudo3D);
            Tag tag_;
            boost::mpi::communicator world_;
            AABB wall_aabb_;
            ADT wall_adt_;
            bool only_aabb_;
            //bool holeless_;
            //AABB wallaabb_;
            //AABB outeraabb_;
            std::unique_ptr<RegularMesh> wallrm_;
            std::unique_ptr<RegularMesh> outerrm_;
            //std::vector<Tag> holebin_;

            //bool do_wallaabb_intersect_outer(const Mesh* dirichlet, const Mesh* farfield, const Mesh* interog);
            bool do_wallaabb_intersect_outer(const Mesh& mesh);
            bool do_wallaabb_intersect_outer_(const mcc& container);
    };
}

#endif

