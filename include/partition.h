#ifndef TAILOR_PARTITION_H
#define TAILOR_PARTITION_H

#include "profiler.h"
#include "request_exchanger.h"
#include "nonresi_exchanger.h"
#include "loadmap.h"
#include "remap.h"
#include "read_mesh.h"
//#include <experimental/filesystem>
//#include <filesystem>

namespace Tailor
{
    class Partition
    {
        friend class Assembler;
        friend class Solver;

        public:

        const boost::mpi::communicator* comm() const;
        void set_comm(boost::mpi::communicator* comm);
        void set_profiler(Profiler* prof);
            Partition(boost::mpi::communicator* comm, LoadEstim, bool pseudo3D, std::string name, Profiler* profiler=nullptr);
            Partition();
            ~Partition();
            void repartition(const std::deque<Mesh>& mesh, bool make_load_balance, RegType regtype, bool merge_bins, bool part_to_commsize_only, LoadEstim, int iter);
            void make(bool merge_bins, bool make_load_balance, int iter);
            SpatialPartitionContainer& spc();
            const SpatialPartitionContainer& spc() const;
            const std::deque<Mesh> mesh() const;
            void read_mesh(const std::vector<std::string>&);
            //const std::deque<Mesh> wall() const;
            //const std::deque<Mesh> dirichlet() const;
            //const std::deque<Mesh> farfield() const;
            //const std::deque<Mesh> interog() const;
            //const std::deque<Mesh> empty() const;
            void rotate(const Tag& mesh, double ang, int axis, const Vector3& pivot);
            void move(const Tag& mesh, const Vector3& v);
            void clear_mesh();
            void set_mesh(std::deque<Mesh>&& mesh);
            std::string name() const;
            const Loadmap* loadmap() const;
            void print_cell_dist(const boost::mpi::communicator* comm, int iter) const;

            template<class Archive> void serialize(Archive & ar, const unsigned int version)
            {
                ar & print_bin_dist_;
                ar & print_mesh_system_size_;
                ar & print_cell_dist_;
                ar & name_;
                ar & pseudo3D_;
                ar & load_estim_type_;
                ar & n_load_balance_;
                ar & spc_;

                assert(print_bin_dist_ == true || print_bin_dist_ == false);
                assert(print_mesh_system_size_ == true || print_mesh_system_size_ == false);
                assert(print_cell_dist_ == true || print_cell_dist_ == false);
                assert(pseudo3D_ == true || pseudo3D_ == false);
            }

        private:

            bool print_bin_dist_;
            bool print_mesh_system_size_;
            bool print_cell_dist_;
            std::string name_;
            bool pseudo3D_;
            LoadEstim load_estim_type_;
            //bool remap_batch_;
            int n_load_balance_;
            boost::mpi::communicator* comm_;
            SpatialPartitionContainer* spc_;
            Loadmap* loadmap_;
            Profiler* profiler_;
            //std::deque<Mesh> wall_;
            //std::deque<Mesh> dirichlet_;
            //std::deque<Mesh> farfield_;
            //std::deque<Mesh> empty_;
            //std::deque<Mesh> interog_;
            std::deque<Mesh> mesh_;

            void make_loadmap_uniproc(const std::deque<Mesh>& mesh);
            void make_loadmap(const std::deque<Mesh>& mesh, bool adaptive, RegType regtype, bool make_load_balance, bool part_to_commsize_only, LoadEstim);
            void remap(bool merge_bins, int iter);
            void connect_cells();
            void read_settings();
            void read_mesh_(BouType, std::string fn);
    };

    std::deque<Mesh> updated_mesh(const Partition& partition, int rank);
    void print_bin_dist(const boost::mpi::communicator* comm, const Partition* partition, int iter);
}

#endif
