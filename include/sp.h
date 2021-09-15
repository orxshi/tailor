#ifndef SP_H
#define	SP_H

#include "regular_mesh.h"
#include <boost/serialization/deque.hpp>
#include <boost/serialization/vector.hpp>
#include "profiler.h"
#include "stencil_walk.h"
#include "bc.h"

namespace Tailor
{
    class SpatialPartition
    {
        friend class DonorInfo;
        friend class DirectCutter;
        friend class DonorSearcher;
        friend class HoleProfiler;
        friend class SpCellMover;
        friend class Disconnecter;
        friend class ArrivalConnecter;
        friend class Solver;
        friend class Assembler;
        friend class DonorRequestExchanger;
        friend class TCRExchanger;
        friend class boost::serialization::access;

        public:

        void set_comm(boost::mpi::communicator* comm);
        void set_profiler(Profiler* prof);

        void set_nei(std::deque<BinRMTag> nei);

        const std::deque<BinRMTag>& nei() const;
        void set_mesh_priority(std::deque<BinRMTag> nei);

        void set_mesh_priority(int mtag, int pri);
        void sort_mesh();

        void determine_non_residents();

        void connect_partition_cells(ArrCon<Nei>& arrival_cell, std::function<bool(const Vector3&, int celltag)> is_resi, Profiler*);
        //void update_ghost_primitives(const ArrCon<Var>& arrival_cell);
        void oga_interpolate(const ArrCon<Var>& arrival_cell);
        void convert_to_fortran() const;
        //void solve(Solver& solver, const Vector3& vel, std::function<void(CellExchanger& cell_exchanger)> exchange_ghosts);
        //void solve(Solver& solver, const Vector3& vel, std::function<void()> exchange_ghosts, bool called_exc_ghosts);
        void init_sod();
        //void init_interior();
        //void init_dirichlet();
        void init_flow();
        void connect_after_exchange(std::function<bool(const Vector3&)> is_resi, Profiler*, std::string);
        void disconnect_orphan_faces();
        void update_prev_cand_donor();
        const std::vector<RegularMesh>& rm() const;
        const Tag& tag() const;
        const AABB& aabb() const;
        const std::deque<Mesh>& mesh() const;
        std::deque<Mesh>& mesh_p();
        const Mesh& mesh(const Tag& t) const;
        const Mesh& mesh(int i) const;
        const RegularMesh& rm(const Tag& mt) const;

        void add_merge_mesh_leave_dups(const SpatialPartition& other_sp);
        void remove_dups(std::string name, Profiler*);
        void remove_ghosts();
        void remove_nonresidents();

        // Load
        //const std::vector<double>& mesh_load() const;
        //double load_solver(double& total_ml, double& total_nwall, double& total_nouter, double& total_npart);
        //double load(LoadEstimType load_estim_type, double& total_ml, double& total_nwall, double& total_nouter, double& total_npart);

        void reset_erase_marks();
        void reset_all_oga_status();
        void update_connectivity(int rank);

        //bool is_resident(const Vector3& cnt, const Outline& outline) const;
        //bool is_resident(const Vector3& cnt) const;
        std::deque<std::vector<Point>> mesh_pts(const Outline& outline) const;
        std::deque<std::vector<Point>> mesh_pts() const;
        bool fully_resident(const MeshCell& mc) const;
        void connect_ghosts_and_neis(bool verbose=false);

        void extend_aabb(const AABB& aabb);

        void set_tag(const Tag& t);
        //void set_comm(boost::mpi::communicator*);

        // Regular map
        void make_regular_maps_for_mesh();
        void make_regular_maps_for_mesh_uniproc(Profiler* profiler);
        void make_regular_maps_for_mesh_serial();

        // Connect cells
        void connect_cells(std::function<bool(const Vector3&)> is_resi, Profiler* profiler, std::string name);

        // Move mesh
        void rotate_mesh(const Tag& _parent_mesh, double angle, int axis, const Vector3& rot_point);
        void move_mesh(const Tag& _parent_mesh, const Vector3& v);

        // Add mesh
        void add_mesh(const Mesh& m);
        void add_mesh(Mesh&& m);
        void add_mesh(const SpatialPartition& other_sp);

        // Merge
        void merge(const SpatialPartition& other_sp);
        void merge_mesh(const SpatialPartition& other_sp);
        void merge_meshes(std::deque<Mesh>& mesh) const;
        void merge_meshes_no_check(std::deque<Mesh>& mesh) const;

        // Setter
        //void set_slave_comm(const boost::mpi::communicator* comm);
        //void set_profiler(Profiler* profiler);
        void set_aabb(const AABB& aabb);

        void remove_dup_cells_and_points();
        //void dummy(LoadEstim load_estim_type, double& total_ml, double& total_nwall, double& total_nouter, double& total_npart);
        void set_neighborship(const Vector3Int& global_nstripe, std::vector<unsigned int>& bin_to_proc);
        void remove_mesh(Tag mt);

        void mark_to_be_erased(const Tag& im, const Tag& ic);
        void erase_marked_cells();

        //std::deque<AABB> mesh_aabb(const Outline& outline) const;
        //std::deque<AABB> mesh_aabb() const;

        void make_cell_adt();
        void make_mesh_aabb();
        const std::deque<ADT>& cell_adt() const;
        const std::deque<AABB>& mesh_aabb() const;
        const ADT& cell_adt(int i) const;
        const AABB& mesh_aabb(int i) const;
        void reset_mesh_connectivity();

        template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            ar & mesh_;
            ar & aabb_;
            ar & tag_;
            ar & nei_;

            assert(!std::isnan(tag_()));
        }

        private:
        std::deque<BinRMTag> nei_;
        std::deque<ADT> cell_adt_;
        std::deque<AABB> mesh_aabb_;
        Profiler* profiler_;
        Tag tag_; // same as bintag;
        boost::mpi::communicator* comm_;
        //boost::mpi::communicator* slave_;
        std::vector<RegularMesh> rm_; 
        std::deque<Mesh> mesh_;
        AABB aabb_;
        //std::vector<double> mesh_load_;

        Mesh& mesh_p(const Tag& t);
        RegularMesh& rm_p(const Tag& mt);
        Tag generate_meshtag() const;
        void make_regular_map_for_mesh(const Tag& mt);
        void make_regular_map_for_mesh_uniproc(const Tag& mt, Profiler* profiler);
        void make_regular_map_for_mesh_serial(const Tag& mt);
        void remove_cell(const Tag& im, const Tag& ic);
        void profstart(std::string fun);
        void profstop(std::string fun);
    };

    void donor_info(int rank, const SpatialPartition& sp, int meshtag);
    double solver_load(const SpatialPartition& sp);
    double field_load(const SpatialPartition& sp);
}

#endif
