#ifndef SPC_H
#define	SPC_H

#include <boost/mpi.hpp>
#include "sp.h"
#include "aerodyncoef.h"
#include "exchanger.h"

namespace Tailor
{
    class SpatialPartitionContainer
    {
        friend class DonorInfo;
        friend class DonorSearcher;
        friend class HoleProfiler;
        friend class CellExchanger;
        friend class Disconnecter;
        friend class ArrivalConnecter;
        friend class Solver;
        friend class Assembler;
        friend class TCRExchanger;

        public:
        void set_comm(boost::mpi::communicator* comm);
        void set_profiler(Profiler* prof);
        int nresi() const;
        void get_coef(const std::vector<AeroCoefPara>& aero_para, int iter) const;
        //bool search_rm(const Vector3& cnt, BinRMTag& brmt, int& rank) const;
        bool search_rm(const Vector3& cnt, Tag& bintag, int& rank, bool verbose=false) const;
        //std::map<BinRMTag, int> search_rm(const AABB& aabb) const;
        //std::map<Tag, int> search_rm(const AABB& aabb, bool verbose=false) const;
        void make_global_adt();
        int mesh_system_size() const;
        void sort_sp_meshes();
        const boost::mpi::communicator* comm() const;
        void remove_ghosts();
        void remove_nonresidents();
        void update_connectivity();
        void reset_all_oga_status();
        bool is_resident(const Vector3& cnt, int celltag, Tag& brmt) const;
        //bool is_resident(const Vector3& cnt, int celltag, Tag& brmt, int& dest_rank) const;
        void connect_partition_cells(ArrCon<Nei>& storage);
        //void update_ghost_primitives(const std::vector<Storage<Var>>& storage);
        //void update_ghost_primitives(const RecvCon<Var>& storage);
        void oga_interpolate(const ArrCon<Var>& stor);
        void convert_to_fortran() const;
        //void solve(Solver& solver, const Vector3& vel, std::function<void(CellExchanger& cell_exchanger)> exchange_ghosts);
        //void solve(Solver& solver, const Vector3& vel, std::function<void()> exchange_ghosts);
        void init_sod();
        //void init_interior();
        //void init_dirichlet();
        void init_flow(FlowInitType flow_init_type);
        void connect_after_exchange(Profiler* profiler, std::string proc);
        void update_prev_cand_donor();
            //void make_outline(std::vector<std::vector<Vector3>> allpts);
            const std::map<int, int>& bintag_proc_map() const;
            //const std::map<int, int>& bintag_proc_map() const;

            //SpatialPartitionContainer(boost::mpi::communicator* comm, Profiler* profiler, bool verbose, const std::map<BinRMTag, int>& bintag_proc_map, int nmesh);
            SpatialPartitionContainer(boost::mpi::communicator* comm, Profiler* profiler, bool verbose, const std::map<int, int>& bintag_proc_map, int nmesh);
            SpatialPartitionContainer(boost::mpi::communicator* comm, Profiler* profiler, bool verbose, int nmesh);
            SpatialPartitionContainer();
            //~SpatialPartitionContainer();

            //void output_load(int iteration);
            
            //void reduce_sp(const std::vector<std::vector<BinRMTag>>& dist, const std::vector<Storage<SpatialPartition>>& storage);
            void reduce_sp(const std::vector<std::vector<BinRMTag>>& dist, const ArrCon<SpatialPartition>& storage, std::string name);

            void make_mesh_aabbs();

            // Connect
            void connect_cells(std::string name);

            // Getter
            const std::vector<SpatialPartition>& sp() const;
            //void info() const;
            const RegularMesh& global_rm() const;
            RegularMesh& global_rm();
            const ADT& global_adt() const;

            // Make rm
            void make_regular_maps();
            void make_regular_maps_uniproc();

            // Move mesh
            void rotate_meshblocks(const Tag& _parent_mesh, double ang, int axis, const Vector3& rot_axis);
            void move_meshblocks(const Tag& _parent_mesh, const Vector3& v);

            // Remap
            //void push_sp(std::vector<SpatialPartition>& sp, const std::vector<bool>& cond);
            //void push_sp(const std::vector<std::vector<BinRMTag>>& dist, const std::vector<Storage<SpatialPartition>>& storage, std::vector<bool>& cond);
            void add_sp(SpatialPartition& other_sp);
            //void add_sp(const std::deque<SpatialPartition>& sp);
            //void add_sp(const std::vector<SpatialPartition>& sp, std::vector<bool>& cond);
            //void add_sp(const std::vector<std::vector<BinRMTag>>& dist, const std::vector<Storage<SpatialPartition>>& storage);
            void add_sp(const std::vector<std::vector<BinRMTag>>& dist, const ArrCon<SpatialPartition>& storage);
            void remap_uniproc(const std::deque<Mesh>& mesh);

            void merge_sp(const std::deque<SpatialPartition>& sp);
            void remove_dup_cells_and_points();
            double dev() const;
            void merge_sp_meshes(std::deque<Mesh>& mesh) const;
            //void print_all_meshes_in_partitions();
            //void print_all_meshes_in_partitions_uniproc();
            bool load();
            void set_global_rm(const RegularMesh& rm);
            std::map<Tag, int> search_rm(const AABB& aabb, const Vector3& cnt, std::pair<Tag, int>& resbin, bool verbose) const;

            //void reduce_sp();
            void reset_mesh_connectivity();

            template<class Archive> void serialize(Archive & ar, const unsigned int version)
            {
                ar & tl_;
                ar & fl_;
                ar & verbose_;
                ar & global_rm_;
                ar & global_adt_;
                ar & dev_;
                ar & sp_;
                ar & bintag_proc_map_;
                ar & nmesh_;

                assert(!isnan(tl_));
                assert(!isnan(fl_));
                assert(!isnan(dev_));
                assert(!isnan(nmesh_));
            }

        private:
            double tl_;
            double fl_;
            ADT global_adt_;
            bool verbose_;
            RegularMesh global_rm_; // global rm / load map rm / bins are empty of cells / only used as a map.
            //double load_;
            double dev_;
            Profiler* profiler_;
            boost::mpi::communicator* comm_;
            std::vector<SpatialPartition> sp_;
            //std::map<BinRMTag, int> bintag_proc_map_;
            std::map<int, int> bintag_proc_map_;
            int nmesh_;
            //std::vector<SpCellMover> sp_cell_mover_;

            //bool master() const;
            void profstart(std::string fun);
            void profstop(std::string fun);
            //bool brmt_rank(int bintag, BinRMTag& brmt, int& rank) const;
    };

    bool duplicate_exist(const SpatialPartitionContainer& spc);
    std::vector<AABB> get_hole_aabb(boost::mpi::communicator& comm, const SpatialPartitionContainer& spc);
    double imbalance(const SpatialPartitionContainer& spc);
    double field_imbalance(const SpatialPartitionContainer& spc);
    double cell_imbalance(const SpatialPartitionContainer& spc, double& max, double& min);
}

#endif
