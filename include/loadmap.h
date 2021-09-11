#ifndef LOADMAP_H
#define	LOADMAP_H

#include <boost/serialization/map.hpp>
#include "regular_mesh.h"
#include "profiler.h"
//#include "spsender.h"
#include "load_calculator.h"
//#include "graph.h"

namespace Tailor
{

    class Loadmap
    {
        int forced_n_refine_; // refines exactly forced_n_refine times.
        int min_forced_n_refine_; // refines minimum min_forced_n_refine times.
        std::string name_;
        RegType reg_type_;
        bool make_load_balance_;
        bool pseudo3D_;
        bool adaptive_;
        //std::shared_ptr<Settings> settings_;
        std::unique_ptr<Graph> graph_;
        //double threshold_;
        bool verbose_;
        bool print_dev_;
        bool print_graph_dist_;
        std::vector<Mesh> arrival_meshes;
        double refine_tol_;
        double final_dev_;
        bool uniproc_;
        Profiler* profiler_;
        //std::vector<DonorInfo> donor_info_;
        bool printmesh_;
        bool dorefine;
        bool printlm_;
        boost::mpi::communicator* comm_;
        //unsigned int nworker_;
        int refine_limit_;
        int nrefine;
        //std::unique_ptr<RegularMesh> rm_;
        RegularMesh* rm_;
        const std::deque<Mesh>* mesh_;
        std::vector<std::vector<BinRMTag>> aug_bin_tag_;
        //std::map<BinRMTag, int> bintag_proc_map_;
        std::map<int, int> bintag_proc_map_;
        //boost::bimap<int, int> mesh_tag_index_map;
        LoadEstim load_estim_type_;

        Tag generate_meshtag() const;
        //void add_mesh(Mesh& m);
        //void add_mesh(const std::deque<Mesh>& mesh);
        //bool resolution_is_enough(BinRMTag& heaviest_bt, int nmake_map, int iteration);
        bool resolution_is_enough(BinRMTag& heaviest_bt, int nmake_map);
        void read_settings();
        //void deviation_in_load(double& dev, BinRMTag& heaviest_bt, bool& collapsed, bool verbose=false);
        void print_lm() const;
        void print_rm(const RegularMesh& rm, const LoadCalculator& lc) const;
        void print_orphan(int iteration) const;
        void print_mesh(int iteration) const;
        void refine(int nmake_map);
        void make_regular_map(RegType, bool adaptive);
        //void make_regular_map_uniproc();
        //void make_regular_map_uniproc_resident();
        void gather_regular_maps();
        AABB max_aabb() const;
        //void remap();
        //void send_sp(const std::vector<MeshTransferInfo>& mti, std::vector<boost::mpi::request>& send_req);
        //void send_sp(const std::vector<MeshTransferInfo>& mti, std::vector<RecvRequest>& recv_req_tup, std::vector<Mesh>& arrival_meshes);
        //void complete_recv(const std::vector<MeshTransferInfo>& mti, std::vector<RecvRequest>& recv_req_tup, std::vector<Mesh>& arrival_meshes);
        //void try_to_recv(std::vector<RecvRequest>& recv_req_tup, std::vector<Mesh>& arrival_meshes);
        //void recv_sp(const std::vector<MeshTransferInfo>& mti, std::vector<Mesh>& arrival_meshes);
        //void add_mesh(std::vector<Mesh>& mesh);
        //void add_arrival_meshes(std::vector<Mesh>& arrival_meshes);
        void reset_rm(int ns);
        //void recv_arrival_meshes(const SpatialPartition& sp);
        void profstart(std::string fun);
        void profstop(std::string fun);
        void calc_load(BinRMTag& heaviest_bt, std::vector<BinRMTag>& sorted_tags, bool verbose);
        BinRMTag get_heaviest_bin(const std::vector<BinRMTag>& sorted_tags, bool verbose);
        void sort_load(std::vector<double>& sorted_load, std::vector<double>& non_sorted_load, std::vector<BinRMTag>& sorted_tags, bool verbose);
        double estimated_deviation(const std::vector<double>& aug_sorted_load);

        public:

        ~Loadmap();
        Loadmap(boost::mpi::communicator* comm, Profiler* profiler, bool verbose, LoadEstim load_estim_type, bool pseudo3D, bool make_load_balance, RegType reg_type, const std::deque<Mesh>*, std::string name);

        void set_comm(boost::mpi::communicator* comm);
        void set_profiler(Profiler* prof);
        void set_refine_tol(double d);
        //const std::map<BinRMTag, int> bintag_proc_map() const;
        const std::map<int, int> bintag_proc_map() const;
        //const size_t nmesh() const;
        const std::vector<std::vector<BinRMTag>>& aug_bin_tag() const;
        //unsigned int nworker() const;
        const RegularMesh& rm() const;
        //const Mesh& mesh(const Tag& tag) const;
        //double final_dev() const;
        //void make_map(const std::vector<std::string>& file_name);
        const std::deque<Mesh>* mesh() const;
        //void merge_donor_info(std::vector<boost::mpi::request>& request);
        //void move_mesh(const std::vector<Vector3>& v);
        //void rotate_mesh(const std::vector<double>& rotation, int axis, const std::vector<Vector3>& rot_point);
        //void read_mesh_all_procs(const std::vector<std::string>& file_name);
        void make(int nmake_map, bool adaptive);
        void get_bin_to_proc();
        void clear();
        //void distribute_mti(std::vector<MeshTransferInfo>& mti) const;
        void print(int iteration);
        //size_t nbin() const;
        double refine_tol() const;
        LoadEstim load_estim_type() const;
        //void distribute_mti_2(std::vector<MeshTransferInfo>& mti) const;
        void remove_cells_from_rm();
        //void add_mesh_from_sp(const std::deque<SpatialPartition>& incoming);
        //void add_arrival_meshes(const std::deque<SpatialPartition>& incoming);
    };

    bool greedy_partition(int nhead, const std::vector<double>& input, std::vector<double>& output, const std::vector<BinRMTag>& input_bin_tag, std::vector<std::vector<BinRMTag>>& output_bin_tag, bool verbose_);
    //void fill_recv_req_tup(int rank, const std::vector<MeshTransferInfo>& mti, std::vector<RecvRequest>& recv_req_tup);
    //void print_rm(const RegularMesh& rm, std::string s, const std::map<BinRMTag, int>& bintag_proc_map_);
    void print_rm(const RegularMesh& rm, std::string s, const std::map<int, int>& bintag_proc_map_);
}

#endif
