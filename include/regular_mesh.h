#ifndef REGULAR_MESH_H
#define REGULAR_MESH_H

#include <map>
#include "core.h"
#include "profiler.h"
#include "rm_index.h"
#include "adt.h"
#include "graph.h"
#include "rm_bin.h"
#include <fstream>
#include <numeric>
#include <boost/serialization/map.hpp>

#define RM_NNEI 9

namespace Tailor
{
    enum class RelativePosition
    {
        undefined = -1,
        swa = 0,
        sca = 1,
        sea = 2,
        cwa = 3,
        cca = 8,
        cea = 4,
        nwa = 5,
        nca = 6,
        nea = 7,

        swf = 9,
        scf = 10,
        sef = 11,
        cwf = 12,
        ccf = 13,
        cef = 14,
        nwf = 15,
        ncf = 16,
        nef = 17,

        swc = 18,
        scc = 19,
        sec = 20,
        cwc = 21,
        ccc = 22,
        cec = 23,
        nwc = 24,
        ncc = 25,
        nec = 26,
    };

    bool in_range(const RegularMeshIndex& index, const RegularMeshIndex& global_min, const RegularMeshIndex& global_max);
    unsigned int mat_to_vec_index(const RegularMeshIndex& index, const vec3<unsigned int>& nstripe);
    RelativePosition get_relative_position(int index, int targetindex, const vec3<unsigned int>& nstripe, size_t maxbinsize, bool& exact);
    int get_binindex(unsigned int index, RelativePosition rp, const vec3<unsigned int>& nstripe, unsigned int binsize);
    vec3<double> llcoor(const RegularMeshIndex& index, const vec3<double>& aabb_min, const vec3<double>& h);

    class RegularMesh
    {
        Tag tag_;
        vec3<int> nstripe_;
        vec3<double> h_;
        AABB aabb_;
        std::vector<Bin> bin_;
        RegularMeshIndex global_min_;
        RegularMeshIndex global_max_;
        std::map<Tag, RegularMesh*> rmtag_address_map_;
        std::map<Tag, Bin*> bintag_address_map_;
        friend class boost::serialization::access;
        friend class HoleMap;

        Bin& bin_p(int row, int col, int depth);
        Bin& bin_p(const RegularMeshIndex& ind);
        Bin& bin_p(const Tag& t);
        //Bin& bin_p(const BinRMTag& bt);
        bool is_resident(const vec3<double>& cnt) const;
        int get_new_bt_p() const;
        int get_new_bt() const;
        int get_new_rmtag() const;
        //void gather_load(std::vector<double>& load_local, int& j, const std::deque<Mesh>& mesh, LoadEstimType load_estim_type);
        void size(size_t& s) const;
        void max_bin_tag(int& s) const;
        BinRMTag index_to_bin(size_t& s, size_t i) const;


        void print_(std::string file_name) const;

        bool get_index_regular(const vec3<double>& cv, std::vector<RegularMeshIndex>& index, bool closest, int rank, int celltag) const;
        RegularMeshIndex get_index_regular_unique(const vec3<double>& cv) const;
        void get_bintag_adaptive_unique_(const vec3<double>& cv, std::vector<BinRMTag>& tag, int rank, int celltag) const;

        public:

        const std::map<Tag, Bin*>& bintag_address_map() const;

        void insert_bin_addresses(RegularMesh*);

        void get_adjacency_from_graph(const Graph& graph);

        RegularMesh();
        RegularMesh(const Tag& rmtag);
        RegularMesh(const RegularMesh& other);
        RegularMesh& operator=(const RegularMesh& other);

        void flood(const Tag& start, const AABB& aabb);

        void move(const vec3<double>& v);
        void rotate(double ang, int axis, const vec3<double>& rot_point);

        // Getters
        const Tag& tag() const;
        size_t size() const;
        int max_bin_tag() const;
        const std::map<Tag, RegularMesh*>& rmtag_address_map() const;
        const Bin& bin(int row, int col, int depth) const;
        const Bin& bin(const RegularMeshIndex& ind) const;
        const Bin& bin(int i) const;
        const Bin& bin(const Tag& t) const;
        //const Bin& bin(const BinRMTag& bt) const;
        RegularMeshIndex global_index(const RegularMeshIndex& i) const;
        const RegularMeshIndex& global_min() const;
        const RegularMeshIndex& global_max() const;
        const vec3<int>& nstripe() const;
        int nstripe(int i) const;
        const std::vector<Bin>& bin() const;
        const vec3<double>& h() const;
        double h(int i) const;
        const AABB& aabb() const;

        // refine
        void refine(const std::vector<Mesh>& meshes, int rank, bool pseudo3D);
        //void refine_adaptive(const std::deque<Mesh>& meshes, const BinRMTag& heaviest_bt, int rank);
        const RegularMesh* refine_adaptive(const std::deque<Mesh>& meshes, const BinRMTag& heaviest_bt, int rank, int& new_rmtag, int& new_bt, bool pseudo3D, RegType regtype);

        // register
        void register_cell(const MeshCell& cell, const Tag& cell_tag, const Tag& mesh_tag, int rank, bool adaptive);
        void register_resident_cell(const MeshCell& cell, const Tag& cell_tag, const Tag& mesh_tag, int rank);
        void register_mesh(const Mesh& mesh, int rank, bool adaptive);
        void register_overlapping_mesh(const Mesh& mesh, int rank, bool adaptive);
        void register_resident_mesh(const Mesh& mesh, int rank);
        //void register_bincells(const std::vector<BinCell>& bincell, const std::deque<Mesh>& meshes, bool pseudo3D, int rank);
        void register_bincells(const std::deque<BinCell>& bincell, const std::deque<Mesh>& meshes, bool pseudo3D, int rank, RegType);

        BinRMTag index_to_bin(size_t i) const;
        //void gather_load(const boost::mpi::communicator& comm, std::vector<double>& load_global, const std::deque<Mesh>& mesh, LoadEstimType load_estim_type);
        //void gather_load(std::vector<double>& load_global, const std::deque<Mesh>& mesh, LoadEstimType load_estim_type);
        void calc_aabb_max();
        bool extend_aabb(const AABB& other_aabb);
        void insert_to_rmtag_address_map(const Tag& t, RegularMesh* rmp);
        void insert_to_bintag_address_map(const Tag& t, Bin* bin);
        void set_rmtag_address_map(const std::map<Tag, RegularMesh*>& map);
        void update_address(std::map<Tag, RegularMesh*>& map, std::map<Tag, Bin*>& binmap);
        void update_address();
        void set_tag(const Tag& t);
        //std::vector<BinRMTag> sort_bin_tags_based_on_load(const boost::mpi::communicator& comm, std::vector<double>& load_global, const std::deque<Mesh>& mesh, LoadEstimType load_estim_type);
        void get_bin_tags(std::vector<BinRMTag>& tag) const;
        void calc_h();
        void clear_bincells();
        //int get_binindex(unsigned int index, RelativePosition rp);
        std::vector<Tag> pneis_of_bin(unsigned int index);
        //void set_bincell(const RegularMeshIndex& ind, const std::vector<BinCell>& bc);
        void set_bincell(const RegularMeshIndex& ind, const std::deque<BinCell>& bc);
        void set_global_min(const RegularMeshIndex& i);
        void set_global_max(const RegularMeshIndex& i);
        void set_global_min(int i, int j, int k);
        void set_global_max(int i, int j, int k);
        void set_nstripe(int x, int y, int z);
        void set_nstripe(const vec3<int>& ns);
        void set_h(const vec3<double>& v);
        void set_h(double v0, double v1, double v2);
        void set_aabb(const AABB& aabb);
        //void set_aabb_min(const vec3<double>& m);
        //void set_aabb_max(const vec3<double>& m);
        vec3<double> centroid(const RegularMeshIndex&) const;
        vec3<double> llcoor(const RegularMeshIndex& index) const;
        void get_bintag_adaptive_unique(const vec3<double>& cv, BinRMTag& tag, int rank, int celltag) const;
        //bool get_bintag_adaptive(const vec3<double>& cv, std::vector<BinRMTag>& tag, bool closest) const;
        void get_bintag_adaptive(const AABB& aabb, std::vector<BinRMTag>& tag, int& nlevel, int celltag) const;
        void get_bintag_adaptive_(const AABB& aabb, std::vector<BinRMTag>& tag, int& nlevel, int celltag) const;
        void print(std::string aux) const;
        void merge(const RegularMesh& other);
        bool point_inside_bin(const RegularMeshIndex& ind, const vec3<double>& p);
        void insert_bins(int mintag, int mesh_load_size, int rank);
        void set_props(const Mesh& mesh, const AABB& forced_min_aabb, int _nstripe=0);
        void info() const;
        bool validate(const RegularMeshIndex& global_ind, RegularMeshIndex& local_ind) const;
        void calc_step_length();
        void clear_cells();
        bool operator<(const RegularMesh& other) const;
        bool operator==(const RegularMesh& r) const;
        template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            assert(!std::isnan(tag_()));

            ar & nstripe_;
            ar & tag_;
            ar & aabb_;
            ar & h_;
            ar & bin_;
            for (auto& i: rmtag_address_map_)
            {
                i.second = nullptr;
            }
            ar & rmtag_address_map_;
            for (auto& i: bintag_address_map_)
            {
                i.second = nullptr;
            }
            ar & bintag_address_map_;
            if (tag_() == 0) {
                update_address();
            }
        }
    };

    //bool get_index_regular(const vec3<double>& cv, std::vector<RegularMeshIndex>& index, const vec3<double>& aabb_min, const vec3<double>& h, const vec3<int>& nstripe, bool closest);
    Tag closest_non_empty_bin(const RegularMeshIndex& ind, const RegularMesh& rm, const Mesh& m);
    void make_adt(const RegularMesh& rm, ADT& adt, int rank);
}

#endif
