#ifndef MESH_H
#define	MESH_H

#include <vector>
#include <numeric>
#include <array>
#include <string>
#include <algorithm>
#include <boost/bimap.hpp>
#include "tag.h"
#include "geom.h"
//#include "meshcell_donorinfo.h"
//#include "spc_donor_info.h"
#include <fstream>
#include <typeinfo>        
#include "nei.h"
#include "var.h"
#include "utility.h"
#include "boost/program_options.hpp"
#include "cantor.h"
#include "freestream.h"
#include "profiler.h"
#include <boost/serialization/list.hpp>
#include "exchanger.h"
#include "adt.h"

namespace Tailor
{
    typedef std::vector<MeshPoint> mpc;
    //typedef std::list<MeshFace> mfc;
    typedef std::vector<int> intvector;
    typedef boost::bimap<int, int> bimap_int;


    enum class StencilWalkResult
    {
        inside_cell = 0,
        inside_hole = 1,
        outside_mesh = 2,
    };

    class Mesh
    {    
        friend class DonorSearcher;
        friend class HoleProfiler;
        friend class BoundaryCondition;
        friend class Solver;
        friend class Assembler;
        friend class Limiter;
        friend class Gradient;
        friend class DonorRequestExchanger;
        friend class TCRExchanger;

        mcc wall_boundaries_;
        mcc symmetry_boundaries_;
        mcc dirichlet_boundaries_;
        mcc empty_boundaries_;
        mcc farfield_boundaries_;
        mcc interog_boundaries_;
        mcc cell_;
        mpc point_;
        //mfc face_;
        Tag tag_;
        //Tag parent_mesh_;
        //bimap_int wall_tag_index_map_;
        //bimap_int dirichlet_tag_index_map_;
        //bimap_int farfield_tag_index_map_;
        //bimap_int empty_tag_index_map_;
        //bimap_int interog_tag_index_map_;
        //AABB hole_aabb_;
        int priority_;
        friend class boost::serialization::access;

        //MeshPoint& point_p(const Tag& ptag);
        MeshPoint* point_p(const Tag& ptag);
        unsigned int point_index(const Tag& t) const;
        void remove_cell_boundary(const Tag& ic);
        void flood_fill(const Tag& start);

        void set_cell(const Tag& t, const MeshCell& c);
        void remove_cell_from_cellhood(const Tag& ic);
        void deparent_cell_from_vertices(const Tag& ic);
        void merge_wall_to_interior(const Mesh& wall_mesh);
        void merge_symmetry_to_interior(const Mesh& symmetry_mesh);
        void merge_dirichlet_to_interior(const Mesh& dirichlet_mesh);
        void merge_farfield_to_interior(const Mesh& farfield_mesh);
        void merge_empty_to_interior(const Mesh& empty_mesh);
        void increase_overlap_thickness_(MeshCell& mc, int& count, int nlayer, const ADT& passive_cell_adt, Mesh& passive_mesh);

        //std::set<Point> bou_raw_point_(const std::vector<MeshCell>* container) const;

        public:

        Mesh();
        //Mesh(const Tag& tag, const Tag& parent_mesh);
        Mesh(const Tag& tag);
        //Mesh(const Tag& tag): Mesh(tag, tag) {};

        void increase_overlap_thickness(int nlayer, const ADT& passive_cell_adt, Mesh& passive_mesh);
        Vector5 uniform_prim(const FlowInit& fs);
        //void update_cell_pnei_addresses();
        //void update_face_parent_addresses();
        //void update_cell_vertex_addresses();
        void reset_to_mid();
        void convert_receptor_to_hole();
        mcc::iterator query_itp(const Tag& ic);
        int priority() const;
        void set_priority(int pri);
        void connect_partition_cells(ArrCon<Nei>& arrival_cell, int rank, std::function<bool(const Vector3&, int celltag)> is_resi, Profiler*);
        void update_ghost_primitives(const ArrCon<Var>& arrival, int rank, double gamma);
        void oga_interpolate(const ArrCon<Var>& arrival, int rank);
        void oga_interpolate(const std::deque<Mesh>& mesh, int rank);
        void equate_RD_to_RDmid();

        void init_flow();
        //void init(const Vector3& vinf_air, const Freestream& fs);
        void init_uniform(const FlowInit&, double gamma);
        void init_gaussian(const FlowInit&, double gamma);
        void init_xsplit(const FlowInit& finit, double gamma);

        void connect_add_bou_to_interior(BouType boutype, int rank);
        void reset_R_checked();
        MeshCell* boundary(BouType btype, const Tag& t);
        const MeshCell* boundary(BouType btype, const Tag& t) const;
        void change_bou(BouType btype, MeshCell& mc, const FaceTag& ftag);

        void set_all_cells_as_interior();
        bool is_gcl_satisfied(int rank) const;
        void calc_mesh_velocities(const Freestream& fs, int rank, double real_time);
        void reset_R();
        void reset_D();
        void set_prim_cons(BouType btype, const Vector5& prim, double gamma);

        void connect_after_exchange(std::function<bool(const Vector3&)> is_resi, int rank, Profiler*, std::string);
        void update_interior_face(MeshFace& mf, MeshFace& cf, MeshCell& mc, MeshCell& cc, MeshFace& gmf);
        //void updateface(MeshFace& mf, MeshFace& cf);
        //void addface(MeshFace& mf, MeshFace& cf, int& facecounter);
        void add_interior_face(MeshFace& mf, MeshFace& cf, MeshCell& mctag, MeshCell& nctag);
        void add_partition_face(MeshFace& mf);
        //void addface_bou(MeshCell& mc, MeshFace& mf, int& facecounter, int rank);
        void addface_bou(MeshFace& mf);
        MeshFace* common_face(MeshCell& mc, const Tag& celltag);
        MeshFace* common_face(MeshCell& mc, const FaceTag& ft);
        const MeshFace* common_face(const MeshCell& mc, const FaceTag& ft) const;
        bool is_potential_parent_cell(const MeshFace& mf, const Tag& celltag);
        void disconnect_orphan_faces();
        void remove_ghosts();
        void remove_nonresidents();

        //void determine_unreached_hole_cells();
        //void determine_unreached_hole_cells_via_aabb(const AABB& cutter_aabb);
        void determine_unreached_hole_cells(const AABB& cutter_hole_aabb);
        void print_donor_info() const;
        void update_prev_cand_donor();

        //const mfc& face() const;

        //size_t mem() const;

        void bbox(Vector3& min_, Vector3& max_) const;

        void reset_erase_marks();

        void add_wall_boundary(const MeshCell& mc);
        void add_wall_boundary(MeshCell&& mc);
        void add_symmetry_boundary(const MeshCell& mc);
        void add_symmetry_boundary(MeshCell&& mc);
        void add_dirichlet_boundary(const MeshCell& mc);
        void add_dirichlet_boundary(MeshCell&& mc);
        void add_farfield_boundary(const MeshCell& mc);
        void add_farfield_boundary(MeshCell&& mc);
        void add_empty_boundary(const MeshCell& mc);
        void add_empty_boundary(MeshCell&& mc);
        void add_interog_boundary(const MeshCell& mc);
        void add_interog_boundary(MeshCell&& mc);

        void update_points_from_cell_vertices(int rank);
        void merge_batch(const Mesh& other_mesh, int rank);
        void remove_merge_duplicate_cells();
        void remove_merge_duplicate_points();
        //void remove_duplicate_cells();
        void remove_duplicate_cells(Profiler* profiler, std::string s);
        void remove_duplicate_points();

        // Getters
        const Tag& tag() const;
        //const Tag& parent_mesh() const;
        //const AABB& hole_aabb() const;
        const bimap_int& dirichlet_tag_index_map() const;
        const mcc& wall_boundaries() const;
        const mcc& symmetry_boundaries() const;
        const mcc& dirichlet_boundaries() const;
        const mcc& farfield_boundaries() const;
        const mcc& empty_boundaries() const;
        const mcc& interog_boundaries() const;
        size_t npartition() const;
        MeshCell* wall_boundary_p(const Tag& t);
        MeshCell* symmetry_boundary_p(const Tag& t);
        MeshCell* dirichlet_boundary_p(const Tag& t);
        MeshCell* farfield_boundary_p(const Tag& t);
        MeshCell* empty_boundary_p(const Tag& t);
        MeshCell* interog_boundary_p(const Tag& t);

        void erase_marked_cells(int rank);
        void mark_to_be_erased(const Tag& ic);
        void prepare_to_remove_cell(MeshCell& mc, int rank);
        void prepare_to_remove_ghost(const Tag& ic);

        void add_interiors_nonsorted_nopoint(mcc& cells);
        void add_points(const mcc& cells, int rank);
        void add_bous(mcc& cells, BouType type);

        size_t calc_nhole_cell() const;
        //const MeshPoint& point(const Tag& ptag) const;
        const MeshPoint* point(const Tag& ptag) const;
        const MeshPoint& point(int) const = delete;
        const MeshCell& cell(const Tag& ctag) const;
        MeshCell& cell_p(const Tag& ctag);
        const MeshCell& cell(int) const = delete;
        //const MeshFace& face(const FaceTag& ctag) const;
        //const MeshFace& face(int) const = delete;
        //MeshFace& face_p(const FaceTag& ctag);
        //MeshFace& face_p(int) = delete;
        const mpc& point() const;
        std::vector<Point> rawpoint() const;
        //std::set<Point> bou_raw_point(BouType btype) const;
        const mcc& cell() const;
        mcc& cell_p();
        const MeshCell* wall_boundary(const Tag& t) const;
        const MeshCell* symmetry_boundary(const Tag& t) const;
        const MeshCell* dirichlet_boundary(const Tag& t) const;
        const MeshCell* farfield_boundary(const Tag& t) const;
        const MeshCell* interog_boundary(const Tag& t) const;
        const MeshCell* empty_boundary(const Tag& t) const;
        bool do_point_exist(const Tag& t) const;

        void move(const Vector3& v);
        void rotate(double angle, const Vector3& axis, const Vector3& rot_point);
        void print_as_vtk (std::string file_name) const;
        void print_as_vtk_geometry(std::string file_name) const;
        void print_wall_as_vtk(std::string file_name) const;
        void print_interog_as_vtk(std::string file_name) const;
        void set_tag(const Tag& t);
        //void set_parent_mesh(const Tag& t);
        void merge(const Mesh& other_mesh);
        void merge_no_check(const Mesh& other_mesh);
        
        // cell and point addition
        void add_interior_cell(const MeshCell& c);
        void add_interior_cell_no_check(const MeshCell& c);
        void add_interior_cell_find(const MeshCell& c);
        void add_interior_cell_nonsorted(const MeshCell& c);
        mcc::iterator add_interior_cell_sorted(const MeshCell& c, bool& exist);
        mcc::iterator add_cell_sorted(const MeshCell& c, bool& exist);
        mcc::iterator add_cell_nonsorted(const MeshCell& c);
        void add_element_nonsorted(const MeshCell& mc);
        void add_element(const MeshCell& mc);
        void add_element(MeshCell&& mc);
        void add_element_no_check(const MeshCell& mc);
        //void insert_cells(const std::vector<MeshCell>& cells);
        void add_point(MeshPoint p, const Tag& t, size_t size=0);
        void add_cell_only(MeshCell c, size_t size=0);
        void add_interior_nonsorted_addpoint(MeshCell mc);

        void opposing_cell_pc1(const MeshCell& mc, const MeshFace& mf, MeshCell*& op_cell, MeshFace*& shared_face);
        void opposing_cell_pc2(const MeshCell& mc, MeshFace& mf, MeshCell*& op_cell, MeshFace*& shared_face);

        // cell and point removal
        void remove_point(Tag ip);
        void remove_cell(const Tag& ic);
        void remove_isolated_points();
        //void remove_merge_dup_pts();
        //void remove_dup_cells();

        // donor search
        void convert_candholes_to_holes();
        void convert_undefined_to_field();
        void convert_undefined_to_nonresident(std::function<bool(const Vector3& cnt)> resi_fn);
        void set_as_hole(const std::set<Tag>& overlaps);
        void set_as_hole(const Tag& t);
        //void set_as_mreceptor(const Tag& t);
        void set_as_field(const Tag& t);
        //void set_status_of_cells(const SPCDonorInfo& spcdi);
        void reset_all_cell_oga_status();
        void set_cell_type(const Tag& my_cell_tag, const MeshCell* other_cell, const Mesh& other_mesh);

        // hole profiling
        void construct_global_aabb(size_t size, size_t nmesh, const std::vector<AABB>& aabbs, int i);
        //void determine_hole_aabb();
        bool flood_fill(const Tag& start, const AABB& aabb, int dummyrank, int sptag, std::ofstream& out);
        //void set_hole_aabb(double minx, double miny, double minz, double maxx, double maxy, double maxz);
        //void set_hole_aabb(const AABB& aabb);

        // connectivity
        void connect_add_bou_to_interior(Mesh& boumesh, BouType boutype, int rank);
        //void connect_cells();
        void connect_cells(std::function<bool(const Vector3&)> is_resi, int rank);
        void connect_cells(std::function<bool(const Vector3&)> is_resi, int rank, Profiler* profiler, std::string name);
        //void connect_add_outer_to_interior(const Mesh& wm, std::string dummyfilename);
        //void connect_add_wall_to_interior(const Mesh& wm);
        void destroy_cell_hood();
        
        // reservations
        void reserve_interior(size_t size);
        void reserve_wall(size_t size);
        void reserve_dirichlet(size_t size);
        void reserve_farfield(size_t size);
        void reserve_empty(size_t size);
        void reserve_interior_only(size_t size);

        //void update_connectivity(int rank);

        // queries
        const MeshCell* query_sorted(const Tag& ic) const;
        const MeshPoint* query_point(const Tag& ic) const;
        const MeshCell* query(const Tag& ic) const;
        const MeshCell* query_bou(const Tag& ic, BouType type) const;

        void check_nei(const MeshCell& mc, const Polyhedron& pol, std::set<Tag>& overlaps, int cuttercelltag, bool verbose) const;
        //void batch_merge(Mesh& m);
        //void set_point_tag_index_map();
        //void set_cell_tag_index_map();
        void shrink_points();
        void shrink_cells();
        void remove_dup_cells_and_points();
        //void make_cell_tags();
        //void merge_using_tags(const Mesh& other_mesh);
        void sort_cells();
        void sort_points();

        void reset_mesh_connectivity();

        void simple_merge(const Mesh& other_mesh, int rank);
        void remove_parent_cells_of_vertices_of_all_cells();
        void set_parent_cell_of_vertices_of_all_cells();
        void remove_parent_cells_of_all_points();
        void add_meshpoint(const MeshCell& c);
        void add_meshpoint(MeshCell& mc, int rank);

        bool operator<(const Mesh& other) const;
        const MeshCell* cell_ptr(const Tag& ctag) const;

        void determine_non_residents(const AABB& aabb);

        template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            ar & tag_;
            //ar & parent_mesh_;
            ar & point_;
            ar & cell_;
            //ar & wall_tag_index_map_;
            //ar & dirichlet_tag_index_map_;
            //ar & farfield_tag_index_map_;
            //ar & empty_tag_index_map_;
            //ar & interog_tag_index_map_;
            ar & wall_boundaries_;
            ar & symmetry_boundaries_;
            ar & dirichlet_boundaries_;
            ar & farfield_boundaries_;
            ar & interog_boundaries_;
            ar & empty_boundaries_;
            //ar & face_;
        }
    };

    std::ifstream& go_to_beg_of_line(std::ifstream& file, int num);
    std::tuple<MeshCell*, MeshCell*> left_and_right_cells(Mesh& mesh, MeshFace& mf, const Tag& mctag);
    std::tuple<const MeshCell*, const MeshCell*> left_and_right_cells(const Mesh& mesh, const MeshFace& mf, const Tag& mctag);
    const MeshCell* opposing_nei(const Mesh& mesh, const MeshFace& mf, const Tag& me);
    //void convert_mesh_to_fortran(const Mesh& mesh);
    bool are_common_faces(const MeshFace& a, const MeshFace& b);
}

#endif
