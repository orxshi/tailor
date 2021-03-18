#ifndef MESHCELL_H
#define MESHCELL_H

//#include "meshpoint.h"
#include "meshface.h"
#include <boost/bimap.hpp>
#include <boost/serialization/utility.hpp>
#include "array.h"
#include "phys.h"
#include "donor.h"

#define NPOINT 8

namespace Tailor
{
    using FinalCandDonor = Array<Tag, 5>;

    class MeshCell;

    //typedef std::vector<MeshCell> mcc;
    typedef std::deque<MeshCell> mcc;
    //typedef std::list<MeshCell> mcc;
    //
    typedef Array<Donor, 10> Donorcon;

    class MeshCell
    {
        using Boundaries = Array<Tag, 6>;
        using MCPoint = Array<MeshPoint, NPOINT>;
        using Pnei = std::vector<Tag>;
        using Face = Array<MeshFace, 6>;

        double sumarea_;
        BouType btype_;
        bool erase_;
        //bool is_wall_;
        //bool is_dirichlet_;
        //bool is_interior_;
        //bool is_partition_;
        //std::vector<Tag> wall_boundary_; // dynamic
        //std::vector<Tag> dirichlet_boundary_; // dynamic
        Boundaries wall_boundary_;
        Boundaries dirichlet_boundary_;
        Boundaries farfield_boundary_;
        Boundaries empty_boundary_;
        Boundaries interog_boundary_;
        Tag interior_boundary_; // only for boundaries.
        Tag first_tag_;
        Tag root_parent_mesh_;
        Tag tag_;
        Tag parent_mesh_;
        Donor donor_;
        //std::pair<Tag, Tag> donor_; // I will also use this for starting seed for field cells.
        //const MeshCell* donor_addr_;
        //std::vector<Tag> pnei_; // principal neighbors. // dynamic
        Pnei pnei_; // principal neighbors. // dynamic
        MCPoint point_;
        //std::vector<MeshFace> face_; //dynamic
        Face face_; //dynamic
        //std::vector<Tag> iface_;
        Polyhedron poly_;
        OGA_cell_type_t OGA_cell_type_;
        //Array<Donor, 5> cand_donor_;
        ///std::vector<Donor> cand_donor_;
        Donorcon cand_donor_;
        //FinalCandDonor cand_donor_mesh_;
        //FinalCandDonor cand_donor_cell_;
        //Array<const MeshCell*, 5> cand_donor_addr_;
        //Array<Donor, 5> prev_cand_donor_;
        //std::vector<Donor> prev_cand_donor_;
        Donorcon prev_cand_donor_;
        //FinalCandDonor prev_cand_donor_mesh_;
        //FinalCandDonor prev_cand_donor_cell_;
        Vector3 vgn_;

        Vector5 prim_;
        Vector5 cons_sp1_;
        Vector5 cons_s_;
        Vector5 cons_n_;
        Vector5 cons_nm1_;
        Vector5 dQ_; 
        Vector5 old_dQ_;
        //double max_eigen_;
        Vector5 R_;
        Matrix5 D_;
        Vector5 R_mid_;
        Matrix5 D_mid_;
        std::vector<double> ls_wx_;
        std::vector<double> ls_wy_;
        std::vector<double> ls_wz_;
        double dtao_;

        void add_cand_donor(const Tag& donor_mesh, const Tag& donor_cell, const MeshCell* addr);
        friend class boost::serialization::access;
        friend class Mesh;
        friend class BoundaryCondition;
        friend class Solver;
        friend class Gradient;
        friend class DonorSearcher;

        public:

        Vector3 vgn() const;
        void mesh_velocity(double dt, const Freestream& fs, const Component& compo);
        void init(const Vector3& vinf_air, const Freestream& fs, const Component& compo);
        std::vector<Tag> boundaries() const;
        //const MeshCell* donor_addr() const;
        void reset_face_tags();
        double prim(int i) const;
        double cons_sp1(int i) const;
        double cons_s(int i) const;
        double cons_n(int i) const;
        double cons_nm1(int i) const;
        const Vector5& prim() const;
        const Vector5& cons_sp1() const;
        const Vector5& cons_s() const;
        const Vector5& cons_n() const;
        const Vector5& cons_nm1() const;
        const Vector5& R() const;
        const Vector5& dQ() const;
        const Vector5& old_dQ() const;
        const Matrix5& D() const;
        double D(int i, int j) const;
        const Matrix5& M() const;
        double R(int i) const;
        void set_prim(const Vector5& other);

        void deparent_neis_from_faces();
        void update_prev_cand_donor();
        //void reset_cand_donor(const Tag& meshtag);
        void add_boundary(const Tag& t, BouType boutype);

        //MeshCell(const Tag& tag, const Tag& parent_mesh): tag_(tag), parent_mesh_(parent_mesh), OGA_cell_type_(OGA_cell_type_t::undefined), is_ghost_(false), boundary_type_(boundary_t::undefined) {}
        MeshCell(const Tag& tag, const Tag& parent_mesh): tag_(tag), parent_mesh_(parent_mesh), OGA_cell_type_(OGA_cell_type_t::undefined) {}
        //MeshCell(): MeshCell(Tag(-1), Tag(-1)) {}
        MeshCell();
        //MeshCell(const Tag& tag, const Tag& parent_mesh, const std::vector<MeshPoint>& point, boundary_t btype=boundary_t::undefined);
        MeshCell(const Tag& tag, const Tag& parent_mesh, const std::vector<MeshPoint>& point, BouType, Shape shape=Shape::undef);
        MeshCell(const std::vector<MeshPoint>& point, BouType btype=BouType::undefined): MeshCell(Tag(-1), Tag(-1), point, btype) {}

        void reset_btype();
        BouType btype() const;
        void set_btype(BouType btype);
        //size_t mem() const;
        void unmark_to_be_erased();
        void mark_to_be_erased();
        size_t npoint() const;
        bool erase() const;
        //void set_partition(bool p);
        //const bool is_partition() const;
        Face& face_p();
        const Face& face() const;
        void make_faces(Shape shape);
        //void set_cand_donor(const std::vector<Donor>& cand_donor);
        void set_cand_donor(const Donorcon& cand_donor);
        //const std::vector<Donor>& cand_donor() const;
        const Donorcon& cand_donor() const;
        //const std::vector<Donor>& prev_cand_donor() const;
        const Donorcon& prev_cand_donor() const;
        const Tag& root_parent_mesh() const;
        const Tag& first_tag() const;
        const Donor& donor() const;
        void fringe_to_field();
        void reset_oga_status();
        const Tag& parent_mesh() const;
        void remove_parent_cells_of_vertices();
        void rotate_points(double ang, const Vector3& rot_axis);
        void rotate_points(double ang, int axis, const Vector3& rot_axis);
        void move_points(const Vector3& final_loc);
        void remove_cand_donor(const Tag& donor_cell, const Tag& donor_mesh);
        void set_donor(const Tag& im, const Tag& ic, const MeshCell* donor_addr);
        void set_oga_cell_type(OGA_cell_type_t t, const Tag& donor_mesh, const Tag& donor_cell, const MeshCell* donor_addr);
        void set_oga_cell_type(OGA_cell_type_t t);
        OGA_cell_type_t oga_cell_type() const;
        void deparent_self_from_vertices();
        void deparent_self_from_faces();
        void remove_all_parent_cells_of_faces();
        void remove_all_nonself_parents_of_faces();
        //void deparent_from_faces(const Tag& t);
        void remove_self_from_neighbors(const Tag& t);
        const Tag& tag() const;
        void set_tag(const Tag& t);
        void set_parent_mesh(const Tag& celltag);
        const Polyhedron& poly() const;
        //const std::vector<MeshPoint>& point() const;
        //const std::array<MeshPoint, NPOINT>& point() const;
        const MCPoint& point() const;
        const MeshPoint& point(int i) const;
        //const std::vector<Tag>& nei() const;
        const Pnei& pnei() const;
        //const Tag& nei(int i) const;
        const Tag& pnei(int i) const;
        void add_vertex(const MeshPoint& p, const Tag& pointtag);
        //void add_pnei(const Tag& celltag);
        void add_pnei(const Tag&);
        std::vector<std::vector<Tag>> face_vertex();
        void set_parent_cell_of_vertices();
        void set_parent_cell_of_faces();
        void remove_neighbor(const Tag& t);
        //void remove_all_non_pnei();
        void remove_all_neighbors();
        std::vector<Point> geom_point() const;
        bool near_boundary() const;
        bool near_wall() const;
        bool near_interog() const;
        bool near_dirichlet() const;
        bool near_farfield() const;
        bool near_empty() const;
        void add_wall_boundary(const Tag& t);
        void add_dirichlet_boundary(const Tag& t);
        void add_farfield_boundary(const Tag& t);
        void add_empty_boundary(const Tag& t);
        const Boundaries& wall_boundary() const;
        const Boundaries& dirichlet_boundary() const;
        const Boundaries& farfield_boundary() const;
        const Boundaries& interog_boundary() const;
        const Boundaries& empty_boundary() const;
        const Tag& interior_boundary() const;
        void set_interior_boundary(const Tag&);
        //bool is_wall() const;
        //bool is_dirichlet() const;
        //bool is_interior() const;
        void set_point_tag(int i, const Tag& t);
        void merge(const MeshCell& other);
        bool operator<(const MeshCell& other) const;
        bool operator==(const MeshCell& other) const;
        void deparent_cell_from_faces(const Tag& tag);
        template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            ar & tag_;
            ar & first_tag_;
            ar & root_parent_mesh_;
            ar & parent_mesh_;
            ar & pnei_;
            ar & point_;
            ar & face_;
            ar & poly_;
            ar & OGA_cell_type_;
            ar & wall_boundary_;
            ar & dirichlet_boundary_;
            ar & farfield_boundary_;
            ar & empty_boundary_;
            ar & interog_boundary_;
            ar & interior_boundary_;
            ar & erase_;
            ar & btype_;
            ar & prim_;
            ar & cons_sp1_;
            ar & cons_s_;
            ar & cons_n_;
            ar & cons_nm1_;
            ar & sumarea_;

            assert(!prim_.isnan());
            assert(!cons_sp1_.isnan());
        }
    };

    int nnon_boundary_face(const MeshCell& mc);
    void minmax(const mcc& cell, Vector3& min, Vector3& max);
}

#endif
