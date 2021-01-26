#ifndef MESHFACE_H
#define	MESHFACE_H

#include "meshpoint.h"
#include "matrix.h"
#include "facetag.h"
#include "boutype.h"
#include "freestream.h"

namespace Tailor
{
    class MeshCell;

    class MeshFace
    {    
        bool R_checked_;
        BouType btype_;
        //Tag tag_;
        //bool is_boundary_;
        std::vector<Tag> parent_cell_;
        std::vector<Tag> mesh_point_;
        Polygon face_;
        FaceTag tag_;
        //MeshFace* faceaddr_;

        vec3<double> vf_;
        double vgn_;
        //Matrix<NVAR, NVAR> ML_;
        //Matrix<NVAR, NVAR> MR_;
        Matrix<NVAR, NVAR> M_;
        double max_eigen_;

        friend class boost::serialization::access;
        friend class Solver;

        public:

        MeshFace();
        //MeshFace(const MeshPoint& mp0, const MeshPoint& mp1);
        MeshFace(std::vector<MeshPoint>& pts);

        bool R_checked() const;
        void set_R_checked(bool b);
        const vec3<double>& vf() const;
        double vgn() const;
        void face_velocity(double dt, const Freestream& fs, const Component& compo);
        //const Matrix<NVAR, NVAR>& ML() const;
        //const Matrix<NVAR, NVAR>& MR() const;
        const Matrix<NVAR, NVAR>& M() const;
        const Tag& left_cell() const;
        const Tag& right_cell() const;

        //void set_faceaddr(MeshFace* addr);
        //MeshFace* faceaddr();
        //const MeshFace* faceaddr() const;
        void reset_btype();
        BouType btype() const;
        void set_btype(BouType);
        //void set_as_boundary();
        bool is_boundary() const;
        //void set_tag(const Tag& t);
        //void set_tag(int mc, int nc, BouType btype);
        void set_tag(const FaceTag& ftag);
        const FaceTag& tag() const;
        void rotate(double angle, int axis, const vec3<double>& rot_point);
        void move(const vec3<double>& v);
        const Polygon& face() const;
        const std::vector<Tag>& parent_cell() const;
        const std::vector<Tag>& mesh_point() const;
        const Tag& mesh_point(int i) const;
        const Tag& parent_cell(int i) const;
        void add_parent_cell(const Tag& celltag);
        void remove_parent_cells();
        void remove_parent_cell(const Tag& ic);
        template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            ar & parent_cell_;
            ar & mesh_point_;
            ar & face_;
            //ar & is_boundary_;
            ar & tag_;
            ar & btype_;
            ar & vgn_;
            ar & vf_;
        }
    };

    FaceTag gen_face_tag(const MeshFace& mf);
    vec3<double> tangent_vector(const vec3<double>& n);
}

#endif
