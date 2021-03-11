#ifndef MESHPOINT_H
#define	MESHPOINT_H

#include <cassert>
#include "tag.h"
#include "geom.h"
#include "array.h"

namespace Tailor
{
    class MeshCell;
    using ParentCell = std::vector<Tag>;

    class MeshPoint
    {    
        bool erase_;
        Tag tag_;
        Tag parent_mesh_;
        Point p_;
        ParentCell parent_cell_;
        friend class boost::serialization::access;

        public:

        MeshPoint(double x, double y, double z, const Tag& tag, const Tag& parent_mesh): tag_(tag), parent_mesh_(parent_mesh), p_(x, y, z), erase_(false) {}
        MeshPoint(double x, double y, double z): MeshPoint(x, y, z, Tag(-1), Tag(-1)) {}
        MeshPoint(const Tag& tag, const Tag& parent_mesh): MeshPoint(0., 0., 0., tag, parent_mesh) {}
        MeshPoint(): MeshPoint(Tag(-1), Tag(-1)) {}

        //size_t mem() const;
        void mark_to_be_erased();
        bool erase() const;
        const Tag& parent_mesh() const;
        void remove_parent_cell(const Tag& ic);
        void rotate_point(double ang, int axis, const Vector3& rot_point);
        void move_point(const Vector3& v);
        //void set_px(double v);
        void remove_from_parent_cells();
        void set_tag(const Tag& t);
        void set_parent_mesh(const Tag& ptag);
        const Tag& tag() const;
        const Point& p() const;
        const ParentCell& parent_cell() const;
        const Tag& parent_cell(int i) const;
        void add_parent_cell(const Tag& celltag);
        void remove_parent_cells();
        void merge(const MeshPoint& other);
        bool operator<(const MeshPoint& other) const;
        bool operator==(const MeshPoint& other) const;
        template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            ar & tag_;
            ar & parent_mesh_;
            ar & p_;
            ar & parent_cell_;
        }
    };
}

#endif
