#include "geom.h"

namespace Tailor
{
    RegularHexahedron::RegularHexahedron()
    {
        set_shape(Shape::hex);
        min_(0) = TAILOR_BIG_POS_NUM;
        min_(1) = TAILOR_BIG_POS_NUM;
        min_(2) = TAILOR_BIG_POS_NUM;

        max_(0) = TAILOR_BIG_NEG_NUM;
        max_(1) = TAILOR_BIG_NEG_NUM;
        max_(2) = TAILOR_BIG_NEG_NUM;
        //assert(!faces().empty());
        //// not faces set here. calling volume would result in error!
    }

    RegularHexahedron::RegularHexahedron(Vector3 min, Vector3 max)
    {
        set_shape(Shape::hex);
        set_bbox(min, max);
        set_vertices_from_bbox();
        set_faces();
        assert(!faces().empty());
    }

    RegularHexahedron::RegularHexahedron(double xmin, double ymin, double zmin, double xmax, double ymax, double zmax)
    {
        set_shape(Shape::hex);
        set_bbox(xmin, ymin, zmin, xmax, ymax, zmax);
        set_vertices_from_bbox();
        set_faces();
        assert(!faces().empty());
    }

    void RegularHexahedron::set_bbox_from_points(const std::vector<Point>& points)
    {
        double xmin = TAILOR_BIG_POS_NUM;
        double ymin = TAILOR_BIG_POS_NUM;
        double zmin = TAILOR_BIG_POS_NUM;

        double xmax = TAILOR_BIG_NEG_NUM;
        double ymax = TAILOR_BIG_NEG_NUM;
        double zmax = TAILOR_BIG_NEG_NUM;

        for (int p=0; p<points.size(); ++p)
        {
            //std::cout << p << std::endl;
            double vx = points[p].r(0);
            double vy = points[p].r(1);
            double vz = points[p].r(2);

            xmin = std::min(xmin, vx);
            ymin = std::min(ymin, vy);
            zmin = std::min(zmin, vz);

            xmax = std::max(xmax, vx);
            ymax = std::max(ymax, vy);
            zmax = std::max(zmax, vz);
        }

        set_min(Vector3(xmin, ymin, zmin));
        set_max(Vector3(xmax, ymax, zmax));
    }

    bool RegularHexahedron::extend(const RegularHexahedron& other)
    {
        bool modified = false;

        double xmin = min(0);
        double ymin = min(1);
        double zmin = min(2);
        double xmax = max(0);
        double ymax = max(1);
        double zmax = max(2);

        if (other.min(0) < xmin)
        {
            xmin = other.min(0);
            modified = true;
        }
        if (other.min(1) < ymin)
        {
            ymin = other.min(1);
            modified = true;
        }
        if (other.min(2) < zmin)
        {
            zmin = other.min(2);
            modified = true;
        }
        if (other.max(0) > xmax)
        {
            xmax = other.max(0);
            modified = true;
        }
        if (other.max(1) > ymax)
        {
            ymax = other.max(1);
            modified = true;
        }
        if (other.max(2) > zmax)
        {
            zmax = other.max(2);
            modified = true;
        }

        //set(xmin, ymin, zmin, xmax, ymax, zmax);
        set_bbox(xmin, ymin, zmin, xmax, ymax, zmax);
        set_vertices_from_bbox();
        set_faces();

        return modified;
    }

    RegularHexahedron::RegularHexahedron(const std::vector<Point>& vertex)
    {
        //std::cout << "setting shape" << std::endl;
        set_shape(Shape::hex);
        //std::cout << "setting bbix from points" << std::endl;
        //std::cout << "vertex.size(): " << vertex.size() << std::endl;
        set_bbox_from_points(vertex);
        //std::cout << "setting vertices from bbox" << std::endl;
        set_vertices_from_bbox();
        //std::cout << "setting faces" << std::endl;
        set_faces();
        //std::cout << "alalidone" << std::endl;
    }

    double RegularHexahedron::volume_of_overlap(const RegularHexahedron& b, RegularHexahedron& aabb)
    {
        double ov = 0.;

        double left = std::max(min(0), b.min(0));
        double right = std::min(max(0), b.max(0));
        double bottom = std::max(min(1), b.min(1));
        double top = std::min(max(1), b.max(1));
        double back = std::max(min(2), b.min(2));
        double front = std::min(max(2), b.max(2));

        if (left >= right) {
            return ov;
        }

        if (bottom >= top) {
            return ov;
        }

        if (back >= front) {
            return ov;
        }

        ov = (top - bottom) * (right - left) * (front - back);

        aabb = AABB(left, bottom, back, right, top, front);

        assert(ov >= 0.);

        return ov;
    }

    RegularHexahedron::RegularHexahedron(const Polyhedron& poly)
    {
        set_shape(Shape::hex);
        set_bbox(poly.min(), poly.max());
        set_vertices_from_bbox();
        set_faces();
    }

    bool RegularHexahedron::do_intersect(const RegularHexahedron& other) const
    {
        for (int i=0; i<TAILOR_N_DIM; ++i)
        {
            double dis1 = min(i) - other.max(i);
            double dis2 = other.min(i) - max(i);

            if (dis1 > TAILOR_ZERO || dis2 > TAILOR_ZERO) {
                return false;
            }
        }        

        return true;
    }

    bool RegularHexahedron::do_contain(const RegularHexahedron& other, bool strict) const
    {
        for (int i=0; i<TAILOR_N_DIM; ++i)
        {
            if (strict)
            {
                if (other.min(i) <= min(i)) {
                    return false;
                }
            }

            if (other.min(i) < min(i)) {
                return false;
            }

            if (strict)
            {
                if (other.max(i) >= max(i)) {
                    return false;
                }
            }

            if (other.max(i) > max(i)) {
                return false;
            }
        }        

        return true;
    }
}
