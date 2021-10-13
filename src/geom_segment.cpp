#include "geom.h"

namespace Tailor
{
    Vector3 Segment::normal() const
    {
        double dx = vertex_[1].r(0) - vertex_[0].r(0);
        double dy = vertex_[1].r(1) - vertex_[0].r(1);
        return Vector3(dy, -dx, 0.);
    }

    void Segment::set_bbox()
    {
        if (vertex_.empty()) {
            return;
        }

        min_(0) = TAILOR_BIG_POS_NUM;
        min_(1) = TAILOR_BIG_POS_NUM;
        min_(2) = TAILOR_BIG_POS_NUM;

        max_(0) = TAILOR_BIG_NEG_NUM;
        max_(1) = TAILOR_BIG_NEG_NUM;
        max_(2) = TAILOR_BIG_NEG_NUM;

        assert(vertex_.size() == 2);

        for (int i=0; i<2; ++i)
        {
            min_(0) = std::min(min_(0), vertex_[i].r(0));
            min_(1) = std::min(min_(1), vertex_[i].r(1));
            min_(2) = std::min(min_(2), vertex_[i].r(2));

            max_(0) = std::max(max_(0), vertex_[i].r(0));
            max_(1) = std::max(max_(1), vertex_[i].r(1));
            max_(2) = std::max(max_(2), vertex_[i].r(2));
        }
    }

    void Segment::rotate_points(double angle, double axis, const Vector3& rot_point)
    {
        RotationMatrix rm;

        for (Point& p: vertex_)
        {
            auto z = p.r();
            rm.rotate(angle, axis, rot_point, z);
            p.set_r(z);
        }
    }

    void Segment::move_points(const Vector3& v)
    {
        for (Point& _p: vertex_)
        {
            _p.set_r(_p.r(0) + v(0), _p.r(1) + v(1), _p.r(2) + v(2));
        }
    }

    //Vector3 Segment::area() const
    //{
        //assert(vertex_.size() == 2);
//
        //double dx = vertex_[1].r(0) - vertex_[0].r(0);
        //double dy = vertex_[1].r(1) - vertex_[0].r(1);
//
        //return Vector3 (dy, -dx);
        ////return Vector3 (-dy, dx);
    //}

    double Segment::len() const
    {
        assert(vertex_.size() == 2);

        //return CGAL::to_double(cgal_segment().squared_length());
        Vector3 v = vertex_[0].r() - vertex_[1].r();

        return v.len();
    }

    /*Vector3 Segment::centroid() const
    {
        Vector3 _centroid;

        //auto p = CGAL::midpoint(cgal_segment().point(0), cgal_segment().point(1));
        //_centroid.set(CGAL::to_double(p.x()), CGAL::to_double(p.y()));

        double cnt[2];
        cnt[0] = 0.5 * (vertex_[0].r(0) + vertex_[1].r(0));
        cnt[1] = 0.5 * (vertex_[0].r(1) + vertex_[1].r(1));

        _centroid.set(cnt[0], cnt[1]);

        return _centroid;
    }*/

    //CGAL_Segment Segment::cgal_segment() const
    //{
        //return CGAL_Segment(CGAL_Point(vertex_[0].r(0), vertex_[0].r(1)), CGAL_Point(vertex_[1].r(0), vertex_[1].r(1)));
    //}

    Segment::Segment(const Segment& other)
    {
        vertex_ = other.vertex();
        set_bbox();
    }

    const Vector3& Segment::min() const
    {
        return min_;
    }
    
    const Vector3& Segment::max() const
    {
        return max_;
    }

    double Segment::min(int i) const
    {
        return min_(i);
    }

    double Segment::max(int i) const
    {
        return max_(i);
    }

    const Point& Segment::vertex(int i) const
    {
        assert(i >= 0);
        assert(i < 2);
        return vertex_[i];
    }

    const Terminal& Segment::vertex() const
    {
        return vertex_;
    }

    bool Segment::do_intersect(const Point& p, int ignoredim) const
    {
        //return CGAL_Segment(CGAL_Point(vertex(0).r(0), vertex(0).r(1)), CGAL_Point(vertex(1).r(0), vertex(1).r(1))).has_on(CGAL_Point(p.r(0), p.r(1)));
        //return CGAL_Segment(vertex(0).cgal_point(), vertex(1).cgal_point()).has_on(CGAL_Point(p.r(0), p.r(1)));

        const Vector3& a = this->vertex_[0].r();
        const Vector3& c = this->vertex_[1].r();
        const Vector3& b = p.r();

        if (orientation(a, b, c, ignoredim) != 0) return false;

        bool all_on_axis[TAILOR_N_DIM] = {false};
        bool b_equal_a = true;
        bool b_equal_c = true;

        for (int i=0; i<TAILOR_N_DIM; ++i)
        {
            if (i == ignoredim) {
                continue;
            }
            //if (b(i) != a(i))
            if (std::abs(b(i) - a(i)) > TAILOR_ZERO)
            {
                b_equal_a  = false;
            }
        }
        for (int i=0; i<TAILOR_N_DIM; ++i)
        {
            if (i == ignoredim) {
                continue;
            }
            //if (b(i) != c(i))
            if (std::abs(b(i) - c(i)) > TAILOR_ZERO)
            {
                b_equal_c  = false;
            }
        }

        if (b_equal_a) {return true;}
        if (b_equal_c) {return true;}

        for (int i=0; i<TAILOR_N_DIM; ++i)
        {
            if (i == ignoredim) {
                continue;
            }
            //if (b(i) == a(i) && b(i) == c(i)) {all_on_axis[i] = true;}
            if (std::abs(b(i) - a(i)) < TAILOR_ZERO && std::abs(b(i) - c(i)) < TAILOR_ZERO) {all_on_axis[i] = true;}
        }

        for (int i=0; i<TAILOR_N_DIM; ++i)
        {
            if (i == ignoredim) {
                continue;
            }
            //if ( !(all_on_axis[i] || ((b(i) < std::max(a(i), c(i))) && (b(i) > std::min(a(i), c(i))))) )
            if ( !(all_on_axis[i] || ((b(i) - std::max(a(i), c(i))) < TAILOR_ZERO && (b(i) - std::min(a(i), c(i)) > TAILOR_ZERO))) )
            {
                return true;
            }
        }

        return false;
    }

    bool Segment::do_intersect(const Segment& s, bool point_on_segment_means_inter, int ignoredim, bool verbose) const
    {
        // Returns true if line segments 'p1q1' and 'p2q2' intersect.
        // If segment touches the other then no intersection.
        // Adapted from http://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/.

        if(this->vertex_.size() != 2)
        {
            std::cout << this->vertex_.size() << std::endl;
        }
        assert(this->vertex_.size() == 2);
        assert(s.vertex_.size() == 2);

        const Point& p1_ = this->vertex_[0];
        const Point& q1_ = this->vertex_[1];
        const Point& p2_ = s.vertex_[0];
        const Point& q2_ = s.vertex_[1];

        const Vector3 p1 = p1_.r();
        const Vector3 q1 = q1_.r();
        const Vector3 p2 = p2_.r();
        const Vector3 q2 = q2_.r();

        //return CGAL::do_intersect(CGAL_Segment(CGAL_Point(p1(0), p1(1)), CGAL_Point(q1(0), q1(1))), CGAL_Segment(CGAL_Point(p2(0), p2(1)), CGAL_Point(q2(0), q2(1))));

        // find orientations.
        int o1 = orientation(p1, q1, p2, ignoredim, verbose);
        int o2 = orientation(p1, q1, q2, ignoredim, verbose);
        int o3 = orientation(p2, q2, p1, ignoredim, verbose);
        int o4 = orientation(p2, q2, q1, ignoredim, verbose);    

        if (verbose)
        {
        std::cout << "p1: " << p1(0) << " " << p1(1) << " " << p1(2) << std::endl;
        std::cout << "q1: " << q1(0) << " " << q1(1) << " " << q1(2) << std::endl;
        std::cout << "p2: " << p2(0) << " " << p2(1) << " " << p2(2) << std::endl;
        std::cout << "q2: " << q2(0) << " " << q2(1) << " " << q2(2) << std::endl;
        std::cout << "o1: " << o1 << std::endl;
        std::cout << "o2: " << o2 << std::endl;
        std::cout << "o3: " << o3 << std::endl;
        std::cout << "o4: " << o4 << std::endl;
        }

        if (o1 != o2 && o3 != o4)
        {
            if (!point_on_segment_means_inter)
            {
                // p2 lies on segment p1q1.     
                if (this->do_intersect(p2_, ignoredim)) {return false;}

                // q2 lies on segment p1q1.
                if (this->do_intersect(q2_, ignoredim)) {return false;}

                // p1 lies on segment p2q2.
                if (s.do_intersect(p1_, ignoredim)) {return false;}

                // q1 lies on segment p2q2.
                if (s.do_intersect(q1_, ignoredim)) {return false;}
            }

            return true;
        }

        return false;
    }

    Segment::Segment(const Point& p0, const Point& p1) 
    {
        //vertex_ = {p0, p1};
        vertex_.push_back(p0);
        vertex_.push_back(p1);
        set_bbox();
        //cgal_segment_ = CGAL_Segment(p0.cgal_point(), p1.cgal_point());
    }
}
