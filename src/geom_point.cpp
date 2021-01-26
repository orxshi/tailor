#include "geom.h"

namespace Tailor
{
    size_t Point::mem() const
    {
        size_t size = 0;

        size += sizeof(double) * 3;

        return size;
    }

    //CGAL_Point Point::cgal_point() const
    //{
        //return CGAL_Point(r_(0), r_(1));
    //}

    bool Point::operator==(const Point& other) const
    {
        return r_ == other.r();
    }

    bool Point::operator<(const Point& other) const
    {
        if (r_(0) < other.r(0))
        {
            return true;
        }
        else if (r_(0) > other.r(0))
        {
            return false;
        }
        else
        {
            if (r_(1) < other.r(1))
            {
                return true;
            }
            else if (r_(1) > other.r(1))
            {
                return false;
            }
            else
            {
                if (r_(2) < other.r(2))
                {
                    return true;
                }
                else if (r_(2) > other.r(2))
                {
                    return false;
                }
            }
        }

        return true;
    }

    Point::Point(const Point& p)
    {
        r_ = p.r_;
        assert(r_(0) == p.r_(0));
    }

    Point& Point::operator=(const Point& p)
    {
        r_ = p.r_;
        assert(r_(0) == p.r_(0));
        return *this;
    }

    const vec3<double>& Point::r() const
    {
        return r_;
    }

    double Point::r(int i) const
    {
        assert(i >= 0);
        assert(i < 3);
        return r_(i);
    }

    void Point::set_r(double x, double y, double z)
    {
        r_.set(x, y, z);
    }

    void Point::set_r(const vec3<double>& _r)
    {
        r_.set(_r(0), _r(1), _r(2));
    }

    void Point::set_r(int i, double v)
    {
        r_.set(i, v);
    }

    Point::Point(double x, double y, double z): r_(x, y, z)
    {
    }

    Point::Point(const vec3<double>& r): r_(r)
    {
    }
}
