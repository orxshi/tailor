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

    const Vector3& Point::r() const
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
        r_(0) = x;
        r_(1) = y;
        r_(2) = z;
    }

    void Point::set_r(const Vector3& r)
    {
        r_ = r;
    }

    void Point::set_r(int i, double v)
    {
        r_(i) = v;
    }

    Point::Point(double x, double y, double z): r_(x, y, z)
    {
    }

    Point::Point(const Vector3& r): r_(r)
    {
    }
}
