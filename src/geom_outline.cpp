#include "geom.h"

namespace Tailor
{
    /*const std::vector<Segment>& Outline::segment() const
    {
        return edge_;
    }

    void Outline::set_tag(const Tag& t)
    {
        tag_ = t;
    }

    const Tag& Outline::tag() const
    {
        return tag_;
    }

    const Polygon_2& Outline::cgal_polygon() const
    {
        return polygon_;
    }

    bool Outline::do_contain(const Polygon_2& pol, const Vector3& v, bool strict) const
    {
        CGAL::Bounded_side bside = pol.bounded_side(CGAL_Point(v(0), v(1), v(2))); 
        if (bside == CGAL::ON_UNBOUNDED_SIDE) { 
            return false;
        }

        if (strict)
        {
            if (bside == CGAL::ON_BOUNDARY) { 
                return false;
            }
        }

        return true;
    }

    bool Outline::do_contain(const Triangle_2& tri, const Vector3& v, bool strict) const
    {
        CGAL::Bounded_side bside = tri.bounded_side(CGAL_Point(v(0), v(1))); 
        if (bside == CGAL::ON_UNBOUNDED_SIDE) { 
            return false;
        }

        if (strict)
        {
            if (bside == CGAL::ON_BOUNDARY) { 
                return false;
            }
        }

        return true;
    }

    void Outline::build_polygon(int dummyrank)
    {
        if (edge_.empty()) {
            return;
        }

        for (const Segment& seg: edge_) {
            polygon_.push_back(CGAL_TPoint(seg.vertex(0).r(0), seg.vertex(0).r(1)));
        }
    }
     
    bool Outline::do_contain(const Vector3& v, bool strict) const
    {
        assert(polygon_.is_simple());
        if (polygon_.is_simple())
        {
            return do_contain(polygon_, v, strict);
        }
        else
        {
            for (auto it=polygon_.vertices_begin(); it!=polygon_.vertices_end(); ++it)
            {
                std::cout << it->x() << " " << it->y() << std::endl;
            }
            assert(polygon_.is_simple());
            for (const auto& tri: tri_)
            {
                if (do_contain(tri, v, strict) == false) {
                    return false;
                }
            }

            return true;
        }
    }

    void Outline::build(const std::vector<std::vector<Vector3>>& point, int dummyrank)
    {
        assert(!point.empty());
        std::vector<Polygon_2> pgn;
        for (const auto& v: point)
        {
            assert(!v.empty());
            std::vector<CGAL_Point> cpts;
            assert(v.size() == 4);
            for (const auto& p: v) {
                cpts.push_back(CGAL_Point(p(0), p(1)));
            }
            pgn.push_back(Polygon_2(cpts.begin(), cpts.end()));
        }

        assert(!pgn.empty());

        std::vector<Polygon_with_holes_2> oi(10);

        CGAL::join(pgn.begin(), pgn.end(), oi.data());

        auto outer_boundary = oi[0].outer_boundary();
        auto cit = outer_boundary.edges_begin();
        auto cit_end = outer_boundary.edges_end();

        for(;cit!=cit_end;++cit)
        {
            edge_.push_back(Segment(Point(CGAL::to_double(cit->source().x()), CGAL::to_double(cit->source().y())), Point(CGAL::to_double(cit->target().x()), CGAL::to_double(cit->target().y()))));
        }

        assert(!edge_.empty());
    }*/
}
