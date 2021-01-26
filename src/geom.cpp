#include "geom.h"

namespace Tailor
{
    int parallel_inter(const Segment& a, const Segment& b, Segment& news)
    {
        // Works only if segments are parallel
        // Returns:
        // 0 if at least one vertex is out of borders
        // 1 if complete overlap
        // 2 if partial overlap
        
        const Segment* small = nullptr;
        const Segment* large = nullptr;

        if (a.len() < b.len())
        {
            small = &a;
            large = &b;
        }
        else
        {
            small = &b;
            large = &a;
        }

        int count = 0;
        bool aa = false;
        bool ab = false;
        bool ba = false;
        bool bb = false;
        if (small->vertex(0) == large->vertex(0)) {
            aa = true;
            ++count;
        }
        if (small->vertex(0) == large->vertex(1)) {
            ab = true;
            ++count;
        }
        if (small->vertex(1) == large->vertex(0)) {
            ba = true;
            ++count;
        }
        if (small->vertex(1) == large->vertex(1)) {
            bb = true;
            ++count;
        }
        
        double minx = std::min(large->vertex(0).r(0), large->vertex(1).r(0));
        double maxx = std::max(large->vertex(0).r(0), large->vertex(1).r(0));
        double miny = std::min(large->vertex(0).r(1), large->vertex(1).r(1));
        double maxy = std::max(large->vertex(0).r(1), large->vertex(1).r(1));

        if (count == 2)
        {
            assert(a.len() == b.len());
            return 1;
        }
        else if (count == 1)
        {
            if (small->vertex(0).r(0) < minx || small->vertex(0).r(1) < miny) {return 0;}
            if (small->vertex(0).r(0) > maxx || small->vertex(0).r(1) > maxy) {return 0;}
            if (small->vertex(1).r(0) < minx || small->vertex(1).r(1) < miny) {return 0;}
            if (small->vertex(1).r(0) > maxx || small->vertex(1).r(1) > maxy) {return 0;}

            if (aa || ab)
            {
                if (aa)
                {
                    news = Segment(small->vertex(1), large->vertex(1));
                }
                if (ab)
                {
                    news = Segment(small->vertex(1), large->vertex(0));
                }
                std::cout << "small: " << small->vertex(0).r(0) << " " << small->vertex(0).r(1) << std::endl;
                std::cout << "small: " << small->vertex(1).r(0) << " " << small->vertex(1).r(1) << std::endl;
                std::cout << "large: " << large->vertex(0).r(0) << " " << large->vertex(0).r(1) << std::endl;
                std::cout << "large: " << large->vertex(1).r(0) << " " << large->vertex(1).r(1) << std::endl;
                return 2;
            }
            else
            {
                if (ba)
                {
                    news = Segment(small->vertex(0), large->vertex(1));
                }
                if (ab)
                {
                    news = Segment(small->vertex(0), large->vertex(0));
                }
                std::cout << "small: " << small->vertex(0).r(0) << " " << small->vertex(0).r(1) << std::endl;
                std::cout << "large: " << large->vertex(0).r(0) << " " << large->vertex(0).r(1) << std::endl;
                std::cout << "large: " << large->vertex(0).r(0) << " " << large->vertex(0).r(1) << std::endl;
                std::cout << "large: " << large->vertex(1).r(0) << " " << large->vertex(1).r(1) << std::endl;
                return 2;
            }

            return 0;
        }
        else
        {
            return 0;
        }
    }

    void compare(std::vector<Segment>& edge)
    {
        for (auto a = edge.begin(); a != edge.end(); ++a)
        {
            for (auto b = std::next(a); b != edge.end(); ++b)
            {
                Segment news;
                int r = parallel_inter(*a, *b, news);

                if (r == 0) {
                    continue;
                }
                else if (r == 1) {
                    edge.erase(b);
                    edge.erase(a);
                    compare(edge);
                    return;
                }
                else if (r == 2) {
                    assert(false);
                    edge.erase(b);
                    edge.erase(a);
                    edge.push_back(news);
                    compare(edge);
                    return;
                }
            }
        }
    }

    Polygon make_outline(const std::vector<Polygon>& pgon)
    {
        std::vector<Segment> edge;
        for (const auto& p: pgon)
        {
            for (const auto& e: p.edge())
            {
                edge.push_back(e);
            }
        }

        std::cout << "before edge size: " << edge.size() << std::endl;

        compare(edge);

        std::cout << "edge size: " << edge.size() << std::endl;
        for (const auto& e: edge)
        {
            std::cout << "e0: " << e.vertex(0).r(0) << " "  << e.vertex(0).r(1) << " " << e.vertex(0).r(2) << std::endl;
            std::cout << "e1: " << e.vertex(1).r(0) << " "  << e.vertex(1).r(1) << " " << e.vertex(1).r(2) << std::endl;
        }

        std::vector<Point> pts;
        for (const auto& e: edge)
        {
            pts.push_back(e.vertex(0));
            pts.push_back(e.vertex(1));
        }

        std::sort(pts.begin(), pts.end());
        auto last = std::unique(pts.begin(), pts.end());
        pts.erase(last, pts.end());

        std::cout << "pts size: " << pts.size() << std::endl;
        for (const auto& pp: pts)
        {
            std::cout << "p: " << pp.r(0) << " "  << pp.r(1) << " " << pp.r(2) << std::endl;
        }
        assert(!pts.empty());
        // sort points based on orientation.
        return Polygon(pts);
    }

    int orientation (const vec3<double>& p, const vec3<double>& q, const vec3<double>& r, int ignoredim, bool verbose)
    {
        // Find orientation of ordered triplet (p, q, r).
        // Return following values.
        // 0 --> p, q and r are colinear
        // 1 --> clockwise
        // 2 --> counterclockwise

        // See 10th slides from following link for derivation of the formula
        // http://www.dcs.gla.ac.uk/~pat/52233/slides/Geometry1x1.pdf
        
        double val;
        if (ignoredim == 0)
        {
            val = (q(2) - p(2)) * (r(1) - q(1)) - (q(1) - p(1)) * (r(2) - q(2));
            if (verbose)
            {
                std::cout << "(q(2) - p(2)): " << (q(2) - p(2)) << std::endl;
                std::cout << "(r(1) - q(1)): " << (r(1) - q(1)) << std::endl;
                std::cout << "(q(1) - p(1)): " << (q(1) - p(1)) << std::endl;
                std::cout << "(r(2) - q(2)): " << (r(2) - q(2)) << std::endl;
                std::cout << "val: " << val << std::endl;
                std::cout << "TAILOR_ZREO: " << TAILOR_ZERO << std::endl;
            }
        }
        else if (ignoredim == 1)
        {
            val = (q(2) - p(2)) * (r(0) - q(0)) - (q(0) - p(0)) * (r(2) - q(2));
        }
        else if (ignoredim == 2)
        {
            val = (q(1) - p(1)) * (r(0) - q(0)) - (q(0) - p(0)) * (r(1) - q(1));
        }
        else
        {
            assert(false);
        }

        //if (val == 0.) return 0;  // colinear
        if (std::abs(val) < TAILOR_ZERO) return 0;  // colinear

        return (val > 0.) ? 1 : 2; // clock or counterclockwise
    }

    vec3<double> normal_plane(const vec3<double>& v0, const vec3<double>& v1, const vec3<double>& v2)
    {
        vec3<double> cr = cross((v0-v1), (v2-v1));

        if (std::abs(cr(0)) <= TAILOR_ZERO)
        {
            cr.set_x(0.);
        }
        if (std::abs(cr(1)) <= TAILOR_ZERO)
        {
            cr.set_y(0.);
        }
        if (std::abs(cr(2)) <= TAILOR_ZERO)
        {
            cr.set_z(0.);
        }

        if (cr.len() == 0.)
        {
            std::cout << "v0: " << v0(0) << " " << v0(1) << " " << v0(2) << std::endl;
            std::cout << "v1: " << v1(0) << " " << v1(1) << " " << v1(2) << std::endl;
            std::cout << "v2: " << v2(0) << " " << v2(1) << " " << v2(2) << std::endl;
            std::cout << "v0-v1: " << (v0-v1)(0) << " " << (v0-v1)(1) << " " << (v0-v1)(2) << std::endl;
            std::cout << "v2-v1: " << (v2-v1)(0) << " " << (v2-v1)(1) << " " << (v2-v1)(2) << std::endl;
            std::cout << "cr: " << cr(0) << " " << cr(1) << " " << cr(2) << std::endl;
        }
        assert(cr.len() != 0.);
        vec3<double> n = cr / cr.len();

        return n;
    }

    bool do_intersect_plane_segment(const vec3<double>& a, const vec3<double>& b, const vec3<double>& c, const Segment& s, vec3<double>& interp)
    {
        // Plane-segment intersection
        // based on http://geomalgorithms.com/a05-_intersect-1.html

        vec3<double> p0 = s.vertex(0).r();
        vec3<double> p1 = s.vertex(1).r();

        //auto w = a - p0;
        auto w = p0 - a;
        auto u = p1 - p0;
        const auto& n = normal_plane(a, b, c);

        // first check if the polygon and the segment are parallel.
        if (std::abs(dotp(n, u)) <= TAILOR_ZERO) {

            //std::cout << "segment and polygron are parallel" << std::endl;
            return false;
        }

        double si = -dotp(n, w) / dotp(n, u);

        /*std::cout << "n: " << n(0) << std::endl;
          std::cout << "n: " << n(1) << std::endl;
          std::cout << "n: " << n(2) << std::endl;

          std::cout << "u: " << u(0) << std::endl;
          std::cout << "u: " << u(1) << std::endl;
          std::cout << "u: " << u(2) << std::endl;

          std::cout << "p0: " << p0(0) << std::endl;
          std::cout << "p0: " << p0(1) << std::endl;
          std::cout << "p0: " << p0(2) << std::endl;

          std::cout << "p1: " << p1(0) << std::endl;
          std::cout << "p1: " << p1(1) << std::endl;
          std::cout << "p1: " << p1(2) << std::endl;

          std::cout << "v0: " << v0(0) << std::endl;
          std::cout << "v0: " << v0(1) << std::endl;
          std::cout << "v0: " << v0(2) << std::endl;

          std::cout << "w: " << w(0) << std::endl;
          std::cout << "w: " << w(1) << std::endl;
          std::cout << "w: " << w(2) << std::endl;

          std::cout << "si " << si << std::endl;*/

        if (si < 0 || si > 1) {
            return false;
        }

        //std::cout << "p0: " << p0(0) << std::endl;
        //std::cout << "p0: " << p0(1) << std::endl;
        //std::cout << "p0: " << p0(2) << std::endl;
        //std::cout << "u: " << u(0) << std::endl;
        //std::cout << "u: " << u(1) << std::endl;
        //std::cout << "u: " << u(2) << std::endl;
        //std::cout << "si: " << si << std::endl;
        interp = p0 + u * si;

        return true;
    }

    

}
