#include "geom.h"

namespace Tailor
{
    const Vector3& Polygon::min() const
    {
        return min_;
    }
    const Vector3& Polygon::max() const
    {
        return max_;
    }
    /*void cut_polygon_same_plane(Polygon& a, const Polygon& b, int ignoredim)
    {
        // cut polygon b from polygon a.
        // polygons must be in same plane.
        
        assert(std::abs(a.vertex(0).r(ignoredim) - a.vertex(1).r(ignoredim)) < TAILOR_ZERO)
        assert(std::abs(a.vertex(0).r(ignoredim) - a.vertex(2).r(ignoredim)) < TAILOR_ZERO)
        assert(std::abs(a.vertex(0).r(ignoredim) - a.vertex(3).r(ignoredim)) < TAILOR_ZERO)

        assert(std::abs(a.vertex(0).r(ignoredim) - b.vertex(0).r(ignoredim)) < TAILOR_ZERO)

        std::vector<double> b_min(TAILOR_N_DIM, TAILOR_BIG_POS_NUM);
        std::vector<double> b_max(TAILOR_N_DIM, TAILOR_BIG_NEG_NUM);

        for (int d=0; d<TAILOR_N_DIM; ++d)
        {
            if (d == ignoredim) {
                continue;
            }
            for (int v=0; v<b.vertex().size(); ++v)
            {
                b_min[d] = std::min(b.vertex(v).r(d), b_min[d]);
                b_max[d] = std::max(b.vertex(v).r(d), b_max[d]);
            }
        }

        for (int d=0; d<TAILOR_N_DIM; ++d)
        {
            if (d == ignoredim) {
                continue;
            }
        }
    }*/
    //bool Polygon::do_intersect_plane(const Polygon& other)
    //{
        // takes in two in-plane polygons.
        //Vector3 interp;
        //bool do_intersect_plane_segment(vertex_[0].r(), vertex_[1].r(), vertex_[2].r(), other.edge(0), interp);
    //}

    Polygon create_with_sweep(const Segment& segment, const Vector3& v)
    {
        assert(segment.vertex().size() == 2);
        std::vector<Point> vtx;
        std::cout << "pushing vertices" << std::endl;
        vtx.reserve(4);
        vtx.push_back(segment.vertex(0));
        vtx.push_back(segment.vertex(1));
        vtx.push_back(Point(vtx[1].r() + v));
        vtx.push_back(Point(vtx[0].r() + v));
        std::cout << "pushed vertices" << std::endl;
        assert(vtx.size() == 4);
        std::cout << "vvvv: " << vtx[0].r(0) << std::endl;
        std::cout << "vvvv: " << vtx[0].r(1) << std::endl;
        std::cout << "vvvv: " << vtx[1].r(0) << std::endl;
        std::cout << "vvvv: " << vtx[1].r(1) << std::endl;
        std::cout << "vvvv: " << vtx[2].r(0) << std::endl;
        std::cout << "vvvv: " << vtx[2].r(1) << std::endl;
        std::cout << "vvvv: " << vtx[3].r(0) << std::endl;
        std::cout << "vvvv: " << vtx[3].r(1) << std::endl;

        return Polygon(vtx);
    }

    /*bool Polygon::do_intersect(const Vector3& v) const
    {
    // any need for this function?
        for (int i=0; i<3; ++i)
        {
            if (r(i) < min_(i) || r(i) > max_(i)) {
                return false;
            }
        }

        Point p_(minx - 5., miny - 5., minz - 5.);
        Point p(r);
        Segment s(p_, p); 

        int counter = 0;
        for (const Edge& f: face_)
        {
            if (f.do_intersect(s)) {
                ++counter;
            }
        }

        if (counter % 2 == 0) {
            return false;
        }

        return true;
    }*/

    const Segment& Polygon::edge(int i) const
    {
        assert(i >= 0);
        assert(i < 2);
        return edge_[i];
    }

    void Polygon::rotate_points(double angle, double axis, const Vector3& rot_point)
    {
        RotationMatrix rm;

        for (Point& p: vertex_)
        {
            auto z = p.r() - rot_point;
            auto newz = rm.rotate(angle, axis, z);
            p.set_r(newz + rot_point);
        }

        set_bbox();

        for (Segment& e: edge_)
        {
            e.rotate_points(angle, axis, rot_point);
        }

        reset_signed_area_ = true;
        reset_normal_ = true;
        reset_centroid_ = true;
    }

    void Polygon::move_points(const Vector3& v)
    {
        for (Point& _p: vertex_)
        {
            _p.set_r(_p.r(0) + v(0), _p.r(1) + v(1), _p.r(2) + v(2));
        }

        set_bbox();

        for (Segment& e: edge_)
        {
            e.move_points(v);
        }

        reset_signed_area_ = true;
        reset_normal_ = true;
        reset_centroid_ = true;
    }

    const Polygon::Vertices& Polygon::vertex() const
    {
        return vertex_;
    }

    const Polygon::Edges& Polygon::edge() const
    {
        return edge_;
    }

    void Polygon::set_bbox()
    {
        assert(!vertex_.empty());

        min_(0) = TAILOR_BIG_POS_NUM;
        min_(1) = TAILOR_BIG_POS_NUM;
        min_(2) = TAILOR_BIG_POS_NUM;
        max_(0) = TAILOR_BIG_NEG_NUM;
        max_(1) = TAILOR_BIG_NEG_NUM;
        max_(2) = TAILOR_BIG_NEG_NUM;

        assert(!vertex_.empty());

        for (int i=0; i<vertex_.size(); ++i)
        {
            min_(0) = std::min(min_(0), vertex_[i].r(0));
            assert(!vertex_.empty());
            min_(1) = std::min(min_(1), vertex_[i].r(1));
            min_(2) = std::min(min_(2), vertex_[i].r(2));

            assert(!vertex_.empty());

            max_(0) = std::max(max_(0), vertex_[i].r(0));
            assert(!vertex_.empty());
            max_(1) = std::max(max_(1), vertex_[i].r(1));
            assert(!vertex_.empty());
            max_(2) = std::max(max_(2), vertex_[i].r(2));

            assert(!vertex_.empty());
        }
    }

    bool Polygon::degenerate() const
    {
        auto coplanar = [&](int j)
        {
            bool cop = true;
            for (int i=0; i<vertex_.size()-1; ++i)
            {
                if (std::abs(vertex_[i].r(j) - vertex_[i+1].r(j)) > TAILOR_ZERO)
                {
                    cop = false;
                    break;
                }
            }
            return cop;
        };

        int ncop = 0;
        for (int i=0; i<TAILOR_N_DIM; ++i)
        {
            if (coplanar(i))
            {
                ++ncop;
                if (ncop == 2)
                {
                    return true;
                }
            }
        }

        return false;
    }

    Polygon::Polygon(): reset_signed_area_(true), reset_normal_(true), reset_centroid_(true), signed_area_(0.)
    {
    }

    Polygon::Polygon(const std::vector<Point>& vertex): reset_signed_area_(true), reset_normal_(true), reset_centroid_(true), signed_area_(0.)
    {
        set_vertices(vertex); // because copy constructor cannot be triggered with init list: vertex_(vertex).
        assert(!vertex_.empty());
        if (!degenerate())
        {
            assert(!vertex_.empty());
            if (vertex_.size() < 3)
            {
                std::cout << "polygon vertex size: " << vertex_.size() << std::endl;
            }
            assert(vertex_.size() >= 3);
            assert(!vertex_.empty());
            set_edge();
            assert(!vertex_.empty());
            set_bbox();
            assert(!vertex_.empty());
            //normal_ = normal();
        }
        assert(!vertex_.empty());
    }

    void Polygon::set_vertices(const std::vector<Point>& vertex)
    {
        vertex_ = vertex;
    }

    double Polygon::normal(int i) const
    {
        //if (normal_) {
            //return (*normal_)(i);
        //}
        //else
        //{
            //normal_ = normal_plane(vertex_[0].r(), vertex_[1].r(), vertex_[2].r());
            //return (*normal_)(i);

            //return normal_plane(vertex_[0].r(), vertex_[1].r(), vertex_[2].r())(i);
            return normal()(i);
        //}
    }

    double Polygon::cross_len() const
    {
        auto v0 = vertex_[0].r();
        auto v1 = vertex_[1].r();
        auto v2 = vertex_[2].r();

        Vector3 cr = cross((v0-v1), (v2-v1));

        if (std::abs(cr(0)) <= TAILOR_ZERO)
        {
            cr(0) = 0.;
        }
        if (std::abs(cr(1)) <= TAILOR_ZERO)
        {
            cr(1) = 0.;
        }
        if (std::abs(cr(2)) <= TAILOR_ZERO)
        {
            cr(2) = 0.;
        }

        return len(cr);
    }

    Vector3 Polygon::normal() const
    {
        if (reset_normal_ == false)
        {
            return normal_;
        }

        //if (normal_) {
            //return *normal_;
        //}

        //normal_ = normal_plane(vertex_[0].r(), vertex_[1].r(), vertex_[2].r());
        //std::cout << "vertex0: " << vertex_[0].r(0) << " " << vertex_[0].r(1) << " " << vertex_[0].r(2) << std::endl;
        //std::cout << "vertex1: " << vertex_[1].r(0) << " " << vertex_[1].r(1) << " " << vertex_[1].r(2) << std::endl;
        //std::cout << "vertex2: " << vertex_[2].r(0) << " " << vertex_[2].r(1) << " " << vertex_[2].r(2) << std::endl;

        reset_normal_ = false;
        normal_ = normal_plane(vertex_[0].r(), vertex_[1].r(), vertex_[2].r());
        return normal_;
        //return normal_plane(vertex_[0].r(), vertex_[1].r(), vertex_[2].r());

        //return *normal_;
    }

    /*Vector3 Polygon::normal() const
    {
        // 3 points or two edges are enough to calculate normal.
        // No need to calculate area first.

        if (normal_) {
            return *normal_;
        }

        if (vertex_.size() < 3)
        {
            std::cout << "vertex size in normal: " << vertex_.size() << std::endl;
        }
        assert(vertex_.size() >= 3);

        Vector3 a = vertex_[0].r();
        Vector3 b = vertex_[1].r();
        Vector3 c = vertex_[2].r();

        auto cr = cross((a-b), (c-b)); 
        if (cr.len() <= TAILOR_ZERO)
        {
            std::cout << "a = " << a(0) << " " << a(1) << " " << a(2) << std::endl;
            std::cout << "b = " << b(0) << " " << b(1) << " " << b(2) << std::endl;
            std::cout << "c = " << c(0) << " " << c(1) << " " << c(2) << std::endl;
            std::cout << "a-b = " << (a-b)(0) << " " << (a-b)(1) << " " << (a-b)(2) << std::endl;
            std::cout << "c-b = " << (c-b)(0) << " " << (c-b)(1) << " " << (c-b)(2) << std::endl;
            std::cout << "cr = " << cr(0) << " " << cr(1) << " " << cr(2) << std::endl;
            std::cout << "cr.len = " << cr.len() << std::endl;
        }
        assert(cr.len() > TAILOR_ZERO);

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

        normal_ = cr / cr.len();

        return *normal_;
    }*/

    void Polygon::set_edge()
    {
        assert(edge_.empty());

        // vertex_ must be ready at this point.
        
        if (vertex_.size() == 2)
        {
            edge_.push_back(Segment(vertex_[0], vertex_[1]));
        }
        else if (vertex_.size() == 3)
        {
            edge_.push_back(Segment(vertex_[0], vertex_[1]));
            edge_.push_back(Segment(vertex_[1], vertex_[2]));
            edge_.push_back(Segment(vertex_[2], vertex_[0]));
        }
        else if (vertex_.size() == 4)
        {
            edge_.push_back(Segment(vertex_[0], vertex_[1]));
            edge_.push_back(Segment(vertex_[1], vertex_[2]));
            edge_.push_back(Segment(vertex_[2], vertex_[3]));
            edge_.push_back(Segment(vertex_[3], vertex_[0]));
        }

        assert(vertex_.size() == 2 || vertex_.size() == 3 || vertex_.size() == 4);
    }

    std::vector<double> Polygon::centroid_ab(int i0, int i1) const
    {
        // Taken from https://en.wikipedia.org/wiki/Centroid#Centroid_of_polygon
        
        std::vector<double> cnt(2);

        cnt[0] = 0.;
        cnt[1] = 0.;

        std::vector<double> m(4);

        for (int i=0; i<vertex_.size()-1; ++i)
        {
            m[0] = vertex_[i].r(i0);
            m[1] = vertex_[i].r(i1);
            m[2] = vertex_[i+1].r(i0);
            m[3] = vertex_[i+1].r(i1);

            double det = m[0] * m[3] - m[1] * m[2];

            cnt[0] += (m[0] + m[2]) * det;
            cnt[1] += (m[1] + m[3]) * det;
        }

        m[0] = vertex_.back().r(i0);
        m[1] = vertex_.back().r(i1);
        m[2] = vertex_[0].r(i0);
        m[3] = vertex_[0].r(i1);

        double det = m[0] * m[3] - m[1] * m[2];

        cnt[0] += (m[0] + m[2]) * det;
        cnt[1] += (m[1] + m[3]) * det;

        assert(!std::isnan(cnt[0]));
        assert(!std::isnan(cnt[1]));

        double oldbef0 = cnt[0];
        double oldbef1 = cnt[1];

        if (vertex_.size() < 3)
        {
            std::cout << "vertex size in centroid ab: " << vertex_.size() << std::endl;
        }
        assert(vertex_.size() >= 3);
        cnt[0] /= (6.*signed_area());       
        cnt[1] /= (6.*signed_area());       

        if (cnt[0] > max_(0))
        {
            std::cout << "cnt[0]: " << cnt[0] << std::endl;
            std::cout << "min_[0]: " << min_(0) << std::endl;
            std::cout << "max_[0]: " << max_(0) << std::endl;
            std::cout << "area: " << signed_area() << std::endl;

            std::cout << "cntbefore[0]: " << oldbef0 << std::endl;
            std::cout << "cntbefore[1]: " << oldbef1 << std::endl;

            std::cout << "r0: " << vertex_[0].r(0) << std::endl;
            std::cout << "r1: " << vertex_[0].r(1) << std::endl;
            std::cout << "r2: " << vertex_[0].r(2) << std::endl;

            std::cout << "r0: " << vertex_[1].r(0) << std::endl;
            std::cout << "r1: " << vertex_[1].r(1) << std::endl;
            std::cout << "r2: " << vertex_[1].r(2) << std::endl;

            std::cout << "r0: " << vertex_[2].r(0) << std::endl;
            std::cout << "r1: " << vertex_[2].r(1) << std::endl;
            std::cout << "r2: " << vertex_[2].r(2) << std::endl;

            std::cout << "r0: " << vertex_[3].r(0) << std::endl;
            std::cout << "r1: " << vertex_[3].r(1) << std::endl;
            std::cout << "r2: " << vertex_[3].r(2) << std::endl;
        }
        assert(cnt[0] <= max_(0));

        assert(signed_area() != 0.);
        assert(!std::isnan(signed_area()));

        assert(!std::isnan(cnt[0]));
        assert(!std::isnan(cnt[1]));

        return cnt;
    }

    Vector3 Polygon::centroid() const
    {
        // https://stackoverflow.com/a/2360507/1128551


        // for now I am using simple arithmetic average.
        // should work good for simple shapes such as triangle and quad.

        if (reset_centroid_ == false)
        {
            return centroid_;
        }

        double cntx = 0.;
        double cnty = 0.;
        double cntz = 0.;
        for (const auto& v: vertex_)
        {
            cntx += v.r(0);
            cnty += v.r(1);
            cntz += v.r(2);
        }

        Vector3 cnt;
        cnt(0) = cntx / vertex_.size();
        cnt(1) = cnty / vertex_.size();
        cnt(2) = cntz / vertex_.size();

        /*std::vector<double> cnt_xy = Polygon::centroid_ab(0, 1);
        std::vector<double> cnt_xz = Polygon::centroid_ab(0, 2);

        Vector3 cnt;
        cnt.set_x(cnt_xy[0]);
        cnt.set_y(cnt_xy[1]);
        cnt.set_z(cnt_xz[1]);*/

        assert(!std::isnan(cnt(0)));
        assert(!std::isnan(cnt(1)));
        assert(!std::isnan(cnt(2)));

        //if (cnt(1) > max_(1))
        //{
            //std::cout << "cnt1: " << cnt(1) << std::endl;
            //std::cout << "maxx: " << max_(1) << std::endl;
        //}

        /*if ((cnt(0) - max_(0)) >= TAILOR_ZERO)
        {
            for (int i=0; i<vertex_.size(); ++i)
            {
                std::cout << "aaa: " << vertex_[i].r(0) << std::endl;
                std::cout << "aaa: " << vertex_[i].r(1) << std::endl;
                std::cout << "aaa: " << vertex_[i].r(2) << std::endl;
            }
            std::cout << "cnt0: " << cnt(0) << std::endl;
            std::cout << "maxx: " << max_(0) << std::endl;
        }*/

        for (int i=0; i<TAILOR_N_DIM; ++i)
        {
            if ((max_(i) - cnt(i)) < 0.)
            {
                if (std::abs((max_(i) - cnt(i))) > TAILOR_ZERO)
                {
                    std::cout << "cnt: " << cnt(i) << std::endl;
                    std::cout << "min: " << min_(i) << std::endl;
                    std::cout << "max: " << max_(i) << std::endl;
                }
                assert(std::abs((max_(i) - cnt(i))) <= TAILOR_ZERO);
                cnt(i) = max_(i);
            }

            if ((min_(i) - cnt(i)) > 0.)
            {
                if (std::abs((cnt(i) - min_(i))) > TAILOR_ZERO)
                {
                    std::cout << "i: " << cnt(i) << std::endl;
                    std::cout << "cnt: " << cnt(i) << std::endl;
                    std::cout << "min: " << min_(i) << std::endl;
                    std::cout << "max: " << max_(i) << std::endl;
                    std::cout << "diffmin: " << cnt(i) - min_(i) << std::endl;
                }
                assert(std::abs((cnt(i) - min_(i))) <= TAILOR_ZERO);
                cnt(i) = min_(i);
            }
        }




        /*if ((cnt(0) - min_(0)) < 0.);
        {
            std::cout << "cnt0: " << cnt(0) << std::endl;
            std::cout << "min0: " << min_(0) << std::endl;
            std::cout << "max0: " << max_(0) << std::endl;
            std::cout << "diffmin0: " << cnt(0) - min_(0) << std::endl;
            for (int i=0; i<vertex_.size(); ++i)
            {
                std::cout << "aaa: " << vertex_[i].r(0) << std::endl;
                std::cout << "aaa: " << vertex_[i].r(1) << std::endl;
                std::cout << "aaa: " << vertex_[i].r(2) << std::endl;
            }
        }*/

        reset_centroid_ = false;
        centroid_ = cnt;

        return centroid_;
    }

    double Polygon::signed_area() const
    {
        if (reset_signed_area_ == false)
        {
            assert(!std::isnan(signed_area_));
            return signed_area_;
        }

        assert(!vertex_.empty());
        int n = vertex_.size();
        double area = 0.;
        double an, ax, ay, az; // abs value of normal and its coords
        int coord;           // coord to ignore: 1=x, 2=y, 3=z
        int i, j, k;         // loop indices

        if (vertex_.size() < 3)
        {
            std::cout << "vertex size in signed area: " << vertex_.size() << std::endl;
        }
        assert(vertex_.size() >= 3);
        auto N = normal();

        // select largest abs coordinate to ignore for projection
        ax = (N(0)>0 ? N(0) : -N(0));    // abs x-coord
        ay = (N(1)>0 ? N(1) : -N(1));    // abs y-coord
        az = (N(2)>0 ? N(2) : -N(2));    // abs z-coord

        coord = 3;                    // ignore z-coord
        if (ax > ay) {
            if (ax > az) coord = 1;   // ignore x-coord
        }
        else if (ay > az) coord = 2;  // ignore y-coord

        // compute area of the 2D projection
        assert(vertex_.size() >= 3);
        switch (coord) {
            case 1:
                for (i=1, j=2, k=0; i<n; i++, j++, k++)
                {
                    assert(i <= 3);
                    //assert(j <= 3);
                    assert(k <= 3);
                    area += (vertex_[i].r(1) * (vertex_[j==n ? 0 : j].r(2) - vertex_[k].r(2)));
                }
                break;
            case 2:
                for (i=1, j=2, k=0; i<n; i++, j++, k++)
                {
                    assert(i <= 3);
                    //assert(j <= 3);
                    assert(k <= 3);
                    area += (vertex_[i].r(2) * (vertex_[j==n ? 0 : j].r(0) - vertex_[k].r(0)));
                }
                break;
            case 3:
                for (i=1, j=2, k=0; i<n; i++, j++, k++)
                {
                    assert(i <= 3);
                    //assert(j <= 3);
                    assert(k <= 3);
                    area += (vertex_[i].r(0) * (vertex_[j==n ? 0 : j].r(1) - vertex_[k].r(1)));
                }
                break;
        }
        switch (coord) {    // wrap-around term
            case 1:
                area += (vertex_[0].r(1) * (vertex_[1].r(2) - vertex_.back().r(2)));
                break;
            case 2:
                area += (vertex_[0].r(2) * (vertex_[1].r(0) - vertex_.back().r(0)));
                break;
            case 3:
                area += (vertex_[0].r(0) * (vertex_[1].r(1) - vertex_.back().r(1)));
                break;
        }

        // scale to get area before projection
        an = sqrt( ax*ax + ay*ay + az*az); // length of normal vector
        switch (coord) {
            case 1:
                area *= (an / (2 * N(0)));
                break;
            case 2:
                area *= (an / (2 * N(1)));
                break;
            case 3:
                area *= (an / (2 * N(2)));
        }

        signed_area_ = -area;
        //return *signed_area_;

        reset_signed_area_ = false;
        //return -area;
        assert(!std::isnan(signed_area_));
        return signed_area_;
    }

    /*double Polygon::signed_area() const
    {
        // Adapted from http://mathworld.wolfram.com/PolygonArea.html.
        //
        if (signed_area_) {
            return *signed_area_;
        }

        assert(vertex_.size() > 0);
        assert(vertex_.size() < 5);

        double sa = 0.;
        //signed_area_(0.);
        //arma::Mat<double>::fixed<2, 2> m;
        std::vector<double> m(4);

        for (int i=0; i<vertex_.size()-1; ++i)
        {
            m[0] = vertex_[i].r(0);
            m[1] = vertex_[i].r(1);
            m[2] = vertex_[i+1].r(0);
            m[3] = vertex_[i+1].r(1);

            double det = m[0] * m[3] - m[1] * m[2];

            //signed_area += arma::det(m);
            sa += det;
        }

        m[0] = vertex_.back().r(0);
        m[1] = vertex_.back().r(1);
        m[2] = vertex_[0].r(0);
        m[3] = vertex_[0].r(1);

        //signed_area += arma::det(m);
        double det = m[0] * m[3] - m[1] * m[2];
        sa += det;
        sa *= 0.5;

        signed_area_ = sa;

        return *signed_area_;
    }*/

    double magcross2D(const Vector3& a, const Vector3& b)
    {
        return std::abs(a(0) * b(1) - a(1) * b(0));
    }

    //Vector3 normalize(const Vector3& a)
    //{
        //double mag = std::sqrt(a(0)*a(0) + a(1)*a(1));
        //Vector3 b;
        //b.set(a(0)/mag, b(0)/mag);
        //return b;
    //}

    const Point& Polygon::vertex(int i) const
    {
        return vertex_[i];
    }

    bool Polygon::do_intersect(const Vector3& r, int ignoredim, bool verbose) const
    {
        // 2D polygon-point inclusion test.

        assert(ignoredim >= 0);
        assert(ignoredim <= TAILOR_N_DIM);

        for (int i=0; i<TAILOR_N_DIM; ++i)
        {
            if (i == ignoredim) {
                continue;
            }
            
            //if (r(i) < min_(i) || r(i) > max_(i)) {
                //return false;
            //}

                if (verbose) {

                    std::cout << "dim: " << i << std::endl;
                    std::cout << "r(i): " << r(i) << std::endl;
                    std::cout << "min_(i): " << min_(i) << std::endl;
                    std::cout << "max_(i): " << max_(i) << std::endl;
                    std::cout << "r(i) < (min_(i) - TAILOR_ZERO): " << (r(i) < (min_(i) - TAILOR_ZERO)) << std::endl;
                    std::cout << "r(i) > (max_(i) + TAILOR_ZERO): " << (r(i) > (max_(i) + TAILOR_ZERO)) << std::endl;
                }

            if ((r(i) < (min_(i) - TAILOR_ZERO)) || (r(i) > (max_(i) + TAILOR_ZERO))) {
                if (verbose) {
                    std::cout << "out of bounds of polygon" << std::endl;
                }
                return false;
            }
        }

        if (verbose) {
            std::cout << "in bounds of polygon" << std::endl;
        }

        //double off = edge_[0].len() / 2.;

        //Vector3 a(min_(0) - off / 2., min_(1) - off, min_(2) - off / 3.);
        auto a = centroid();
        //a.set(ignoredim, 0.);
        a(ignoredim) = min_(ignoredim);
        for (int i=0; i<TAILOR_N_DIM; ++i)
        {
            if (i == ignoredim) {
                continue;
            }

            a(i) = a(i) - (max_(i) - min_(i));
            break;
        }

        Point p_(a);
        Point p(r);
        Segment s(p_, p);

        if (verbose) {
        std::cout << "segm vertex: " << s.vertex(0).r(0) << " " << s.vertex(0).r(1) << " " << s.vertex(0).r(2) << std::endl;
        std::cout << "segm vertex: " << s.vertex(1).r(0) << " " << s.vertex(1).r(1) << " " << s.vertex(1).r(2) << std::endl;
        }

        int counter = 0;
        for (const Segment& f: edge_)
        {
            //Vector3 t0(f.vertex(0).r());
            //Vector3 t1(f.vertex(1).r());
            //t0.set(ignoredim, 0.);
            //t1.set(ignoredim, 0.);
            //Point pp0(t0);
            //Point pp1(t1);
            //Segment ff(pp0, pp1);
            if (verbose)
            {
                std::cout << "ff: " << f.vertex(0).r(0) << " " << f.vertex(0).r(1) << " " << f.vertex(0).r(2) << std::endl;
                std::cout << "ff: " << f.vertex(1).r(0) << " " << f.vertex(1).r(1) << " " << f.vertex(1).r(2) << std::endl;
            }
            //if (ff.do_intersect(s, true, ignoredim)) {
            if (f.do_intersect(s, true, ignoredim, verbose)) {
                if (verbose) 
                {
                    std::cout << "inter: yes" <<  std::endl;
                }
                ++counter;
            }
            else
            {
                if (verbose) 
                {
                    std::cout << "inter: no" <<  std::endl;
                }
            }
        }

        if (verbose)
        {
            std::cout << "counter: " << counter << std::endl;
        }

        if (counter == 0 || counter % 2 == 0) {
            return false;
        }

        return true;
    }

    bool Polygon::do_plane_intersect(const Segment& s, Vector3& interp, bool verbose) const
    {
        const Vector3& p0 = s.vertex(0).r();
        const Vector3& p1 = s.vertex(1).r();

        auto w = p0 - vertex_[0].r();
        auto u = p1 - p0;
        const auto& n = normal();

        if (verbose)
        {
            //if (n(0) == 1 && n(1) == 0 && n(2) == 0)
            {
            std::cout << "p0: " << p0(0) << " " << p0(1) << " " << p0(2) << std::endl;
            std::cout << "p1: " << p1(0) << " " << p1(1) << " " << p1(2) << std::endl;
            std::cout << "v0: " << vertex_[0].r(0) << " " << vertex_[0].r(1) << " " << vertex_[0].r(2) << std::endl;
            std::cout << "w: " << w(0) << " " << w(1) << " " << w(2) << std::endl;
            std::cout << "u: " << u(0) << " " << u(1) << " " << u(2) << std::endl;
            std::cout << "n: " << n(0) << " " << n(1) << " " << n(2) << std::endl;
            std::cout << "dot: " << std::abs(dot(n, u)) << std::endl;
            }
        }

        // first check if the polygon and the segment are parallel.
        if (std::abs(dot(n, u)) <= TAILOR_ZERO) {

            return false;
        }

        double si = -dot(n, w) / dot(n, u);

        if (verbose)
        {
            //if (n(0) == 1 && n(1) == 0 && n(2) == 0)
            {
            std::cout << "-dot(n, w): " << -dot(n, w) << std::endl;
            std::cout << "dot(n, u): " << dot(n, u) << std::endl;
            std::cout << "si: " << si << std::endl;
            }
        }

        //if (si < 0 || si > 1) {
            //return false;
        //}

        if ((si < -TAILOR_ZERO) || (si > (1 + TAILOR_ZERO))) {
            return false;
        }

        interp = p0 + u * si;

        if (verbose)
        {
            //if (n(0) == 1 && n(1) == 0 && n(2) == 0)
            {
            std::cout << "interp(0)" << interp(0) << std::endl;
            std::cout << "interp(1)" << interp(1) << std::endl;
            std::cout << "interp(2)" << interp(2) << std::endl;
            std::cout << "cond0: " << (interp(0) == p1(0)) << std::endl;
            std::cout << "cond1: " << (interp(1) == p1(1)) << std::endl;
            std::cout << "cond2: " << (interp(2) == p1(2)) << std::endl;
            std::cout << "diff: " << (interp(0) - vertex_[0].r(0)) - (w(0) + si * u(0)) << std::endl;
            std::cout << "diff2: " << interp(0) - p1(0) << std::endl;
            std::cout << "eq: " << (interp(0) == p1(0)) << std::endl;
            std::cout << "eq2: " << si - 1.0 << std::endl;
            }
        }

        if (std::abs(si - 0.) < TAILOR_ZERO || std::abs(si - 1.) < TAILOR_ZERO)
        {
            int eqdim = -1;
            for (int i=0; i<TAILOR_N_DIM; ++i)
            {
                eqdim = i;
                for (auto it = vertex_.begin(); it != std::prev(vertex_.end()); ++it)
                {
                    if (it->r(i) != std::next(it)->r(i))
                    {
                        eqdim = -1;
                    }
                }
                
                //for (const auto& v: vertex_)
                //{
                    //if (v.r(i) != v.r(i))
                    //{
                        //eqdim = -1;
                    //}
                //}

                if (eqdim != -1) {
                    break;
                }
            }

            if (verbose)
            {
                std::cout << "eqdim: " << eqdim << std::endl;
            }
             
            if (eqdim != -1)
            {
                if (verbose)
                {
                    //std::cout << "condd: " << (std::abs(interp(eqdim) - p1(eqdim)) > TAILOR_ZERO) << std::endl;
                    std::cout << "condd: " << (std::abs(interp(eqdim) - p1(eqdim)) > 0.) << std::endl;
                }
                //if (std::abs(interp(eqdim) - p1(eqdim)) > TAILOR_ZERO)
                if (std::abs(interp(eqdim) - p1(eqdim)) > 0.)
                {
                    if (p0(eqdim) < vertex_[0].r(eqdim))
                    {
                        if (p1(eqdim) < vertex_[0].r(eqdim))
                        {
                            if (verbose)
                            {
                                std::cout << "condd11" << std::endl;
                            }
                            return false;
                        }
                        else
                        {
                            if (verbose)
                            {
                                std::cout << "condd12" << std::endl;
                            }
                            return true;
                        }
                    }
                    else if (p0(eqdim) > vertex_[0].r(eqdim))
                    {
                        if (p1(eqdim) < vertex_[0].r(eqdim))
                        {
                            if (verbose)
                            {
                                std::cout << "condd21" << std::endl;
                            }
                            return true;
                        }
                        else
                        {
                            if (verbose)
                            {
                                std::cout << "condd22" << std::endl;
                            }
                            return false;
                        }
                    }
                    else
                    {
                        //for (const auto& v: vertex_)
                        //{
                        //    std::cout << "r(0): " << v.r(0) << std::endl;
                        //    std::cout << "r(1): " << v.r(1) << std::endl;
                        //    std::cout << "r(2): " << v.r(2) << std::endl;
                        //}
                        //std::cout << "si: " << si << std::endl;
                        //std::cout << "p0(0): " << p0(0) << std::endl;
                        //std::cout << "p0(1): " << p0(1) << std::endl;
                        //std::cout << "p0(2): " << p0(2) << std::endl;
                        //std::cout << "p1(0): " << p1(0) << std::endl;
                        //std::cout << "p1(1): " << p1(1) << std::endl;
                        //std::cout << "p1(2): " << p1(2) << std::endl;
                        //std::cout << "interp(0): " << p1(0) << std::endl;
                        //std::cout << "interp(1): " << p1(1) << std::endl;
                        //std::cout << "interp(2): " << p1(2) << std::endl;
                        //std::cout << "eqdim: " << eqdim << std::endl;
                        //std::cout << "vertex_[0].r(eqdim): " << vertex_[0].r(eqdim) << std::endl;
                        //assert(false);
                    }
                    
                    /*if (p1(eqdim) < vertex_[0].r(eqdim))
                    {
                        if (verbose)
                        {
                            std::cout << "condd1" << std::endl;
                        }
                        return false;
                    }
                    else if (p1(eqdim) > vertex_[0].r(eqdim))
                    {
                        if (verbose)
                        {
                            std::cout << "condd2" << std::endl;
                        }
                        return true;
                    }
                    else
                    {
                        assert(false);
                    }*/
                }
            }

            /*for (const auto& t: s.vertex())
            {
                bool found = true;
                for (int i=0; i<TAILOR_N_DIM; ++i)
                {
                    if (std::abs(t.r(i) - interp(i)) > TAILOR_ZERO)
                    {
                        found = false;
                        break;
                    }
                }

                if (found)
                {
                    std::vector<Point> pts(vertex().begin(), vertex().end());

                    std::vector<double> nn = {std::abs(n(0)), std::abs(n(1)), std::abs(n(2))};
                    auto it = std::max_element(nn.begin(), nn.end());
                    int domdim = std::distance(nn.begin(), it);

                    auto aabb = AABB(pts);
                    for (int i=0; i<TAILOR_N_DIM; ++i)
                    {
                        if (i == domdim)
                        {
                            if (i == 0)
                            //if (t.r(i) - aabb.min(i) < -TAILOR_ZERO || aabb.max(i) - t.r(i) < -TAILOR_ZERO)
                            {
                                std::cout << "i: " << i << std::endl;
                                std::cout << "aabb.min(i): " << aabb.min(i) << std::endl;
                                std::cout << "aabb.max(i): " << aabb.max(i) << std::endl;
                                std::cout << "p0(i): " << p0(i) << std::endl;
                                std::cout << "u(i): " << u(i) << std::endl;
                                std::cout << "sa(i): " << s.vertex(0).r(i) << std::endl;
                                std::cout << "sb(i): " << s.vertex(1).r(i) << std::endl;
                                //std::cout << "interp(i): " << interp(i) << std::endl;
                                std::cout << "si: " << si << std::endl;
                                std::cout << "cond1: " << (si < 0.) << std::endl;
                                std::cout << "diffmin: " << t.r(i) - aabb.min(i) << std::endl;
                                std::cout << "diffmax: " << aabb.max(i) - t.r(i)<< std::endl;
                                std::cout << "v(0): " << vertex_[0].r(i) << std::endl;
                                std::cout << "v(1): " << vertex_[1].r(i) << std::endl;
                                std::cout << "v(2): " << vertex_[2].r(i) << std::endl;
                                std::cout << "v(3): " << vertex_[3].r(i) << std::endl;
                            }
                            assert(t.r(i) - aabb.min(i) > -TAILOR_ZERO);
                            assert(aabb.max(i) - t.r(i) > -TAILOR_ZERO);
                        }
                    }
                }
            }*/
        }


        return true;
    }

    bool Polygon::do_intersect(const Segment& s, bool& just_on_face, bool verbose) const
    {
        Vector3 interp;
        //bool inter = do_intersect_plane_segment(vertex_[0].r(), vertex_[1].r(), vertex_[2].r(), s, interp);
        bool inter = do_plane_intersect(s, interp, verbose);

        if (!inter) {
            if (verbose) {
                std::cout << "not intersecting plane" << std::endl;
            }
            return false;
        }

        if (verbose) {
            std::cout << "intersecting plane" << std::endl;
        }

        // this check may make si==1 check redundant.
        if (std::abs(interp(0) - s.vertex(0).r(0)) < TAILOR_ZERO)
        {
            if (std::abs(interp(1) - s.vertex(0).r(1)) < TAILOR_ZERO)
            {
                if (std::abs(interp(2) - s.vertex(0).r(2)) < TAILOR_ZERO)
                {
                    just_on_face = true;
                    if (verbose) {
                        std::cout << "just on face" << std::endl;
                    }
                    return true;
                }
            }
        }
        if (std::abs(interp(0) - s.vertex(1).r(0)) < TAILOR_ZERO)
        {
            if (std::abs(interp(1) - s.vertex(1).r(1)) < TAILOR_ZERO)
            {
                if (std::abs(interp(2) - s.vertex(1).r(2)) < TAILOR_ZERO)
                {
                    just_on_face = true;
                    if (verbose) {
                        std::cout << "just on face" << std::endl;
                    }
                    return true;
                }
            }
        }

        const auto& n = normal();
        std::vector<double> nn = {std::abs(n(0)), std::abs(n(1)), std::abs(n(2))};
        auto it = std::max_element(nn.begin(), nn.end());
        int ignoredim = std::distance(nn.begin(), it);

        //interp.set(ignoredim, 0.);
        interp(ignoredim) = min_(ignoredim);
        return do_intersect(interp, ignoredim, verbose);
    }

    /*bool Polygon::do_intersect(const Segment& s) const
    {
        //https://stackoverflow.com/questions/5666222/3d-line-plane-intersection
        //https://developer.blender.org/diffusion/B/browse/master/source/blender/blenlib/intern/math_geom.c;6856ea06425357921bd9a7e5089c2f2adbf1413a$1352

        Vector3 p0 = s.vertex(0).r();
        Vector3 p1 = s.vertex(1).r();
        Vector3 v0 = vertex_[0].r();

        auto n = normal();

        auto u = p1 - p0;
        auto h = p0 - v0;
        auto dot = dotp(n, u);

        if (std::abs(dot) > TAILOR_ZERO)
        {
            double lambda = -dotp(n, h) / dot;
            if (lambda < 0. || lambda > 1.)
            {
                return false;
            }
            else
            {
                return true;
            }
        }
        else
        {
            return false;
        }
    }*/

    /*bool Polygon::do_intersect(const Segment& s) const
    {
        // https://math.stackexchange.com/questions/47594/plane-intersecting-line-segment

        Vector3 p0 = s.vertex(0).r();
        Vector3 p1 = s.vertex(1).r();
        Vector3 v0 = vertex_[0].r();
        if (vertex_.size() < 3)
        {
            std::cout << "vertex size in signed area: " << vertex_.size() << std::endl;
        }
        auto n = normal();

        auto u = p0 - v0;
        auto v = p1 - v0;

        double uu = dotp(n, u);
        double vv = dotp(n, v);

        double res = uu * vv;

        std::cout << "p0: " << p0(0) << std::endl;
        std::cout << "p0: " << p0(1) << std::endl;
        std::cout << "p0: " << p0(2) << std::endl;

        std::cout << "p1: " << p1(0) << std::endl;
        std::cout << "p1: " << p1(1) << std::endl;
        std::cout << "p1: " << p1(2) << std::endl;

        std::cout << "n: " << n(0) << std::endl;
        std::cout << "n: " << n(1) << std::endl;
        std::cout << "n: " << n(2) << std::endl;

        std::cout << "u: " << u(0) << std::endl;
        std::cout << "u: " << u(1) << std::endl;
        std::cout << "u: " << u(2) << std::endl;

        std::cout << "v: " << v(0) << std::endl;
        std::cout << "v: " << v(1) << std::endl;
        std::cout << "v: " << v(2) << std::endl;

        std::cout << "uu: " << uu << std::endl;
        std::cout << "vv: " << vv << std::endl;

        std::cout << "res: " << res << std::endl;

        if (res <= TAILOR_ZERO) {
            return true;
        }

        return false;
    }*/

    /*bool Polygon::do_intersect(const Segment& s) const
    {
        // based on http://geomalgorithms.com/a05-_intersect-1.html

        Vector3 p0 = s.vertex(0).r();
        Vector3 p1 = s.vertex(1).r();
        Vector3 v0 = vertex_[0].r();

        auto w = v0 - p0;
        auto u = p1 - p0;
        auto n = normal();

        // first check if the polygon and the segment are parallel.
        if (std::abs(dotp(n, u)) <= TAILOR_ZERO) {

            std::cout << "segment and polygron are parallel" << std::endl;
            return false;
        }

        std::cout << "segment and polygron are not parallel" << std::endl;

        double si = -dotp(n, w) / dotp(n, u);

            std::cout << "n: " << n(0) << std::endl;
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

            std::cout << "si " << si << std::endl;

        if (si < 0 || si > 1) {
            return false;
        }

        return true;
    }*/

    //std::vector<CGAL_Point> Polygon::cgalpoint() const
    //{
        //std::vector<CGAL_Point> pts;
        //for (const Point& p: vertex_)
        //{
            //pts.push_back(CGAL_Point(p.r(0), p.r(1)));
        //}
//
        //return pts;
    //}

    //Polygon_2 Polygon::cgalpolygon() const
    //{
        //auto pts = this->cgalpoint();
        //return Polygon_2(pts.begin(), pts.end());
    //}

}
