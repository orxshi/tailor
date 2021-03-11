#include "geom.h"

namespace Tailor
{
    Polyhedron create_with_sweep(const Polygon& polygon, const Vector3& v)
    {
        if (polygon.vertex().size() == 3)
        {
            std::vector<Point> vtx;
            vtx.reserve(6);
            vtx.push_back(polygon.vertex(0));
            vtx.push_back(polygon.vertex(1));
            vtx.push_back(polygon.vertex(2));
            vtx.push_back(Point(vtx[0].r() + v));
            vtx.push_back(Point(vtx[1].r() + v));
            vtx.push_back(Point(vtx[2].r() + v));

            return Polyhedron(vtx, Shape::pri);
        }
        else if (polygon.vertex().size() == 4)
        {
            std::vector<Point> vtx;
            vtx.reserve(8);
            vtx.push_back(polygon.vertex(0));
            vtx.push_back(polygon.vertex(1));
            vtx.push_back(polygon.vertex(2));
            vtx.push_back(polygon.vertex(3));
            vtx.push_back(Point(vtx[0].r() + v));
            vtx.push_back(Point(vtx[1].r() + v));
            vtx.push_back(Point(vtx[2].r() + v));
            vtx.push_back(Point(vtx[3].r() + v));

            return Polyhedron(vtx, Shape::hex);
        }
        else
        {
            assert(false);
        }
    }

    size_t Polyhedron::mem() const
    {
        size_t size = 0;
        size += vertex_.mem();
        size += face_.mem();
        size += sizeof(double);
        size += sizeof(int);
        size += sizeof(double) * 6;

        return size;
    }

    Polyhedron::Polyhedron(): reset_centroid_(true), reset_volume_(true), volume_(0.)
    {
    }

    Polyhedron::Polyhedron(const std::vector<Point>& vertex, Shape shape): shape_(shape), reset_centroid_(true), reset_volume_(true), volume_(0.)
    {
        if (vertex.size() > 8)
        {
            std::cout << "vertex size: " << vertex.size() << std::endl;
            std::cout << "shape: " << static_cast<int>(shape) << std::endl;
        }
        assert(vertex.size() <= 8);
        if (shape_ == Shape::tet)
        {
            assert(vertex.size() == 4);
        }
        set_vertices(vertex); // because copy constructor cannot be triggered with init list: vertex_(vertex).
        set_faces();
        set_bbox();
        assert(!face_.empty());
    }

    Shape Polyhedron::shape() const
    {
        return shape_;
    }

    void Polyhedron::set_min(Vector3 min)
    {
        min_ = min;
    }

    void Polyhedron::set_max(Vector3 max)
    {
        max_ = max;
    }

    void Polyhedron::set_vertices_from_bbox()
    {
        //assert(vertex_.size() == 0);

        vertex_.reserve(8);

        vertex_[0].set_r(max(0), max(1), min(2));
        vertex_[1].set_r(min(0), max(1), min(2));
        vertex_[2].set_r(min(0), min(1), min(2));
        vertex_[3].set_r(max(0), min(1), min(2));

        vertex_[4].set_r(max(0), max(1), max(2));
        vertex_[5].set_r(min(0), max(1), max(2));
        vertex_[6].set_r(min(0), min(1), max(2));
        vertex_[7].set_r(max(0), min(1), max(2));
    }

    void Polyhedron::set_bbox(Vector3 min, Vector3 max)
    {
        min_ = min;
        max_ = max;
    }

    void Polyhedron::set_bbox(double xmin, double ymin, double zmin, double xmax, double ymax, double zmax)
    {
        min_(0) = xmin;
        min_(1) = ymin;
        min_(2) = zmin;

        max_(0) = xmax;
        max_(1) = ymax;
        max_(2) = zmax;
    }

    void Polyhedron::set_shape(Shape shape)
    {
        shape_ = shape;
    }

    void Polyhedron::set_bbox()
    {
        min_(0) = TAILOR_BIG_POS_NUM;
        min_(1) = TAILOR_BIG_POS_NUM;
        min_(2) = TAILOR_BIG_POS_NUM;

        max_(0) = TAILOR_BIG_NEG_NUM;
        max_(1) = TAILOR_BIG_NEG_NUM;
        max_(2) = TAILOR_BIG_NEG_NUM;

        for (int i=0; i<vertex_.size(); ++i)
        {
            min_(0) = std::min(min_(0), vertex_[i].r(0));
            min_(1) = std::min(min_(1), vertex_[i].r(1));
            min_(2) = std::min(min_(2), vertex_[i].r(2));

            max_(0) = std::max(max_(0), vertex_[i].r(0));
            max_(1) = std::max(max_(1), vertex_[i].r(1));
            max_(2) = std::max(max_(2), vertex_[i].r(2));
        }
    }

    void Polyhedron::set_vertices(const std::vector<Point>& vertex)
    {
        if (vertex.size() > 8)
        {
            std::cout << "vertex size: " << vertex.size() << std::endl;
        }
        assert(vertex.size() <= 8);
        vertex_ = vertex;
        if (shape_ == Shape::tet)
        {
            assert(vertex.size() == 4);
            assert(vertex_.size() == 4);
        }
    }

    const Polygon& Polyhedron::face(int i) const
    {
        assert(i >= 0);
        assert(i < 6);
        return face_[i];
    }

    void Polyhedron::rotate_points(double angle, double axis, const Vector3& rot_point)
    {
        RotationMatrix rm;

        for (Point& p: vertex_)
        {
            auto z = p.r() - rot_point;
            auto newz = rm.rotate(angle, axis, z);
            p.set_r(newz + rot_point);
        }

        set_bbox();
        for (Polygon& f: face_)
        {
            f.rotate_points(angle, axis, rot_point);
        }

        reset_centroid_ = true;
    }

    const Vector3& Polyhedron::min() const
    {
        //std::cout << "returning hedron min" << std::endl;
        //std::cout << "min0: " << min_(0) << std::endl;
        //std::cout << "min1: " << min_(1) << std::endl;
        //std::cout << "min2: " << min_(2) << std::endl;
        return min_;
    }
    
    const Vector3& Polyhedron::max() const
    {
        return max_;
    }

    double Polyhedron::min(int i) const
    {
        assert(i < 3);
        assert(i >= 0);
        //std::cout << "returning hedron min i" << std::endl;
        //std::cout << "minn0: " << min_(0) << std::endl;
        //std::cout << "minn1: " << min_(1) << std::endl;
        //std::cout << "minn2: " << min_(2) << std::endl;
        return min_(i);
    }

    double Polyhedron::max(int i) const
    {
        assert(i < 3);
        assert(i >= 0);
        return max_(i);
    }

    void Polyhedron::move_points(const Vector3& v)
    {
        for (Point& _p: vertex_)
        {
            _p.set_r(_p.r(0) + v(0), _p.r(1) + v(1), _p.r(2) + v(2));
        }

        set_bbox();

        for (Polygon& f: face_)
        {
            f.move_points(v);
        }

        reset_centroid_ = true;
    }

    const Polyhedron::Vertices& Polyhedron::vertices() const
    {
        return vertex_;
    }

    const Polyhedron::Faces& Polyhedron::faces() const
    {
        return face_;
    }

    bool Polyhedron::face_cross_len() const
    {
        for (const auto& ff: face_)
        {
            if (ff.cross_len() == 0.)
            {
                assert(!ff.vertex().empty());
                //for (const auto& v: ff.vertex())
                //{
                    //std::cout << "problem: " << v.r(0) << " " << v.r(1) << " " << v.r(2) << std::endl;
                //}
                return true;
            }
        }

        return false;
    }

    //double Polyhedron::volume_triface(std::vector<Polygon>::const_iterator begin, std::vector<Polygon>::const_iterator end) const
    double Polyhedron::volume_triface(const Polygon* begin, const Polygon* end) const
    {
        // no need to call directly.
        // use Polyhedron::volume() directly.

        double sum = 0.;

        //std::cout << "****************" << std::endl;
        //std::cout << "nface: " << end - begin << std::endl;
        //std::cout << "----------------" << std::endl;

        for (auto f = begin; f != end; ++f)
        {
            if (f->vertex().size() < 3)
            {
                std::cout << "vertex size in hedron volume triface: " << f->vertex().size() << std::endl;
            }
            assert(f->vertex().size() >= 3);
            assert(!std::isnan(f->signed_area()));
            assert(f->signed_area() != 0.);

            assert(!std::isnan(f->centroid()(0)));
            assert(!std::isnan(f->centroid()(1)));
            assert(!std::isnan(f->centroid()(2)));

            assert(!std::isnan(f->normal()(0)));
            assert(!std::isnan(f->normal()(1)));
            assert(!std::isnan(f->normal()(2)));

            assert(f->signed_area() != 0.);

            //assert(f->centroid()(0) != 0.);
            //assert(f->centroid()(1) != 0.);
            //assert(f->centroid()(2) != 0.);

            //assert(f->normal()(0) != 0.);
            //assert(f->normal()(1) != 0.);
            //assert(f->normal()(2) != 0.);

            //if (dotp(f->centroid(), f->normal()) == 0.)
            //{
                //std::cout << "cnt 0: " << f->centroid()(0) << std::endl;
                //std::cout << "cnt 1: " << f->centroid()(1) << std::endl;
                //std::cout << "cnt 2: " << f->centroid()(2) << std::endl;

                //std::cout << "norm 0: " << f->normal()(0) << std::endl;
                //std::cout << "norm 1: " << f->normal()(1) << std::endl;
                //std::cout << "norm 2: " << f->normal()(2) << std::endl;

                //std::cout << "area: " << std::abs(f->signed_area()) << std::endl;
            //}
            //assert(dotp(f->centroid(), f->normal()) != 0.);
            
            sum += dot(f->centroid(), f->normal()) * std::abs(f->signed_area());

            //std::cout << "tsum: " << dotp(f->centroid(), f->normal()) * std::abs(f->signed_area()) << std::endl;
            //std::cout << "----------------" << std::endl;
        }

        assert(!std::isnan(sum));
        //assert(sum != 0.);

        double volume = sum / 3.;

        return volume;
    }

    Vector3 Polyhedron::centroid() const
    {
        // this assumes that only vertices have weight.

        if (reset_centroid_ == false)
        {
            return centroid_;
        }

        double cnt_x = 0.;
        double cnt_y = 0.;
        double cnt_z = 0.;

        for (const auto& p: vertex_)
        {
            cnt_x += p.r(0);
            cnt_y += p.r(1);
            cnt_z += p.r(2);
        }

        cnt_x /= vertex_.size();
        cnt_y /= vertex_.size();
        cnt_z /= vertex_.size();

        Vector3 cnt;
        cnt(0) = cnt_x;
        cnt(1) = cnt_y;
        cnt(2) = cnt_z;

        assert((cnt(0) - min_(0)) > -TAILOR_ZERO);
        assert((cnt(1) - min_(1)) > -TAILOR_ZERO);
        assert((cnt(2) - min_(2)) > -TAILOR_ZERO);
        assert((cnt(0) - max_(0)) < TAILOR_ZERO);
        assert((cnt(1) - max_(1)) < TAILOR_ZERO);
        assert((cnt(2) - max_(2)) < TAILOR_ZERO);

        //if (cnt(1) > max_(1))
        //{
            //std::cout << "cnt(0): " << cnt(0) << std::endl;
            //std::cout << "cnt(1): " << cnt(1) << std::endl;
            //std::cout << "cnt(2): " << cnt(2) << std::endl;
            //std::cout << "min_(0): " << min_(0) << std::endl;
            //std::cout << "min_(1): " << min_(1) << std::endl;
            //std::cout << "min_(2): " << min_(2) << std::endl;
            //std::cout << "max_(0): " << max_(0) << std::endl;
            //std::cout << "max_(1): " << max_(1) << std::endl;
            //std::cout << "max_(2): " << max_(2) << std::endl;
        //}

        reset_centroid_ = false;
        centroid_ = cnt;

        return centroid_;
    }

    /*Vector3 Polyhedron::centroid() const
    {
        // based on http://wwwf.imperial.ac.uk/~rn/centroid.pdf

        double cnt_x = 0.;
        double cnt_y = 0.;
        double cnt_z = 0.;

        for (auto f = face_.begin(); f != face_.end(); ++f)
        //for (auto f = begin; f != end; ++f)
        {
            assert(!std::isnan(f->vertex(0).r(0)));
            assert(!std::isnan(f->vertex(0).r(1)));
            assert(!std::isnan(f->vertex(0).r(2)));

            assert(!std::isnan(f->vertex(1).r(0)));
            assert(!std::isnan(f->vertex(1).r(1)));
            assert(!std::isnan(f->vertex(1).r(2)));

            assert(!std::isnan(f->vertex(2).r(0)));
            assert(!std::isnan(f->vertex(2).r(1)));
            assert(!std::isnan(f->vertex(2).r(2)));

            double cx = (f->vertex(0).r(0) + f->vertex(1).r(0) + f->vertex(2).r(0)) / 3.;
            double cy = (f->vertex(0).r(1) + f->vertex(1).r(1) + f->vertex(2).r(1)) / 3.;
            double cz = (f->vertex(0).r(2) + f->vertex(1).r(2) + f->vertex(2).r(2)) / 3.;

            assert(!std::isnan(cx));
            assert(!std::isnan(cy));
            assert(!std::isnan(cz));

            Vector3 n = f->normal();
            assert(!std::isnan(n(0)));
            assert(!std::isnan(n(1)));
            assert(!std::isnan(n(2)));

            cnt_x += (cx * cx) * n(0);
            cnt_y += (cy * cy) * n(1);
            cnt_z += (cz * cz) * n(2);
        }

        Vector3 cnt;
        cnt.set_x(cnt_x);
        cnt.set_y(cnt_y);
        cnt.set_z(cnt_z);

        assert(!std::isnan(cnt_x));
        assert(!std::isnan(cnt_y));
        assert(!std::isnan(cnt_z));

        assert(volume() != 0.);
        assert(!std::isnan(volume()));
        cnt = cnt / (2. * volume());

        if (cnt(0) < min_(0))
        {
            std::cout << "cnt(0): " << cnt(0) << std::endl;
            std::cout << "cnt(1): " << cnt(1) << std::endl;
            std::cout << "cnt(2): " << cnt(2) << std::endl;

            std::cout << "min(0): " << min_(0) << std::endl;
            std::cout << "min(1): " << min_(1) << std::endl;
            std::cout << "min(2): " << min_(2) << std::endl;

            std::cout << "max(0): " << max_(0) << std::endl;
            std::cout << "max(1): " << max_(1) << std::endl;
            std::cout << "max(2): " << max_(2) << std::endl;

            std::cout << "vol: " << volume() << std::endl;

            std::cout << "shape: " << static_cast<int>(shape_) << std::endl;

            for (const Point& p: vertex_)
            {
                std::cout << "r: " << p.r(0) << std::endl;
                std::cout << "r: " << p.r(1) << std::endl;
                std::cout << "r: " << p.r(2) << std::endl;
                std::cout << std::endl;
            }

            for (const Polygon& f: face_)
            {
                for (const Point& p: f.vertex())
                {
                    std::cout << "vertex: " << p.r(0) << " " << p.r(1) << " " << p.r(2) << std::endl;
                }
                std::cout << "norm: " << f.normal()(0) << " " << f.normal()(1) << " " << f.normal()(2) << std::endl;
                std::cout << "signedarea: " << f.signed_area() << std::endl;
            }
        }

        assert(cnt(0) >= min_(0));
        assert(cnt(1) >= min_(1));
        assert(cnt(2) >= min_(2));

        assert(cnt(0) <= max_(0));
        assert(cnt(1) <= max_(1));
        assert(cnt(2) <= max_(2));

        return cnt;
    }*/

    bool Polyhedron::degenerate() const
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

        for (int i=0; i<TAILOR_N_DIM; ++i)
        {
            if (coplanar(i))
            {
                return true;
            }
        }

        return false;
    }

    double Polyhedron::volume() const
    {
        // based on https://cfd.ku.edu/papers/1999-AIAAJ.pdf

        //if (volume_) {
            //return *volume_;
        //}
        //
        if (reset_volume_ == false)
        {
            return volume_;
        }

        double volume;

        if (shape_ == Shape::tet)
        {
            if (degenerate())
            {
                return 0.;
            }
            for (const auto& ff: face_)
            {
                if (ff.vertex().size() < 3)
                {
                    std::cout << "vertex size in hedron volume triface: " << ff.vertex().size() << std::endl;
                }
                assert(ff.vertex().size() >= 3);
            }

            //volume = volume_triface(face_.begin(), face_.end());
            volume = volume_triface(face_.pbegin(), face_.pend());
        }
        else if (shape_ == Shape::hex)
        {
            if (degenerate())
            {
                return 0.;
            }
            // break each face of hex into 2 triangles.
            // add triangles to container.

            std::array<Polygon, 12> faces;

            assert(!vertex_.empty());
            assert(!face_.empty());

            int i = 0;
            for (const Polygon& f: face_)
            {
                faces[i]   = Polygon(std::vector<Point>{f.vertex(0), f.vertex(1), f.vertex(2)});
                assert(faces[i].vertex().size() >= 3);
                faces[i+1] = Polygon(std::vector<Point>{f.vertex(0), f.vertex(2), f.vertex(3)});
                assert(faces[i+1].vertex().size() >= 3);

                i = i + 2;
            }

            for (const auto& ff: faces)
            {
                if (ff.vertex().size() < 3)
                {
                    std::cout << "vertex size in hedron volume triface: " << ff.vertex().size() << std::endl;
                }
                assert(ff.vertex().size() >= 3);
            }

            volume = volume_triface(faces.begin(), faces.end());
        }
        else if (shape_ == Shape::pri)
        {
            if (degenerate())
            {
                return 0.;
            }
            // break each quad face of prism into 2 triangles: 2 tri + 2 * 3 tri = 8 tri
            // add triangles to container.

            std::array<Polygon, 8> faces;

            int i = 0;
            for (const Polygon& f: face_)
            {
                if (f.vertex().size() == 3)
                {
                    faces[i] = f;
                    ++i;
                }
                else
                {
                    faces[i]   = Polygon(std::vector<Point>{f.vertex(0), f.vertex(1), f.vertex(2)});
                    faces[i+1] = Polygon(std::vector<Point>{f.vertex(0), f.vertex(2), f.vertex(3)});
                    i = i + 2;
                }
            }

            for (const auto& ff: faces)
            {
                if (ff.vertex().size() < 3)
                {
                    std::cout << "vertex size in hedron volume triface: " << ff.vertex().size() << std::endl;
                }
                assert(ff.vertex().size() >= 3);
            }

            volume = volume_triface(faces.begin(), faces.end());
        }

        reset_volume_ = false;
        volume_ = std::abs(volume);

        return volume_;
    }

    void Polyhedron::set_faces()
    {
        if (vertex_.empty()) {
            return;
        }

        face_.clear();

        if (shape_ == Shape::quad)
        {
            if (vertex_.size() != 4)
            {
                std::cout << "vertex size: " << vertex_.size() << std::endl;
            }
            face_.push_back(Polygon(std::vector<Point>{vertex_[0], vertex_[1], vertex_[2], vertex_[3]}));
        }
        else if (shape_ == Shape::tri)
        {
            if (vertex_.size() != 3)
            {
                std::cout << "vertex size: " << vertex_.size() << std::endl;
            }
            face_.push_back(Polygon(std::vector<Point>{vertex_[0], vertex_[1], vertex_[2]}));
        }
        else if (shape_ == Shape::tet)
        {
            if (vertex_.size() != 4)
            {
                std::cout << "vertex size: " << vertex_.size() << std::endl;
            }
            assert(vertex_.size() == 4);
            face_.push_back(Polygon(std::vector<Point>{vertex_[1], vertex_[2], vertex_[0]})); // back face
            face_.push_back(Polygon(std::vector<Point>{vertex_[0], vertex_[3], vertex_[1]})); // front left
            face_.push_back(Polygon(std::vector<Point>{vertex_[0], vertex_[2], vertex_[3]})); // front right
            face_.push_back(Polygon(std::vector<Point>{vertex_[2], vertex_[1], vertex_[3]})); // bottom
        }
        else if (shape_ == Shape::pri)
        {
            assert(vertex_.size() == 6);
            //face_.push_back(Polygon(std::vector<Point>{vertex_[1], vertex_[0], vertex_[2]})); // back face
            face_.push_back(Polygon(std::vector<Point>{vertex_[0], vertex_[1], vertex_[2]})); // back face
            face_.push_back(Polygon(std::vector<Point>{vertex_[5], vertex_[4], vertex_[3]})); // front face
            face_.push_back(Polygon(std::vector<Point>{vertex_[1], vertex_[0], vertex_[3], vertex_[4]})); // left face
            face_.push_back(Polygon(std::vector<Point>{vertex_[3], vertex_[0], vertex_[2], vertex_[5]})); // oblique face
            face_.push_back(Polygon(std::vector<Point>{vertex_[2], vertex_[1], vertex_[4], vertex_[5]})); // bottom face
        }
        else if (shape_ == Shape::hex)
        {
            assert(vertex_.size() == 8);
            face_.push_back(Polygon(std::vector<Point>{vertex_[0], vertex_[1], vertex_[2], vertex_[3]})); // back face
            face_.push_back(Polygon(std::vector<Point>{vertex_[6], vertex_[5], vertex_[4], vertex_[7]})); // front face
            face_.push_back(Polygon(std::vector<Point>{vertex_[2], vertex_[1], vertex_[5], vertex_[6]})); // left face
            face_.push_back(Polygon(std::vector<Point>{vertex_[4], vertex_[0], vertex_[3], vertex_[7]})); // right face
            face_.push_back(Polygon(std::vector<Point>{vertex_[1], vertex_[0], vertex_[4], vertex_[5]})); // top face
            face_.push_back(Polygon(std::vector<Point>{vertex_[3], vertex_[2], vertex_[6], vertex_[7]})); // bottom face
        }
        else 
        {
            assert(false);
        }
        assert(!face_.empty());
    }

    /*void erase_rubbing_faces(std::vector<Polygon>& container, int ignoredim)
    {
        std::vector<bool> to_be_erased(container.size(), false);

        for (auto face = container.begin(); face != container.end();)
        {
            assert(face.vertex().size() == 4);
            assert(std::abs(face.vertex(0).r(ignoredim) - face.vertex(1).r(ignoredim)) < TAILOR_ZERO);
            assert(std::abs(face.vertex(0).r(ignoredim) - face.vertex(2).r(ignoredim)) < TAILOR_ZERO);
            assert(std::abs(face.vertex(0).r(ignoredim) - face.vertex(3).r(ignoredim)) < TAILOR_ZERO);

            for (auto comp = std::next(face); comp != container.end();)
            {
                // ignore non-in-plane faces.
                if (std::abs(face.vertex(0).r(ignoredim) - comp.vertex(0).r(ignoredim)) > TAILOR_ZERO) {
                    continue;
                }

                bool inter1 = face.edge(0).do_intersect(comp.edge(0), false, ignoredim);
                if (!inter1) {
                    continue;
                }
                bool inter2 = face.edge(3).do_intersect(comp.edge(3), false, ignoredim);
                if (!inter2) {
                    continue;
                }

                if (std::abs(face.signed_area()) < std::abs(comp.signed_area()))
                {
                    to_be_erased[face - face.begin()] = true;
                }
                else {
                    to_be_erased[comp - face.begin()] = true;
                }
                
                // subtract smaller faces from bigger faces.
            }
        }

        container.erase(std::remove_if(container.begin(), container.end(), [&](const auto& c){return to_be_erased[&c - container.data()] == true;}), container.end());
    }*/

    /*void ignore_contained_aabbs(std::vector<AABB>& aabb)
    {
        for (auto a = aabb.begin(); a != aabb.end();)
        {
            if (to_be_ignored[a - aabb.begin()]) {
                continue;
            }

            for (auto b = std::next(a); b != aabb.end();)
            {
                bool breakall = false;
                if (a.do_contain(b))
                {
                    to_be_ignored[b - aabb.begin()] = true;
                    break;
                }
                if (b.do_contain(a))
                {
                    to_be_ignored[a - aabb.begin()] = true;
                    break;
                }
            }
        }

        aabb.erase(std::remove_if(aabb.begin(), aabb.end(), [&](const auto& c){return to_be_ignored[&c - aabb.data()] == true;}), aabb.end());
    }*/

    /*void make_outline(std::vector<AABB>& aabb)
    {
        ignore_contained_aabbs(aabb);

        for (const auto& b: aabb)
        {
            xz.push_back(b.face[0]);
            xz.push_back(b.face[1]);
            yz.push_back(b.face[2]);
            yz.push_back(b.face[3]);
            xy.push_back(b.face[4]);
            xy.push_back(b.face[5]);
        }

        erase_rubbing_faces_(xz, 1);
        erase_rubbing_faces_(yz, 0);
        erase_rubbing_faces_(xy, 2);
    }*/

    bool Polyhedron::do_intersect(const Vector3& r, bool verbose) const
    {
        for (int i=0; i<TAILOR_N_DIM; ++i)
        {
            if (min_(i) - r(i) > TAILOR_ZERO) {
                return false;
            }

            if (r(i) - max_(i) > TAILOR_ZERO) {
                return false;
            }

            //if (r(i) < min_(i) || r(i) > max_(i)) {
                //return false;
            //}
        }

        //double off = std::pow(volume(), 1/3.);
        auto a = centroid();
        if (verbose)
        {
            std::cout << "reset_centroid_: " << reset_centroid_ << std::endl;
            std::cout << "a: " << a(0) << " " << a(1) << " " << a(2) << std::endl;
            //std::cout << "centroid: " << centroid_(0) << " " << centroid_(1) << " " << centroid_(2) << std::endl;
            std::cout << "min: " << min_(0) << " " << min_(1) << " " << min_(2) << std::endl;
            std::cout << "max: " << max_(0) << " " << max_(1) << " " << max_(2) << std::endl;
        }
        a(0) = a(0) - (max_(0) - min_(0));

        //Point p_(min_(0) - off, min_(1) - off/2., min_(2) - off/3.);
        Point p_(a);
        Point p(r);
        Segment s(p_, p); 

        int counter = 0;
        assert(!face_.empty());
        for (const Polygon& f: face_)
        {
            if (verbose) {
                std::cout << "checking face" << std::endl;
            }

            //const auto& n = f.normal();
            //auto cnt = f.centroid();
            //Point p_(cnt(0) + off * n(0), cnt(1) + off * n(1), cnt(2) + off * n(2));
            //Segment s(p_, p); 

            if (verbose)
            {
                std::cout << "seg vertex: " << s.vertex(0).r(0) << " " << s.vertex(0).r(1) << " " << s.vertex(0).r(2) << std::endl;
                std::cout << "seg vertex: " << s.vertex(1).r(0) << " " << s.vertex(1).r(1) << " " << s.vertex(1).r(2) << std::endl;
                //std::cout << "normal: " << n(0) << " " << n(1) << " " << n(2) << std::endl;
                //std::cout << "off: " << off << std::endl;
            }

            /*bool point_out_aabb = false;
            for (int i=0; i<TAILOR_N_DIM; ++i)
            {
                if ((r(i) < (f.min()(i) - TAILOR_ZERO)) || (r(i) > (f.max()(i) + TAILOR_ZERO))) {
                    point_out_aabb = true;
                    if (verbose) {
                        std::cout << "point_out_aabb: " << point_out_aabb  << std::endl;
                    }
                    break;
                }
            }
            
            if (point_out_aabb) {
                continue;
            }*/

            if (verbose)
            {
                for (const Point& pp: f.vertex())
                {
                    std::cout << "face vertex: " << pp.r(0) << " " << pp.r(1) << " " << pp.r(2) << std::endl;
                }
            }

            bool just_on_face;
            if (f.do_intersect(s, just_on_face, verbose))
            {
                if (just_on_face)
                {
                    if (verbose)
                    {
                        std::cout << "just on face. so the point is inside hedron." << std::endl;
                    }
                    return true;
                }
                if (verbose)
                {
                    std::cout << "yesinter" << std::endl;
                }
                ++counter;
            }
            else
            {
                if (verbose)
                {
                    std::cout << "nointer" << std::endl;
                }
            }
        }

        if (verbose)
        {
            std::cout << "ninter: " << counter << std::endl;
        }

        if (counter == 0 || counter % 2 == 0) {
            return false;
        }

        return true;
    }

    bool Polyhedron::do_intersect(const Segment& s) const
    {
        for (const Polygon& f: face_)
        {
            bool just_on_face;
            if (f.do_intersect(s, just_on_face))
            {
                return true;
            }
        }

        for (int i=0; i<TAILOR_N_DIM; ++i)
        {
            if (s.min(i) < min_(i)) {
                return false;
            }

            if (s.max(i) > max_(i)) {
                return false;
            }
        }        

        return true;
    }
}
