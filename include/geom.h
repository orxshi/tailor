#ifndef GEOMETRIC_SHAPE_H
#define	GEOMETRIC_SHAPE_H

#include <vector>
#include <algorithm>
#include "core.h"
//#include "vec3.h"
#include "array.h"
#include "tag.h"
#include "rotationmatrix.h"

//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
//#include <CGAL/Polygon_2.h>

//typedef CGAL::Exact_predicates_exact_constructions_kernel CGAL_kernel;
//typedef CGAL::Polygon_2<CGAL_kernel> Polygon_2;
//typedef CGAL_kernel::Point_2 CGAL_Point;

namespace Tailor
{
    class Point
    {
        Vector3 r_;
        friend class boost::serialization::access;

        public:

        Point(): Point(0., 0., 0.) {};
        Point(double x, double y, double z);
        Point(const Vector3& r);
        Point(const Point& p);

        size_t mem() const;
        const Vector3& r() const;
        double r(int i) const;
        void set_r(double x, double y, double z);
        void set_r(const Vector3& _r);
        void set_r(int i, double v);
        bool operator==(const Point& other) const;
        bool operator<(const Point& other) const;
        Point& operator=(const Point& p);
        template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            ar & r_;
        }
    };

    //using Terminal = Array<Point, 2>; 
    using Terminal = std::vector<Point>;

    class Segment
    {
        public:

        Segment() = default;
        Segment(const Point&, const Point&);
        Segment(const Segment& other);

        Vector3 normal() const;
        //Vector3 centroid() const;
        const Terminal& vertex() const;
        const Point& vertex(int i) const;
        bool do_intersect(const Point&, int ignoredim) const;
        bool do_intersect(const Segment& s, bool point_on_segment_means_inter, int ignoredim, bool verbose=false) const;
        void rotate_points(double ang, const Vector3& axis, const Vector3& rot_axis);
        const Vector3& min() const;
        const Vector3& max() const;
        double min(int i) const;
        double max(int i) const;
        void move_points(const Vector3& v);
        double len() const;
        //Vector3 area() const;
        template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            ar & vertex_;
            ar & min_;
            ar & max_;
        }

        private:
        
        Terminal vertex_;
        Vector3 min_;
        Vector3 max_;

        void set_bbox();

        friend class boost::serialization::access;
    };

    class Polygon
    {
        using Vertices = Array<Point, 4>; 
        using Edges = Array<Segment, 4>; 

        public:

        Polygon();
        Polygon(const std::vector<Point>& vertex);

        double cross_len() const;
        const Vector3& min() const;
        const Vector3& max() const;
        bool degenerate() const;
        void set_vertices(const std::vector<Point>& vertex);
        void set_edge();
        void rotate_points(double ang, const Vector3& axis, const Vector3& rot_axis);
        void move_points(const Vector3& v);
        const Vertices& vertex() const;
        const Point& vertex(int i) const;
        const Edges& edge() const;
        const Segment& edge(int i) const;
        double signed_area() const;
        Vector3 centroid() const;
        bool do_intersect(const Segment& s, bool& just_on_face, bool verbose=false) const;
        bool do_intersect(const Vector3& r, int ignoredim, bool verbose=false) const;
        bool do_plane_intersect(const Segment& s, Vector3& interp, bool verbose) const;
        //bool do_intersect(const Point&) const;
        //bool do_intersect(const Vector3& r) const;
        //bool do_intersect(const Polygon&) const;
        double longest_edge_length() const;
        Vector3 normal() const;
        double normal(int i) const;
        template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            ar & vertex_;
            ar & edge_;
            ar & min_;
            ar & max_;
            ar & reset_signed_area_;
            ar & reset_normal_;
            ar & reset_centroid_;
            ar & signed_area_;
            ar & normal_;
            ar & centroid_;
            //ar & nvertex_;

            assert(!std::isnan(signed_area_));
        }

        private:

        Vector3 min_;
        Vector3 max_;
        Vertices vertex_;
        Edges edge_;
        mutable bool reset_centroid_;
        mutable bool reset_signed_area_;
        mutable bool reset_normal_;
        mutable double signed_area_;
        mutable Vector3 centroid_;
        mutable Vector3 normal_;
        //mutable boost::optional<double> signed_area_;
        //int nvertex_;
        //mutable boost::optional<Vector3> normal_; // don't use this directly even for internal usage. use normal().
        void set_bbox();
        std::vector<double> centroid_ab(int i0, int i1) const;

        friend class boost::serialization::access;
    };

    /*class Face
    {
        using Vertices = Array<Point, 8>; 

        public:

        void set_vertices(const std::vector<Point>& vertex);

        private:

        Vertices vertex_;
    };*/

    enum Shape
    {
        tet,
        hex,
        pri,
        quad,
        tri,
        undef,
    };

    class Polyhedron
    {
        using Vertices = Array<Point, 8>; 
        using Faces = Array<Polygon, 6>; 

        public:

        Polyhedron();
        Polyhedron(const std::vector<Point>& vertex, Shape shape);

        bool face_cross_len() const;
        bool degenerate() const;
        Shape shape() const;
        void rotate_points(double ang, const Vector3& axis, const Vector3& rot_axis);
        void move_points(const Vector3& v);
        const Vertices& vertices() const;
        const Faces& faces() const;
        Vector3 centroid() const;
        bool do_intersect(const Vector3& r, bool verbose=false) const;
        bool do_intersect(const Segment&) const;
        double volume() const;
        const Vector3& min() const;
        const Vector3& max() const;
        double min(int i) const;
        double max(int i) const;
        const Polygon& face(int i) const;
        void set_min(Vector3 min);
        void set_max(Vector3 max);
        void set_bbox(Vector3 min, Vector3 max);
        void set_bbox(double xmin, double ymin, double zmin, double xmax, double ymax, double zmax);
        void set_shape(Shape shape);
        void set_faces();
        void set_vertices_from_bbox();
        size_t mem() const;

        template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            ar & vertex_;
            ar & face_;
            ar & min_;
            ar & max_;
            ar & shape_;
            ar & reset_centroid_;
            ar & reset_volume_;
            ar & centroid_;
            ar & volume_;

            assert(!std::isnan(volume_));
        }

        private:

        mutable bool reset_centroid_;
        mutable bool reset_volume_;
        mutable Vector3 centroid_;
        mutable double volume_;
        Vertices vertex_;
        Faces face_;
        //mutable boost::optional<double> volume_;
        Shape shape_;
        //int nvertex_;

        protected:

        Vector3 min_; 
        Vector3 max_; 

        //double volume_triface(std::array<Polygon>::const_iterator begin, std::array<Polygon>::const_iterator end) const;
        double volume_triface(const Polygon* begin, const Polygon* end) const;
        void set_bbox();
        void set_vertices(const std::vector<Point>& vertex);
        //Vector3 centroid(std::array<Polygon>::const_iterator begin, std::array<Polygon>::const_iterator end) const;

        friend class boost::serialization::access;
    };

    class RegularHexahedron: public Polyhedron
    {
        // Opposing faces are of equal size.
        // Dimensions can be different.

        public:

        RegularHexahedron();
        RegularHexahedron(Vector3 min, Vector3 max);
        RegularHexahedron(double minx, double miny, double maxx, double maxy, double zmin, double zmax);
        RegularHexahedron(const std::vector<Point>& points);
        RegularHexahedron(const Polyhedron&);

        double volume_of_overlap(const RegularHexahedron& a, RegularHexahedron& aabb);
        bool extend(const RegularHexahedron& other);
        void extend(double aabb_margin);
        using Polyhedron::do_intersect;
        bool do_intersect(const RegularHexahedron&) const;
        bool do_contain(const RegularHexahedron& other, bool strict) const;

        template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            ar & boost::serialization::base_object<Polyhedron>(*this);
        }

        private:

        void set_bbox_from_points(const std::vector<Point>& points);
        friend class boost::serialization::access;
    };

    using AABB = RegularHexahedron;

    class Outline
    {
        public:

            /*void build(const std::vector<std::vector<Vector3>>& point, int dummyrank);
            //bool contain(Vector3 p, bool strict) const;
            const Tag& tag() const;
            void set_tag(const Tag& t);
            const std::vector<Segment>& segment() const;
            //bool do_intersect(const std::vector<Vector3>& vertices) const;
            bool do_contain(const Vector3& v, bool strict) const;
            //Polygon_2 polygon() const;
            void build_polygon(int dummyrank);*/

            template<class Archive> void serialize(Archive & ar, const unsigned int version)
            {
                ar & tag_;
                ar & edge_;
                //ar & polygon_;
            }
            
        private:
            friend class boost::serialization::access;
            Tag tag_;
            std::vector<Segment> edge_;

    };

    int orientation (const Vector3& p, const Vector3& q, const Vector3& r, int ignoredim, bool verbose=false);
    Polyhedron create_with_sweep(const Polygon& polygon, const Vector3& v);
    Polygon create_with_sweep(const Segment& segment, const Vector3& v);
    Vector3 normal_plane(const Vector3& v0, const Vector3& v1, const Vector3& v2);
    bool do_intersect_plane_segment(const Vector3& a, const Vector3& b, const Vector3& c, const Segment& s, Vector3& interp);
    Polygon make_outline(const std::vector<Polygon>& pgon);
}

#endif
