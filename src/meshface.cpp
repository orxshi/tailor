#include <meshface.h>

namespace Tailor
{
    MeshFace::MeshFace(): R_checked_(false), btype_(BouType::undefined), vgn_(0.)
    {
    }

    Vector3 tangent_vector(const Vector3& n)
    {
        Vector3 t = n;

        for (int i=0; i<TAILOR_N_DIM; ++i)
        {
            if (std::abs(n(i)) > TAILOR_ZERO)
            {
                int j = 0;
                for (; j<TAILOR_N_DIM; ++j)
                {
                    if (j != i) {
                        break;
                    }
                }

                int k = 0;
                for (; k<TAILOR_N_DIM; ++k)
                {
                    if (k != i && k != j) {
                        break;
                    }
                }

                t(j) = 1.;
                t(k) = 1.;
                t(i) = (-t(j) * n(j) - t(k) * n(k)) / n(i);

                t = normalize(t);
                return t;
            }
        }

        return t;
    }
    bool MeshFace::R_checked() const
    {
        return R_checked_;
    }

    void MeshFace::set_R_checked(bool b)
    {
        R_checked_ = b;
    }

    void MeshFace::reset_btype()
    {
        if (!is_boundary())
        {
            //if (btype_ != BouType::interior) {
                btype_ = BouType::undefined;
            //}
        }
    }
    const Tag& MeshFace::left_cell() const
    {
        if (parent_cell().size() != 2)
        {
            std::cout << "btype: " << static_cast<int>(btype_) << std::endl;
        }
        assert(parent_cell().size() == 2);
        return parent_cell(1);
    }

    const Tag& MeshFace::right_cell() const
    {
        assert(!parent_cell().empty());
        assert(parent_cell().size() <= 2);
        return parent_cell(0);
    }

    const Matrix5& MeshFace::M() const
    {
        return M_;
    }

    //const Matrix<NVAR, NVAR>& MeshFace::ML() const
    //{
    //    return ML_;
    //}

    //const Matrix<NVAR, NVAR>& MeshFace::MR() const
    //{
    //    return MR_;
    //}

    double sum_face_velocities(const Polygon& polygon, double dt, const Vector3& vel)
    {
        assert(!(vel(0) == 0. && vel(1) == 0.));
        double sum = 0.;

        Vector3 d = vel * dt;
        for (const Segment& seg: polygon.edge())
        {
            auto normal = seg.normal();
            Polygon swept = create_with_sweep(seg, d);
            auto vb = normal * (swept.signed_area() / dt / seg.len());
            double vbn = vb(0)*normal(0) + vb(1)*normal(1);
            sum += vbn;
        }

        return sum;
    }

    double sum_face_velocities2(const Polyhedron& hedron, double dt, const Vector3& vel)
    {
        assert(!(vel(0) == 0. && vel(1) == 0. && vel(2) == 0.));
        double sum = 0.;
        //Vector3 sum;
        //sum.set_x(0.);
        //sum.set_y(0.);
        //sum.set_z(0.);

        std::cout << "dt: " << dt << std::endl;

        Vector3 d = vel * dt;
        for (const Polygon& face: hedron.faces())
        {
            Polyhedron swept = create_with_sweep(face, d);
            //auto vgn = face.normal() * (swept.volume() / dt) * face.signed_area();
            //auto vgn = face.normal() * swept.volume() / face.signed_area() / dt;
            //auto vgn = swept.volume() * face.signed_area() / dt;
            auto vgn = swept.volume() / dt;
            sum = sum + vgn;

            std::cout << "swept volume: " << swept.volume() << std::endl;
            std::cout << "face area: " << face.signed_area() << std::endl;
            std::cout << "normal(0): " << face.normal(0) << std::endl;
            std::cout << "normal(1): " << face.normal(1) << std::endl;
            std::cout << "normal(2): " << face.normal(2) << std::endl;
            std::cout << "vgn: " << vgn << std::endl;
            //std::cout << "vgn(0): " << vgn(0) << std::endl;
            //std::cout << "vgn(1): " << vgn(1) << std::endl;
            //std::cout << "vgn(2): " << vgn(2) << std::endl;
        }


        return sum;
    }

    double sum_face_velocities(const Polyhedron& hedron, double dt, const Vector3& vel)
    {
        assert(!(vel(0) == 0. && vel(1) == 0. && vel(2) == 0.));
        double sum = 0.;
        //Vector3 sum;
        //sum.set_x(0.);
        //sum.set_y(0.);
        //sum.set_z(0.);

        std::cout << "dt: " << dt << std::endl;

        Vector3 d = vel * dt;
        for (const Polygon& face: hedron.faces())
        {
            Polyhedron swept = create_with_sweep(face, d);
            //auto vgn = face.normal() * (swept.volume() / dt) * face.signed_area();
            //auto vgn = face.normal() * swept.volume() / dt;
            //auto vgn = swept.volume() / dt;
            double vgn;
            if (dot(vel, face.normal()) >= 0.)
            {
                vgn = swept.volume() / dt;
            }
            else
            {
                vgn = swept.volume() * -1. / dt;
            }

            sum = sum + vgn;

            std::cout << "swept volume: " << swept.volume() << std::endl;
            //std::cout << "face area: " << face.signed_area() << std::endl;
            std::cout << "normal(0): " << face.normal(0) << std::endl;
            std::cout << "normal(1): " << face.normal(1) << std::endl;
            std::cout << "normal(2): " << face.normal(2) << std::endl;
            std::cout << "vgn: " << vgn << std::endl;
            //std::cout << "vgn(0): " << vgn(0) << std::endl;
            //std::cout << "vgn(1): " << vgn(1) << std::endl;
            //std::cout << "vgn(2): " << vgn(2) << std::endl;
        }

        //std::cout << "sum(0): " << sum(0) << std::endl;
        //std::cout << "sum(1): " << sum(1) << std::endl;
        //std::cout << "sum(2): " << sum(2) << std::endl;

        return sum;
    }

    //Vector3 MeshFace::vgn() const
    //{
        //return vgn_;
    //}

    double MeshFace::vgn() const
    {
        return vgn_;
    }
    const Vector3& MeshFace::vf() const
    {
        return vf_;
    }

    void MeshFace::face_velocity(const Freestream& fs, const Component& compo, double real_time)
    {
        // https://www.lehman.edu/faculty/anchordoqui/chapter06.pdf




        /*if (vel(0) == 0. && vel(1) == 0. && vel(2) == 0.)
        {
            vgn_ = 0.;
            return;
        }

        Polyhedron swept_hedron = create_with_sweep(face(), vel * dt);
        vgn_ = swept_hedron.volume() / face().signed_area() / dt;

        assert(face().signed_area() != 0.);
        assert(dt != 0.);

        if (dotp(vel, face().normal()) < 0.)
        {
            vgn_ *= -1.;
        }

        assert(!std::isnan(vgn_));*/

        //double pinf = fs.pinf_;
        //double rhoinf = fs.rhoinf_;
        //double machfoil = fs.machfoil_;
        //double aoa_foil_x = fs.aoa_foil_x_;
        //double aoa_foil_z = fs.aoa_foil_z_;
        double rotation = compo.rotation_;
        double oscillation = compo.oscillation;
        double rotaxis = compo.rotaxis_;
        double rpm = compo.rpm_;
        auto pivot = compo.pivot_;
        //double mach = compo.mach_;
        double dirx = compo.dirx_;
        double dirz = compo.dirz_;

        //double cinf = std::sqrt(fs.gamma_ * pinf / rhoinf);

        Vector3 vel(0., 0., 0.);
        if (rotation)
        {
            if (!oscillation)
            {
                // rpm_foil is positive if ccw.
                double om = rpm * 2. * PI / 60.; // rad/s
                Vector3 omega(0., 0., om);

                if (rotaxis == 0)
                {
                    assert(false);
                }
                else if (rotaxis == 1)
                {
                    assert(false);
                }
                else if (rotaxis == 2)
                {
                    auto cnt = face_.centroid();
                    auto r = cnt - pivot;
                    r(2) = 0.;

                    vel = cross(omega, r);
                    if (vel(2) != 0.)
                    {
                        std::cout << "vel(2): " << vel(2) << std::endl;
                    }
                    assert(vel(2) == 0.);

                    //std::cout << cnt(0) << " " << cnt(1) << " " << cnt(2) << " " << vel(0) << " " << vel(1) << " " << vel(2) << std::endl;

                    //auto rvec = Vector3(cnt(0), cnt(1), pivot(2));
                    //r = (rvec - pivot).len();
                    //alpha = std::atan2(cnt(0), cnt(1));
                    //thetax = -std::sin(alpha);
                    //thetay = std::cos(alpha);
                    //vel = Vector3(
                    //        omega * r * thetax,
                    //        omega * r * thetay,
                    //        0.);
                }
            }
            else
            {
                double dirz = compo.dirz_;
                double reduced_freq = compo.reduced_freq;

                double u = compo.u;
                double v = compo.v;
                double w = compo.w;

                Freestream fs;
                fs.read();

                //Vector3 U(u, v, w);
                //double u_ref = U.len();

                double u_ref = fs.u_;

                double chord = compo.chord;
                double aoa_mean_deg = compo.aoa_mean_deg;
                double aoa_o_deg = compo.aoa_o_deg;
                double aoa_mean = deg_to_rad(aoa_mean_deg);
                double aoa_o = deg_to_rad(aoa_o_deg);

                double omega = 2. * reduced_freq * u_ref / chord; // angular frequency.
                
                double rad_vel_z = aoa_o * omega * std::cos(omega * real_time); // z-component of angular velocity.

                Vector3 rad_vel(0., 0., -rad_vel_z);

                Vector3 pivot = compo.pivot_;
                auto rotaxis = compo.rotaxis_;

                assert(rotaxis == 2);

                if (rotaxis == 2)
                {
                    auto cnt = face_.centroid();
                    auto r = cnt - pivot;
                    r(2) = 0.;
                    
                    //vel = cross(rad_vel, project);
                    if (btype_ == BouType::wall) // TODO just to test.
                    //if (btype_ == BouType::wall || btype_ == BouType::interior)
                    //if (btype_ == BouType::wall || btype_ == BouType::interior && btype_ == BouType::farfield && btype_ == BouType::empty) // TODO just to test.
                    //if (btype_ != BouType::partition) // TODO just to test.
                    {
                        vel = cross(rad_vel, r);
                    }

                    //std::cout << omega << " " << real_time << " " << omega * real_time << " " << std::sin(omega * real_time) << " " << std::cos(omega * real_time) << " " << rad_vel(0) << " " << rad_vel(1) << " " << rad_vel(2) << " " << vel(0) << " " << vel(1) << " " << vel(2) << " " << r(0) << " " << r(1) << " " << r(2) << std::endl;
                    //std::cout << "omega: " << omega << std::endl;
                    //std::cout << "real_time: " << real_time << std::endl;
                    //std::cout << "rad_vel_z: " << rad_vel_z << std::endl;
                    //std::cout << "vel(0): " << vel(0) << std::endl;
                    //std::cout << "vel(1): " << vel(1) << std::endl;
                    //std::cout << "vel(2): " << vel(2) << std::endl;
                    //std::cout << "r(0): " << r(0) << std::endl;
                    //std::cout << "r(1): " << r(1) << std::endl;
                    //std::cout << "r(2): " << r(2) << std::endl;
                    //std::cout << "u: " << u << std::endl;
                    //std::cout << "v: " << v << std::endl;
                    //std::cout << "w: " << w << std::endl;
                    //assert(false);
                }
            }
        }
        //else
        {
            //vel = Vector3(0., 0., 0.);
            //vel = Vector3(
                    //machfoil * cinf * std::cos(deg_to_rad(aoa_foil_x)),
                    //machfoil * cinf * std::cos(deg_to_rad(90. - aoa_foil_x)),
                    //machfoil * cinf * std::cos(deg_to_rad(aoa_foil_z)));

            //if (compo.mach_ != 0.)
            //{
            //    vel = Vector3(
            //            mach * cinf * std::cos(deg_to_rad(dirx)),
            //            mach * cinf * std::cos(deg_to_rad(90. - dirx)),
            //            mach * cinf * std::cos(deg_to_rad(dirz))
            //            );
            //}
            //else
            {
                vel += Vector3(
                        compo.u,
                        compo.v,
                        compo.w
                        );

                //if (btype_ == BouType::farfield)
                //{
                    //std::cout << "vel(0): " << vel(0) << std::endl;
                    //std::cout << "vel(1): " << vel(1) << std::endl;
                    //std::cout << "vel(2): " << vel(2) << std::endl;
                //}
            }
        }

        //https://core.ac.uk/download/pdf/81977729.pdf

        auto n = face_.normal();

        /*
        if (vel.len() > TAILOR_ZERO)
        {
            Vector3 d = vel * dt;
            double swvol = 0.;
            if (std::abs(dotp(n, d)) > TAILOR_ZERO)
            {
                Polyhedron swept = create_with_sweep(face_, d);
                //std::cout << "dt: " << dt << std::endl;
                //std::cout << "vel: " << vel(0) << " " << vel(1) << " " << vel(2) << std::endl;
                //std::cout << "d: " << d(0) << " " << d(1) << " " << d(2) << std::endl;
                //for (const auto& vtx: swept.vertices())
                //{
                    //std::cout << vtx.r(0) << " " << vtx.r(1) << " " << vtx.r(2) << std::endl;
                //}
                swvol = swept.volume();
            }
            //vgn_ = normalize(vel) * swept.volume() / std::abs(face_.signed_area()) / dt;
            //vgn_ = n * swept.volume() / std::abs(face_.signed_area()) / dt;
            //vgn_ = swept.volume() / std::abs(face_.signed_area()) / dt;
            vgn_ = swvol / std::abs(face_.signed_area()) / dt;
            if (dotp(n, vel) < 0.)
            {
                vgn_ = vgn_ * -1;
            }
        }*/

        vgn_ = dot(vel, n);
        vf_ = vel;
        //if (vf_(2) != 0.)
        //{
            //std::cout << "vf_(2): " << vf_(2) << std::endl;
        //}
        //assert(vf_(2) == 0.);

        //vgn_ = face_.normal() * dotp(vgn_, face_.normal());

        //vgn_ = dotp(vel, face().normal());
        //vgn_ = vel;
        //assert(!std::isnan(vgn_));
        //if (rpm != 0.)
        //{
        //std::cout << "cnt: " << face_.centroid()(0) << " " << face_.centroid()(1) << " " << face_.centroid()(2) << std::endl;
        //std::cout << "pivot: " << pivot(0) << " " << pivot(1) << " " << pivot(2) << std::endl;
        //std::cout << "rpm: " << rpm << std::endl;
        //std::cout << "alpha: " << alpha << std::endl;
        //std::cout << "rotaxis: " << rotaxis << std::endl;
        //assert(false);
        //}
        //if (std::isnan(vgn_))
        //{
        //    std::cout << "vel(0): " << vel(0) << std::endl;
        //    std::cout << "vel(1): " << vel(1) << std::endl;
        //    std::cout << "vel(2): " << vel(2) << std::endl;
        //    std::cout << "omega: " << omega << std::endl;
        //    std::cout << "r: " << r << std::endl;
        //    std::cout << "thetax" << thetax << std::endl;
        //    std::cout << "thetay" << thetay << std::endl;
        //    std::cout << "alpha: " << alpha << std::endl;
        //    std::cout << "rotation: " << rotation << std::endl;
        //    std::cout << "rpm: " << rpm << std::endl;
        //}
        //assert(!std::isnan(vgn_));

        //std::cout << "vol: " << swept_hedron.volume() << std::endl;
        //for (const Point& vtx: face_.vertex())
        //{
        //std::cout << "p: " << vtx.r(0) << std::endl;
        //std::cout << "p: " << vtx.r(1) << std::endl;
        //std::cout << "p: " << vtx.r(2) << std::endl;
        //}
    }

    //const Tag& MeshFace::tag() const
    //{
        //return tag_;
    //}

    //void MeshFace::set_as_boundary()
    //{
        //is_boundary_ = true;
    //}

    void MeshFace::set_btype(BouType type)
    {
        btype_ = type;
    }
    
    bool MeshFace::is_boundary() const
    {
        if (btype_ == BouType::wall || btype_ == BouType::dirichlet || btype_ == BouType::farfield || btype_ == BouType::empty || btype_ == BouType::interog || btype_ == BouType::symmetry) {
            return true;
        }

        return false;
    }

    const Polygon& MeshFace::face() const
    {
        return face_;
    }

    void MeshFace::move(const Vector3& v)
    {
        //Point p0, p1;
        //p0.set_r(segment_.vertex(0).r(0) + v(0), segment_.vertex(0).r(1) + v(1));
        //p1.set_r(segment_.vertex(1).r(0) + v(0), segment_.vertex(1).r(1) + v(1));
        //segment_ = Segment(p0, p1);
        face_.move_points(v);
    }

    void MeshFace::rotate(double angle, int axis, const Vector3& rot_point)
    {
        face_.rotate_points(angle, axis, rot_point);
    }

    const Tag& MeshFace::mesh_point(int i) const
    {
        assert(i >= 0);
        assert(i < mesh_point_.size());

        return mesh_point_[i];
    }

    const std::vector<Tag>& MeshFace::mesh_point() const
    {
        return mesh_point_;
    }

    //void MeshFace::set_tag(const Tag& t)
    //void MeshFace::set_tag(int mc, int nc, BouType btype)
    //{
        //tag_ = t;
        //tag_.set(mc, nc, btype);
    //}

    void MeshFace::set_tag(const FaceTag& ftag)
    {
        tag_ = ftag;
    }

    BouType MeshFace::btype() const
    {
        return btype_;
    }

    //void MeshFace::set_faceaddr(MeshFace* addr)
    //{
    //    if (addr != nullptr)
    //    {
    //        assert(tag_ == addr->tag());
    //    } 
    //    faceaddr_ = addr;
    //}

    //MeshFace* MeshFace::faceaddr()
    //{
    //    return faceaddr_;
    //}

    //const MeshFace* MeshFace::faceaddr() const
    //{
    //    return faceaddr_;
    //}

    MeshFace::MeshFace(std::vector<MeshPoint>& pts): btype_(BouType::undefined)
    {
        vgn_ = 0.;
        vf_ = 0.;

        assert(!pts.empty());
        std::vector<Point> rawpts;
        rawpts.reserve(pts.size());

        for (const MeshPoint& mp: pts)
        {
            mesh_point_.push_back(mp.tag());
            rawpts.push_back(mp.p());
        }

        assert(!rawpts.empty());
        for (const auto& a: rawpts)
        {
            assert(!std::isnan(a.r(0)));
            assert(!std::isnan(a.r(1)));
            assert(!std::isnan(a.r(2)));
        }
        face_ = Polygon(rawpts);
        assert(!face_.vertex().empty());
    }

    void MeshFace::remove_parent_cell(const Tag& ic)
    {
        assert(ic.isvalid());
        parent_cell_.erase(std::remove_if(parent_cell_.begin(), parent_cell_.end(), [&](const auto& a){return a == ic;}), parent_cell_.end());
    }

    /*const Tag& MeshFace::parent_mesh() const
    {
        return parent_mesh_;
    }*/

    void MeshFace::remove_parent_cells()
    {
        parent_cell_.clear();
    }

    /*void MeshFace::set_tag(const Tag& t)
    {
        tag_ = t;
    }
    void MeshFace::set_parent_mesh(const Tag& ptag)
    {
        parent_mesh_ = ptag;
    }*/

    const std::vector<Tag>& MeshFace::parent_cell() const
    {
        return parent_cell_;
    }

    const Tag& MeshFace::parent_cell(int i) const
    {
        return parent_cell_[i];
    }

    const FaceTag& MeshFace::tag() const
    {
        return tag_;
    }

    void MeshFace::add_parent_cell(const Tag& celltag)
    {
        assert(celltag.isvalid());
        parent_cell_.push_back(celltag);

        //auto it = std::unique(parent_cell_.begin(), parent_cell_.end());
        //if (it != parent_cell_.end())
        //{
        //    for (const auto& pc: parent_cell_)
        //    {
        //        std::cout << "pc: " << pc.tag()() << std::endl;
        //    }
        //}
        //assert(it == parent_cell_.end());
    }
}
