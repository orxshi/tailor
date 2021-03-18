#include <meshcell.h>

namespace Tailor
{
    Vector3 MeshCell::vgn() const
    {
        return vgn_;
    }

    void MeshCell::mesh_velocity(double dt, const Freestream& fs, const Component& compo)
    {
        // https://www.lehman.edu/faculty/anchordoqui/chapter06.pdf

        double pinf = fs.pinf_;
        double rhoinf = fs.rhoinf_;
        double rotation = compo.rotation_;
        double rotaxis = compo.rotaxis_;
        double rpm = compo.rpm_;
        auto pivot = compo.pivot_;
        double mach = compo.mach_;
        double dirx = compo.dirx_;
        double dirz = compo.dirz_;

        double cinf = std::sqrt(fs.gamma_ * pinf / rhoinf);

        Vector3 vel;
        double omega, alpha, thetax, thetay, r;
        if (rotation)
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
                auto cnt = poly_.centroid();
                auto r = cnt - pivot;

                vel = cross(omega, r);
            }
        }
        else
        {
            vel = Vector3(
                    mach * cinf * std::cos(deg_to_rad(dirx)),
                    mach * cinf * std::cos(deg_to_rad(90. - dirx)),
                    mach * cinf * std::cos(deg_to_rad(dirz))
                    );
        }

        vgn_ = vel;
    }
    void MeshCell::init(const Vector3& vinf_air, const Freestream& fs, const Component& compo)
    {
        double rotation = compo.rotation_;
        double rotaxis = compo.rotaxis_;
        double rpm = compo.rpm_;
        auto pivot = compo.pivot_;
        double mach = compo.mach_;
        double dirx = compo.dirx_;
        double dirz = compo.dirz_;

        double cinf = std::sqrt(fs.gamma_ * fs.pinf_ / fs.rhoinf_);

        Vector3 vel(0.,0.,0.);

        //if (rotation)
        //{
        //    double om = rpm * 2. * PI / 60.; // rad/s
        //    Vector3 omega(0., 0., om);

        //    if (rotaxis == 0)
        //    {
        //        assert(false);
        //    }
        //    else if (rotaxis == 1)
        //    {
        //        assert(false);
        //    }
        //    else if (rotaxis == 2)
        //    {
        //        auto cnt = poly_.centroid();
        //        auto r = cnt - pivot;
        //        //auto rvec = Vector3(cnt(0), cnt(1), pivot(2));
        //        //double r = (rvec - pivot).len();
        //        //double alpha = std::atan2(cnt(0), cnt(1));
        //        //double thetax = -std::sin(alpha);
        //        //double thetay = std::cos(alpha);
        //        //vel = Vector3(
        //                //omega * r * thetax,
        //                //omega * r * thetay,
        //                //0.);
        //        vel = cross(omega, r);
        //    }
        //}
        //else
        //{

            //vel = Vector3(
                    //mach * cinf * std::cos(deg_to_rad(dirx)),
                    //mach * cinf * std::cos(deg_to_rad(90 - dirx)),
                    //mach * cinf * std::cos(deg_to_rad(dirz))
                    //);
        //}

        assert(!vel.isnan());

        Vector3 vinf = vinf_air - vel;

        Vector5 prim;
        prim_(0) = fs.rhoinf_;
        prim_(1) = vinf(0);
        prim_(2) = vinf(1);
        prim_(3) = vinf(2);
        prim_(4) = fs.pinf_;

        assert(!prim_.isnan());

        cons_sp1_ = prim_to_cons(prim_, fs.gamma_);
    }
    void minmax(const mcc& cell, Vector3& min, Vector3& max)
    {
        double xmin = TAILOR_BIG_POS_NUM;
        double ymin = TAILOR_BIG_POS_NUM;
        double zmin = TAILOR_BIG_POS_NUM;

        double xmax = TAILOR_BIG_NEG_NUM;
        double ymax = TAILOR_BIG_NEG_NUM;
        double zmax = TAILOR_BIG_NEG_NUM;

        for (const auto& mc: cell)
        {
            for (const auto& p: mc.poly().vertices())
            {
                double vx = p.r(0);
                double vy = p.r(1);
                double vz = p.r(2);

                xmin = std::min(xmin, vx);
                ymin = std::min(ymin, vy);
                zmin = std::min(zmin, vz);

                xmax = std::max(xmax, vx);
                ymax = std::max(ymax, vy);
                zmax = std::max(zmax, vz);
            }
        }

        min = Vector3(xmin, ymin, zmin);
        max = Vector3(xmax, ymax, zmax);
    }

    void MeshCell::reset_btype()
    {
        btype_ = BouType::undefined;
        for (auto& mf: face_)
        {
            mf.reset_btype();
        }
    }

    //const MeshCell* MeshCell::donor_addr() const
    //{
        //return donor_addr_;
    //}

    void MeshCell::reset_face_tags()
    {
        for (MeshFace& mf: face_)
        {
            mf.set_tag(FaceTag());
        }
    }

    void MeshCell::set_prim(const Vector5& other)
    {
        prim_ = other;
    }

    double MeshCell::D(int i, int j) const
    {
        return D(i,j);
    }

    const Vector5& MeshCell::R() const
    {
        return R_;
    }

    const Vector5& MeshCell::dQ() const
    {
        return dQ_;
    }

    const Vector5& MeshCell::old_dQ() const
    {
        return old_dQ_;
    }

    const Matrix5& MeshCell::D() const
    {
        return D_;
    }

    double MeshCell::R(int i) const
    {
        return R_(i);
    }

    const Vector5& MeshCell::prim() const
    {
        return prim_;
    }

    const Vector5& MeshCell::cons_sp1() const
    {
        return cons_sp1_;
    }

    const Vector5& MeshCell::cons_s() const
    {
        return cons_s_;
    }

    const Vector5& MeshCell::cons_n() const
    {
        return cons_n_;
    }

    const Vector5& MeshCell::cons_nm1() const
    {
        return cons_nm1_;
    }

    double MeshCell::prim(int i) const
    {
        return prim_(i);
    }

    double MeshCell::cons_sp1(int i) const
    {
        return cons_sp1_(i);
    }

    double MeshCell::cons_s(int i) const
    {
        return cons_s_(i);
    }

    double MeshCell::cons_n(int i) const
    {
        return cons_n_(i);
    }

    double MeshCell::cons_nm1(int i) const
    {
        return cons_nm1_(i);
    }

    void MeshCell::set_btype(BouType btype)
    {
        btype_ = btype;
    }

    BouType MeshCell::btype() const
    {
        return btype_;
    }

    void MeshCell::add_boundary(const Tag& t, BouType boutype)
    {
        if (boutype == BouType::wall) {
            wall_boundary_.push_back(t);
        }
        else if (boutype == BouType::dirichlet) {
            dirichlet_boundary_.push_back(t);
        }
        else if (boutype == BouType::farfield) {
            farfield_boundary_.push_back(t);
        }
        else if (boutype == BouType::empty) {
            empty_boundary_.push_back(t);
        }
        else if (boutype == BouType::interog) {
            interog_boundary_.push_back(t);
        }
        else {
            assert(false);
        }
    }
    void MeshCell::mark_to_be_erased()
    {
        erase_ = true;
    }

    void MeshCell::unmark_to_be_erased()
    {
        erase_ = false;
    }

    bool MeshCell::erase() const
    {
        return erase_;
    }

    bool MeshCell::operator<(const MeshCell& other) const
    {
        return tag_ < other.tag();
    }

    bool MeshCell::operator==(const MeshCell& other) const
    {
        return tag_ == other.tag();
    }

    void MeshCell::merge(const MeshCell& other)
    {
    }

    size_t MeshCell::npoint() const
    {
        return point_.size();
    }

    const Tag& MeshCell::interior_boundary() const
    {
        return interior_boundary_;
    }

    /*bool MeshCell::is_interior() const
      {
      return is_interior_;
      }

      bool MeshCell::is_wall() const
      {
      return is_wall_;
      }

      bool MeshCell::is_dirichlet() const
      {
      return is_dirichlet_;
      }*/

    bool MeshCell::near_boundary() const
    {
        if (near_wall() || near_dirichlet() || near_farfield() || near_interog())
        {
            return true;
        }

        return false;
    }

    bool MeshCell::near_empty() const
    {
        if (!empty_boundary_.empty())
        {
            return true;
        }

        return false;
    }

    bool MeshCell::near_interog() const
    {
        if (!interog_boundary_.empty())
        {
            return true;
        }

        return false;
    }

    bool MeshCell::near_wall() const
    {
        if (!wall_boundary_.empty())
        {
            return true;
        }

        return false;
    }

    bool MeshCell::near_dirichlet() const
    {
        if (!dirichlet_boundary_.empty())
        {
            return true;
        }

        return false;
    }

    bool MeshCell::near_farfield() const
    {
        if (!farfield_boundary_.empty())
        {
            return true;
        }

        return false;
    }

    //const std::vector<Tag>& MeshCell::wall_boundary() const
    const MeshCell::Boundaries& MeshCell::wall_boundary() const
    {
        return wall_boundary_;
    }

    //const std::vector<Tag>& MeshCell::dirichlet_boundary() const
    const MeshCell::Boundaries& MeshCell::dirichlet_boundary() const
    {
        return dirichlet_boundary_;
    }

    const MeshCell::Boundaries& MeshCell::farfield_boundary() const
    {
        return farfield_boundary_;
    }

    const MeshCell::Boundaries& MeshCell::interog_boundary() const
    {
        return interog_boundary_;
    }

    const MeshCell::Boundaries& MeshCell::empty_boundary() const
    {
        return empty_boundary_;
    }

    void MeshCell::add_wall_boundary(const Tag& t)
    {
        wall_boundary_.push_back(t);
    }

    void MeshCell::add_dirichlet_boundary(const Tag& t)
    {
        dirichlet_boundary_.push_back(t);
    }

    void MeshCell::add_farfield_boundary(const Tag& t)
    {
        farfield_boundary_.push_back(t);
    }

    void MeshCell::add_empty_boundary(const Tag& t)
    {
        empty_boundary_.push_back(t);
    }

    void MeshCell::set_interior_boundary(const Tag& t)
    {
        interior_boundary_ = t;
    }

    //void MeshCell::set_partition(bool p)
    //{
    //is_partition_ = p;
    //}
    //const bool MeshCell::is_partition() const
    //{
    //return is_partition_;
    //}

    //const std::vector<MeshFace>& MeshCell::face() const
    const MeshCell::Face& MeshCell::face() const
    {
        return face_;
    }

    //std::vector<MeshFace>& MeshCell::face_p()
    MeshCell::Face& MeshCell::face_p()
    {
        return face_;
    }

    void MeshCell::make_faces(Shape shape)
    {
        assert(face_.empty());
        assert(!point_.empty());

        if (shape == Shape::quad)
        {
            std::vector<MeshPoint> pts0;
            pts0 = {point_[0], point_[1], point_[2], point_[3]};
            face_.push_back(MeshFace(pts0));
            assert(!face_.back().face().vertex().empty());
        }
        else if (shape == Shape::tri)
        {
            std::vector<MeshPoint> pts0;
            pts0 = {point_[0], point_[1], point_[2]};
            face_.push_back(MeshFace(pts0));
            assert(!face_.back().face().vertex().empty());
        }
        else if (shape == Shape::tet)
        {
            std::vector<MeshPoint> pts0, pts1, pts2, pts3;
            pts0 = {point_[1], point_[2], point_[0]};
            pts1 = {point_[0], point_[3], point_[1]};
            pts2 = {point_[0], point_[2], point_[3]};
            pts3 = {point_[2], point_[1], point_[3]};
            face_.push_back(MeshFace(pts0));
            assert(!face_.back().face().vertex().empty());
            face_.push_back(MeshFace(pts1));
            assert(!face_.back().face().vertex().empty());
            face_.push_back(MeshFace(pts2));
            assert(!face_.back().face().vertex().empty());
            face_.push_back(MeshFace(pts3));
            assert(!face_.back().face().vertex().empty());
        }
        else if (shape == Shape::pri)
        {
            std::vector<MeshPoint> pts0, pts1, pts2, pts3, pts4;
            //pts0 = {point_[1], point_[0], point_[2]};
            pts0 = {point_[0], point_[1], point_[2]};
            pts1 = {point_[5], point_[4], point_[3]};
            pts2 = {point_[1], point_[0], point_[3], point_[4]};
            pts3 = {point_[3], point_[0], point_[2], point_[5]};
            pts4 = {point_[2], point_[1], point_[4], point_[5]};
            face_.push_back(MeshFace(pts0));
            assert(!face_.back().face().vertex().empty());
            face_.push_back(MeshFace(pts1));
            assert(!face_.back().face().vertex().empty());
            face_.push_back(MeshFace(pts2));
            assert(!face_.back().face().vertex().empty());
            face_.push_back(MeshFace(pts3));
            assert(!face_.back().face().vertex().empty());
            face_.push_back(MeshFace(pts4));
            assert(!face_.back().face().vertex().empty());
        }
        else if (shape == Shape::hex)
        {
            std::vector<MeshPoint> pts0, pts1, pts2, pts3, pts4, pts5;
            pts0 = {point_[0], point_[1], point_[2], point_[3]};
            pts1 = {point_[6], point_[5], point_[4], point_[7]};
            pts2 = {point_[2], point_[1], point_[5], point_[6]};
            pts3 = {point_[4], point_[0], point_[3], point_[7]};
            pts4 = {point_[1], point_[0], point_[4], point_[5]};
            pts5 = {point_[3], point_[2], point_[6], point_[7]};
            face_.push_back(MeshFace(pts0));
            assert(!face_.back().face().vertex().empty());
            face_.push_back(MeshFace(pts1));
            assert(!face_.back().face().vertex().empty());
            face_.push_back(MeshFace(pts2));
            assert(!face_.back().face().vertex().empty());
            face_.push_back(MeshFace(pts3));
            assert(!face_.back().face().vertex().empty());
            face_.push_back(MeshFace(pts4));
            assert(!face_.back().face().vertex().empty());
            face_.push_back(MeshFace(pts5));
            assert(!face_.back().face().vertex().empty());
        }
        else
        {
            assert(false);
        }
    }

    //void MeshCell::set_cand_donor(const std::vector<Donor>& cand_donor)
    void MeshCell::set_cand_donor(const Donorcon& cand_donor)
    {
        cand_donor_ = cand_donor;
    }

    /*const FinalCandDonor& MeshCell::cand_donor_mesh() const
    {
        return cand_donor_mesh_;
    }
    const FinalCandDonor& MeshCell::cand_donor_cell() const
    {
        return cand_donor_cell_;
    }
    const FinalCandDonor& MeshCell::prev_cand_donor_mesh() const
    {
        return prev_cand_donor_mesh_;
    }
    const FinalCandDonor& MeshCell::prev_cand_donor_cell() const
    {
        return prev_cand_donor_cell_;
    }*/

    /*const Array<Donor, 5>& MeshCell::cand_donor() const
    {
        return cand_donor_;
    }

    const Array<Donor, 5>& MeshCell::prev_cand_donor() const
    {
        return prev_cand_donor_;
    }*/

    //void MeshCell::set_boundary_type(boundary_t t)
    //{
    //boundary_type_ = t;
    //}

    //boundary_t MeshCell::boundary_type() const
    //{
    //return boundary_type_;
    //}

    const Tag& MeshCell::root_parent_mesh() const
    {
        return root_parent_mesh_;
    }

    const Tag& MeshCell::first_tag() const
    {
        return first_tag_;
    }

    const Donor& MeshCell::donor() const
    {
        return donor_;
    }

    //const std::vector<Donor>& MeshCell::cand_donor() const
    const Donorcon& MeshCell::cand_donor() const
    {
        return cand_donor_;
    }

    //const std::vector<Donor>& MeshCell::prev_cand_donor() const
    const Donorcon& MeshCell::prev_cand_donor() const
    {
        return prev_cand_donor_;
    }

    /*const Tag& MeshCell::donor_mesh() const
      {
      return donor_mesh_;
      }

      const Tag& MeshCell::donor_cell() const
      {
      return donor_cell_;
      }*/

    std::vector<Point> MeshCell::geom_point() const
    {
        std::vector<Point> pts;

        for (const MeshPoint& mp: point_)
            pts.push_back(mp.p());

        return pts;
    }

    //bool MeshCell::is_ghost() const
    //{
    //assert(boundary_type_ != boundary_t::undefined);
    //return is_ghost_;
    //}

    void MeshCell::update_prev_cand_donor()
    {
        prev_cand_donor_ = cand_donor_;
    }

    /*void MeshCell::reset_cand_donor(const Tag& meshtag)
      {
      assert(OGA_cell_type_ == OGA_cell_type_t::field);

      bool anyleft = true;

      while(anyleft)
      {
      for (int i=0; i<cand_donor_mesh_.size(); ++i)
      {
      if (cand_donor_mesh_[i] == meshtag)
      {
      cand_donor_mesh_.erase(cand_donor_mesh_[i]);
      cand_donor_cell_.erase(cand_donor_cell_[i]);
      break;
      }
      }

      anyleft = false;
      for (int i=0; i<cand_donor_mesh_.size(); ++i)
      {
      if (cand_donor_mesh_[i] == meshtag)
      {
      anyleft = true;
      }
      }
      }
      }*/

    void MeshCell::reset_oga_status()
    {
        cand_donor_.clear();
        set_oga_cell_type(OGA_cell_type_t::undefined); 
        donor_.mesh_tag_ = Tag();
        donor_.cell_tag_ = Tag();
    }

    const Tag& MeshCell::parent_mesh() const
    {
        return parent_mesh_;
    }

    void MeshCell::remove_all_neighbors()
    {
        //nei_.clear();
        pnei_.clear();
    }
    /*void MeshCell::remove_all_non_pnei()
      {
      nei_.clear();
      }*/

    void MeshCell::remove_parent_cells_of_vertices()
    {
        for (MeshPoint& _p: point_)
        {
            _p.remove_parent_cells();
        }
    }

    void MeshCell::set_point_tag(int i, const Tag& t)
    {
        assert(i >= 0);
        assert(i < point_.size());

        point_[i].set_tag(t);
    }

    /*const bool MeshCell::is_resident() const
      {
      return residency_;
      }*/

    void MeshCell::rotate_points(double ang, int axis, const Vector3& rot_axis)
    {
        for (MeshFace& mf: face_)
        {
            mf.rotate(ang, axis, rot_axis);
        }
        for (MeshPoint& point: point_)
        {
            point.rotate_point(ang, axis, rot_axis);
        }

        poly_.rotate_points(ang, axis, rot_axis);
        //polytope_->rotate_points(ang, rot_axis);
    }

    void MeshCell::move_points(const Vector3& v)
    {
        for (MeshFace& mf: face_)
        {
            mf.move(v);
        }
        for (MeshPoint& point: point_)
        {
            point.move_point(v);
        }

        poly_.move_points(v);
        //polytope_->move_points(v);
    }

    void MeshCell::set_donor(const Tag& im, const Tag& ic, const MeshCell* donor)
    {
        donor_.mesh_tag_ = im;
        donor_.cell_tag_ = ic;
        donor_.addr_ = donor;
        assert(im.isvalid());
        //assert(im() == 0 || im() == 1);
    }

    void MeshCell::remove_cand_donor(const Tag& donor_cell, const Tag& donor_mesh)
    {
        auto it = std::find_if(cand_donor_.begin(), cand_donor_.end(), [&](const auto& d){return (d.mesh_tag_ == donor_mesh && d.cell_tag_ == donor_cell);});
        assert(it != cand_donor_.end());
        cand_donor_.erase(it);
    }

    //void MeshCell::set_oga_cell_type(OGA_cell_type_t t, const Tag& donor_mesh, const Tag& donor_cell, const MeshCell* donor_addr)
    //void MeshCell::set_oga_cell_type(OGA_cell_type_t t)
    //{
        //if (t == OGA_cell_type_t::receptor || t == OGA_cell_type_t::mandat_receptor)
        //{
            //assert(donor_mesh() == 0 || donor_mesh() == 1);
        //}
        //set_oga_cell_type(t);
        //set_donor(donor_mesh, donor_cell, donor_addr);
        //add_cand_donor(donor_mesh, donor_cell, donor_addr); // do this seperately.
    //}

    void MeshCell::set_oga_cell_type(OGA_cell_type_t t)
    {
        OGA_cell_type_ = t;
    }

    OGA_cell_type_t MeshCell::oga_cell_type() const
    {
        return OGA_cell_type_;
    }

    void MeshCell::deparent_self_from_vertices()
    {
        for (MeshPoint& p: point_)
            p.remove_parent_cell(tag_);
    }

    //void MeshCell::deparent_from_faces(const Tag& t)
    //{
    //for (MeshFace& mf: face_)
    //{
    //mf.remove_parent_cell(t);
    //}
    //}

    void MeshCell::deparent_self_from_faces()
    {
        for (MeshFace& mf: face_)
        {
            mf.remove_parent_cell(tag_);
        }
    }

    void MeshCell::deparent_neis_from_faces()
    {
        // removes non-self parent cells of faces.

        bool goon = true;
        while(goon)
        {
            goon = false;
            for (MeshFace& mf: face_)
            {
                if (!mf.is_boundary())
                {
                    for (auto pc = mf.parent_cell().begin(); pc != mf.parent_cell().end(); ++pc)
                    {
                        if (*pc != tag_)
                        {
                            mf.remove_parent_cell(*pc);
                            goon = true;
                            break;
                        }
                    }
                }
            }
        }
    }

    void MeshCell::remove_all_nonself_parents_of_faces()
    {
        for (auto& mf: face_)
        {
            if (mf.is_boundary()) {continue;}

            for (const auto& pc: mf.parent_cell())
            {
                if (pc != tag_)
                {
                    mf.remove_parent_cell(pc);
                }
            }
        }
    }

    void MeshCell::remove_all_parent_cells_of_faces()
    {
        for (auto& mf: face_)
        {
            mf.remove_parent_cells();
        }
    }

    void MeshCell::deparent_cell_from_faces(const Tag& tag)
    {
        for (MeshFace& mf: face_)
        {
            if (mf.is_boundary()) {
                continue;
            }
            mf.remove_parent_cell(tag);
        }
    }

    /*void MeshCell::remove_self_from_neighbors(const Tag& t)
      {
      assert(t() == -1);
      nei_.erase(std::remove(nei_.begin(), nei_.end(), t), nei_.end());
      }*/

    void MeshCell::remove_neighbor(const Tag& t)
    {
        assert(t.isvalid());
        pnei_.erase(std::remove_if(pnei_.begin(), pnei_.end(), [&](const auto& a){return a == t;}), pnei_.end());
    }

    //MPI_Datatype mpi_dtype_meshcell = MPI_DATATYPE_NULL;

    /*void MeshCell::make_layout()
      {
    // members.
    // 1: int parent_mesh;
    // 1: int donor_mesh;
    // 1: int donor_cell;
    // 1: int partition;
    // 2: Array<int, max_nnei> nei;
    // 3: Array<int, max_npoint> ipoint;
    // 4: Array<int, max_nbface> bface;
    // 5: Array<int, max_niface> iface;
    // 6: Array<int, max_npface> pface;
    // 7: Array<MeshPoint, max_npoint> point;
    // 8: Polygon polygon;
    // 9: OGA_cell_type_t OGA_cell_type;
    // 10: Array<int, max_ncandidate_donor> candidate_donor_mesh;
    // 10: Array<int, max_ncandidate_donor> candidate_donor_cell;
    int flag;
    MPI_Initialized(&flag);
    if (mpi_dtype_meshcell == MPI_DATATYPE_NULL && flag == 1)
    {
    int nblock = 10;
    // block count.
    int block_count[nblock];
    block_count[0] = 4;
    block_count[1] = 1;
    block_count[2] = 1;
    block_count[3] = 1;
    block_count[4] = 1;
    block_count[5] = 1;
    block_count[6] = 1;
    block_count[7] = 1;
    block_count[8] = 1;
    block_count[9] = 2;
    // address.
    MPI_Aint addr[nblock];
    MPI_Get_address(&parent_mesh, &addr[0]);
    MPI_Get_address(&nei, &addr[1]);
    MPI_Get_address(&ipoint, &addr[2]);
    MPI_Get_address(&bface, &addr[3]);
    MPI_Get_address(&iface, &addr[4]);
    MPI_Get_address(&pface, &addr[5]);
    MPI_Get_address(&point, &addr[6]);
    MPI_Get_address(&polygon, &addr[7]);
    MPI_Get_address(&OGA_cell_type, &addr[8]);
    MPI_Get_address(&candidate_donor_mesh, &addr[9]);
    // offset.
    MPI_Aint offset[nblock];
    offset[0] = 0;
    offset[1] = offset[0] + addr[1] - addr[0];
    offset[2] = offset[1] + addr[2] - addr[1];
    offset[3] = offset[2] + addr[3] - addr[2];
    offset[4] = offset[3] + addr[4] - addr[3];
    offset[5] = offset[4] + addr[5] - addr[4];
    offset[6] = offset[5] + addr[6] - addr[5];
    offset[7] = offset[6] + addr[7] - addr[6];
    offset[8] = offset[7] + addr[8] - addr[7];
    offset[9] = offset[8] + addr[9] - addr[8];
    // block type.
    MPI_Datatype mdt0;
    MPI_Datatype mdt1;
    MPI_Datatype mdt2;
    MPI_Datatype mdt3;
    MPI_Datatype mdt4;
    MPI_Datatype mdt5;
    MPI_Datatype mdt6;
    nei.make_layout(mdt0, MPI_INT);
    ipoint.make_layout(mdt1, MPI_INT);
    bface.make_layout(mdt2, MPI_INT);
    iface.make_layout(mdt3, MPI_INT);
    pface.make_layout(mdt4, MPI_INT);
    point.make_layout(mdt5, mpi_dtype_meshpoint);
    candidate_donor_mesh.make_layout(mdt6, MPI_INT);

    MPI_Datatype block_type[nblock];
    block_type[0] = MPI_INT;
    block_type[1] = mdt0;
    block_type[2] = mdt1;
    block_type[3] = mdt2;
    block_type[4] = mdt3;
    block_type[5] = mdt4;
    block_type[6] = mdt5;
    block_type[7] = mpi_dtype_polygon;
    block_type[8] = MPI_INT;
    block_type[9] = mdt6;
    // define and commit type.
    MPI_Type_create_struct(nblock, block_count, offset, block_type, &mpi_dtype_meshcell);
    MPI_Type_commit(&mpi_dtype_meshcell);
    // free MPI datatypes.
    MPI_Type_free(&mdt0);
    MPI_Type_free(&mdt1);
    MPI_Type_free(&mdt2);
    MPI_Type_free(&mdt3);
    MPI_Type_free(&mdt4);
    MPI_Type_free(&mdt5);
    MPI_Type_free(&mdt6);
}
}*/

std::vector<std::vector<Tag>> MeshCell::face_vertex()
{
    std::vector<std::vector<Tag>> face_vertex;

    if (point_.size() == 2)
    {
        face_vertex.resize(1);

        face_vertex[0].push_back(point_[0].tag());
        face_vertex[0].push_back(point_[1].tag());
    }
    else if (point_.size() == 3)
    {
        face_vertex.resize(3);

        face_vertex[0].push_back(point_[0].tag());
        face_vertex[0].push_back(point_[1].tag());

        face_vertex[1].push_back(point_[1].tag());
        face_vertex[1].push_back(point_[2].tag());

        face_vertex[2].push_back(point_[2].tag());
        face_vertex[2].push_back(point_[0].tag());
    }
    else if (point_.size() == 4)
    {
        face_vertex.resize(4);

        face_vertex[0].push_back(point_[0].tag());
        face_vertex[0].push_back(point_[1].tag());

        face_vertex[1].push_back(point_[1].tag());
        face_vertex[1].push_back(point_[2].tag());

        face_vertex[2].push_back(point_[2].tag());
        face_vertex[2].push_back(point_[3].tag());

        face_vertex[3].push_back(point_[3].tag());
        face_vertex[3].push_back(point_[0].tag());
    }
    else
    {
        std::cout << "invalid point size in face_vertex()" << std::endl;
        exit(0);
    }

    return face_vertex;
}

/*void MeshCell::set_ipoint(int i, const Tag& ip)
  {
  ipoint_[i] = ip;
  }
  void MeshCell::set_ipoint(const Tag& t)
  {
  std::fill(ipoint_.begin(), ipoint_.end(), t);
  }*/
void MeshCell::set_tag(const Tag& t)
{
    assert(t.isvalid());
    if (!first_tag_.isvalid())
        first_tag_ = t;
    tag_ = t;
}
void MeshCell::set_parent_mesh(const Tag& celltag)
{
    assert(celltag.isvalid());
    if (!root_parent_mesh_.isvalid())
        root_parent_mesh_ = celltag;
    parent_mesh_ = celltag;
}
const Polyhedron& MeshCell::poly() const
{
    return poly_;
}
/*const Polytope& MeshCell::polytope() const
  {
  return *polytope_;
  }*/
/*std::shared_ptr<Polytope> MeshCell::polytope() const
  {
  return polytope_;
  }*/
//const std::vector<MeshPoint>& MeshCell::point() const
//const std::array<MeshPoint, NPOINT>& MeshCell::point() const
//{
//return point_;
//}
const MeshCell::MCPoint& MeshCell::point() const
{
    return point_;
}
const MeshPoint& MeshCell::point(int i) const
{
    return point_[i];
}
/*const std::vector<Tag>& MeshCell::ipoint() const
  {
  return ipoint_;
  }*/
/*const std::vector<Tag>& MeshCell::nei() const
  {
  return nei_;
  }*/
const MeshCell::Pnei& MeshCell::pnei() const
{
    return pnei_;
}
/*const Tag& MeshCell::ipoint(int i) const
  {
  return ipoint_[i];
  }*/
/*const Tag& MeshCell::nei(int i) const
  {
  return nei_[i];
  }*/
const Tag& MeshCell::pnei(int i) const
{
    return pnei_[i];
}
const Tag& MeshCell::tag() const
{
    return tag_;
}

MeshCell::MeshCell(): sumarea_(0.), erase_(false), OGA_cell_type_(OGA_cell_type_t::undefined), dtao_(0.)
{
    prim_     = 0.;
    cons_sp1_ = 0.;
    cons_s_   = 0.;
    cons_n_   = 0.;
    cons_nm1_ = 0.;
    dQ_       = 0.;
    old_dQ_   = 0.;
    R_        = 0.;
    R_mid_    = 0.;
    D_        = 0.;
    D_mid_    = 0.;
}

MeshCell::MeshCell(const Tag& tag, const Tag& parent_mesh, const std::vector<MeshPoint>& point, BouType btype, Shape shape): tag_(tag), parent_mesh_(parent_mesh), OGA_cell_type_(OGA_cell_type_t::undefined), erase_(false), btype_(btype), sumarea_(0.)
{
    vgn_ = 0;

    assert(btype != BouType::undefined);
    assert(tag.isvalid());

    //if (btype == boundary_t::wall) {
    //is_wall_ = true;
    //}
    //else if (btype == boundary_t::dirichlet) {
    //is_dirichlet_ = true;
    //}
    //else {
    //is_interior_ = true;
    //}

    point_ = point;
    //assert(point.size() <= NPOINT);
    //std::copy(point.begin(), point.end(), point_.begin());
    //this->point_ = point;
    make_faces(shape);

    std::vector<Point> pts;
    for (const MeshPoint& mp: point)
    {
        pts.push_back(mp.p());
    }

    if (shape == Shape::tet)
    {
        assert(point.size() == 4);
        assert(pts.size() == 4);
    }

    if (pts.size() > 8)
    {
        std::cout << "pts size: " << pts.size() << std::endl;
    }
    assert(pts.size() <= 8);
    poly_ = Polyhedron(pts, shape);
    //poly_.set_vertex(pts);
    //poly_.set_edge();

    first_tag_ = tag_;
    root_parent_mesh_ = parent_mesh_;

    prim_     = 0.;
    cons_sp1_ = 0.;
    cons_s_   = 0.;
    cons_n_   = 0.;
    cons_nm1_ = 0.;
    dQ_       = 0.;
    old_dQ_   = 0.;
    R_        = 0.;
    R_mid_    = 0.;

    D_ = 0.;
    D_mid_ = 0.;
}

void MeshCell::add_vertex(const MeshPoint& p, const Tag& pointtag)
{
    point_.add(p);
    //assert(npoint_ < NPOINT);
    //point_[npoint_] = p;
    //++npoint_;
    //point_.push_back(p);
}

/*void MeshCell::add_nei(const Tag& celltag)
  {
  nei_.push_back(celltag);
  }*/

void MeshCell::add_pnei(const Tag& celltag)
{
    pnei_.push_back(celltag);

    //auto it = std::unique(pnei_.begin(), pnei_.end());
    //assert(it == pnei_.end());
}

void MeshCell::add_cand_donor(const Tag& donor_mesh, const Tag& donor_cell, const MeshCell* addr)
{
    cand_donor_.push_back(Donor(donor_mesh, donor_cell, addr));
    //cand_donor_mesh_.add(donor_mesh);
    //cand_donor_cell_.add(donor_cell);
    //cand_donor_addr_.add(donor);
    /*assert(donor_mesh.isvalid());
      assert(donor_cell.isvalid());

      int count_cell = cand_cell_tag_index_map.left.count(donor_cell());
      int count_mesh = cand_mesh_tag_index_map.left.count(donor_mesh());

      if ((count_cell == 0 && count_mesh == 0) || (count_cell == 0 && count_mesh != 0) || (count_cell != 0 && count_mesh == 0))
      {
      candidate_donor_mesh.push_back(donor_mesh);
      candidate_donor_cell.push_back(donor_cell);
    // insert to bimap.
    cand_cell_tag_index_map.insert(boost::bimap<int, int>::value_type(donor_cell(), candidate_donor_cell.size() - 1));
    cand_mesh_tag_index_map.insert(boost::bimap<int, int>::value_type(donor_mesh(), candidate_donor_mesh.size() - 1));

    return;
    }

    //std::cout << "count_mesh = " << count_mesh << std::endl;
    //std::cout << "count_cell = " << count_cell << std::endl;
    //assert(false);*/
}

void MeshCell::set_parent_cell_of_vertices()
{
    for (MeshPoint& _p: point_)
    {
        assert(tag_.isvalid());
        _p.add_parent_cell(tag_);
    }
}

void MeshCell::set_parent_cell_of_faces()
{
    for (auto& mf: face_)
    {
        mf.add_parent_cell(tag_);
    }
}

/*size_t MeshCell::mem() const
{
    size_t size = 0;

    size += sizeof(bool) * 5;
    size += wall_boundary_.mem();
    size += dirichlet_boundary_.mem();
    size += farfield_boundary_.mem();
    size += empty_boundary_.mem();
    size += sizeof(int) * 5;
    size += sizeof(int) * 2;
    size += pnei_.mem();
    size += point_.mem();
    size += face_.mem();
    size += poly_.mem();
    size += sizeof(int);
    //size += cand_donor_mesh_.mem();
    //size += cand_donor_cell_.mem();

    return size;
}*/

int nnon_boundary_face(const MeshCell& mc)
{
    int count = 0;
    for (const MeshFace& mf: mc.face())
    {
        if (!mf.is_boundary())
        {
            ++count;
        }
    }
    return count;
}
}


