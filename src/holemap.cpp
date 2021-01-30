#include "holemap.h"

namespace Tailor
{
    //bool HoleMap::holeless() const
    //{
        //assert(holeless_ == false);
        //return holeless_;
    //}

    //const AABB& HoleMap::wallaabb() const
    //{
        //return wallaabb_;
    //}

    const Tag& HoleMap::tag() const
    {
        return tag_;
    }

    bool HoleMap::do_wallaabb_intersect_outer(const Mesh& mesh)
    {
        if (do_wallaabb_intersect_outer_(mesh.dirichlet_boundaries())) {
            std::cout << "mesh: " << mesh.tag()() << std::endl;
            std::cout << "wall min(0): " << wall_aabb_.min(0) << std::endl;
            std::cout << "wall min(1): " << wall_aabb_.min(1) << std::endl;
            std::cout << "wall min(2): " << wall_aabb_.min(2) << std::endl;
            std::cout << "wall max(0): " << wall_aabb_.max(0) << std::endl;
            std::cout << "wall max(1): " << wall_aabb_.max(1) << std::endl;
            std::cout << "wall max(2): " << wall_aabb_.max(2) << std::endl;
            assert(false);
            return true;
        }
        if (do_wallaabb_intersect_outer_(mesh.farfield_boundaries())) {
            std::cout << "mesh: " << mesh.tag()() << std::endl;
            std::cout << "wall min(0): " << wall_aabb_.min(0) << std::endl;
            std::cout << "wall min(1): " << wall_aabb_.min(1) << std::endl;
            std::cout << "wall min(2): " << wall_aabb_.min(2) << std::endl;
            std::cout << "wall max(0): " << wall_aabb_.max(0) << std::endl;
            std::cout << "wall max(1): " << wall_aabb_.max(1) << std::endl;
            std::cout << "wall max(2): " << wall_aabb_.max(2) << std::endl;
            assert(false);
            return true;
        }
        if (do_wallaabb_intersect_outer_(mesh.interog_boundaries())) {
            //if (world_.rank() == 4)
            {
            std::cout << "mesh: " << mesh.tag()() << std::endl;
            std::cout << "wall min(0): " << wall_aabb_.min(0) << std::endl;
            std::cout << "wall min(1): " << wall_aabb_.min(1) << std::endl;
            std::cout << "wall min(2): " << wall_aabb_.min(2) << std::endl;
            std::cout << "wall max(0): " << wall_aabb_.max(0) << std::endl;
            std::cout << "wall max(1): " << wall_aabb_.max(1) << std::endl;
            std::cout << "wall max(2): " << wall_aabb_.max(2) << std::endl;
                mesh.print_as_vtk_geometry("io.vtk");
                mesh.print_wall_as_vtk("wall.vtk");
                mesh.print_interog_as_vtk("interog.vtk");
            }
            assert(false);
            return true;
        }

        return false;
    }

    bool HoleMap::do_wallaabb_intersect_outer_(const mcc& container)
    {
        if (!container.empty())
        {
            for (const MeshCell& mc: container)
            {
                //if (wallaabb_.do_intersect(AABB(mc.poly()))) // hope this handles horizontal AABB.
                //if (wall_aabb_.do_intersect(AABB(mc.poly()))) // hope this handles horizontal AABB.
                if (wall_aabb_.do_intersect(mc.poly())) // hope this handles horizontal AABB.
                {
                    //std::cout << "mc.aabb.min(0): " << AABB(mc.poly()).min(0) << std::endl;
                    //std::cout << "mc.aabb.min(1): " << AABB(mc.poly()).min(1) << std::endl;
                    //std::cout << "mc.aabb.min(2): " << AABB(mc.poly()).min(2) << std::endl;
                    //std::cout << "mc.aabb.max(0): " << AABB(mc.poly()).max(0) << std::endl;
                    //std::cout << "mc.aabb.max(1): " << AABB(mc.poly()).max(1) << std::endl;
                    //std::cout << "mc.aabb.max(2): " << AABB(mc.poly()).max(2) << std::endl;
                    return true;
                }
            }
        }

        return false;
    }

    /*bool HoleMap::do_wallaabb_intersect_outer(const Mesh* dirichlet, const Mesh* farfield, const Mesh* interog)
    {
        if (dirichlet != nullptr)
        {
            for (const MeshCell& mc: dirichlet->cell())
            {
                if (wallaabb_.do_intersect(AABB(mc.poly()))) // hope this handles horizontal AABB.
                {
                    return true;
                }
            }
        }

        if (farfield != nullptr)
        {
            for (const MeshCell& mc: farfield->cell())
            {
                if (wallaabb_.do_intersect(AABB(mc.poly()))) // hope this handles horizontal AABB.
                {
                    return true;
                }
            }
        }

        if (interog != nullptr)
        {
            for (const MeshCell& mc: interog->cell())
            {
                if (wallaabb_.do_intersect(AABB(mc.poly()))) // hope this handles horizontal AABB.
                {
                    return true;
                }
            }
        }

        return false;
    }*/

    bool HoleMap::is_inside_holebin(const vec3<double>& querypoint) const
    {
        //if (holeless_)
        //{
            //return false;
        //}

        //if (wallaabb_.do_intersect(querypoint) == false)
        if (wall_aabb_.do_intersect(querypoint) == false)
        {
            return false;
        }

        if (only_aabb_)
        {
            return true;
        }

        //std::vector<BinRMTag> tag;
        //BinRMTag tag;
        //wallrm_->get_bintag_adaptive(querypoint, tag, false);
        //wallrm_->get_bintag_adaptive_unique(querypoint, tag, -1, -1);

        auto adtpoint = ADTPoint(querypoint);
        std::vector<int> res = wall_adt_.search(adtpoint);


        //if (!tag.isvalid())
        if (res.empty())
        //if (tag.empty())
        {
            return false;
        }

        assert(res.size() == 1);

        //for (const auto& t: tag)
        {
            //if (wallrm_->bin(BinRMTag(Tag(res[0]), Tag(0))).load() == -2) {
            if (wallrm_->bin(Tag(res[0])).load() == -2) {
            //if (wallrm_->bin(tag).load() == -2) {
                return true;
            }
            //auto it = std::find_if(holebin_.begin(), holebin_.end(), [&](const Tag& hbt){return hbt == t.bintag();});
            //if (it != holebin_.end()) {
                //return true;
            //}
        }

        return false;
    }

    //void HoleMap::make(const Mesh& mesh, bool pseudo3D)
    void HoleMap::make(const Mesh* mesh, bool pseudo3D)
    {
        //std::vector<bool> gholeless;
        //boost::mpi::all_gather(world_, holeless_, gholeless);
        //holeless_ = all_of(gholeless.begin(), gholeless.end(), [&](bool b){return b == false;});

        //bool outer_inter = do_wallaabb_intersect_outer(*mesh);
        //only_aabb_ = !outer_inter;

        //if (only_aabb_ == false)
        //if (outer_inter)
        {
            assert(wallrm_ == nullptr);
            wallrm_ = std::make_unique<RegularMesh>();
            wallrm_->set_tag(Tag(0));
            wallrm_->set_nstripe(2, 2, 2);
            wallrm_->set_aabb(wall_aabb_);
            wallrm_->calc_step_length();
            wallrm_->insert_bins(0, 0, world_.rank());
            wallrm_->insert_bin_addresses(wallrm_.get());

            assert(outerrm_ == nullptr);
            outerrm_ = std::make_unique<RegularMesh>();
            outerrm_->set_tag(Tag(0));
            outerrm_->set_nstripe(2, 2, 2);
            outerrm_->set_aabb(wall_aabb_);
            outerrm_->calc_step_length();
            outerrm_->insert_bins(0, 0, world_.rank());
            outerrm_->insert_bin_addresses(outerrm_.get());

            if (mesh != nullptr)
            {
                for (const MeshCell& mc: mesh->wall_boundaries()) {
                    wallrm_->register_resident_cell(mc, mc.tag(), mesh->tag(), world_.rank());
                }
                for (const MeshCell& mc: mesh->dirichlet_boundaries()) {
                    outerrm_->register_resident_cell(mc, mc.tag(), mesh->tag(), world_.rank());
                }
                for (const MeshCell& mc: mesh->farfield_boundaries()) {
                    outerrm_->register_resident_cell(mc, mc.tag(), mesh->tag(), world_.rank());
                }
                for (const MeshCell& mc: mesh->empty_boundaries()) {
                    outerrm_->register_resident_cell(mc, mc.tag(), mesh->tag(), world_.rank());
                }
                for (const MeshCell& mc: mesh->interog_boundaries()) {
                    outerrm_->register_resident_cell(mc, mc.tag(), mesh->tag(), world_.rank());
                }
            }

            bool refine = true;
            int nrefine = 1;
            while(refine)
            {
                refine = false;

                {
                    std::vector<RegularMesh> grm;
                    boost::mpi::all_gather(world_, *wallrm_, grm);
                    //wallrm_->update_address();
                    for (int i = 0; i<grm.size(); ++i)
                    {
                        if (i == world_.rank()) {
                            continue;
                        }

                        wallrm_->merge(grm[i]);
                    }
                }
                {
                    std::vector<RegularMesh> grm;
                    boost::mpi::all_gather(world_, *outerrm_, grm);
                    //outerrm_->update_address();
                    for (int i = 0; i<grm.size(); ++i)
                    {
                        if (i == world_.rank()) {
                            continue;
                        }

                        outerrm_->merge(grm[i]);
                    }
                }

                for (const Bin& wb: wallrm_->bin())
                {
                    if (wb.cell().empty()) {
                        continue;
                    }

                    if (outerrm_->bin(wb.tag()).cell().empty() == false)
                    {
                        refine = true;
                        break;
                    }
                }

                if (refine)
                {
                    wallrm_.reset();
                    outerrm_.reset();

                    ++nrefine;
                    int ns = 2. * nrefine;

                    wallrm_ = std::make_unique<RegularMesh>();
                    wallrm_->set_tag(Tag(0));
                    wallrm_->set_nstripe(ns, ns, ns);
                    wallrm_->set_aabb(wall_aabb_);
                    wallrm_->calc_step_length();
                    wallrm_->insert_bins(0, 0, world_.rank());
                    wallrm_->insert_bin_addresses(wallrm_.get());

                    outerrm_ = std::make_unique<RegularMesh>();
                    outerrm_->set_tag(Tag(0));
                    outerrm_->set_nstripe(ns, ns, ns);
                    outerrm_->set_aabb(wall_aabb_);
                    outerrm_->calc_step_length();
                    outerrm_->insert_bins(0, 0, world_.rank());
                    outerrm_->insert_bin_addresses(outerrm_.get());
                }
                else {
                    break;
                }
            }

            for (Bin& wb: wallrm_->bin_)
            {
                if (wb.cell().empty() == false) {
                    wb.set_load(-2);
                }
            }

            //if (world_.rank() == 0)
            //{
                //if (mesh != nullptr && mesh->tag()() == 5)
                //{
                    //wallrm_->print("hm");
                //}
            //}

            wallrm_->clear_cells();
            outerrm_.reset();

            for (Bin& b: wallrm_->bin_)
            {
                assert(b.tag().isvalid());
                assert(b.tag()() < wallrm_->bin().size());
                if (b.load() == -2)
                {
                    wallrm_->flood(b.tag(), wall_aabb_);
                }
            }

            for (Bin& b: wallrm_->bin_)
            {
                if (b.load() == -1)
                {
                    b.set_load(-2);
                }
            }

            assert(wallrm_ != nullptr);
            make_adt(*wallrm_, wall_adt_, -1);
        }
    }

    /*void HoleMap::make(const Mesh* wall, const Mesh* dirichlet, const Mesh* farfield, const Mesh* interog, bool pseudo3D)
    {
        if (wall->cell().empty())
        {
            holeless_ = true;
            return;
        }
        else
        {
            holeless_= false;
        }

        bool outer_inter = do_wallaabb_intersect_outer(dirichlet, farfield, interog);
        only_aabb_ = !outer_inter;

        if (only_aabb_ == false)
        {
            // TODO beware you have not insertted bins!
            wallrm_ = std::make_unique<RegularMesh>();
            wallrm_->set_tag(0);
            wallrm_->set_nstripe(2, 2, 2);
            wallrm_->set_aabb(wallaabb_);
            wallrm_->calc_step_length();

            outerrm_ = std::make_unique<RegularMesh>();
            outerrm_->set_tag(0);
            outerrm_->set_nstripe(2, 2, 2);
            outerrm_->set_aabb(outeraabb_);
            outerrm_->calc_step_length();

            wallrm_->register_mesh(*wall, world_.rank(), false);

            std::vector<Mesh> outervector;

            if (dirichlet != nullptr) {
                outerrm_->register_overlapping_mesh(*dirichlet, world_.rank(), false);
                outervector.push_back(*dirichlet);
            }
            if (farfield != nullptr) {
                outerrm_->register_overlapping_mesh(*farfield, world_.rank(), false);
                outervector.push_back(*farfield);
            }
            if (interog != nullptr) {
                outerrm_->register_overlapping_mesh(*interog, world_.rank(), false);
                outervector.push_back(*interog);
            }

            std::vector<Mesh> wallvector = {*wall};

            bool refine = true;
            while(refine)
            {
                refine = false;
                for (const Bin& wb: wallrm_->bin())
                {
                    if (wb.cell().empty()) {
                        continue;
                    }

                    if (outerrm_->bin(wb.tag()).cell().empty() == false)
                    {
                        refine = true;
                        break;
                    }
                }

                if (refine)
                {
                    wallrm_->refine(wallvector, world_.rank(), pseudo3D);
                    outerrm_->refine(outervector, world_.rank(), pseudo3D);
                }
                else {
                    break;
                }
            }

            for (Bin& wb: wallrm_->bin_)
            {
                if (wb.cell().empty() == false) {
                    wb.set_load(-2);
                    //holebin_.push_back(wb.tag());
                }
            }

            wallrm_->clear_cells();
            outerrm_.reset();

            for (Bin& b: wallrm_->bin_)
            {
                if (b.load() == -2)
                {
                    wallrm_->flood(b.tag()(), wallaabb_);
                }
            }

            for (Bin& b: wallrm_->bin_)
            {
                if (b.load() == -1)
                {
                    b.set_load(-2);
                }
            }
        }
    }*/

    HoleMap::HoleMap(const MPI_Comm& comm, const Mesh* mesh, bool pseudo3D, const std::vector<AABB>& hole_aabb): world_(comm, boost::mpi::comm_attach)
    {
        if (mesh != nullptr)
        {
            tag_ = mesh->tag();
        }
        /*{
            std::vector<Point> pts;
            std::set<Point> s;

            s = mesh.bou_raw_point(BouType::wall);
            pts.insert(pts.end(), s.begin(), s.end());

            wallaabb_ = AABB(pts);
        }
        {
            std::vector<Point> pts;
            std::set<Point> s;

            s = mesh.bou_raw_point(BouType::dirichlet);
            pts.insert(pts.end(), s.begin(), s.end());

            s = mesh.bou_raw_point(BouType::farfield);
            pts.insert(pts.end(), s.begin(), s.end());

            s = mesh.bou_raw_point(BouType::empty);
            pts.insert(pts.end(), s.begin(), s.end());

            s = mesh.bou_raw_point(BouType::interog);
            pts.insert(pts.end(), s.begin(), s.end());

            outeraabb_ = AABB(pts);
        }*/

        wall_aabb_ = hole_aabb[tag_()];
        make(mesh, pseudo3D);
    }

    /*HoleMap::HoleMap(const MPI_Comm& comm, const Mesh* wall, const Mesh* dirichlet, const Mesh* farfield, const Mesh* interog, bool pseudo3D): world_(comm, boost::mpi::comm_attach), wallaabb_(wall->rawpoint()), tag_(wall->tag())
    {
        std::vector<Point> pts;
        if (dirichlet != nullptr)
        {
            auto dirichlet_raw = dirichlet->rawpoint();
            pts.insert(pts.end(), dirichlet_raw.begin(), dirichlet_raw.end());
        }
        if (farfield != nullptr)
        {  
            auto farfield_raw = farfield->rawpoint();
            pts.insert(pts.end(), farfield_raw.begin(), farfield_raw.end());
        }
        if (interog != nullptr)
        {  
            auto interog_raw = interog->rawpoint();
            pts.insert(pts.end(), interog_raw.begin(), interog_raw.end());
        }

        outeraabb_ = AABB(pts);
        make(wall, dirichlet, farfield, interog, pseudo3D);
    }*/
}
