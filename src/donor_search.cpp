#include "donor_search.h"

namespace Tailor
{
    DonorSearcher::DonorSearcher(boost::mpi::communicator* comm, SpatialPartitionContainer* spc, Profiler* profiler, DonorSearchAlgo dsa, bool mergebins, bool verbose, bool print_info): comm_(comm), profiler_(profiler), verbose_(verbose), uniproc_(false), dsa_(dsa), mergebins_(mergebins), print_info_(print_info)
    {
        if (comm_->size() == 0) {
            uniproc_ = true;
        }
        spc_ = spc;

        read_priority();
    }

    void DonorSearcher::increase_overlap_thickness(int nlayer)
    {
        // Increase volume of region of mandatory receptors.
        // Only applicable to mreceptors near interog boundary.

        if (nlayer == 1) {
            return;
        }

        assert(false);

        assert(spc_->sp_.size() == 1);
        auto& sp = spc_->sp_.front();

        for (int i=0; i<sp.mesh_.size(); ++i)
        {
            for (int j=0; j<sp.mesh_.size(); ++j)
            {
                if (i == j) {
                    continue;
                }

                Mesh& active_mesh = sp.mesh_[i];
                Mesh& passive_mesh = sp.mesh_[j];
                const ADT& passive_cell_adt = sp.cell_adt(j);

                active_mesh.increase_overlap_thickness(nlayer, passive_cell_adt, passive_mesh);
            }
        }
    }

    void DonorSearcher::convert_receptor_to_hole()
    {
        for (auto& sp: spc_->sp_)
        {
            for (Mesh& m: sp.mesh_)
            {
                m.convert_receptor_to_hole();
            }
        }
    }

    void DonorSearcher::read_priority()
    {
        std::ifstream in;
        in.open("priority");

        if (!in.is_open())
        {
            return;
        }

        int meshtag = -1;
        int pri = -1;

        while (in >> meshtag >> pri)
        {
            for (auto& sp: spc_->sp_) {
                sp.set_mesh_priority(meshtag, pri);
            }
        }

        in.close();
    }

    void DonorSearcher::donor_search(int nassemble, bool pseudo3D)
    {
        //std::cout << comm_->rank() << " start ds" << std::endl;
        assert(!spc_->sp_.empty());

        if (profiler_ != nullptr) {profiler_->start("asm-ds-reset-oga");}
        for (SpatialPartition& _sp: spc_->sp_) {
            reset_oga_status(_sp);
        }
        if (profiler_ != nullptr) {profiler_->stop("asm-ds-reset-oga");}

        if (profiler_ != nullptr) {profiler_->start("asm-ds-make-adt");}
        for (SpatialPartition& _sp: spc_->sp_)
        {
            if (_sp.mesh().size() <= 1) {
                break;
            }
            _sp.make_cell_adt();
        }
        if (profiler_ != nullptr) {profiler_->stop("asm-ds-make-adt");}

        if (profiler_ != nullptr) {profiler_->bstart("asm-ds-holeaabb");}
        auto hole_aabb = get_hole_aabb(*comm_, *spc_);
        if (profiler_ != nullptr) {profiler_->bstop("asm-ds-holeaabb");}

        //std::cout << comm_->rank() << " before holemap" << std::endl;

        if (profiler_ != nullptr) {profiler_->start("asm-ds-hole-map");}
        std::deque<HoleMap> holemap;
        for (SpatialPartition& _sp: spc_->sp_)
        {
            //if (_sp.mesh().size() <= 1) {
                //break;
            //}

            const Mesh* msh = nullptr;

            //for (const Mesh& m: _sp.mesh()) {
            for (int i=0; i<spc_->mesh_system_size(); ++i) {
                if (i < _sp.mesh().size())
                {
                    msh = &_sp.mesh()[i];
                }
                //holemap.push_back(HoleMap(*comm_, m, pseudo3D, hole_aabb));
                holemap.push_back(HoleMap(*comm_, msh, pseudo3D, hole_aabb));
                //std::cout << "min(0): " << i << " " << hole_aabb[i].min(0) << std::endl;
                //std::cout << "min(1): " << i << " " << hole_aabb[i].min(1) << std::endl;
                //std::cout << "max(0): " << i << " " << hole_aabb[i].max(0) << std::endl;
                //std::cout << "max(1): " << i << " " << hole_aabb[i].max(1) << std::endl;
            }
        }
        //std::cout << comm_->rank() << " after holemap" << std::endl;
        if (profiler_ != nullptr) {profiler_->stop("asm-ds-hole-map");}

        if (profiler_ != nullptr) {profiler_->start("asm-ds-mesh-aabb");}
        for (SpatialPartition& _sp: spc_->sp_) {
            if (_sp.mesh().size() <= 1) {
                break;
            }
            _sp.make_mesh_aabb();
        }
        if (profiler_ != nullptr) {profiler_->stop("asm-ds-mesh-aabb");}

        int nmesh = 0;
        int nsp = 0;
        int nwall = 0;
        int nsymmetry = 0;
        int ndirichlet = 0;
        int nfarfield = 0;
        int nempty = 0;
        int ninterog = 0;
        int nsw = 0;
        int adthit = 0;
        //std::vector<double> ds;
        int ncell = 0;
        int totaliter = 0;
        int nseed = 0;
        int adtsearch = 0;
        int swsearch = 0;
        int min_iter = TAILOR_BIG_POS_NUM;
        int max_iter = TAILOR_BIG_NEG_NUM;

        if (profiler_ != nullptr) {profiler_->start("asm-ds-nonresi");}
        for (SpatialPartition& _sp: spc_->sp_) {
            determine_nonresidents(_sp);
        }
        if (profiler_ != nullptr) {profiler_->stop("asm-ds-nonresi");}

        if (profiler_ != nullptr) {profiler_->start("asm-ds-mixed");}
        for (SpatialPartition& _sp: spc_->sp_) {
            if (_sp.mesh().size() <= 1) {
                break;
            }
            donor_search_mixed(_sp, holemap, adthit, nassemble, adtsearch, swsearch, min_iter, max_iter, hole_aabb);
        }
        if (profiler_ != nullptr) {profiler_->stop("asm-ds-mixed");}

        //if (profiler_ != nullptr) {profiler_->start("ds-mrecep");}
        //if (!master() || uniproc_)
        //{
            for (SpatialPartition& _sp: spc_->sp_) {
                donor_search_adt_mreceptor(_sp, holemap, hole_aabb);
            }
        //}
        //if (profiler_ != nullptr) {profiler_->stop("ds-mrecep");}

        //if (profiler_ != nullptr) {profiler_->start("handle-wall");}
        //for (auto& sp: spc_->sp_) {
            //handle_wall_contained_meshes(sp, holemap);
        //}
        //if (profiler_ != nullptr) {profiler_->stop("handle-wall");}

        //if (profiler_ != nullptr) {profiler_->start("convert");}
        //for (auto& sp: spc_->sp_)
        //{
            //convert_undefined_to_nonresident(sp);
            //convert_undefined_to_field(sp);
        //}
        //if (profiler_ != nullptr) {profiler_->stop("convert");}

        nsp += spc_->sp().size();
        for (const auto& sp: spc_->sp_)
        {
            nmesh += sp.mesh().size();
            //auto myoutline = std::find_if(spc_->outline().begin(), spc_->outline().end(), [&](const Outline& ol){return ol.tag()() == world_.rank();});
            //assert(myoutline != spc_->outline().end());
            for (const Mesh& m: sp.mesh_)
            {
                ncell += m.cell().size();
                //for (const MeshCell& mc: m.cell()) // remove later to gain speed.
                //{
                //if (myoutline->do_contain(mc.polygon().centroid(), false)) {
                //++ncell;
                //}
                //}
                nwall += m.wall_boundaries().size();
                nsymmetry += m.symmetry_boundaries().size();
                ndirichlet += m.dirichlet_boundaries().size();
                nfarfield += m.farfield_boundaries().size();
                nempty += m.empty_boundaries().size();
                ninterog += m.interog_boundaries().size();
            }
        }

        if (print_info_)
        {
            std::ofstream out;
            std::string fn = "ds-";
            fn.append(std::to_string(comm_->rank()));
            fn.append("-iter-");
            fn.append(std::to_string(nassemble));
            fn.append(".csv");
            out.open(fn);
            out << "nsp: " << nsp << std::endl;
            out << "nmesh: " << nmesh << std::endl;
            out << "ncell: " << ncell << std::endl;
            out << "nwall: " << nwall << std::endl;
            out << "nsymmetry: " << nsymmetry << std::endl;
            out << "ndirichlet: " << ndirichlet << std::endl;
            out << "nfarfield: " << nfarfield << std::endl;
            out << "nempty: " << nempty << std::endl;
            out << "ninterog: " << ninterog << std::endl;
            out << "nsw: " << nsw << std::endl;
            out << "adthit: " << adthit << std::endl;
            out << "spsize: " << spc_->sp_.size() << std::endl;
            out << "totaliter: " << totaliter << std::endl;
            out << "nseed: " << totaliter << std::endl;
            out << "adtsearch: " << adtsearch << std::endl;
            out << "swsearch: " << swsearch << std::endl;
            out << "min_iter: " << min_iter << std::endl;
            out << "max_iter: " << max_iter << std::endl;
            //int nc = 0;
            //for (const auto& sp: spc_->sp_)
            //{
            //    //out << "sptag: " << sp.tag()() << std::endl;
            //    //out << "meshsize: " << sp.mesh_.size() << std::endl;
            //    for (const Mesh& m: sp.mesh_)
            //    {
            //        if (sp.mesh_.size() != 1)
            //        {
            //            nc += m.cell().size();
            //        }
            //        //out << "cellsize: " << m.cell().size() << std::endl;
            //    }
            //}
            //out << "totalncell: " << nc << std::endl;
            ////double dsall = 0.;
            //for (int i=0; i<ds.size(); ++i)
            //{
            //    //dsall += ds[i];
            //    out << "ds: " << ds[i] << std::endl;
            //}
            ////out << "ds: " << ds << std::endl;
            out.close();
        }
    }

    void set_type_no_overlap(Mesh& active_mesh, const HoleMap* passivehm, const MeshCell& mc)
    {
        const Vector3& target = mc.poly().centroid();
        const Tag& ct = mc.tag();
        auto cell_type = mc.oga_cell_type();

        if (passivehm != nullptr && passivehm->is_inside_holebin(target))
        {
            active_mesh.set_as_hole(ct);
        }
        else
        {
            if (cell_type != OGA_cell_type_t::receptor && cell_type != OGA_cell_type_t::mandat_receptor)
            {
                active_mesh.set_as_field(ct);
            }
        }
    }


    void ds_adt_per_cell(Mesh& active_mesh, const Mesh& passive_mesh, const ADT& passive_cell_adt, const HoleMap* passivehm, const MeshCell& mc, const SpatialPartition& sp, int& adthit)
    {
        ++adthit;
        const Tag& ct = mc.tag();
        const Vector3& target = mc.poly().centroid();
        ADTPoint targetadt(target, mc.tag()());
        std::vector<int> res = passive_cell_adt.search(targetadt);

        if (res.empty())
        {
            set_type_no_overlap(active_mesh, passivehm, mc);
            return;
        }


        bool found = false;
        for (int r: res)
        {
            const MeshCell& passive_cell = passive_mesh.cell(Tag(r));

        //if (active_mesh.tag()() == 3) {
            //if (mc.tag()() == 4488) {
                //std::cout << "type: " << static_cast<int>(mc.oga_cell_type()) << std::endl;
                //std::cout << "res size: " << res.size() << std::endl;
                //std::cout << "inter: " << passive_cell.poly().do_intersect(target) << std::endl;
            //}
            //assert(active_mesh.cell(Tag(4488)).oga_cell_type() != OGA_cell_type_t::field);
        //}
            if (passive_cell.poly().do_intersect(target))
            {
                active_mesh.set_cell_type(ct, &passive_cell, passive_mesh);
                found = true;
                //return;
            }
        }

        if (!found) {
            set_type_no_overlap(active_mesh, passivehm, mc);
        }

        //if (active_mesh.oga_cell_type() == OGA_cell_type_t::receptor || active_mesh.oga_cell_type() == OGA_cell_type_t::mandat_receptor)
        //{
            //assert(passive_cell.oga_cell_type() == OGA_cell_type_t::field);
        //}
    }

    void ds_adt(Mesh& active_mesh, const Mesh& passive_mesh, const ADT& passive_cell_adt, const HoleMap* passivehm, const SpatialPartition& sp, int rank, int& adthit)
    {
        for (const MeshCell& mc: active_mesh.cell())
        {
            //if (mc.oga_cell_type() == OGA_cell_type_t::hole || !sp.is_resident(mc.poly().centroid()))
            if (mc.oga_cell_type() == OGA_cell_type_t::hole || mc.oga_cell_type() == OGA_cell_type_t::non_resident || mc.oga_cell_type() == OGA_cell_type_t::ghost)
            {
                continue;
            }

            ds_adt_per_cell(active_mesh, passive_mesh, passive_cell_adt, passivehm, mc, sp, adthit);

            //if (active_mesh.tag()() == 0)
            //{
                //if (mc.tag()() == 50)
                //{
                    //std::cout << "oga type: " << static_cast<int>(mc.oga_cell_type()) << std::endl;
                    //assert(false);
                //}
            //}
        }
    }

    ADT make_cell_adt(const Mesh& mesh)
    {
        //std::vector<ADTPoint> adt_pts;
        std::deque<ADTPoint> adt_pts;
        //adt_pts.reserve(mesh.cell().size());
        for (const MeshCell& mc: mesh.cell())
        {
            std::vector<Point> pts;
            pts.reserve(mc.point().size());
            for (auto& p: mc.point()) {
                pts.push_back(p.p());
            }
            adt_pts.push_back(ADTPoint(pts.begin(), pts.end(), mc.tag()()));
        }
        if (adt_pts.empty())
        {
            return ADT();
        }
        return ADT(adt_pts);
    }

    ADT make_wall_adt(const Mesh& mesh)
    {
        //std::vector<ADTPoint> adt_pts;
        std::deque<ADTPoint> adt_pts;
        //adt_pts.reserve(mesh.wall_boundaries().size());
        for (const MeshCell& mc: mesh.wall_boundaries())
        {
            //adt_pts.push_back(ADTPoint(mc.point().raw_points(), mc.tag()()));
            std::vector<Point> pts;
            pts.reserve(mc.point().size());
            for (auto& p: mc.point()) {
                pts.push_back(p.p());
            }
            adt_pts.push_back(ADTPoint(pts.begin(), pts.end(), mc.tag()()));
        }
        if (adt_pts.empty())
        {
            return ADT();
        }
        return ADT(adt_pts);
    }

    ADT make_symmetry_adt(const Mesh& mesh)
    {
        //std::vector<ADTPoint> adt_pts;
        std::deque<ADTPoint> adt_pts;
        //adt_pts.reserve(mesh.symmetry_boundaries().size());
        for (const MeshCell& mc: mesh.symmetry_boundaries())
        {
            //adt_pts.push_back(ADTPoint(mc.point().raw_points(), mc.tag()()));
            std::vector<Point> pts;
            pts.reserve(mc.point().size());
            for (auto& p: mc.point()) {
                pts.push_back(p.p());
            }
            adt_pts.push_back(ADTPoint(pts.begin(), pts.end(), mc.tag()()));
        }
        if (adt_pts.empty())
        {
            return ADT();
        }
        return ADT(adt_pts);
    }

    ADT make_dirichlet_adt(const Mesh& mesh)
    {
        //std::vector<ADTPoint> adt_pts;
        std::deque<ADTPoint> adt_pts;
        //adt_pts.reserve(mesh.dirichlet_boundaries().size());
        for (const MeshCell& mc: mesh.dirichlet_boundaries())
        {
            //adt_pts.push_back(ADTPoint(mc.point().raw_points(), mc.tag()()));
            std::vector<Point> pts;
            pts.reserve(mc.point().size());
            for (auto& p: mc.point()) {
                pts.push_back(p.p());
            }
            adt_pts.push_back(ADTPoint(pts.begin(), pts.end(), mc.tag()()));
        }
        if (adt_pts.empty())
        {
            return ADT();
        }
        return ADT(adt_pts);
    }

    const MeshCell* get_starting_cell(const Mesh& mesh, const RegularMesh& rm, const Bin& bin, HoleCutType hole_cut_type)
    {
        // A mesh cell in bin will be tried to be found.
        // If the bin is empty, neighbor bins will be checked out.
        // If even neighbors are empty, any cell will be returned.
        // Returned cell is used as starting cell for stencil walk algorithm.
        // For implicit hole_cut_type, starting cell is a must.
        // For direct cut, absense of starting cell means that the cell is in hole, so nullptr will be returned.

        Tag tag;
        if (bin.cell().empty())
        {
            if (hole_cut_type == HoleCutType::implicit)
            {
                tag = closest_non_empty_bin(bin.index(), rm, mesh);

                // worst case: get any cell in mesh.
                if (!tag.isvalid())
                {
                    for (const Bin& b: rm.bin())
                    {
                        if (!b.cell().empty())
                        {
                            tag = b.cell().front().cell();
                        }
                    }
                }
            }
        }
        else
        {
            for (const BinCell& bc: bin.cell())
            {
                tag = bc.cell();
                break;
            }
        }

        if (!tag.isvalid())
        {
            return nullptr;
        }
        return &(mesh.cell(tag));
    }

    const MeshCell* get_prev_seed(const Mesh& passive_mesh, const MeshCell& mc)
    {
        for (const Donor& donor: mc.prev_cand_donor())
        {
            if (donor.mesh_tag_ != passive_mesh.tag()) {continue;}
            const MeshCell* starting_cell = passive_mesh.query(donor.cell_tag_);
            if (starting_cell == nullptr) {
                continue;
            }
            return starting_cell;
        }

        return nullptr;
    }

    /*void DonorSearcher::apply_stencil_walk_to_bincells(const Bin& active_bin, const RegularMesh& active_rm, const RegularMesh& passive_rm, Mesh& active_mesh, const Mesh& passive_mesh, HoleCutType hole_cut_type, size_t& nsw, int dummyrank, const std::deque<Mesh>& wall, const std::deque<Mesh>& dirichlet, const HoleMap& holemap, bool use_past_seeds, int& totaliter, int& nseed)
    {
        int dummy;
        std::vector<BinRMTag> tag;
        passive_rm.get_bintag_adaptive(active_bin.aabb(), tag, dummy);
        for (const BinRMTag& brmt: tag)
        {
            const Bin& passive_bin = passive_rm.bin(brmt.bintag());

            for (const BinCell& bc: active_bin.cell())
            {
                Tag ct = bc.cell();
                Tag mt = bc.mesh();
                assert(ct.isvalid());
                assert(mt.isvalid());
                assert(active_mesh.query(ct));
                MeshCell& mc = active_mesh.cell_p(ct);
                const Point& target = mc.poly().centroid();
                auto cell_type = mc.oga_cell_type();
                //if (cell_type == OGA_cell_type_t::hole ||
                //cell_type == OGA_cell_type_t::receptor ||
                //cell_type == OGA_cell_type_t::mandat_receptor) {
                //continue;
                //}
                if (cell_type == OGA_cell_type_t::hole) {
                    continue;
                }

                bool is_resident = active_bin.is_resident(target.r());
                if (!is_resident) {
                    continue;
                }

                const MeshCell* starting_cell = nullptr;

                bool prevseedfound = false;
                if (use_past_seeds)
                {
                    for (int i=0; i<mc.prev_cand_donor_mesh().size(); ++i)
                    {
                        if (mc.prev_cand_donor_mesh()[i] != passive_mesh.tag()) {continue;}
                        if (passive_mesh.query(mc.prev_cand_donor_cell()[i]) == nullptr) {continue;}
                        starting_cell = &passive_mesh.cell(mc.prev_cand_donor_cell()[i]);
                        break;
                        //mc.reset_cand_donor(mc.prev_cand_donor_mesh()[i]);
                    }
                }

                if (starting_cell == nullptr) {
                    starting_cell = get_starting_cell(passive_mesh, passive_rm, passive_bin, hole_cut_type);
                }
                else
                {
                    ++nseed;
                    prevseedfound = true;
                }

                // starting cell is not mandatory for direct cut.
                if (starting_cell == nullptr)
                {
                    if (hole_cut_type_ == HoleCutType::holemap)
                    {
                        if (holemap.is_inside_holebin(target.r()))
                        {
                            active_mesh.set_as_hole(ct);
                        }
                        else
                        {
                            if (cell_type != OGA_cell_type_t::receptor && cell_type != OGA_cell_type_t::mandat_receptor)
                            {
                                active_mesh.set_as_field(ct);
                            }
                        }
                    }
                    else if (hole_cut_type_ == HoleCutType::direct)
                    {
                        if (cell_type != OGA_cell_type_t::receptor && cell_type != OGA_cell_type_t::mandat_receptor)
                        {
                            active_mesh.set_as_field(ct);
                        }
                    }
                    continue;
                }

                assert(passive_mesh.query(starting_cell->tag()));

                ++nsw;

                const MeshCell* closest_cell;
                int iter = 0;
                StencilWalkResult result = stencil_walk(passive_mesh, target, *starting_cell, closest_cell, dummyrank, iter);
                totaliter += iter;
                if (iter > 10 && prevseedfound)
                {
                    if (settings_->dsspat_ == "rm")
                    {
                        std::cout << "rank: " << dummyrank << std::endl;
                        std::cout << "mc: " << mc.tag()() << std::endl;
                        std::cout << "cnt 0: " << mc.poly().centroid()(0) << std::endl;
                        std::cout << "cnt 1: " << mc.poly().centroid()(1) << std::endl;
                        std::cout << "cnt 2: " << mc.poly().centroid()(2) << std::endl;
                        std::cout << "starting cell: " << starting_cell->tag()() << std::endl;
                        std::string fn = "passive";
                        fn.append(std::to_string(dummyrank));
                        fn.append(".vtk");
                        passive_mesh.print_as_vtk(fn);
                        for (int i=0; i<mc.prev_cand_donor_mesh().size(); ++i)
                        {
                            std::cout << "prev cand mesh: " << mc.prev_cand_donor_mesh()[i]() << std::endl;
                            std::cout << "prev cand cell: " << mc.prev_cand_donor_cell()[i]() << std::endl;
                        }
                        assert(false);
                    }
                }

                if (result != StencilWalkResult::inside_cell)
                {
                    if (hole_cut_type_ == HoleCutType::holemap)
                    {
                        if (holemap.is_inside_holebin(target.r()))
                        {
                            active_mesh.set_as_hole(ct);
                        }
                        else
                        {
                            if (cell_type != OGA_cell_type_t::receptor && cell_type != OGA_cell_type_t::mandat_receptor)
                            {
                                active_mesh.set_as_field(ct);
                            }
                        }
                    }
                    else if (hole_cut_type_ == HoleCutType::direct)
                    {
                        if (cell_type != OGA_cell_type_t::receptor && cell_type != OGA_cell_type_t::mandat_receptor)
                        {
                            active_mesh.set_as_field(ct);
                        }
                    }
                }
                else
                {
                    active_mesh.set_cell_type(ct, closest_cell);
                }
            }
        }
    }*/

    void ds_stencil_walk(Mesh& active_mesh, const Mesh& passive_mesh, const MeshCell& mc, const MeshCell& starting_cell, int dummyrank, const HoleMap* passivehm, int& min_iter, int& max_iter)
    {
        const Tag& ct = mc.tag();
        const MeshCell* closest_cell;
        int iter = 0;
        const Vector3& target = mc.poly().centroid();
        StencilWalkResult result = stencil_walk(passive_mesh, target, starting_cell, closest_cell, dummyrank, iter, false);
        min_iter = std::min(min_iter, iter);
        max_iter = std::max(max_iter, iter);

        if (result != StencilWalkResult::inside_cell)
        {
            set_type_no_overlap(active_mesh, passivehm, mc);
            return;
        }

        active_mesh.set_cell_type(ct, closest_cell, passive_mesh);
    }

    /*std::deque<ADT> DonorSearcher::make_adt(SpatialPartition& sp)
      {
      std::deque<ADT> cell_adt;

      if (master() && !uniproc_)
      {return cell_adt;}

      for (Mesh& m: sp.mesh_) {
      cell_adt.push_back(make_cell_adt(m));
      }

      return cell_adt;
      }

      std::deque<AABB> DonorSearcher::make_aabb(SpatialPartition& sp)
      {
      std::deque<AABB> aabbs;
      for (Mesh& m: sp.mesh_) {
      aabbs.push_back(AABB(m.point()));
      }

      return aabbs;
      }*/

    void DonorSearcher::handle_wall_contained_meshes(SpatialPartition& sp, const std::deque<HoleMap>& holemap, const std::vector<AABB>& hole_aabb)
    {
        auto aabbs = sp.mesh_aabb();
        for (const auto& hm: holemap)
        {
            auto msh = std::find_if(sp.mesh().begin(), sp.mesh().end(), [&](const Mesh& m){return m.tag() == hm.tag();});
            if (msh != sp.mesh().end()) {
                continue;
            }

            for (int i=0; i<sp.mesh_.size(); ++i)
            {
                Mesh& activemesh = sp.mesh_[i];

                //if (hm.holeless()) {
                    //continue;
                //}

                //if (hm.wallaabb().do_contain(aabbs[i]))
                if (hole_aabb[activemesh.tag()()].do_contain(aabbs[i], false))
                {
                    for (MeshCell& mc: activemesh.cell_)
                    {
                        activemesh.set_as_hole(mc.tag());
                    }
                }
            }
        }
    }

    void DonorSearcher::donor_search_adt_mreceptor(SpatialPartition& sp, const std::deque<HoleMap>& holemap, const std::vector<AABB>& hole_aabb)
    {
        assert(!sp.mesh_.empty());
        for (int i=0; i<sp.mesh_.size(); ++i)
        {
            assert(!sp.mesh_[i].cell().empty());
        }

        //auto cell_adt = make_adt(sp);
        //auto aabbs = make_aabb(sp);

        for (int i=0; i<sp.mesh_.size(); ++i)
        {
            Mesh& active_mesh = sp.mesh_[i];
            const AABB& active_aabb = sp.mesh_aabb(i);

            for (int j=0; j<sp.mesh_.size(); ++j)
            {
                if (i == j) {
                    continue;
                }

                const Mesh& passive_mesh = sp.mesh_[j];
                const AABB& passive_aabb = sp.mesh_aabb(j);
                const ADT& passive_cell_adt = sp.cell_adt(j);
                auto iterpassivehm = std::find_if(holemap.begin(), holemap.end(), [&](const HoleMap& hm){return hm.tag() == sp.mesh_[j].tag();});
                const HoleMap* passivehm = nullptr;
                if (iterpassivehm != holemap.end())
                {
                    passivehm = &(*iterpassivehm);
                }
                //assert(active_mesh.parent_mesh() != passive_mesh.parent_mesh());

                // if active and passive meshes do not overlap by any means then continue.
                if (!active_aabb.do_intersect(passive_aabb))
                {
                    //if (hole_cut_type_ == HoleCutType::holemap)
                    {
                        //if (passivehm != nullptr && !passivehm->holeless() && !active_aabb.do_intersect(passivehm->wallaabb())) {
                        //if (passivehm != nullptr && !active_aabb.do_intersect(passivehm->wallaabb())) {
                        if (passivehm != nullptr && !active_aabb.do_intersect(hole_aabb[passivehm->tag()()])) {
                            continue;
                        }
                    }
                    //else if (hole_cut_type_ == HoleCutType::direct)
                    //{
                        //if (!active_aabb.do_intersect(passive_mesh.hole_aabb())) {
                        //continue; // if active_aabb intersects hole of passive, direct cut before this function already marked hole cells.
                        //}
                    //}
                }

                for (MeshCell& mc: active_mesh.cell_)
                {
                    auto cell_type = mc.oga_cell_type();

                    if (cell_type == OGA_cell_type_t::hole)
                    {
                        for (const auto& ct: mc.pnei())
                        {
                            const MeshCell& nei = active_mesh.cell(ct);
                            auto neitype = nei.oga_cell_type();
                            if (neitype != OGA_cell_type_t::hole && neitype != OGA_cell_type_t::mandat_receptor)
                            {
                                const Point& target = nei.poly().centroid();

                                ADTPoint targetadt(target.r(), nei.tag()());
                                std::vector<int> res = passive_cell_adt.search(targetadt);

                                for (int r: res)
                                {
                                    const MeshCell& passive_cell = passive_mesh.cell(Tag(r));
                                    if (passive_cell.poly().do_intersect(target.r()))
                                    {
                                        //if (passive_cell.oga_cell_type() != OGA_cell_type_t::field) {
                                            //continue;
                                        //}

                                        mc.set_oga_cell_type(OGA_cell_type_t::mandat_receptor);
                                        mc.add_cand_donor(passive_cell.parent_mesh(), passive_cell.tag(), &passive_cell);
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }

                //if (world_.rank() == 3 && active_mesh.tag()() == 0) {
                //std::cout << "zype: " << static_cast<int>(sp.mesh(Tag(0)).cell(Tag(27420)).oga_cell_type()) << std::endl;
                //}
            }
        }
    }

    void getdiff(double start, std::vector<double>& ds, int rank, int i)
    {
        double stop = MPI_Wtime();
        double elapsed = stop - start;

        if (i >= ds.size())
        {
            ds.push_back(elapsed);
        }
        else
        {
            ds[i] += elapsed;
        }
    }

    bool meshes_overlap(const AABB& active_aabb, const AABB& passive_aabb)
    {
        if (!active_aabb.do_intersect(passive_aabb))
        {
            return false;
        }

        return true;
    }

    bool mesh_contained_by_hole(const HoleMap* passivehm, const AABB& active_aabb, Mesh& active_mesh, const std::vector<AABB>& hole_aabb)
    {
        //if (passivehm != nullptr && !passivehm->holeless() && passivehm->wallaabb().do_contain(active_aabb))
        //if (passivehm != nullptr && passivehm->wallaabb().do_contain(active_aabb))
        if (passivehm != nullptr && hole_aabb[passivehm->tag()()].do_contain(active_aabb, false))
        {
            for (const MeshCell& mc: active_mesh.cell())
            {
                if (mc.oga_cell_type() == OGA_cell_type_t::non_resident || mc.oga_cell_type() == OGA_cell_type_t::ghost) {
                    continue;
                }
                active_mesh.set_as_hole(mc.tag());
            }

            return true;
        }

        return false;
    }

    void DonorSearcher::donor_search_mixed(SpatialPartition& sp, const std::deque<HoleMap>& holemap, int& adthit, int nassemble, int& adtsearch, int& swsearch, int& min_iter, int& max_iter, const std::vector<AABB>& hole_aabb)
    {
        if (sp.mesh_.size() == 1) {
            return;
        }

        for (int i=0; i<sp.mesh_.size(); ++i)
        {
            assert(!sp.mesh_[i].cell().empty());
        }

        for (int i=0; i<sp.mesh_.size(); ++i)
        {
            for (int j=0; j<sp.mesh_.size(); ++j)
            {
                if (i == j) {
                    continue;
                }

                Mesh& active_mesh = sp.mesh_[i];
                const Mesh& passive_mesh = sp.mesh_[j];
                const AABB& active_aabb = sp.mesh_aabb(i);
                const AABB& passive_aabb = sp.mesh_aabb(j);
                auto iterpassivehm = std::find_if(holemap.begin(), holemap.end(), [&](const HoleMap& hm){return hm.tag() == sp.mesh_[j].tag();});
                const HoleMap* passivehm = nullptr;
                if (iterpassivehm != holemap.end())
                {
                    passivehm = &(*iterpassivehm);
                }
                //assert(passivehm != holemap.end());
                //assert(active_mesh.parent_mesh() != passive_mesh.parent_mesh());
                const ADT& passive_cell_adt = sp.cell_adt(j);

                if (!meshes_overlap(active_aabb, passive_aabb)) {
                    continue;
                }

                //if (passivehm != holemap.end())
                if (passivehm != nullptr)
                {
                    if (mesh_contained_by_hole(passivehm, active_aabb, active_mesh, hole_aabb)) {
                        continue;
                    }
                }

                if (nassemble == 0)
                {
                    ds_adt(active_mesh, passive_mesh, passive_cell_adt, passivehm, sp, comm_->rank(), adthit);
                    continue;
                }

                for (const MeshCell& mc: active_mesh.cell())
                {
                    //if (mc.oga_cell_type() == OGA_cell_type_t::hole || !sp.is_resident(mc.poly().centroid()))
                    if (mc.oga_cell_type() == OGA_cell_type_t::hole || mc.oga_cell_type() == OGA_cell_type_t::non_resident || mc.oga_cell_type() == OGA_cell_type_t::ghost)
                    {
                        continue;
                    }

                    const MeshCell* starting_cell = nullptr;
                    if (dsa_ == DonorSearchAlgo::mixed)
                    {
                        starting_cell = get_prev_seed(passive_mesh, mc);
                    }

                    if (starting_cell == nullptr || starting_cell->oga_cell_type() == OGA_cell_type_t::ghost)
                    {
                        //std::cout << "rank: " << comm_->rank() << std::endl;
                        //std::cout << "m: " << mc.parent_mesh()() << std::endl;
                        //std::cout << "mc: " << mc.tag()() << std::endl;
                        //assert(false);
                        ds_adt_per_cell(active_mesh, passive_mesh, passive_cell_adt, passivehm, mc, sp, adthit);
                        ++adtsearch;
                        continue;
                    }

                    const MeshCell* closest_cell;
                    int iter = 0;
                    ds_stencil_walk(active_mesh, passive_mesh, mc, *starting_cell, comm_->rank(), passivehm, min_iter, max_iter);
                    ++swsearch;
                }
            }
        }
        //for (const auto& m: sp.mesh())
        //{
        //    if (m.tag()() == 3)
        //    {
        //        for (const auto& mc: m.cell())
        //        {
        //            if (mc.tag()() == 4488) {
        //                std::cout << "rank: " << comm_->rank() << std::endl;
        //                std::cout << "type: " << static_cast<int>(mc.oga_cell_type()) << std::endl;
        //            }
        //        }
        //    }
        //}
    }

    /*void DonorSearcher::donor_search_adt(SpatialPartition& sp, const std::deque<HoleMap>& holemap, size_t& adtmiss, size_t& adthit, std::vector<double>& ds)
    {
        if (master() && !uniproc_)
        {return;}

        assert(!sp.mesh_.empty());
        if (sp.mesh_.size() == 1) {
            return;
        }

        for (int i=0; i<sp.mesh_.size(); ++i)
        {
            assert(!sp.mesh_[i].cell().empty());
        }

        for (int i=0; i<sp.mesh_.size(); ++i)
        {

            for (int j=0; j<sp.mesh_.size(); ++j)
            {
                if (i == j) {
                    continue;
                }
                
                double start = MPI_Wtime();

                Mesh& active_mesh = sp.mesh_[i];
                const Mesh& passive_mesh = sp.mesh_[j];
                const AABB& active_aabb = sp.mesh_aabb(i);
                const AABB& passive_aabb = sp.mesh_aabb(j);
                auto passivehm = std::find_if(holemap.begin(), holemap.end(), [&](const HoleMap& hm){return hm.tag() == sp.mesh_[j].tag();});
                assert(passivehm != holemap.end());
                assert(active_mesh.parent_mesh() != passive_mesh.parent_mesh());
                const ADT& passive_cell_adt = sp.cell_adt(j);

                // if active and passive meshes do not overlap by any means then continue.
                if (!active_aabb.do_intersect(passive_aabb))
                {
                    if (hole_cut_type_ == HoleCutType::holemap)
                    {
                        if (!passivehm->holeless() && !active_aabb.do_intersect(passivehm->wallaabb())) {
                            getdiff(start, ds, world_.rank(), i);
                            continue;
                        }
                    }
                    else if (hole_cut_type_ == HoleCutType::direct)
                    {
                        //if (!active_aabb.do_intersect(passive_mesh.hole_aabb())) {
                            getdiff(start, ds, world_.rank(), i);
                            continue; // if active_aabb intersects hole of passive, direct cut before this function already marked hole cells.
                        //}
                    }
                }

                for (const MeshCell& mc: active_mesh.cell())
                {
                    const Tag& ct = mc.tag();
                    const Point& target = mc.poly().centroid();
                    auto cell_type = mc.oga_cell_type();

                    if (cell_type == OGA_cell_type_t::hole) {
                        continue;
                    }

                    //if (!myoutline->do_contain(mc.polygon().centroid(), false)) { // here is the problem. time costly for some procs than others.
                        //continue;
                    //}

                    ADTPoint targetadt(target.r(), mc.tag()());
                    std::vector<int> res = passive_cell_adt.search(targetadt);

                    assert(target.r(0) >= sp.mesh_aabb(i).min(0));
                    assert(target.r(1) >= sp.mesh_aabb(i).min(1));
                    assert(target.r(2) >= sp.mesh_aabb(i).min(2));

                    assert(target.r(0) <= sp.mesh_aabb(i).max(0));
                    assert(target.r(1) <= sp.mesh_aabb(i).max(1));
                    assert(target.r(2) <= sp.mesh_aabb(i).max(2));

                    if (res.empty())
                    {
                        if (hole_cut_type_ == HoleCutType::holemap)
                        {
                            if (passivehm->is_inside_holebin(target.r()))
                            {
                                active_mesh.set_as_hole(ct);
                            }
                            else
                            {
                                if (cell_type != OGA_cell_type_t::receptor && cell_type != OGA_cell_type_t::mandat_receptor)
                                {
                                    active_mesh.set_as_field(ct);
                                }
                            }
                        }
                        else if (hole_cut_type_ == HoleCutType::direct)
                        {
                            // because hole cell is already continued. if res is empty then active cell must be field.
                            if (cell_type != OGA_cell_type_t::receptor && cell_type != OGA_cell_type_t::mandat_receptor)
                            {
                                active_mesh.set_as_field(ct);
                            }
                        }

                        ++adtmiss;
                    }
                    else {
                        ++adthit;
                    }

                    for (int r: res)
                    {
                        const MeshCell& passive_cell = passive_mesh.cell(Tag(r));
                        if (passive_cell.poly().do_intersect(target.r()));
                        {
                            active_mesh.set_cell_type(ct, &passive_cell);
                            break;
                        }
                    }
                }

                getdiff(start, ds, world_.rank(), i);
            }
        }
    }*/

    /*void DonorSearcher::donor_search_rm(SpatialPartition& sp, const std::deque<Mesh>& wall, const std::deque<Mesh>& dirichlet, const std::deque<HoleMap>& holemap, size_t& nsw, std::vector<double>& ds, int& totaliter, int& nseed)
    {
        if (master() && !uniproc_)
        {return;}

        assert(!sp.mesh_.empty());
        if (sp.mesh_.size() == 1) {
            return;
        }

        for (int i=0; i<sp.mesh_.size(); ++i)
        {
            assert(!sp.mesh_[i].cell().empty());
        }

        for (int i=0; i<sp.mesh_.size(); ++i)
        {
            for (int j=0; j<sp.mesh_.size(); ++j)
            {
                if (i == j) {
                    continue;
                }

                auto start = MPI_Wtime();

                Mesh& active_mesh = sp.mesh_[i];
                const Mesh& passive_mesh = sp.mesh_[j];
                const AABB& active_aabb = sp.mesh_aabb(i);
                const AABB& passive_aabb = sp.mesh_aabb(j);
                auto passivehm = std::find_if(holemap.begin(), holemap.end(), [&](const HoleMap& hm){return hm.tag() == sp.mesh_[j].tag();});
                assert(passivehm != holemap.end());
                assert(active_mesh.parent_mesh() != passive_mesh.parent_mesh());
                const RegularMesh& active_rm = sp.rm_[i];
                const RegularMesh& passive_rm = sp.rm_[j];

                if (settings_->dsspat_ == "rm")
                {
                    if (!passivehm->holeless() && passivehm->wallaabb().do_contain(active_rm.aabb()))
                    {
                        for (const Bin& active_bin: active_rm.bin())
                        {
                            for (const BinCell& bc: active_bin.cell())
                            {
                                Tag ct = bc.cell();
                                active_mesh.set_as_hole(ct);
                            }
                        }
                        getdiff(start, ds, world_.rank(), i);
                        continue;
                    }
                }
                
                // if active and passive meshes do not overlap by any means then continue.
                if (!active_aabb.do_intersect(passive_aabb))
                {
                    if (hole_cut_type_ == HoleCutType::holemap)
                    {
                        if (!passivehm->holeless() && !active_aabb.do_intersect(passivehm->wallaabb())) {
                            getdiff(start, ds, world_.rank(), i);
                            continue;
                        }
                    }
                    else if (hole_cut_type_ == HoleCutType::direct)
                    {
                        //if (!active_aabb.do_intersect(passive_mesh.hole_aabb())) {
                            getdiff(start, ds, world_.rank(), i);
                            continue; // if active_aabb intersects hole of passive, direct cut before this function already marked hole cells.
                        //}
                    }
                }

                for (const Bin& active_bin: active_rm.bin())
                {
                    if (hole_cut_type_ == HoleCutType::holemap)
                    {
                        if (!passivehm->holeless() && passivehm->wallaabb().do_contain(active_bin.aabb()))
                        {
                            for (const BinCell& bc: active_bin.cell())
                            {
                                Tag ct = bc.cell();
                                active_mesh.set_as_hole(ct);
                            }
                            continue;
                        }
                    }
                    else if (hole_cut_type_ == HoleCutType::direct)
                    {
                        if (active_bin.aabb().do_intersect(passive_mesh.hole_aabb())) {
                            continue; // since cells inside the bin should already be marked as hole cells in direct function.
                        }
                    }
                    if (!active_bin.aabb().do_intersect(passive_rm.aabb())) {
                        continue;
                    }

                    apply_stencil_walk_to_bincells(active_bin, active_rm, passive_rm, active_mesh, passive_mesh, hole_cut_type_, nsw, world_.rank(), wall, dirichlet, *passivehm, settings_->use_past_seeds_, totaliter, nseed);
                }

                getdiff(start, ds, world_.rank(), i);
            }
        }
    }*/

    void DonorSearcher::convert_undefined_to_field()
    {
        for (auto& sp: spc_->sp_)
        {
            for (Mesh& m: sp.mesh_)
            {
                m.convert_undefined_to_field();
            }
        }
    }

    void DonorSearcher::determine_nonresidents(SpatialPartition& sp)
    {
        if (uniproc_) {
            return;
        }

        for (Mesh& m: sp.mesh_)
        {
            for (MeshCell& mc: m.cell_)
            {
                if (mc.oga_cell_type() == OGA_cell_type_t::ghost) {
                    continue;
                }

                if (mergebins_)
                {
                    int celltag = -1;
                    Tag brmt;

                    if (!spc_->is_resident(mc.poly().centroid(), celltag, brmt))
                    {
                        mc.set_oga_cell_type(OGA_cell_type_t::non_resident);
                    }
                }
                else
                {
                    assert(false);
                    //if (!sp.is_resident(mc.poly().centroid()))
                    {
                        mc.set_oga_cell_type(OGA_cell_type_t::non_resident);
                    }
                }

            }
        }
    }

    void DonorSearcher::convert_undefined_to_nonresident(SpatialPartition& sp)
    {
        assert(false);
        if (uniproc_) {
            return;
        }

        //std::vector<Outline>::const_iterator myoutline;

        //if (mergebins_)
        {
            //assert(false);
            //if (!spc_->outline().empty()) // this is another way of checking for merge_bins option.
            //{
                //myoutline = std::find_if(spc_->outline().begin(), spc_->outline().end(), [&](const Outline& ol){return ol.tag()() == world_.rank();});
                //assert(myoutline != spc_->outline().end());
            //}
        }

        for (Mesh& m: sp.mesh_)
        {
            for (MeshCell& mc: m.cell_)
            {
                //if (mergebins_)
                {
                    //assert(false);
                    //if (!spc_->outline().empty())
                    //{
                        //if (!myoutline->do_contain(mc.poly().centroid(), false))
                        //{
                            //mc.set_oga_cell_type(OGA_cell_type_t::non_resident);
                        //}
                    //}
                }
                //else
                {
                    //if (!sp.is_resident(mc.poly().centroid()))
                    {
                        mc.set_oga_cell_type(OGA_cell_type_t::non_resident);
                    }
                }
            }
        }
    }

    void DonorSearcher::reset_oga_status(SpatialPartition& sp)
    {
        for (Mesh& m: sp.mesh_) {
            m.reset_all_cell_oga_status();
        }
    }

    void DonorSearcher::check_donor_validity()
    {
        for (auto& sp: spc_->sp_)
        {
            for (Mesh& m: sp.mesh_)
            {
                for (MeshCell& mc: m.cell_)
                {
                    auto type = mc.oga_cell_type();
                    if (type == OGA_cell_type_t::receptor || type == OGA_cell_type_t::mandat_receptor)
                    {
                            auto mesh_it = std::find_if(sp.mesh_.begin(), sp.mesh_.end(), [&](const Mesh& m){return m.tag() == mc.donor().mesh_tag_;});
                            assert(mesh_it != sp.mesh_.end());
                            assert(mesh_it->tag() != m.tag());

                            auto cit = mesh_it->query(mc.donor().cell_tag_);
                            assert(cit != nullptr);

                            if (cit->oga_cell_type() != OGA_cell_type_t::field)
                            {
                                std::cout << "oga: " << static_cast<int>(cit->oga_cell_type()) << std::endl;
                            }
                            assert(cit->oga_cell_type() == OGA_cell_type_t::field);
                    }
                }
            }
        }
    }

    void DonorSearcher::convert_orphan_to_field(const ArrCon<DonorInfo2>& arrival)
    {
        assert(false);
        // if all donors of receptor are mandat receptor, instead of designating the cell as oprhan, just set to field.

        for (auto& sp: spc_->sp_)
        {
            for (Mesh& m: sp.mesh_)
            {
                for (MeshCell& mc: m.cell_)
                {
                    auto type = mc.oga_cell_type();
                    if (type == OGA_cell_type_t::receptor)
                    {
                        assert(!mc.cand_donor().empty());

                        for (auto& cand_donor: mc.cand_donor())
                        {
                            const auto& cdm = cand_donor.mesh_tag_;
                            const auto& cdc = cand_donor.cell_tag_;

                            auto mesh_it = std::find_if(sp.mesh_.begin(), sp.mesh_.end(), [&](const Mesh& m){return m.tag() == cdm;});
                            assert(mesh_it != sp.mesh_.end());
                            assert(mesh_it->tag() != m.tag());

                            auto cit = mesh_it->query(Tag(cdc));

                            OGA_cell_type_t type = OGA_cell_type_t::undefined;
                            assert(cit != nullptr);

                            //if (!sp.is_resident(cit->poly().centroid()))
                            {
                                //for (const auto& receiver: nonresi)
                                {
                                    //if (receiver.tag() == sp.tag()())
                                    {
                                        auto nr = std::find_if(arrival.begin(), arrival.end(), [&](const DonorInfo2& n){return (n.meshtag_ == cdm() && n.celltag_ == cdc());});
                                        assert(nr != arrival.end());
                                        type = nr->oga_cell_type_;
                                        //break;
                                    }
                                }
                            }
                            //else
                            //{
                            //    type = cit->oga_cell_type();
                            //}

                            assert(type != OGA_cell_type_t::undefined);

                            if (type == OGA_cell_type_t::mandat_receptor)
                            {
                                mc.set_oga_cell_type(OGA_cell_type_t::field);
                                //if (mc.parent_mesh()() == 1)
                                //{
                                    //assert(mc.tag()() != 174);
                                //}
                            }
                            break;
                        }
                    }
                }
            }
        }
    }

    void DonorSearcher::receptor_to_field(const ArrCon<DonorInfo2>& arrival)
    {
        for (auto& sp: spc_->sp_)
        {
            for (Mesh& m: sp.mesh_)
            {
                for (MeshCell& mc: m.cell_)
                {
                    auto type = mc.oga_cell_type();
                    if (type == OGA_cell_type_t::receptor)
                    {
                        assert(!mc.cand_donor().empty());

                        for (auto cand_donor = mc.cand_donor_.begin(); cand_donor != mc.cand_donor_.end();)
                        {
                            const auto& cdm = cand_donor->mesh_tag_;
                            const auto& cdc = cand_donor->cell_tag_;

                            auto mesh_it = std::find_if(sp.mesh_.begin(), sp.mesh_.end(), [&](const Mesh& m){return m.tag() == cdm;});
                            assert(mesh_it != sp.mesh_.end());
                            assert(mesh_it->tag() != m.tag());

                            auto cit = mesh_it->query(Tag(cdc));

                            OGA_cell_type_t type = OGA_cell_type_t::undefined;
                            assert(cit != nullptr);

                            //if (!sp.is_resident(cit->poly().centroid()))

                            int celltag;
                            Tag brmt;

                            if (!spc_->is_resident(cit->poly().centroid(), celltag, brmt))
                            {
                                //for (const auto& receiver: nonresi)
                                {
                                    //if (receiver.tag() == sp.tag()())
                                    {
                                        auto nr = std::find_if(arrival.begin(), arrival.end(), [&](const DonorInfo2& n){return (n.meshtag_ == cdm() && n.celltag_ == cdc());});
                                        if (nr == arrival.end())
                                        {
                                            std::cout << "rank: " << comm_->rank() << std::endl;
                                            std::cout << "cdm: " << cdm() << std::endl;
                                            std::cout << "cdc: " << cdc() << std::endl;
                                            std::cout << "arrival size: " << arrival.size() << std::endl;
                                            std::cout << "sp aabb min(0)" << sp.aabb().min(0) << std::endl;
                                            std::cout << "sp aabb min(1)" << sp.aabb().min(1) << std::endl;
                                            std::cout << "sp aabb min(2)" << sp.aabb().min(2) << std::endl;
                                            std::cout << "sp aabb max(0)" << sp.aabb().max(0) << std::endl;
                                            std::cout << "sp aabb max(1)" << sp.aabb().max(1) << std::endl;
                                            std::cout << "sp aabb max(2)" << sp.aabb().max(2) << std::endl;
                                            std::cout << "cnt(0)" << cit->poly().centroid()(0) << std::endl;
                                            std::cout << "cnt(1)" << cit->poly().centroid()(1) << std::endl;
                                            std::cout << "cnt(2)" << cit->poly().centroid()(2) << std::endl;
                                        }
                                        assert(nr != arrival.end());
                                        type = nr->oga_cell_type_;
                                        //break;
                                    }
                                }
                            }
                            else
                            {
                                type = cit->oga_cell_type();
                            }

                            assert(type != OGA_cell_type_t::undefined);

                            //if (type == OGA_cell_type_t::mandat_receptor)
                            if (type == OGA_cell_type_t::receptor || type == OGA_cell_type_t::mandat_receptor || type == OGA_cell_type_t::orphan)
                            {
                                mc.cand_donor_.erase(cand_donor);
                            }
                            else
                            {
                                ++cand_donor;
                            }
                        }

                        if (mc.cand_donor().empty())
                        {
                            mc.set_oga_cell_type(OGA_cell_type_t::field);
                        }
                    }
                }
            }
        }
    }

    /*bool no_valid_cand(const MeshCell& mc, const SpatialPartition& sp, const std::vector<Receiver<DonorInfo2>>& nonresi)
    {
        for (auto cand_donor = mc.cand_donor().begin(); cand_donor != mc.cand_donor().end(); ++cand_donor)
        {
            const auto& cdm = cand_donor->mesh_tag_;
            const auto& cdc = cand_donor->cell_tag_;

            auto mesh_it = std::find_if(sp.mesh_.begin(), sp.mesh_.end(), [&](const Mesh& m){return m.tag() == cdm;});
            assert(mesh_it != sp.mesh_.end());
            assert(mesh_it->tag() != m.tag());

            auto cit = mesh_it->query(Tag(cdc));

            OGA_cell_type_t type = OGA_cell_type_t::undefined;
            assert(cit != nullptr);

            if (cit->oga_cell_type() == OGA_cell_type_t::non_resident)
            {
                assert(!sp.is_resident(cit->poly().centroid()));
                for (const auto& receiver: nonresi)
                {
                    if (receiver.tag() == sp.tag()())
                    {
                        auto nr = std::find_if(receiver.arrival().begin(), receiver.arrival().end(), [&](const auto& n){return (n.meshtag_ == cdm() && n.celltag_ == cdc());});
                        assert(nr != receiver.arrival().end());
                        type = nr->oga_cell_type_;
                        break;
                    }
                }
            }
            else
            {
                type = cit->oga_cell_type();
            }

            assert(type != OGA_cell_type_t::undefined);

            if (type != OGA_cell_type_t::receptor && type != OGA_cell_type_t::mandat_receptor && type != OGA_cell_type_t::orphan)
            {
                return false;
            }
        }

        return true;
    }*/

    /*void force_to_field(const MeshCell& mc, SpatialPartition& sp, std::vector<const Nei*>& nei, const std::vector<Receiver<DonorInfo2>>& nonresi)
    {
        for (auto cand_donor = mc.cand_donor().begin(); cand_donor != mc.cand_donor().end();)
        {
            const auto& cdm = cand_donor->mesh_tag_;
            const auto& cdc = cand_donor->cell_tag_;

            auto mesh_it = std::find_if(sp.mesh_.begin(), sp.mesh_.end(), [&](const Mesh& m){return m.tag() == cdm;});
            assert(mesh_it != sp.mesh_.end());
            assert(mesh_it->tag() != m.tag());

            auto cit = mesh_it->query(Tag(cdc));

            OGA_cell_type_t type = OGA_cell_type_t::undefined;
            assert(cit != nullptr);

            if (cit->oga_cell_type() == OGA_cell_type_t::non_resident)
            {
                assert(!sp.is_resident(cit->poly().centroid()));
                for (const auto& receiver: nonresi)
                {
                    if (receiver.tag() == sp.tag()())
                    {
                        auto nr = std::find_if(receiver.arrival().begin(), receiver.arrival().end(), [&](const auto& n){return (n.meshtag_ == cdm() && n.celltag_ == cdc());});
                        assert(nr != receiver.arrival().end());
                        nei.push_back(&(*nr));
                        break;
                    }
                }
                break;
            }
            else
            {
                cit->set_oga_cell_type(OGA_cell_type_t::field);
                break;
            }
        }
    }*/

    void DonorSearcher::handle_donor_conflict(const ArrCon<DonorInfo2>& arrival)
    {
        TCRExchanger tcr_exc(spc_, comm_, &arrival);
        tcr_exc.exchange(false, "tcr", profiler_);

        if (profiler_ != nullptr) {profiler_->start("asm-ds-handle");}
        for (auto& sp: spc_->sp_)
        {
            //for (const auto& receiver: tcr_exc.receiver())
            {
                //if (receiver.tag() == sp.tag()())
                {
                    //for (const auto& n: receiver.arrival())
                    for (const auto& n: tcr_exc.arrival())
                    {
                        //auto mesh = std::find_if(sp.mesh_.begin(), sp.mesh_.end(), [&](const auto& m){return (m.tag() == n.cell_.parent_mesh());});
                        auto mesh = std::find_if(sp.mesh_.begin(), sp.mesh_.end(), [&](const Mesh& m){return (m.tag()() == n.mesh_cell_.first);});
                        if (mesh == sp.mesh_.end())
                        {
                            std::cout << "rank: " << comm_->rank() << std::endl;
                            std::cout << "sp: " << sp.tag()() << std::endl;
                            std::cout << "sp mesh size" << sp.mesh().size() << std::endl;
                            std::cout << "n mesh tag: " << n.mesh_cell_.first << std::endl;
                            std::cout << "n cell tag: " << n.mesh_cell_.second << std::endl;
                            std::cout << "n source rank: " << n.source_rank_ << std::endl;
                            std::cout << "n source tag: " << n.source_tag_ << std::endl;
                        }
                        assert(mesh != sp.mesh_.end());
                        mesh->set_as_field(Tag(n.mesh_cell_.second));
                    }
                }
            }
        }
        if (profiler_ != nullptr) {profiler_->stop("asm-ds-handle");}
    }

    void DonorSearcher::determine_orphan(const ArrCon<DonorInfo2>& arrival)
    {
        // handle orphans.
        // remove invalid candidate donors. this is why this is called before determine_best_donor() to avoid comparing invalid candidate donors.

        for (auto& sp: spc_->sp_)
        {
            for (Mesh& m: sp.mesh_)
            {
                for (MeshCell& mc: m.cell_)
                {
                    auto type = mc.oga_cell_type();
                    //if (type == OGA_cell_type_t::receptor || type == OGA_cell_type_t::mandat_receptor)
                    if (type == OGA_cell_type_t::mandat_receptor)
                    {
                        if (type == OGA_cell_type_t::receptor)
                        {
                            assert(!mc.cand_donor().empty());
                        }

                        for (auto cand_donor = mc.cand_donor_.begin(); cand_donor != mc.cand_donor_.end();)
                        {
                            const auto& cdm = cand_donor->mesh_tag_;
                            const auto& cdc = cand_donor->cell_tag_;

                            auto mesh_it = std::find_if(sp.mesh_.begin(), sp.mesh_.end(), [&](const Mesh& m){return m.tag() == cdm;});
                            assert(mesh_it != sp.mesh_.end());
                            assert(mesh_it->tag() != m.tag());

                            auto cit = mesh_it->query(Tag(cdc));

                            OGA_cell_type_t type = OGA_cell_type_t::undefined;
                            assert(cit != nullptr);

                            if (cit->oga_cell_type() == OGA_cell_type_t::non_resident)
                            {
                                //assert(!sp.is_resident(cit->poly().centroid()));
                                //for (const auto& receiver: nonresi)
                                {
                                    //if (receiver.tag() == sp.tag()())
                                    {
                                        auto nr = std::find_if(arrival.begin(), arrival.end(), [&](const DonorInfo2& n){return (n.meshtag_ == cdm() && n.celltag_ == cdc());});
                                        assert(nr != arrival.end());
                                        type = nr->oga_cell_type_;
                                        //break;
                                    }
                                }
                            }
                            else
                            {
                                type = cit->oga_cell_type();
                            }

                            assert(type != OGA_cell_type_t::undefined);

                            if (type == OGA_cell_type_t::receptor || type == OGA_cell_type_t::mandat_receptor || type == OGA_cell_type_t::orphan)
                            {
                                //if (comm_->rank() == 1 && mc.parent_mesh()() == 1 && mc.tag()() == 25433)
                                //{
                                    //std::cout << "cand mesh: " << cdm() << std::endl;
                                    //std::cout << "cand cell: " << cdc() << std::endl;
                                    //std::cout << "cand type: " << static_cast<int>(type) << std::endl;
                                //}
                                mc.cand_donor_.erase(cand_donor);
                            }
                            else
                            {
                                ++cand_donor;
                            }
                        }

                        if (mc.cand_donor().empty())
                        {
                            mc.set_oga_cell_type(OGA_cell_type_t::orphan);
                            //assert(false);
                            //mc.set_oga_cell_type(OGA_cell_type_t::field);
                        }
                    }
                }
            }
        }
    }

    void DonorSearcher::determine_best_donor(const ArrCon<DonorInfo2>& arrival)
    {
        for (auto& sp: spc_->sp_)
        {
            for (Mesh& m: sp.mesh_)
            {
                for (MeshCell& mc: m.cell_)
                {
                    auto type = mc.oga_cell_type();

                    if (type != OGA_cell_type_t::receptor && type != OGA_cell_type_t::mandat_receptor) {
                        continue;
                    }

                    assert(!mc.cand_donor().empty());

                    double current_area = TAILOR_BIG_POS_NUM;
                    Tag cdm_best;
                    Tag cdc_best;
                    const MeshCell* donor_best = nullptr;
                    for (auto& cand_donor: mc.cand_donor())
                    {
                        const auto& cdm = cand_donor.mesh_tag_;
                        const auto& cdc = cand_donor.cell_tag_;

                        auto mesh_it = std::find_if(sp.mesh_.begin(), sp.mesh_.end(), [&](const Mesh& m){return m.tag() == cdm;});
                        assert(mesh_it != sp.mesh_.end());
                        assert(mesh_it->tag() != m.tag());

                        auto cit = mesh_it->query(cdc);
                        assert(cit != nullptr);

                        double sa = 0.;
                        OGA_cell_type_t donor_type(OGA_cell_type_t::undefined);

                        int celltag;
                        Tag brmt;
                        //if (!sp.is_resident(cit->poly().centroid()))
                        if (!spc_->is_resident(cit->poly().centroid(), celltag, brmt))
                        {
                            //for (const auto& receiver: nonresi)
                            {
                                //auto nr = std::find_if(receiver.arrival().begin(), receiver.arrival().end(), [&](const auto& n){return (n.meshtag_ == cdm() && n.celltag_ == cdc());});
                                auto nr = std::find_if(arrival.begin(), arrival.end(), [&](const DonorInfo2& n){return (n.meshtag_ == cdm() && n.celltag_ == cdc());});
                                if (nr != arrival.end())
                                {
                                    donor_type = nr->oga_cell_type_;
                                    //sa = nr->volume_;
                                    sa = cit->poly().volume();
                                    //break;
                                }
                            }
                        }
                        else
                        {
                            donor_type = cit->oga_cell_type();
                            sa = cit->poly().volume();
                        }

                        if (donor_type == OGA_cell_type_t::receptor || donor_type == OGA_cell_type_t::mandat_receptor || donor_type == OGA_cell_type_t::orphan) {
                            continue;
                        }

                        if (sa < current_area)
                        {
                            current_area = sa; 
                            cdm_best = cdm;
                            cdc_best = cdc;
                            donor_best = &(*cit);
                        }
                    }

                    if (cdm_best.isvalid() && cdc_best.isvalid())
                    {
                        assert(donor_best != nullptr);
                        mc.set_donor(cdm_best, cdc_best, donor_best);
                        //if (donor_best->oga_cell_type() != OGA_cell_type_t::field)
                        //{
                            //std::cout << static_cast<int>(donor_best->oga_cell_type()) << std::endl;
                        //}
                        //assert(donor_best->oga_cell_type() == OGA_cell_type_t::field);
                    }
                    else
                    {
                        mc.set_oga_cell_type(OGA_cell_type_t::orphan);
                    }

                    if (mc.oga_cell_type() == OGA_cell_type_t::receptor || mc.oga_cell_type() == OGA_cell_type_t::mandat_receptor)
                    {
                        assert(mc.donor().addr_ != nullptr);
                    }
                }
            }
        }

        for (auto& sp: spc_->sp_)
        {
            for (Mesh& m: sp.mesh_)
            {
                for (MeshCell& mc: m.cell_)
                {
                    auto type = mc.oga_cell_type();

                    if (type != OGA_cell_type_t::receptor && type != OGA_cell_type_t::mandat_receptor) {
                        continue;
                    }

                    assert(mc.donor().addr_ != nullptr);
                    assert(mc.donor().addr_->oga_cell_type() != OGA_cell_type_t::receptor && mc.donor().addr_->oga_cell_type() != OGA_cell_type_t::mandat_receptor);
                }
            }
        }
    }
}

