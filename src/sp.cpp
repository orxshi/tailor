#include "sp.h"

namespace Tailor
{
    void SpatialPartition::set_comm(boost::mpi::communicator* comm)
    {
        assert(comm != nullptr);
        comm_ = comm;
    }

    void SpatialPartition::set_profiler(Profiler* prof)
    {
        assert(prof != nullptr);
        profiler_ = prof;
    }
    void SpatialPartition::reset_mesh_connectivity()
    {
        for (auto& m: mesh_)
        {
            m.reset_mesh_connectivity();
        }
    }
    std::vector<double> mesh_load(const SpatialPartition& sp)
    {
        std::vector<double> ml(sp.mesh().size(), 0.);

        for (int i=0; i<sp.mesh().size(); ++i)
        {
            const auto& m = sp.mesh()[i];

            ml[i] += m.cell().size();
        }

        return ml;
    }
    std::vector<double> mesh_field_load(const SpatialPartition& sp)
    {
        std::vector<double> ml(sp.mesh().size(), 0.);

        for (int i=0; i<sp.mesh().size(); ++i)
        {
            int count = 0;
            const auto& m = sp.mesh()[i];
            for (const auto& mc: m.cell())
            {
                if (mc.oga_cell_type() == OGA_cell_type_t::field)
                {
                    ++count;
                }
            }

            ml[i] += count;
        }

        return ml;
    }

    double solver_load(const SpatialPartition& sp)
    {
        auto ml = mesh_load(sp);
        return std::accumulate(ml.begin(), ml.end(), 0.);
    }

    double field_load(const SpatialPartition& sp)
    {
        auto ml = mesh_field_load(sp);
        return std::accumulate(ml.begin(), ml.end(), 0.);
    }

    void SpatialPartition::set_nei(std::deque<BinRMTag> nei)
    {
        nei_ = nei;
    }

    const std::deque<BinRMTag>& SpatialPartition::nei() const
    {
        return nei_;
    }

    void SpatialPartition::set_mesh_priority(int mtag, int pri)
    {
        assert(mtag != -1);
        assert(pri != -1);

        for (auto& m: mesh_)
        {
            if (mtag == m.tag()()) {
                m.set_priority(pri);
            }
        }
    }

    void SpatialPartition::sort_mesh()
    {
        // sort meshes in acsending order based on tags.
        std::sort(mesh_.begin(), mesh_.end(), [&](auto& a, auto& b){return a.tag()() < b.tag()();});
    }

    void SpatialPartition::update_connectivity(int rank)
    {
        for (auto& mesh: mesh_)
        {
            //mesh.update_connectivity(rank);
        }
    }
    void SpatialPartition::reset_all_oga_status()
    {
        for (auto& mesh: mesh_)
        {
            mesh.reset_all_cell_oga_status();
        }
    }

    void SpatialPartition::determine_non_residents()
    {
        for (Mesh& m: mesh_)
        {
            m.determine_non_residents(aabb_);
        }
    }

    void SpatialPartition::oga_interpolate(const ArrCon<Var>& arrival)
    {
        for (Mesh& mesh: mesh_)
        {
            mesh.oga_interpolate(arrival, -1);
        }
    }

    //void SpatialPartition::update_ghost_primitives(const ArrCon<Var>& arrival)
    //{
        //for (Mesh& mesh: mesh_)
        //{
            //mesh.update_ghost_primitives(arrival, comm_->rank());
        //}
    //}

    //void SpatialPartition::connect_partition_cells(const std::vector<MeshCell>& arrival_cell)
    void SpatialPartition::connect_partition_cells(ArrCon<Nei>& arrival_cell, std::function<bool(const vec3<double>&, int celltag)> is_resi, Profiler* profiler)
    {
        //assert(!arrival_cell.empty());
        //if (arrival_cell.empty()) {
            //return;
        //}

        for (Mesh& mesh: mesh_)
        {
            /*if (arrival_cell.empty())
            {
                for (const MeshCell& mc: mesh.cell())
                {
                    if (mc.btype() == BouType::partition && is_resident(mc.poly().centroid()))
                    {
                        std::cout << "mc: " << mc.tag()() << std::endl;
                        std::string fn = "ppp";
                        fn.append(std::to_string(comm_->rank()));
                        fn.append(".vtk");
                        mesh.print_as_vtk(fn);
                        assert(false);
                    }
                }
            }*/
            mesh.connect_partition_cells(arrival_cell, comm_->rank(), is_resi, profiler);
        }
    }

    void donor_info(int rank, const SpatialPartition& sp, int meshtag)
    {
        std::ofstream out;
        std::string fn = "mesh-donor-info-";
        fn.append(std::to_string(rank));
        fn.append("-");
        fn.append(std::to_string(sp.tag()()));
        fn.append(".csv");
        out.open(fn);

        for (int i=0; i<sp.mesh().size(); ++i)
        {
            const Mesh& m = sp.mesh()[i];
            //if (m.tag()() != meshtag) {continue;}
            for (const MeshCell& mc: m.cell())
            {
                out << mc.tag()();
                out << " ";
                out << m.tag()();
                out << " ";
                out << mc.oga_cell_type();
                out << "\n";
            }
        }

        out.close();
    }

    void SpatialPartition::convert_to_fortran() const
    {
        for (const Mesh& m: mesh_)
        {
            //convert_mesh_to_fortran(m);
        }
    }

    void SpatialPartition::init_sod()
    {
        Freestream fs;
        fs.read();

        vararray L, R;
        L[0] = 1.;
        L[1] = 0.;
        L[2] = 0.;
        L[3] = 0.;
        L[4] = 1.;

        R[0] = 0.125;
        R[1] = 0.;
        R[2] = 0.;
        R[3] = 0.;
        R[4] = 0.1;

        for (auto& mesh: mesh_)
        {
            mesh.init_sod(L, R, fs);
        }
    }

    void SpatialPartition::init()
    {
        //namespace po = boost::program_options;
        //po::options_description op;

        //po::options_description desc{"compo"};
        //desc.add_options()
        //    ("rotation", po::value<bool>()->default_value(false), "")
        //    ("rotaxis", po::value<int>()->default_value(0), "")
        //    ("rpm", po::value<int>()->default_value(0), "")
        //    ("pivotx", po::value<double>()->default_value(0), "")
        //    ("pivoty", po::value<double>()->default_value(0), "")
        //    ("pivotz", po::value<double>()->default_value(0), "")
        //    ;

        //op.add(desc);
        //std::string s = "compo_";
        //s.append(std::to_string(tag_()));
        //s.append(".ini");
        //std::ifstream settings_file(s);

        //boost::program_options::variables_map vm;
        //po::store(po::parse_config_file(settings_file, op, true), vm);
        //po::notify(vm);

        //bool rotation = vm["rotation"].as<bool>();
        //int rotaxis = vm["rotaxis"].as<int>();
        //int rpm = vm["rpm"].as<int>();
        //double pivotx = vm["pivotx"].as<double>();
        //double pivoty = vm["pivoty"].as<double>();
        //double pivotz = vm["pivotz"].as<double>();
        //auto pivot = vec3<double>(pivotx, pivoty, pivotz);




        Freestream fs;
        fs.read();
        

        double cinf = std::sqrt(fs.gamma_ * fs.pinf_ / fs.rhoinf_);
        vec3<double> vinf_air;

        if (fs.velair_ != 0.)
        {
            vinf_air = vec3<double>(
                    fs.velair_ * std::cos(deg_to_rad(fs.aoa_air_x_)),
                    fs.velair_ * std::cos(deg_to_rad(90. - fs.aoa_air_x_)),
                    fs.velair_ * std::cos(deg_to_rad(fs.aoa_air_z_)));
        }
        else
        {
            vinf_air = vec3<double>(
                    fs.machair_ * cinf * std::cos(deg_to_rad(fs.aoa_air_x_)),
                    fs.machair_ * cinf * std::cos(deg_to_rad(90. - fs.aoa_air_x_)),
                    fs.machair_ * cinf * std::cos(deg_to_rad(fs.aoa_air_z_)));
        }

        //vec3<double> vinf_foil(
                //fs.machfoil_ * cinf * std::cos(deg_to_rad(fs.aoa_foil_x_)),
                //fs.machfoil_ * cinf * std::cos(deg_to_rad(90. - fs.aoa_foil_x_)),
                //fs.machfoil_ * cinf * std::cos(deg_to_rad(fs.aoa_foil_z_)));
        //vec3<double> vinf = vinf_air;
        //vinf.set(0, 0, 0); // TODO test

        for (int i=0; i<mesh_.size(); ++i)
        {
            //vinf_air.set(0.);
            mesh_[i].init(vinf_air, fs);

            //if (std::abs(prim[3]) > TAILOR_ZERO)
            //{
            //    std::cout << "prim[3]: " << prim[3] << std::endl;
            //    std::cout << "fs.aoa_air_z_: " << fs.aoa_air_z_ << std::endl;
            //    std::cout << "machir" << fs.machair_ << std::endl;
            //    std::cout << "cinf" << cinf << std::endl;
            //}
            //assert(std::abs(prim[3]) < TAILOR_ZERO);
            //mesh_[i].set_prim_cons(BouType::dirichlet, prim, fs.gamma_);
            //mesh_[i].set_prim_cons(BouType::farfield, prim);
            //mesh_[i].set_prim_cons(BouType::interior, prim, fs.gamma_);
        }
    }

    /*void SpatialPartition::init_interior()
    {
        std::ifstream in;
        in.open("init_dirichlet.bc");
        if (!in.is_open()) {
            return;
        }

        for (int i=0; i<mesh_.size(); ++i)
        {
            go_to_beg_of_line(in, i*NVAR);
            double rhoinf, pinf, mach, aoa;
            in >> rhoinf;
            in >> pinf;
            in >> mach;
            in >> aoa;
            vararray prim;
            prim[0] = rhoinf;
            prim[1] = mach * std::cos(deg_to_rad(aoa));
            prim[2] = mach * std::sin(deg_to_rad(aoa));
            prim[3] = 0.;
            prim[4] = pinf;
            mesh_[i].set_prim_cons(BouType::interior, prim);
        }

        in.close();
    }/

    void SpatialPartition::init_sod()
    {
        for (int i=0; i<mesh_.size(); ++i)
        {
            BoundaryCondition::init_sod(mesh_[i]);
        }
    }

    /*void SpatialPartition::init_interior()
    {
        std::ifstream in;
        in.open("init_interior.bc");
        assert(in.is_open());

        for (int i=0; i<mesh_.size(); ++i)
        {
            go_to_beg_of_line(in, i*NVAR);
            vararray prim;
            for (int j=0; j<NVAR; ++j)
            {
                in >> prim[j];
            }
            mesh_[i].set_prim_cons(boundary_t("interior"), prim);
        }

        in.close();
    }*/

    //void SpatialPartition::solve(Solver& solver, const vec3<double>& vel, std::function<void(CellExchanger& cell_exchanger)> exchange_ghosts)
    //void SpatialPartition::solve(Solver& solver, const vec3<double>& vel, std::function<void()> exchange_ghosts, bool called_exc_ghosts)
    //{
        //for (int i=0; i<mesh_.size(); ++i)
        //{
            //if (called_exc_ghosts == false)
            //{
                //if (i != 0) {
                    //called_exc_ghosts = true;
                //}
            //}
            //solver.solve(mesh_[i], vel, exchange_ghosts, comm_->rank(), called_exc_ghosts);
        //}
    //}

    void SpatialPartition::connect_after_exchange(std::function<bool(const vec3<double>&)> is_resi, Profiler* profiler, std::string proc)
    {
        for (Mesh& m: mesh_)
        {
            m.connect_after_exchange(is_resi, comm_->rank(), profiler, proc);
        }
    }

    void SpatialPartition::update_prev_cand_donor()
    {
        for (Mesh& m: mesh_)
        {
            m.update_prev_cand_donor();
        }
    }
    const std::deque<ADT>& SpatialPartition::cell_adt() const
    {
        return cell_adt_;
    }

    const std::deque<AABB>& SpatialPartition::mesh_aabb() const
    {
        return mesh_aabb_;
    }

    const AABB& SpatialPartition::mesh_aabb(int i) const
    {
        return mesh_aabb_[i];
    }

    const ADT& SpatialPartition::cell_adt(int i) const
    {
        return cell_adt_[i];
    }

    void SpatialPartition::make_mesh_aabb()
    {
        mesh_aabb_.clear();
        for (const Mesh& m: mesh_)
        {
            mesh_aabb_.push_back(AABB(m.rawpoint()));
        }
    }

    void SpatialPartition::make_cell_adt()
    {
        cell_adt_.clear();
        for (const Mesh& m: mesh_)
        {
            //std::vector<ADTPoint> adt_pts;
            std::deque<ADTPoint> adt_pts;
            //adt_pts.reserve(m.cell().size());
            for (const MeshCell& mc: m.cell())
            {
                std::vector<Point> pts;
                pts.reserve(mc.point().size());
                for (auto it=mc.point().begin(); it != mc.point().end(); ++it) {
                    pts.push_back(it->p());
                }
                adt_pts.push_back(ADTPoint(pts.begin(), pts.end(), mc.tag()()));
            }
            //assert(std::unique(adt_pts.begin(), adt_pts.end(), [&](const auto& a, const auto& b){return a.idx() == b.idx();}) == adt_pts.end());

            if (adt_pts.empty())
            {
                cell_adt_.push_back(ADT());
            }
            else
            {
                cell_adt_.push_back(ADT(adt_pts));
            }
        }
    }

    void SpatialPartition::extend_aabb(const AABB& aabb)
    {
        aabb_.extend(aabb);
    }

    void SpatialPartition::profstart(std::string fun)
    {
        if (profiler_ == nullptr) {
            return;
        }

        profiler_->start(fun);
    }

    void SpatialPartition::profstop(std::string fun)
    {
        if (profiler_ == nullptr) {
            return;
        }

        profiler_->stop(fun);
    }

    //void SpatialPartition::set_profiler(Profiler* profiler)
    //{
        //profiler_ = profiler;
    //}

    void SpatialPartition::merge_meshes_no_check(std::deque<Mesh>& mesh) const
    {
        for (const Mesh& m: mesh_)
        {
            auto mit = std::find_if(mesh.begin(), mesh.end(), [&](const Mesh& tmpm){return tmpm.tag() == m.tag();});
            
            if (mit != mesh.end())
            {
                //mit->reserve_interior(m.cell().size()); // wrong. should be current size + incoming size.
                mit->reserve_interior(mit->cell().size() + m.cell().size());
                for (const MeshCell& mc: m.cell())
                {
                    mit->add_element_no_check(mc);
                }
                for (const MeshCell& mc: m.wall_boundaries())
                {
                    int count = std::count_if(mit->wall_boundaries().begin(), mit->wall_boundaries().end(), [&](const MeshCell& tmc){return tmc.tag() == mc.tag();});
                    if (count == 0)
                    {
                        mit->add_element(mc);
                    }
                }
                for (const MeshCell& mc: m.dirichlet_boundaries())
                {
                    int count = std::count_if(mit->dirichlet_boundaries().begin(), mit->dirichlet_boundaries().end(), [&](const MeshCell& tmc){return tmc.tag() == mc.tag();});
                    if (count == 0)
                    {
                        mit->add_element(mc);
                    }
                }
                for (const MeshCell& mc: m.farfield_boundaries())
                {
                    int count = std::count_if(mit->farfield_boundaries().begin(), mit->farfield_boundaries().end(), [&](const MeshCell& tmc){return tmc.tag() == mc.tag();});
                    if (count == 0)
                    {
                        mit->add_element(mc);
                    }
                }
                for (const MeshCell& mc: m.empty_boundaries())
                {
                    int count = std::count_if(mit->empty_boundaries().begin(), mit->empty_boundaries().end(), [&](const MeshCell& tmc){return tmc.tag() == mc.tag();});
                    if (count == 0)
                    {
                        mit->add_element(mc);
                    }
                }
                for (const MeshCell& mc: m.interog_boundaries())
                {
                    int count = std::count_if(mit->interog_boundaries().begin(), mit->interog_boundaries().end(), [&](const MeshCell& tmc){return tmc.tag() == mc.tag();});
                    if (count == 0)
                    {
                        mit->add_element(mc);
                    }
                }
            }
            else
            {
                mesh.push_front(m);
            }
        }
    }

    void SpatialPartition::merge_meshes(std::deque<Mesh>& mesh) const
    {
        assert(!mesh_.empty());
        for (const Mesh& m: mesh_)
        {
            auto mit = std::find_if(mesh.begin(), mesh.end(), [&](const Mesh& tmpm){return tmpm.tag() == m.tag();});
            
            if (mit != mesh.end())
            {
                for (const MeshCell& mc: m.cell())
                {
                    int count = std::count_if(mit->cell().begin(), mit->cell().end(), [&](const MeshCell& tmc){return tmc.tag() == mc.tag();});
                    if (count == 0)
                    {
                        mit->add_element_nonsorted(mc);
                    }
                }
                for (const MeshCell& mc: m.wall_boundaries())
                {
                    int count = std::count_if(mit->wall_boundaries().begin(), mit->wall_boundaries().end(), [&](const MeshCell& tmc){return tmc.tag() == mc.tag();});
                    if (count == 0)
                    {
                        mit->add_element_nonsorted(mc);
                    }
                }
                for (const MeshCell& mc: m.dirichlet_boundaries())
                {
                    int count = std::count_if(mit->dirichlet_boundaries().begin(), mit->dirichlet_boundaries().end(), [&](const MeshCell& tmc){return tmc.tag() == mc.tag();});
                    if (count == 0)
                    {
                        mit->add_element_nonsorted(mc);
                    }
                }
                for (const MeshCell& mc: m.farfield_boundaries())
                {
                    int count = std::count_if(mit->farfield_boundaries().begin(), mit->farfield_boundaries().end(), [&](const MeshCell& tmc){return tmc.tag() == mc.tag();});
                    if (count == 0)
                    {
                        mit->add_element_nonsorted(mc);
                    }
                }
                for (const MeshCell& mc: m.empty_boundaries())
                {
                    int count = std::count_if(mit->empty_boundaries().begin(), mit->empty_boundaries().end(), [&](const MeshCell& tmc){return tmc.tag() == mc.tag();});
                    if (count == 0)
                    {
                        mit->add_element_nonsorted(mc);
                    }
                }
                for (const MeshCell& mc: m.interog_boundaries())
                {
                    int count = std::count_if(mit->interog_boundaries().begin(), mit->interog_boundaries().end(), [&](const MeshCell& tmc){return tmc.tag() == mc.tag();});
                    if (count == 0)
                    {
                        mit->add_element_nonsorted(mc);
                    }
                }
            }
            else
            {
                mesh.push_front(m);
            }
        }
    }

    //void SpatialPartition::set_slave_comm(const boost::mpi::communicator* comm)
    //{
        //slave_ = comm;
    //}

    void SpatialPartition::remove_dup_cells_and_points()
    {
        for (Mesh& m: mesh_)
        {
            m.remove_dup_cells_and_points();
        }
    }

    void SpatialPartition::remove_cell(const Tag& im, const Tag& ic)
    {
        mesh_p(im).remove_cell(ic);
    }

    void SpatialPartition::reset_erase_marks()
    {
        for (Mesh& m: mesh_)
        {
            m.reset_erase_marks();
        }
    }

    void SpatialPartition::mark_to_be_erased(const Tag& im, const Tag& ic)
    {
        Mesh& m = mesh_p(im);
        m.mark_to_be_erased(ic);
    }

    void SpatialPartition::remove_ghosts()
    {
        for (Mesh& m: mesh_)
        {
            m.remove_ghosts();
        }
    }

    void SpatialPartition::remove_nonresidents()
    {
        for (Mesh& m: mesh_)
        {
            m.remove_nonresidents();
        }
    }

    void SpatialPartition::erase_marked_cells()
    {
        for (Mesh& m: mesh_)
        {
            m.erase_marked_cells(comm_->rank());
        }
    }

    void SpatialPartition::disconnect_orphan_faces()
    {
        for (Mesh& m: mesh_)
        {
            m.disconnect_orphan_faces();
        }
    }

    void SpatialPartition::add_mesh(const SpatialPartition& other_sp)
    {
        assert(tag_ == other_sp.tag());

        for (const Mesh& om: other_sp.mesh())
        {
            auto msh = std::find_if(mesh_.begin(), mesh_.end(), [&](const Mesh& m){return m.tag() == om.tag();});
            if (msh == mesh_.end())
            {
                add_mesh(om);
            }
        }
    }

    void SpatialPartition::add_merge_mesh_leave_dups(const SpatialPartition& other_sp)
    {
        //assert(tag_ == other_sp.tag());

        for (const Mesh& om: other_sp.mesh())
        {
            auto msh = std::find_if(mesh_.begin(), mesh_.end(), [&](const Mesh& m){return m.tag() == om.tag();});
            if (msh != mesh_.end())
            {
                //msh->merge_no_check(om);
                //std::cout << "merging - " << world_.rank() << std::endl;
                msh->merge_batch(om, comm_->rank());
                //msh->merge(om);
                //std::cout << "merged - " << world_.rank() << std::endl;
                //msh->remove_dup_cells();
                //msh->shrink_points();
                //msh->set_point_tag_index_map();
                //msh->merge_using_tags(om);
            }
            else
            {
                //std::cout << "adding - " << world_.rank() << std::endl;
                add_mesh(om);
                //std::cout << "added - " << world_.rank() << std::endl;
            }
        }
    }

    void SpatialPartition::remove_dups(std::string s, Profiler* profiler)
    {
        for (Mesh& msh: mesh_)
        {
            if (profiler != nullptr) {profiler->start(s + "-rem-dup-cells");}
            msh.remove_duplicate_cells(profiler, s);
            if (profiler != nullptr) {profiler->stop(s + "-rem-dup-cells");}

            if (profiler != nullptr) {profiler->start(s + "-rem-dup-pts");}
            msh.remove_duplicate_points();
            if (profiler != nullptr) {profiler->stop(s + "-rem-dup-pts");}

            if (profiler != nullptr) {profiler->start(s + "-rem-rem-parent");}
            msh.remove_parent_cells_of_vertices_of_all_cells();
            if (profiler != nullptr) {profiler->stop(s + "-rem-rem-parent");}

            if (profiler != nullptr) {profiler->start(s + "-rem-set-parent");}
            msh.set_parent_cell_of_vertices_of_all_cells();
            if (profiler != nullptr) {profiler->stop(s + "-rem-set-parent");}

            if (profiler != nullptr) {profiler->start(s + "-rem-update-pts");}
            msh.update_points_from_cell_vertices(comm_->rank());
            if (profiler != nullptr) {profiler->stop(s + "-rem-update-pts");}
        }
    }

    void SpatialPartition::merge_mesh(const SpatialPartition& other_sp)
    {
        assert(tag_ == other_sp.tag());

        for (const Mesh& om: other_sp.mesh())
        {
            auto msh = std::find_if(mesh_.begin(), mesh_.end(), [&](const Mesh& m){return m.tag() == om.tag();});
            if (msh != mesh_.end())
            {
                //msh->merge_no_check(om);
                msh->merge(om);
                msh->remove_merge_duplicate_cells();
                msh->shrink_points();
                //msh->set_point_tag_index_map();
                //msh->merge_using_tags(om);
            }
        }
    }

    void SpatialPartition::merge(const SpatialPartition& other_sp)
    {
        assert(tag_ == other_sp.tag());

        for (const Mesh& om: other_sp.mesh())
        {
            auto msh = std::find_if(mesh_.begin(), mesh_.end(), [&](const Mesh& m){return m.tag() == om.tag();});
            if (msh == mesh_.end())
            {
                add_mesh(om);
            }
            else
            {
                msh->remove_parent_cells_of_vertices_of_all_cells();
                msh->set_parent_cell_of_vertices_of_all_cells();
                msh->simple_merge(om, comm_->rank());
                msh->add_points(om.cell(), comm_->rank());
                //for (const MeshCell& mc: om.cell()) {
                    //msh->add_meshpoint(mc, world_.rank());
                //}

                //msh->shrink_points();
            }
        }
    }

    /*void SpatialPartition::print_result(int iteration) const
    {
        if (master()) return;

        size_t norphan = 0.;
        for (const Mesh& m: mesh_)
        {
            for (const MeshCell& mc: m.cell())
            {
                if (mc.oga_cell_type() == OGA_cell_type_t::orphan)
                    ++norphan;
            }
        }

        std::fstream out;
        out.open("orphan-vs-iter.csv", std::fstream::out | std::fstream::app);
        if (iteration == 0)
            out << "iter,orphan\n";
        out << iteration << "," << norphan << "\n";
        out.close();

        //if (printmesh_)
        {
            int i=0;
            for (const Mesh& m: mesh_)
            {
                std::string s = "sp_";
                s.append(std::to_string(tag_()));
                s.append("_mesh_");
                s.append(std::to_string(i));
                s.append("_iter_");
                s.append(std::to_string(iteration));
                s.append(".vtk");
                m.print_as_vtk(s);
                ++i;
            }
        }
    }*/

    std::deque<std::vector<Point>> SpatialPartition::mesh_pts() const
    {
        std::deque<std::vector<Point>> pts(mesh_.size());

        int i = 0;
        for (const Mesh& m: mesh_)
        {
            for (const MeshCell& mc: m.cell())
            {
                if (!is_resident(mc.poly().centroid())) {
                    continue;
                }

                std::vector<Point>& pt = pts[i];

                for (const MeshPoint& mp: mc.point())
                {
                    pt.push_back(mp.p());
                }
            }

            ++i;
        }

        return pts;
    }

    std::deque<std::vector<Point>> SpatialPartition::mesh_pts(const Outline& outline) const
    {
        assert(false);
        std::deque<std::vector<Point>> pts(mesh_.size());

        int i = 0;
        for (const Mesh& m: mesh_)
        {
            for (const MeshCell& mc: m.cell())
            {
                //if (!is_resident(mc.poly().centroid(), outline)) {
                    //continue;
                ////////}

                std::vector<Point>& pt = pts[i];

                for (const MeshPoint& mp: mc.point())
                {
                    pt.push_back(mp.p());
                }
            }

            ++i;
        }

        return pts;
    }

    /*std::deque<AABB> SpatialPartition::mesh_aabb() const
    {
        std::deque<AABB> aabb(mesh_.size());

        for (AABB& ab: aabb)
        {
            vec3<double> minn(TAILOR_BIG_POS_NUM, TAILOR_BIG_POS_NUM, TAILOR_BIG_POS_NUM);
            vec3<double> maxx(TAILOR_BIG_NEG_NUM, TAILOR_BIG_NEG_NUM, TAILOR_BIG_NEG_NUM);
            ab.set_bbox(minn, maxx);
        }

        int i = 0;
        for (const Mesh& m: mesh_)
        {
            AABB& ab = aabb[i];

            for (const MeshCell& mc: m.cell())
            {
                if (!is_resident(mc.polygon().centroid())) {
                    continue;
                }

                vec3<double> min_, max_;
                min_.set(ab.min(0), ab.min(1), ab.min(2));
                max_.set(ab.max(0), ab.max(1), ab.max(2));

                for (const MeshPoint& mp: mc.point())
                {
                    double vx = mp.p().r(0);
                    double vy = mp.p().r(1);
                    double vz = mp.p().r(2);

                    min_.set_x(std::min(min_(0), vx));
                    min_.set_y(std::min(min_(1), vy));
                    min_.set_z(std::min(min_(2), vz));
                    max_.set_x(std::max(max_(0), vx));
                    max_.set_y(std::max(max_(1), vy));
                    max_.set_z(std::max(max_(2), vz));
                }

                ab.set_bbox(min_, max_);
            }
            ++i;
        }

        return aabb;
    }*/

    /*std::deque<AABB> SpatialPartition::mesh_aabb(const Outline& outline) const
    {
        assert(false);

        std::deque<AABB> aabb(mesh_.size());

        for (AABB& ab: aabb)
        {
            vec3<double> minn(TAILOR_BIG_POS_NUM, TAILOR_BIG_POS_NUM, TAILOR_BIG_POS_NUM);
            vec3<double> maxx(TAILOR_BIG_NEG_NUM, TAILOR_BIG_NEG_NUM, TAILOR_BIG_NEG_NUM);
            ab.set_bbox(minn, maxx);
        }

        int i = 0;
        for (const Mesh& m: mesh_)
        {
            AABB& ab = aabb[i];

            for (const MeshCell& mc: m.cell())
            {
                //if (!is_resident(mc.polygon().centroid(), outline)) {
                    //continue;
                //}

                vec3<double> min_, max_;
                min_.set(ab.min(0), ab.min(1), ab.min(2));
                max_.set(ab.max(0), ab.max(1), ab.max(2));

                for (const MeshPoint& mp: mc.point())
                {
                    double vx = mp.p().r(0);
                    double vy = mp.p().r(1);
                    double vz = mp.p().r(2);

                    min_.set_x(std::min(min_(0), vx));
                    min_.set_y(std::min(min_(1), vy));
                    min_.set_z(std::min(min_(2), vz));
                    max_.set_x(std::max(max_(0), vx));
                    max_.set_y(std::max(max_(1), vy));
                    max_.set_z(std::max(max_(2), vz));
                }

                ab.set_bbox(min_, max_);
            }
            ++i;
        }

        return aabb;
    }*/

    //void SpatialPartition::set_comm(boost::mpi::communicator* comm)
    //{
        //comm_ = comm;
    //}

    void SpatialPartition::set_aabb(const AABB& aabb)
    {
        aabb_ = aabb;
    }
}
