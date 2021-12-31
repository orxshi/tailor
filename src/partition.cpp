#include "partition.h"

namespace Tailor
{ 
    const boost::mpi::communicator* Partition::comm() const
    {
        return comm_;
    }
    void Partition::set_comm(boost::mpi::communicator* comm)
    {
        assert(comm != nullptr);
        comm_ = comm;

        if (spc_ != nullptr)
        {
            spc_->set_comm(comm_);
        }
        if (loadmap_ != nullptr)
        {
            loadmap_->set_comm(comm_);
        }
    }

    void Partition::set_profiler(Profiler* prof)
    {
        assert(prof != nullptr);
        profiler_ = prof;

        if (spc_ != nullptr)
        {
            spc_->set_profiler(profiler_);
        }
        if (loadmap_ != nullptr)
        {
            spc_->set_profiler(profiler_);
        }
    }

    Partition::Partition(): print_bin_dist_(false), print_mesh_system_size_(false), print_cell_dist_(false), name_(""), pseudo3D_(false), load_estim_type_(LoadEstim::solver), n_load_balance_(0), comm_(nullptr), spc_(nullptr), loadmap_(nullptr), profiler_(nullptr)
    {
    }

    Partition::Partition(boost::mpi::communicator* comm, LoadEstim load_estim_type, bool pseudo3D, std::string name, Profiler* profiler): comm_(comm), n_load_balance_(0), load_estim_type_(load_estim_type), pseudo3D_(pseudo3D), loadmap_(nullptr), spc_(nullptr), profiler_(profiler), name_(name), print_bin_dist_(false), print_mesh_system_size_(false), print_cell_dist_(false)
    {
        read_settings();
    }

    //void Partition::rotate(const Tag& mesh, double ang, int axis, const Vector3& pivot)
    void Partition::rotate(const Tag& mesh, double ang, const Vector3& axis, const Vector3& pivot)
    {
        spc_->rotate_meshblocks(mesh, ang, axis, pivot);
    }

    void Partition::move(const Tag& mesh, const Vector3& v)
    {
        spc_->move_meshblocks(mesh, v);
    }

    void Partition::make(bool merge_bins, bool make_load_balance, int iter)
    {
        if (comm_->size() == 1)
        {
            make_loadmap_uniproc(mesh_);
        }
        else
        {
            if (load_estim_type_ == LoadEstim::solver)
            {
                repartition(mesh_, make_load_balance, RegType::aabb, merge_bins, false, load_estim_type_, iter);
                assert(!spc_->sp().empty());
            }
            else if (load_estim_type_ == LoadEstim::area || load_estim_type_ == LoadEstim::hybrid || load_estim_type_ == LoadEstim::type1)
            {
                repartition(mesh_, true, RegType::centroid, merge_bins, true, LoadEstim::solver, iter);

                if (make_load_balance)
                {
                    auto s = name_;
                    s.append("-updated-mesh");
                    if (profiler_ != nullptr) {profiler_->start(s);}
                    auto m = updated_mesh(*this, comm_->rank());
                    if (profiler_ != nullptr) {profiler_->stop(s);}

                    repartition(m, make_load_balance, RegType::aabb, merge_bins, false, load_estim_type_, iter);
                }
            }
            else
            {
                assert(false);
            }
        }

        connect_cells();

            //std::cout << "uuuuuuuuuuuuuuuuuuuu" << std::endl;
        
            //for (const auto& m: mesh_)
            //{
            //    std::string fn = "sp-";
            //    fn.append(std::to_string(comm_->rank()));
            //    fn.append("-");
            //    fn.append(std::to_string(m.tag()()));
            //    fn.append(".vtk");
            //    m.print_as_vtk_geometry(fn);
            //}

            //comm_->barrier();
            //assert(false);
    }

    SpatialPartitionContainer& Partition::spc()
    {
        return *spc_;
    }

    const SpatialPartitionContainer& Partition::spc() const
    {
        return *spc_;
    }

    void Partition::set_mesh(std::deque<Mesh>&& mesh)
    {
        auto s = name_;
        s.append("-set-mesh");
        if (profiler_ != nullptr) {profiler_->start(s);}
        mesh_ = mesh;
        if (profiler_ != nullptr) {profiler_->stop(s);}
    }

    void Partition::clear_mesh()
    {
        auto s = name_;
        s.append("-clear-mesh");
        if (profiler_ != nullptr) {profiler_->start(s);}
        mesh_.clear();
        //mesh_.shrink_to_fit();
        std::deque<Mesh>().swap(mesh_);
        if (profiler_ != nullptr) {profiler_->stop(s);}
    }

    std::deque<Mesh> updated_mesh(const Partition& partition, int rank)
    {
        //int nmesh = 0;
        //int ncell = 0;
        //for (const auto& sp: partition.spc().sp())
        //{
            //nmesh += sp.mesh().size();
            //for (const auto& m: sp.mesh())
            //{
                //ncell += m.cell().size();
            //}
        //}
        //std::cout << rank << " part begin updatedmesh " << partition.spc().sp().size() << " " << nmesh << " " << ncell << std::endl;

        // If you don't update mesh, you are taking risk of duplicated cells in different processors. After exchange of cells, a cell might reside in different processors. Initially, with graph partitioning, cells were non-duplicated.

        // update from spc after exchange.

        //for (const auto& sp: partition.spc().sp())
        //{
            //for (const auto& m: sp.mesh())
            //{
                //for (const auto& mc: m.cell())
                //{
                    //assert(std::count(m.cell().begin(), m.cell().end(), mc) == 1);
                //}
            //}
        //}

        std::deque<Mesh> mesh;
        partition.spc().merge_sp_meshes(mesh);

        //for (const auto& sp: partition.spc().sp())
        //{
            //for (const auto& m: sp.mesh())
            //{
                //for (const auto& mc: m.cell())
                //{
                    //assert(std::count(m.cell().begin(), m.cell().end(), mc) == 1);
                //}
            //}
        //}

        //std::cout << rank << " part end updatedmesh" << std::endl;

        return mesh;
    }

    void Partition::make_loadmap_uniproc(const std::deque<Mesh>& mesh)
    {
        spc_ = new SpatialPartitionContainer(comm_, profiler_, false, mesh.size());
        spc_->remap_uniproc(mesh);
        //connect_cells();
    }

    void Partition::repartition(const std::deque<Mesh>& mesh, bool make_load_balance, RegType regtype, bool merge_bins, bool part_to_commsize_only, LoadEstim let, int iter)
    {
            {
                std::string mems = "bef-lm";
                mem_usage(comm_, mems);
            }

        make_loadmap(mesh, make_load_balance, regtype, make_load_balance, part_to_commsize_only, let);

            {
                std::string mems = "aft-lm";
                mem_usage(comm_, mems);
            }

        remap(merge_bins, iter);

            {
                std::string mems = "aft-remap";
                mem_usage(comm_, mems);
            }

        assert(!spc_->sp().empty());
    }

    void Partition::connect_cells()
    {
        auto s = name_;
        s.append("-connect-cells");
        if (profiler_ != nullptr) {profiler_->start(s);}
        spc_->connect_cells(name_);
        if (profiler_ != nullptr) {profiler_->stop(s);}
    }

    void Partition::make_loadmap(const std::deque<Mesh>& mesh, bool adaptive, RegType reg_type, bool make_load_balance, bool part_to_commsize_only, LoadEstim let)
    {
        if (loadmap_ == nullptr) {
            loadmap_ = new Loadmap(comm_, profiler_, false, let, pseudo3D_, make_load_balance, reg_type, &mesh, name_);
        }
        //if (part_to_commsize_only)
        //{
            //loadmap_->set_refine_tol(100.);
        //}
        loadmap_->make(n_load_balance_, adaptive);
        ++n_load_balance_;
    }

    void Partition::read_settings()
    {
        namespace po = boost::program_options;

        po::options_description all_options;

        po::options_description desc{"Partition"};
        desc.add_options()
            ("partition.print-bin-dist", po::value<bool>()->default_value(false), "")
            ("partition.print-cell-dist", po::value<bool>()->default_value(false), "")
            ("partition.print-mesh-system-size", po::value<bool>()->default_value(false), "")
            ;

        all_options.add(desc);

        std::ifstream settings_file("settings.ini");

        boost::program_options::variables_map vm;
        po::store(po::parse_config_file(settings_file, all_options, true), vm);
        po::notify(vm);

        print_bin_dist_ = vm["partition.print-bin-dist"].as<bool>();
        print_cell_dist_ = vm["partition.print-cell-dist"].as<bool>();
        print_mesh_system_size_ = vm["partition.print-mesh-system-size"].as<bool>();
    }

    void Partition::remap(bool merge_bins, int iter)
    {
        auto s = name_;
        s.append("-del-spc");
        if (profiler_ != nullptr) {profiler_->start(s);}
            {
                std::string mems = "bef-del-spc";
                mem_usage(comm_, mems);
            }
        delete spc_;
            {
                std::string mems = "aft-del-spc";
                mem_usage(comm_, mems);
            }
        spc_ = new SpatialPartitionContainer(comm_, profiler_, false, loadmap_->bintag_proc_map(), mesh_.size());
            {
                std::string mems = "aft-ctr-spc";
                mem_usage(comm_, mems);
            }
        if (profiler_ != nullptr) {profiler_->stop(s);}

        assert(mesh_.front().tag()() == 0);

        if (print_bin_dist_) {
            print_bin_dist(comm_, this, iter);
        }

        Remapper remapper(spc_, comm_, name_);
        remapper.remap(*loadmap_, merge_bins, iter, name_, profiler_);
        //assert(!spc_->sp().empty());
            {
                std::string mems = "aft-rmp";
                mem_usage(comm_, mems);
            }

        auto s1 = name_;
        s1.append("-del-lm");
        if (profiler_ != nullptr) {profiler_->start(s1);}
        delete loadmap_;
            {
                std::string mems = "aft-del-lm";
                mem_usage(comm_, mems);
            }
        loadmap_ = nullptr;
        if (profiler_ != nullptr) {profiler_->stop(s1);}
    }

    const std::deque<Mesh> Partition::mesh() const
    {
        return mesh_;
    }

    bool read_msh(Mesh& mesh, std::string file_name, std::string type, int rank, int comm_size)
    {
        assert(file_name != "");

        if (type != "interior")
        {
            file_name.append("_");
            file_name.append(type);
            if (comm_size != 1)
            {
                file_name.append("_");
                file_name.append(std::to_string(rank+1)); // +1 is due to gmsh convention.
            }
        }
        else
        {
            if (comm_size != 1)
            {
                file_name.append("_");
                file_name.append(std::to_string(rank+1)); // +1 is due to gmsh convention.
            }
        }
        file_name.append(".msh");
        //std::cout << type << std::endl;

        //if (!std::experimental::filesystem::exists(file_name)) {
        //if (!std::filesystem::exists("helloworld.txt")) {
        std::ifstream infile(file_name);
        if (!infile.good()) {
            return false;
        }

        assert(mesh.tag().isvalid());
        read_mesh_GMSH(mesh, file_name);
        //if (type == "wall") {
            //std::cout << file_name << std::endl;
            //std::cout << "size: "  << mesh.wall_boundaries().size() << std::endl;
        //}
        return true;
    }

    void Partition::read_mesh_(BouType btype, std::string fn)
    {
        //std::deque<Mesh>* bcontainer;
        std::string s;

        if (btype == BouType::wall)
        {
            //bcontainer = &wall_;
            s = "wall";
        }
        else if (btype == BouType::symmetry)
        {
            //bcontainer = &wall_;
            s = "symmetry";
        }
        else if (btype == BouType::dirichlet)
        {
            //bcontainer = &dirichlet_;
            s = "dirichlet";
        }
        else if (btype == BouType::farfield)
        {
            //bcontainer = &farfield_;
            s = "farfield";
        }
        else if (btype == BouType::empty)
        {
            //bcontainer = &empty_;
            s = "empty";
        }
        else if (btype == BouType::interog)
        {
            //bcontainer = &interog_;
            s = "interog";
        }

        for (int i=0; i<comm_->size(); ++i)
        {
            Mesh bm;
            bm.set_tag(mesh_.back().tag());

            read_msh(bm, fn, s, i, comm_->size());
            auto s1 = name_;
            s1.append("-read-mesh");
            if (profiler_ != nullptr) {profiler_->start(s1);}
            mesh_.back().connect_add_bou_to_interior(bm, btype, comm_->rank());
            if (profiler_ != nullptr) {profiler_->stop(s1);}
        }


        //if (read_msh(bm, fn, s, comm_->rank(), comm_->size()))
        //{
            //bcontainer->push_back(std::move(bm));
            //mesh_.back().connect_add_bou_to_interior(bcontainer->back(), btype, comm_->rank());
            //mesh_.back().connect_add_bou_to_interior(bm, btype, comm_->rank());
        //}

        for (const auto& mc: mesh_.back().cell())
        {
            for (const Tag& bc: mc.symmetry_boundary())
            {
                assert(mesh_.back().query_bou(bc, BouType::symmetry) != nullptr);
            }
        }
        for (const auto& mc: mesh_.back().cell())
        {
            for (const Tag& bc: mc.dirichlet_boundary())
            {
                assert(mesh_.back().query_bou(bc, BouType::dirichlet) != nullptr);
            }
        }

        for (const auto& mc: mesh_.back().cell())
        {
            for (const Tag& bc: mc.farfield_boundary())
            {
                assert(mesh_.back().query_bou(bc, BouType::farfield) != nullptr);
            }
        }

        for (const auto& mc: mesh_.back().cell())
        {
            for (const Tag& bc: mc.interog_boundary())
            {
                assert(mesh_.back().query_bou(bc, BouType::interog) != nullptr);
            }
        }
    }

    void Partition::read_mesh(const std::vector<std::string>& file_name)
    {
        auto s1 = name_;
        s1.append("-read-mesh");
        if (profiler_ != nullptr) {profiler_->start(s1);}
        for (int i=0; i<file_name.size(); ++i)
        {
            auto fn = file_name[i];

            Mesh m;
            m.set_tag(Tag(i));
            read_msh(m, fn, "interior", comm_->rank(), comm_->size());
            if (m.cell().empty())
            {
                std::cout << "fn: " << fn << std::endl;
            }
            assert(!m.cell().empty());

            m.connect_add_bou_to_interior(BouType::wall, comm_->rank());
            m.connect_add_bou_to_interior(BouType::symmetry, comm_->rank());
            m.connect_add_bou_to_interior(BouType::dirichlet, comm_->rank());
            m.connect_add_bou_to_interior(BouType::farfield, comm_->rank());
            m.connect_add_bou_to_interior(BouType::interog, comm_->rank());
            m.connect_add_bou_to_interior(BouType::empty, comm_->rank());

            /*
            bool allempty = true;

            if (!m.cell().empty()) {
                allempty = false;
            }
            if (!m.wall_boundaries().empty()) {
                allempty = false;
            }
            if (!m.farfield_boundaries().empty()) {
                allempty = false;
            }
            if (!m.interog_boundaries().empty()) {
                allempty = false;
            }
            if (!m.empty_boundaries().empty()) {
                allempty = false;
            }
            if (!m.dirichlet_boundaries().empty()) {
                allempty = false;
            }
            assert(allempty == false);*/

            //comm_->barrier(); // delete later.

            mesh_.push_back(std::move(m));
            //for (const auto& m: mesh_)
            //{
            //    std::string fn = "readmeshrior-";
            //    fn.append(std::to_string(comm_->rank()));
            //    fn.append("-");
            //    fn.append(std::to_string(m.tag()()));
            //    fn.append(".vtk");
            //    m.print_as_vtk_geometry(fn);
            //}


            //read_mesh_(BouType::wall, fn);
            //read_mesh_(BouType::dirichlet, fn);
            //read_mesh_(BouType::farfield, fn);
            //read_mesh_(BouType::interog, fn);
            //read_mesh_(BouType::empty, fn);

            //std::cout << comm_->rank() << " " << mesh_.back().cell().size() << " " << mesh_.back().point().size() << std::endl;
            //std::cout << comm_->rank() << " " << wall_.back().cell().size() << " " << wall_.back().point().size() << std::endl;

            //mesh_.back().update_cell_vertex_addresses();

            //for (const auto& m: mesh_)
            //{
            //    std::string fn = "pre-";
            //    fn.append(std::to_string(comm_->rank()));
            //    fn.append("-");
            //    fn.append(std::to_string(m.tag()()));
            //    fn.append(".vtk");
            //    m.print_as_vtk(fn);
            //}
        }
        if (profiler_ != nullptr) {profiler_->stop(s1);}

        if (print_mesh_system_size_)
        {
            for (int i=0; i<file_name.size(); ++i)
            {
                int total_cell_size = 0;
                int local_cell_size = mesh_[i].cell().size();
                //std::cout << file_name[i] << " "  << comm_->rank() << " " << local_cell_size << std::endl;
                boost::mpi::reduce(*comm_, local_cell_size, total_cell_size, std::plus<int>(), 0);
                if (comm_->rank() == 0)
                {
                    std::ofstream out;
                    out.open("cell-size.dat", std::ios_base::app);
                    out << i << " " << total_cell_size << std::endl;
                    out.close();
                }
            }
        }
    }

    Partition::~Partition()
    {
        if (spc_ != nullptr) {
            delete spc_;
            spc_ = nullptr;
        }
        if (loadmap_ != nullptr) {
            delete loadmap_;
            loadmap_ = nullptr;
        }
    }

    /*const std::deque<Mesh> Partition::wall() const
    {
        return wall_;
    }
    const std::deque<Mesh> Partition::dirichlet() const
    {
        return dirichlet_;
    }
    const std::deque<Mesh> Partition::farfield() const
    {
        return farfield_;
    }
    const std::deque<Mesh> Partition::interog() const
    {
        return interog_;
    }
    const std::deque<Mesh> Partition::empty() const
    {
        return empty_;
    }*/

    std::string Partition::name() const
    {
        return name_;
    }

    void Partition::print_cell_dist(const boost::mpi::communicator* comm, int iter) const
    {
        if  (print_cell_dist_)
        {
            std::ofstream out;
            std::string s = name_;
            s.append("-ncell-");
            s.append(std::to_string(iter));
            s.append(".dat");
            out.open(s, std::ios_base::app);
            int nmesh = 0;
            int ncell = 0;
            for (const auto& sp: spc_->sp())
            {
                nmesh += sp.mesh().size();
                for (const auto& m: sp.mesh())
                {
                    ncell += m.cell().size();
                }
            }
            out << comm->rank() << " " << nmesh << " " << ncell << std::endl;
            out.close();
        }
    }

    void print_bin_dist(const boost::mpi::communicator* comm, const Partition* partition, int iter)
    {
        if (comm->rank() == 0)
        {
            std::ofstream out;
            std::string s = partition->name();
            s.append("-bin-dist.dat");
            out.open(s);
            auto abt = partition->loadmap()->aug_bin_tag();

            for (int i=0; i<abt.size(); ++i)
            {
                out << i << " | " << abt[i].size() << " | ";
                for (auto ii: abt[i]) {
                    out << "(" << ii.bintag()() << ","; 
                }
                out << std::endl;
            }
            out.close();
        }
    }

    const Loadmap* Partition::loadmap() const
    {
        return loadmap_;
    }
}
