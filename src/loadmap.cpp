#include "loadmap.h"

namespace Tailor
{
    void Loadmap::set_comm(boost::mpi::communicator* comm)
    {
        assert(comm != nullptr);
        comm_ = comm;
    }

    void Loadmap::set_profiler(Profiler* prof)
    {
        assert(prof != nullptr);
        profiler_ = prof;
    }

    void Loadmap::set_refine_tol(double d)
    {
        refine_tol_ = d;
    }

    void Loadmap::remove_cells_from_rm()
    {
        rm_->clear_cells();
    }

    Loadmap::Loadmap(boost::mpi::communicator* comm, Profiler* profiler, bool verbose, LoadEstim load_estim_type, bool pseudo3D, bool make_load_balance, RegType reg_type, const std::deque<Mesh>* mesh, std::string name): profiler_(profiler), comm_(comm), verbose_(verbose), load_estim_type_(load_estim_type), pseudo3D_(pseudo3D), make_load_balance_(make_load_balance), reg_type_(reg_type), mesh_(mesh), rm_(nullptr), name_(name)
    {
        read_settings();
        nrefine = 0;

        if (refine_limit_ == 0) {
            dorefine = false;
        }

        //threshold_ = 10.;
    }

    void Loadmap::reset_rm(int ns)
    {
        assert(rm_ == nullptr);
        //rm_.reset();
        delete rm_;
        //rm_ = std::make_unique<RegularMesh>();
        rm_ = new RegularMesh;
        rm_->set_tag(Tag(0));
        if (pseudo3D_) {
            rm_->set_nstripe(ns, ns, 1);
        }
        else {
            rm_->set_nstripe(ns, ns, ns);
        }
    }

    void Loadmap::profstart(std::string fun)
    {
        if (profiler_ == nullptr) {
            return;
        }

        profiler_->start(fun);
    }

    void Loadmap::profstop(std::string fun)
    {
        if (profiler_ == nullptr) {
            return;
        }

        profiler_->stop(fun);
    }

    //void Loadmap::recv_arrival_meshes(const SpatialPartition& sp)
    //{
        //arrival_meshes.insert(arrival_meshes.end(), sp.mesh().begin(), sp.mesh().end());
    //}

    void Loadmap::print(int iteration)
    {
        print_orphan(iteration);
        print_mesh(iteration);
        print_lm();
    }

    //size_t Loadmap::nbin() const
    //{
        //return rm_->size();
    //}

    /*void Loadmap::distribute_mti_2(std::vector<MeshTransferInfo>& mti) const
    {
        // prepare
        std::vector<int> mti_proc;
        std::vector<int> dest_rank(comm_->size());

        //if (master())
        {
            assert(!aug_bin_tag_.empty());
            for (int p=0; p<comm_->size(); ++p)
            {
                //assert(!aug_bin_tag_[p].empty());
                for (const BinRMTag& t: aug_bin_tag_[p])
                {
                    const Bin& bin = rm_->bin(t);
                    //assert(!bin.cell().empty());

                    assert(bin.parent_rm().isvalid());
                    bin.prepare_transfer_info(p+1, mti, mti_proc, dest_rank, rm_->size());
                }
            }

            assert(!mti.empty());
            //assert(mti_proc.empty());
        }

        // make sure mesh tags are unique
        for (const MeshTransferInfo& mti_: mti)
        {
            for (const MeshTransferInfo::Dest& dest: mti_.dest())
            {
                std::vector<Tag> temp;
                for (const auto& m: dest.mesh_)
                {
                    temp.push_back(m.tag_);
                }
                std::sort(temp.begin(), temp.end());
                auto it = std::adjacent_find(temp.begin(), temp.end());
                assert(it == temp.end());
            }
        }

        //boost::mpi::broadcast(*comm_, mti_proc, MASTER);
        //boost::mpi::broadcast(*comm_, dest_rank, MASTER);

        // send/recv mti to sources.
        for (int i=0; i<mti_proc.size(); ++i)
        {
            //if (master())
            {
                comm_->send(mti_proc[i], 0, mti[i]);
            }
            //else
            {
                if (comm_->rank() == mti_proc[i])
                {
                    MeshTransferInfo mti_temp;
                    comm_->recv(MASTER, 0, mti_temp);
                    mti.push_back(mti_temp);
                }
            }
        }

        std::vector<boost::mpi::request> req;
        
        // isend mti to dest.
        //if (master())
        {
            {
                for (const auto& mti_: mti)
                {
                    for (int i=0; i<mti_.rank().size(); ++i)
                    //for (int d: mti_.rank())
                    //for (const MeshTransferInfo::Dest& dest: mti_.dest())
                    {
                        int d = mti_.rank()[i];
                        auto dst = mti_.dest()[i];
                        //if (dest.rank_ == mti_.source()) continue;
                        if (d == mti_.source()) continue;
                        //std::cout << "I am " << mti_.source() << " sending to " << d << " ncell " << dst.ncell() << std::endl;
                        //req.push_back(world.isend(dest.rank_, 0, mti_));
                        req.push_back(comm_->isend(d, 0, mti_));
                    }
                }
            }
        }

        //std::cout << "receiving - " << world.rank() << std::endl;

        //world.barrier();

        // irecv mti for dest.
        //if (!master())
        {
            for (int i=0; i<dest_rank.size(); ++i)
            {
                if (comm_->rank() == i)
                {
                    for (int j=0; j<dest_rank[i]; ++j)
                    {
                        MeshTransferInfo mti_temp;
                        //boost::mpi::status = world.recv(MASTER, 0, mti_temp);
                        comm_->recv(MASTER, 0, mti_temp);
                        //std::cout << "I am " << world.rank() << " receiving from " << status.source() << " bin tag " << mti_temp << std::endl;
                        //std::cout << "I am " << world.rank() << " receiving from master" << std::endl;
                        mti.push_back(mti_temp);
                    }
                }
            }
        }

        boost::mpi::wait_all(req.begin(), req.end());
    }*/

    void Loadmap::clear()
    {
        unsigned int nsx = rm_->nstripe(0);
        unsigned int nsy = rm_->nstripe(1);
        unsigned int nsz = rm_->nstripe(2);
        //donor_info_.clear();
        nrefine = 0;
        //rm_.reset();
        delete rm_;
        //rm_ = std::make_unique<RegularMesh>();
        rm_ = new RegularMesh;
        //rm_ = RegularMesh();
        rm_->set_nstripe(nsx, nsy, nsz);
        rm_->set_tag(Tag(0));
        //mesh_.clear();
        //mesh_tag_index_map.clear();
        aug_bin_tag_.clear();
        bintag_proc_map_.clear();
    }

    /*void Loadmap::move_mesh(const std::vector<Vector3>& v)
    {
        //for (int i=0; i<mesh_.size(); ++i)
        int i = 0;
        for (Mesh& m: mesh_)
        {
            //mesh_[i].move(v[i]);
            m.move(v[i]);
            ++i;
        }
    }*/

    /*void Loadmap::rotate_mesh(const std::vector<double>& rotation, int axis, const std::vector<Vector3>& rot_point)
    {
        int i = 0;
        for (Mesh& m: mesh_)
        {
            //for (int i=0; i<mesh_.size(); ++i)
            //mesh_[i].rotate(rotation[i], rot_axis[i]);
            m.rotate(rotation[i], axis, rot_point[i]);
            ++i;
        }
    }*/

    void Loadmap::print_orphan(int iteration) const
    {
        //if (!master()) return;

        size_t norphan = 0.;
        for (const Mesh& m: *mesh_)
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
    }

    void Loadmap::print_mesh(int iteration) const
    {
        //if (master()) return;

        if (printmesh_)
        {
            //std::cout << "printing lm mesh" << std::endl;
            for (const Mesh& m: *mesh_)
            {
                std::string s = name_;
                s.append("-lm-rank-");
                s.append(std::to_string(comm_->rank()));
                s.append("-mesh-");
                s.append(std::to_string(m.tag()()));
                s.append("-iter-");
                s.append(std::to_string(iteration));
                s.append(".vtk");
                m.print_as_vtk(s);
            }
        }
    }

    bool Loadmap::resolution_is_enough(BinRMTag& heaviest_bt, int nmake_map)
    {
        std::ofstream out;
        if (print_dev_)
        {
            if (comm_->rank() == 0)
            {
                std::string fname = name_;
                fname.append("-dev-vs-refine-");
                fname.append(std::to_string(nmake_map));
                fname.append(".csv");
                out.open(fname, std::fstream::app);
                if (nrefine == 0) {
                    out << "ref,dev\n";
                }
            }
        }

        bool ire = true;
        bool collapse = true;

        std::vector<double> non_sorted_load;
        std::vector<double> sorted_load;
        std::vector<BinRMTag> sorted_tags;
        sort_load(sorted_load, non_sorted_load, sorted_tags, verbose_);

        //if (master())
        {
            heaviest_bt = get_heaviest_bin(sorted_tags, verbose_);

            if (nrefine == 0)
            {
                //graph_ = std::make_unique<Graph>(sorted_load[0], sorted_load[1], sorted_load[2], sorted_load[3], sorted_load[4], sorted_load[5], sorted_load[6], sorted_load[7]);
                std::vector<int> vertex_id;
                std::vector<int> weight(non_sorted_load.begin(), non_sorted_load.end());
                for (const Bin& b: rm_->bin()) {
                    vertex_id.push_back(b.tag()());
                }
                if (pseudo3D_) {
                    graph_ = std::make_unique<Graph>(vertex_id, weight, 4);
                }
                else {
                    assert(std::all_of(weight.begin(), weight.end(), [&](auto d){return d == 0.;}) == false);
                    graph_ = std::make_unique<Graph>(vertex_id, weight, 8);
                }
            }

            if (sorted_load.size() < comm_->size()) {
                collapse = false;
                ire = false;
            }

            if (make_load_balance_ == false)
            {
                if (sorted_load.size() >= comm_->size())
                {
                    collapse = true;
                    ire = true;
                }
            }

            if (collapse)
            {
                //std::ofstream out;
                //if (comm_->rank() == 0)
                //{
                    //std::string fname = "heaviest";
                    //fname.append(std::to_string(nrefine));
                    //fname.append(".dat");
                    //out.open(fname, std::fstream::app);
                    //out << heaviest_bt.rmtag()() << " " << heaviest_bt.bintag()() << " " << sorted_load.front() << std::endl;
                //}

                graph_->connect(comm_->rank());
                //graph_->print();
                double dev;
                bool balance = graph_->partition(comm_->size(), dev, aug_bin_tag_, nrefine, comm_->rank());
                assert(!aug_bin_tag_.empty());

                collapse = true;

                if (balance == false) {
                    collapse = false;
                    ire = false;
                }


                if (collapse)
                {
                    ire = true;
                    if (print_dev_)
                    {
                        if (comm_->rank() == 0) {
                            out << nrefine << "," << dev << "\n";
                        }
                    }
                }

                if (dev > refine_tol_) {
                    collapse = false;
                    ire = false;
                }

                if (make_load_balance_ == false)
                {
                    if (sorted_load.size() >= comm_->size())
                    {
                        collapse = true;
                        ire = true;
                    }
                }
            }
        }

        if (print_dev_)
        {
            if (comm_->rank() == 0) {
                out.close();
            }
        }

        return ire;
    }

    /*bool Loadmap::resolution_is_enough(BinRMTag& heaviest_bt, int nmake_map, int iteration)
    {
        bool ire;

        std::vector<double> sorted_load;
        std::vector<BinRMTag> sorted_tags;
        sort_load(sorted_load, sorted_tags, verbose_);
        if (master())
        {
            heaviest_bt = get_heaviest_bin(sorted_tags, verbose_);

            std::vector<double> aug_sorted_load;
            aug_bin_tag_.clear();
            assert(!sorted_tags.empty());
            bool collapsed = greedy_partition(comm_->size(), sorted_load, aug_sorted_load, sorted_tags, aug_bin_tag_, verbose_);
            double dev;
            if (!aug_sorted_load.empty()) {
                dev = estimated_deviation(aug_sorted_load);
            }
            final_dev_ = dev;

            std::ofstream out;
            std::string fname = "dev-vs-refine-";
            fname.append(std::to_string(nmake_map));
            fname.append(".csv");
            out.open(fname, std::fstream::app);
            if (nrefine == 0)
                out << "ref,dev\n";

            if (!collapsed)
            {
                ire = false;
            }
            else
            {
                ire = true;
                out << nrefine << "," << dev << "\n";
                out.close();
            }

            out.close();

            if (verbose_) {
                std::cout << "lm master outputs in reso - " << world.rank() << std::endl;
            }

            if (dev > refine_tol_)
                ire = false;
        }

        if (!uniproc_)
        {
            broadcast(world, ire, MASTER);
            broadcast(world, heaviest_bt, MASTER);
        }

        return ire;
    }*/

    Tag Loadmap::generate_meshtag() const
    {
        int i = TAILOR_BIG_NEG_NUM;
        Tag t;

        if (mesh_->empty())
        {
            t.set(0);
            return t;
        }
        for (const Mesh& m: *mesh_)
        {
            i = std::max(i, m.tag()()); 
        }

        t.set(i+1);

        return t;
    }

    void Loadmap::read_settings()
    {
        namespace po = boost::program_options;

        po::options_description all_options;

        po::options_description desc{"Loadmap options"};
        desc.add_options()
            ("loadmap.refine", po::value<bool>()->default_value(true), "Can refine load-map")
            ("loadmap.printlm", po::value<bool>()->default_value(false), "Print load-map")
            ("loadmap.print-lm-mesh", po::value<bool>()->default_value(false), "Print meshes in load-map")
            ("loadmap.refine-tol", po::value<double>()->default_value(10.), "Refinement tolerance")
            ("loadmap.refine-limit", po::value<int>()->default_value(1200), "Maximum number of refinement")
            ("loadmap.print-dev", po::value<bool>()->default_value(false), "")
            ("loadmap.print-graph-dist", po::value<bool>()->default_value(false), "")
            ;

        all_options.add(desc);

        std::ifstream settings_file("settings.ini");

        boost::program_options::variables_map vm;
        po::store(po::parse_config_file(settings_file, all_options, true), vm);
        po::notify(vm);

        dorefine = vm["loadmap.refine"].as<bool>();
        printlm_ = vm["loadmap.printlm"].as<bool>();
        printmesh_ = vm["loadmap.print-lm-mesh"].as<bool>();
        refine_tol_ = vm["loadmap.refine-tol"].as<double>();
        refine_limit_ = vm["loadmap.refine-limit"].as<int>();
        print_dev_ = vm["loadmap.print-dev"].as<bool>();
        print_graph_dist_ = vm["loadmap.print-graph-dist"].as<bool>();
    }

    AABB Loadmap::max_aabb() const
    {
        AABB aabb;

        for (const Mesh& m: *mesh_) {
            aabb.extend(AABB(m.rawpoint()));
        }

        return aabb;
    }

    /*void Loadmap::make_regular_map_uniproc_resident()
    {
        Vector3 global_steplength_;
        Vector3 global_aabb_min_;
        
        AABB aabb = max_aabb();

        rm_->set_aabb(aabb);
        rm_->calc_step_length();

        global_aabb_min_ = aabb.min();
        global_steplength_ = rm_->h();

        //rm_->set_aabb_min(global_aabb_min_);

        rm_->insert_bins(0, mesh_->size(), comm_->rank());

        for (const Mesh& mb: *mesh_) {
            rm_->register_resident_mesh(mb, true, comm_->rank());
        }
    }*/

    /*void Loadmap::make_regular_map_uniproc()
    {
        Vector3 global_steplength_;
        Vector3 global_aabb_min_;
        
        AABB aabb = max_aabb();

        rm_->set_aabb(aabb);
        rm_->calc_step_length();

        global_aabb_min_ = aabb.min();
        global_steplength_ = rm_->h();

        //rm_->set_aabb_min(global_aabb_min_);

        rm_->insert_bins(0, mesh_->size(), comm_->rank());

        for (const Mesh& mb: *mesh_) {
            rm_->register_mesh(mb, true, comm_->rank(), true);
        }
    }*/

    void Loadmap::make_regular_map(RegType regtype, bool adaptive)
    {
        std::string s = name_;
        s.append("-register");
        if (profiler_ != nullptr) {profiler_->start(s);}

        if (adaptive) {
            reset_rm(2);
        }
        else {
            reset_rm(std::ceil(std::sqrt(comm_->size())));
        }

        std::vector<AABB> aabbs;

        AABB aabb = max_aabb();
        if (profiler_ != nullptr) {profiler_->stop(s);}

        boost::mpi::all_gather(*comm_, aabb, aabbs);

        if (profiler_ != nullptr) {profiler_->start(s);}
        for (const AABB& ab: aabbs) {
            aabb.extend(ab);
        }

        rm_->set_aabb(AABB(aabb.min(), aabb.max()));
        rm_->calc_h();

        rm_->insert_bins(0, mesh_->size(), comm_->rank());
        rm_->insert_bin_addresses(rm_);
        assert(!rm_->bintag_address_map().empty());

        assert(!mesh_->empty());

        if (regtype == RegType::aabb)
        {
            for (const Mesh& mb: *mesh_) {
                rm_->register_mesh(mb, comm_->rank(), true);
            }
        }
        else if (regtype == RegType::centroid)
        {
            for (const Mesh& mb: *mesh_) {
                rm_->register_resident_mesh(mb, comm_->rank());
            }
        }
        if (profiler_ != nullptr) {profiler_->stop(s);}
    }

    /*void Loadmap::add_mesh(const std::deque<Mesh>& mesh)
    {
        for (const Mesh& m: mesh)
        {
            mesh_.push_back(m);
            mesh_tag_index_map.insert(boost::bimap<int, int>::value_type(mesh_.back().tag()(), mesh_.size() - 1));
        }
    }

    void Loadmap::add_mesh(Mesh& m)
    {
        if (m.tag()() == -1)
        {
            Tag t = generate_meshtag();
            m.set_tag(t);
        }
        assert(m.parent_mesh()() != -1);
        mesh_.push_back(m);
        mesh_tag_index_map.insert(boost::bimap<int, int>::value_type(m.tag()(), mesh_.size() - 1));
    }*/

    //void Loadmap::generate_aug_bin_tag()
    //{
        //if (!master()) {
            //return;
        //}

        //aug_bin_tag_.resize(comm_->size());

        //BinRMTag t(Tag(0), Tag(0));
        //int n = -1;
        //std::generate(aug_bin_tag_.begin(), aug_bin_tag_.end(),
                //[&n, &t]()
                //{
                //t.set_bin_tag(Tag(++n));
                //return t;
                //});
    //}

    void Loadmap::refine(int nmake_map)
    {
        nrefine = 0;
        BinRMTag heaviest_bt;
        bool ire = resolution_is_enough(heaviest_bt, nmake_map);

        if (dorefine)
        {
            while (!ire)
            {
                ++nrefine;
                int new_rmtag;
                int new_bt;
                for (int i=0; i<TAILOR_ZERO; ++i) {
                    assert(rm_->aabb().min(i) != rm_->aabb().max(i));
                }
                const RegularMesh* inited_rm = rm_->refine_adaptive(*mesh_, heaviest_bt, comm_->rank(), new_rmtag, new_bt, pseudo3D_, reg_type_);
                for (int i=0; i<TAILOR_ZERO; ++i) {
                    assert(inited_rm->aabb().min(i) != inited_rm->aabb().max(i));
                }

                std::vector<double> loads;
                std::vector<double> loads_r;
                LoadCalculator lc(load_estim_type_, comm_);
                //lc.gather_load(*rm_->bin(heaviest_bt.bintag()).rm(), loads, loads_r, *mesh_);
                lc.gather_load(*rm_->bin(heaviest_bt.bintag()).rm(), loads, *mesh_);
                //rm_->bin(heaviest_bt.bintag()).rm()->update_address();

                std::vector<int> tag;
                for (const Bin& b: inited_rm->bin())
                {
                    tag.push_back(b.tag()());
                }
                std::vector<int> loads_int(loads.begin(), loads.end());
                graph_->refine(heaviest_bt.bintag()(), heaviest_bt.rmtag()(), new_rmtag, tag, loads_int);

                ire = resolution_is_enough(heaviest_bt, nmake_map);

                if (nrefine == refine_limit_) {
                    break;
                }
            }
        }

        for (const auto& a: aug_bin_tag_)
        {
            for (const auto& brmt: a)
            {
                assert(brmt.isvalid());
            }
        }
    }

    void Loadmap::make(int nmake_map, bool adaptive)
    {
        make_regular_map(reg_type_, adaptive);

        std::string s = name_;
        s.append("-refine");
        if (profiler_ != nullptr) {profiler_->bstart(s);}
        if (adaptive) {
            refine(nmake_map);
        }
        else
        {
            aug_bin_tag_.resize(comm_->size());
            std::vector<BinRMTag> tag;
            rm_->get_bin_tags(tag);
            int j = 0;
            for (int i=0; i<tag.size(); ++i)
            {
                if (j >= aug_bin_tag_.size()) {j = 0;}
                aug_bin_tag_[j].push_back(tag[i]);
                ++j;
            }
        }

        get_bin_to_proc();

        if (graph_ != nullptr) {
            rm_->get_adjacency_from_graph(*graph_);
        }

        if (profiler_ != nullptr) {profiler_->bstop(s);}

        if (print_graph_dist_)
        {
            if (comm_->rank() == 0) {
                if (graph_ != nullptr) {
                    graph_->print_to_file();
                }
            }
        }

        if (profiler_ != nullptr) {profiler_->start(s);}
        graph_.reset();
        if (profiler_ != nullptr) {profiler_->stop(s);}
    }

    Loadmap::~Loadmap()
    {
        if (rm_ != nullptr) {
            delete rm_;
            rm_ = nullptr;
        }
    }

    void Loadmap::gather_regular_maps()
    {
        std::vector<RegularMesh> other_rm;

        assert(rm_ != nullptr);
        boost::mpi::all_gather(*comm_, *rm_, other_rm);

        //if (!master()) return;

        assert(!other_rm.empty());
        assert(other_rm.size() == comm_->size());

        for (int i=0; i<other_rm.size(); ++i)
        {
            if (i == comm_->rank()) {
                continue;
            }
            rm_->merge(other_rm[i]);
        }
    }

    void Loadmap::sort_load(std::vector<double>& sorted_load, std::vector<double>& non_sorted_load, std::vector<BinRMTag>& sorted_tags, bool verbose)
    {
        LoadCalculator lc(load_estim_type_, comm_);
        //std::vector<double> load_global_r;
        //sorted_tags = lc.sorted_bin_tags(*rm_, sorted_load, non_sorted_load, load_global_r, *mesh_);
        sorted_tags = lc.sorted_bin_tags(*rm_, sorted_load, non_sorted_load, *mesh_);

        //if (master()) {
            std::sort(sorted_load.begin(), sorted_load.end(), std::greater<double>());
        //}
    }

    BinRMTag Loadmap::get_heaviest_bin(const std::vector<BinRMTag>& sorted_tags, bool verbose)
    {
        return sorted_tags.front();

        if (verbose_) {
            std::cout << "lm sorted - " << comm_->rank() << std::endl;
        }
    }

    double Loadmap::estimated_deviation(const std::vector<double>& aug_sorted_load)
    {
        assert(!aug_sorted_load.empty());
        auto result = std::minmax_element(aug_sorted_load.begin(), aug_sorted_load.end());
        double minload = *result.first;
        double maxload = *result.second;
        double dev = (maxload - minload) * 100. / maxload;
        return dev;
    }

    /*void Loadmap::deviation_in_load(double& dev, BinRMTag& heaviest_bt, bool& collapsed, bool verbose)
    {
        if (verbose_) {
            std::cout << "lm calculating dev-in-load - " << world.rank() << std::endl;
        }
        verbose = false;

        std::vector<double> sorted_load;
        LoadCalculator lc(load_estim_type_, world, settings_);
        auto sorted_tags = lc.sorted_bin_tags(*rm_, sorted_load, mesh_);
        //std::vector<BinRMTag> sorted_tags = rm_->sort_bin_tags_based_on_load(world, sorted_load, mesh_, load_estim_type_);
        if (verbose_) {
            std::cout << "lm sorted tags - " << world.rank() << std::endl;
        }
        if (master())
        {
            assert(!sorted_load.empty());
            assert(!sorted_tags.empty());
        }
        if (!master()) return;
        heaviest_bt = sorted_tags.front();
        std::sort(sorted_load.begin(), sorted_load.end(), std::greater<double>());

        if (verbose_) {
            std::cout << "lm sorted - " << world.rank() << std::endl;
        }

        std::vector<double> aug_sorted_load;
        aug_bin_tag_.clear();
        assert(!sorted_tags.empty());
        if (uniproc_) {
            collapsed = greedy_partition(world.size(), sorted_load, aug_sorted_load, sorted_tags, aug_bin_tag_, verbose_);
        }
        else {
            collapsed = greedy_partition(comm_->size(), sorted_load, aug_sorted_load, sorted_tags, aug_bin_tag_, verbose_);
        }
        if (verbose_) {
            std::cout << "lm greedy partitioned - " << world.rank() << std::endl;
        }
        //assert(!aug_bin_tag_.empty());

        if (collapsed)
        {
            auto result = std::minmax_element(aug_sorted_load.begin(), aug_sorted_load.end());
            double minload = *result.first;
            double maxload = *result.second;
            //dev = static_cast<double>(maxload - minload) * 100. / static_cast<double>(maxload);
            dev = (maxload - minload) * 100. / maxload;

            //for (double ll: aug_sorted_load)
            //{
                //std::cout << "lm load: " << ll << std::endl;
            //}

            if (master())
            {
                for (int i=0; i<sorted_load.size(); ++i)
                {
                    std::cout << "z = " << i << "  " << sorted_load[i] << std::endl;
                }
                for (int i=0; i<aug_sorted_load.size(); ++i)
                {
                    std::cout << "i = " << i << "  " << aug_sorted_load[i] << std::endl;
                    for (int j=0; j<aug_bin_tag_[i].size(); ++j)
                    {
                        auto aaa = aug_bin_tag_[i][j];
                        std::cout << "j = " << j << "   " << aaa.bintag()() << "   " << rm_->bin(aaa).load_basic() << std::endl;
                    }
                }
            }
        }
    }*/

    void Loadmap::get_bin_to_proc()
    {
        //if (master())
        {
            assert(!aug_bin_tag_.empty());

            for (int j=0; j<aug_bin_tag_.size(); ++j)
            {
                //std::cout << "proc - " << j+1 << std::endl;
                for(const BinRMTag& i: aug_bin_tag_[j])
                {
                    //bintag_proc_map_.insert(std::pair<BinRMTag, int>(i, j));
                    bintag_proc_map_.insert(std::pair<int, int>(i.bintag()(), j));
                    //std::cout << i.rmtag()() << " " << i.bintag()() << std::endl;
                }
            }
        }

        /*for (const auto& a: bintag_proc_map_)
        {
            if (comm_->rank() == 3)
            {
                std::cout << "bintag: " << a.first.bintag()() << std::endl;
                std::cout << "proc: " << a.second << std::endl;
                assert(false);
            }
        }*/

        //if (!uniproc_) {
            //broadcast(*comm_, bintag_proc_map_, MASTER);
        //}
    }

    void Loadmap::print_lm() const
    {
        if (!printlm_) {return;}
        
        LoadCalculator lc(load_estim_type_, comm_);
        print_rm(*rm_, lc);
    }

    /*void print_rm(const RegularMesh& rm, std::string ss, const std::map<BinRMTag, int>& bintag_proc_map_)
    {
        std::string s = ss;
        s.append(std::to_string(rm.tag()()));
        s.append(".vtk");
        print(s, bintag_proc_map_, *rm_);

        for (const Bin& b: rm.bin())
        {
            if (b.rm() == nullptr) continue;
            print_rm(*b.rm(), ss , bintag_proc_map_);
        }
    }*/

    void Loadmap::print_rm(const RegularMesh& rm, const LoadCalculator& lc) const
    {
        //if (!master()) return;

        std::string s = name_;
        s.append("-lm-");
        s.append(std::to_string(rm.tag()()));
        s.append(".vtk");
        print_regmesh(s, bintag_proc_map_, lc, rm, *mesh_);

        std::cout << "bin sizee: " << rm_->bin().size() << std::endl;

        for (const Bin& b: rm.bin())
        {
            if (b.rm() == nullptr) continue;
            print_rm(*b.rm(), lc);
        }
    }

    bool spatial_partition(int nhead, const std::vector<double>& input, std::vector<double>& output, const std::vector<BinRMTag>& input_bin_tag, std::vector<std::vector<BinRMTag>>& output_bin_tag, bool verbose_)
    {
        if (input.size() < nhead)
            return false;

        output.reserve(nhead);
        output_bin_tag.resize(nhead);

        for (int i=0; i<nhead; ++i)
        {
            output.push_back(input[i]);
            output_bin_tag[i].push_back(input_bin_tag[i]);
        }

        if (input.size() == nhead)
            return true;

        for (auto i=input.begin()+nhead; i!=input.end(); ++i)
        {
            int dist = std::distance(input.begin(), i);
            double current = TAILOR_BIG_POS_NUM;
            int j = -1;
            for (int k=0; k<nhead; ++k)
            {
                if (output[k] < current)
                {
                    current = output[k];
                    j = k;
                }
            }
            if (j == -1)
            {
                if (verbose_) {
                    std::cout << "output.size = " << output.size() << std::endl;
                    std::cout << "current = " << current << std::endl;
                    std::cout << "nhead = " << nhead << std::endl;
                }
            }
            assert(j != -1);
            if (*i != 0)
            {
                output[j] += *i;
                output_bin_tag[j].push_back(input_bin_tag[dist]);
            }
        }

        int k = 0;
        for (auto i=input.rbegin(); i!=input.rend()-nhead; ++i)
        {
            if (*i != 0) continue;
            int dist = std::distance(input.begin(), i.base()) - 1;

            output[k] += *i;
            output_bin_tag[k].push_back(input_bin_tag[dist]);

            ++k;
            if (k >= nhead) k = 0;
        }

        return true;
    }

    bool greedy_partition(int nhead, const std::vector<double>& input, std::vector<double>& output, const std::vector<BinRMTag>& input_bin_tag, std::vector<std::vector<BinRMTag>>& output_bin_tag, bool verbose_)
    {
        if (input.size() < nhead)
            return false;

        output.reserve(nhead);
        output_bin_tag.resize(nhead);

        for (int i=0; i<nhead; ++i)
        {
            output.push_back(input[i]);
            output_bin_tag[i].push_back(input_bin_tag[i]);
        }

        if (input.size() == nhead)
            return true;

        for (auto i=input.begin()+nhead; i!=input.end(); ++i)
        {
            int dist = std::distance(input.begin(), i);
            double current = TAILOR_BIG_POS_NUM;
            int j = -1;
            for (int k=0; k<nhead; ++k)
            {
                if (output[k] < current)
                {
                    current = output[k];
                    j = k;
                }
            }
            if (j == -1)
            {
                if (verbose_) {
                    std::cout << "output.size = " << output.size() << std::endl;
                    std::cout << "current = " << current << std::endl;
                    std::cout << "nhead = " << nhead << std::endl;
                }
            }
            assert(j != -1);
            if (*i != 0)
            {
                output[j] += *i;
                output_bin_tag[j].push_back(input_bin_tag[dist]);
            }
        }

        int k = 0;
        for (auto i=input.rbegin(); i!=input.rend()-nhead; ++i)
        {
            if (*i != 0) continue;
            int dist = std::distance(input.begin(), i.base()) - 1;

            output[k] += *i;
            output_bin_tag[k].push_back(input_bin_tag[dist]);

            ++k;
            if (k >= nhead) k = 0;
        }

        return true;
    }
}
