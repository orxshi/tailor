#include "assembler.h"

namespace Tailor
{
    void Assembler::move(const Tag& mesh, const Vector3& v)
    {
        partition_->move(mesh, v);
    }

    void Assembler::set_comm(boost::mpi::communicator* comm)
    {
        assert(comm != nullptr);
        comm_ = comm;

        if (partition_ != nullptr)
        {
            partition_->set_comm(comm_);
        }
    }

    void Assembler::set_profiler(Profiler* prof)
    {
        profiler_ = prof;

        if (partition_ != nullptr)
        {
            partition_->set_profiler(profiler_);
        }
    }
    Assembler::Assembler():
        comm_(nullptr),
        profiler_(nullptr),
        rebalance_thres_(0.),
        print_ds_info_(false),
        print_vtk_(false),
        print_pre_vtk_(false),
        print_repart_info_(false),
        print_imbalance_(false),
        mergebins_(false),
        donor_search_algo_(DonorSearchAlgo::alladt),
        nassemble_(0),
        load_estim_type_(LoadEstim::solver),
        can_rebalance_(true),
        force_rebalance_(false),
        make_load_balance_(true),
        pseudo3D_(false),
        print_map_(false),
        partition_(nullptr)
    {
    }

    Assembler::Assembler(boost::mpi::communicator* comm, Profiler* profiler, const std::vector<std::string>& filename):
        comm_(comm),
        profiler_(profiler),
        rebalance_thres_(0.),
        print_ds_info_(false),
        print_vtk_(false),
        print_pre_vtk_(false),
        print_repart_info_(false),
        print_imbalance_(false),
        mergebins_(false),
        donor_search_algo_(DonorSearchAlgo::alladt),
        nassemble_(0),
        load_estim_type_(LoadEstim::solver),
        can_rebalance_(true),
        force_rebalance_(false),
        make_load_balance_(true),
        pseudo3D_(false),
        print_map_(false),
        partition_(nullptr)
    {
        read_settings();

        partition_ = new Partition(comm_, load_estim_type_, pseudo3D_, "asm", profiler_);
        partition_->read_mesh(filename);

        if (print_pre_vtk_)
        {
            for (const auto& m: partition_->mesh())
            {
                std::string fn = "asm-pre-";
                fn.append(std::to_string(comm_->rank()));
                fn.append("-");
                fn.append(std::to_string(m.tag()()));
                fn.append(".vtk");
                m.print_as_vtk_geometry(fn);
            }
        }

        partition_->make(mergebins_, make_load_balance_, nassemble_);

        partition_->clear_mesh();
        assert(!partition_->spc().sp().empty());
    }

    int Assembler::nassemble() const
    {
        return nassemble_;
    }

    void Assembler::rotate(const Tag& mesh, double ang, int axis, const Vector3& pivot)
    {
        partition_->rotate(mesh, ang, axis, pivot);
    }

    void Assembler::read_settings()
    {
        namespace po = boost::program_options;

        po::options_description all_options;

        po::options_description desc{"Assembler options"};
        desc.add_options()
            ("assembler.rebalance-thres", po::value<double>()->default_value(60.), "")
            ("general.pseudo3D", po::value<bool>()->default_value(false), "2.5D simulation with 1 layer in depth")
            ("assembler.load-estim", po::value<int>()->default_value(2), "Load estimation method")
            ("assembler.make-load-balance", po::value<bool>()->default_value(true), "Make load balance (for area or hybrid load estimation)")
            ("assembler.donor-search-algo", po::value<int>()->default_value(0), "Use adt for donor search in all iterations or use adt in first iteration and whenever seed is not found.")
            ("assembler.can-rebalance", po::value<bool>()->default_value(true), "Make load rebalance if needed.")
            ("assembler.force-rebalance", po::value<bool>()->default_value(false), "Force load rebalance even if relabalance is not needed.")
            ("assembler.print-map", po::value<bool>()->default_value(false), "Print map.")
            ("assembler.merge-bins", po::value<bool>()->default_value(true), "")
            ("assembler.print-ds-info", po::value<bool>()->default_value(false), "")
            ("assembler.print-vtk", po::value<bool>()->default_value(false), "")
            ("assembler.print-pre-vtk", po::value<bool>()->default_value(false), "")
            ("assembler.print-repart-info", po::value<bool>()->default_value(false), "")
            ("assembler.print-imbalance", po::value<bool>()->default_value(false), "")
            ;

        all_options.add(desc);

        std::ifstream settings_file("settings.ini");

        boost::program_options::variables_map vm;
        po::store(po::parse_config_file(settings_file, all_options, true), vm);
        po::notify(vm);

        pseudo3D_ = vm["general.pseudo3D"].as<bool>();
        load_estim_type_ = static_cast<LoadEstim>(vm["assembler.load-estim"].as<int>());
        can_rebalance_ = vm["assembler.can-rebalance"].as<bool>();
        force_rebalance_ = vm["assembler.force-rebalance"].as<bool>();
        make_load_balance_ = vm["assembler.make-load-balance"].as<bool>();
        donor_search_algo_ = static_cast<DonorSearchAlgo>(vm["assembler.donor-search-algo"].as<int>());
        print_map_ = vm["assembler.print-map"].as<bool>();
        mergebins_ = vm["assembler.merge-bins"].as<bool>();
        print_ds_info_ = vm["assembler.print-ds-info"].as<bool>();
        rebalance_thres_ = vm["assembler.rebalance-thres"].as<double>();
        print_vtk_ = vm["assembler.print-vtk"].as<bool>();
        print_pre_vtk_ = vm["assembler.print-pre-vtk"].as<bool>();
        print_repart_info_ = vm["assembler.print-repart-info"].as<bool>();
        print_imbalance_ = vm["assembler.print-imbalance"].as<bool>();
    }

    bool Assembler::repartition()
    {
        bool repart = false;

        if (nassemble_ < 1) {
            return repart;
        }

        auto aa = [&]()
        {
            if (profiler_ != nullptr) {profiler_->start("asm-rem-nonresi");}
            partition_->spc_->remove_nonresidents();
            if (profiler_ != nullptr) {profiler_->stop("asm-rem-nonresi");}

            if (profiler_ != nullptr) {profiler_->start("asm-reset-status");}
            reset_oga_status();
            if (profiler_ != nullptr) {profiler_->stop("asm-reset-status");}

            if (profiler_ != nullptr) {profiler_->start("asm-reset-connect");}
            partition_->spc_->reset_mesh_connectivity();
            if (profiler_ != nullptr) {profiler_->stop("asm-reset-connect");}

            if (mergebins_)
            {
                partition_->set_mesh(std::move(partition_->spc_->sp_.front().mesh_));
                partition_->repartition(partition_->mesh(), make_load_balance_, RegType::aabb, mergebins_, false, load_estim_type_, nassemble_);
            }
            else
            {
                assert(false);
            }

            partition_->clear_mesh();
            partition_->connect_cells();

            repart = true;

            if (print_repart_info_)
            {
                if (comm_->rank() == 0)
                {
                    std::ofstream out;
                    std::string s = "asm";
                    s.append("-repart.dat");
                    out.open(s, std::ios_base::app);
                    out << nassemble_ << std::endl;
                    out.close();
                }
            }
        };

        double cratio = imbalance(partition_->spc());

        if (print_imbalance_)
        {
            if (comm_->rank() == 0)
            {
                std::ofstream out;
                std::string s = "asm";
                s.append("-imbalance.dat");
                out.open(s, std::ios_base::app);
                out << nassemble_ << " " << cratio << std::endl;
                out.close();
            }
        }

        if (can_rebalance_ || force_rebalance_)
        {
            if (force_rebalance_)
            {
                aa();
            }
            else
            {
                assert(false);
                //if (tl > rebalance_thres_)
                {
                    aa();
                }
            }
        }

        return repart;
    }

    Assembler::~Assembler()
    {
        delete partition_;
    }

    void Assembler::assemble()
    {
        assert(!partition_->spc().sp().empty());

        partition_->print_cell_dist(comm_, nassemble_);

        donor_search();
        ++nassemble_;

        if (print_map_)
        {
            partition_->spc().global_rm().print("asm");
        }

        if (print_vtk_)
        {
            for (const auto& sp: partition_->spc().sp())
            {
                for (const auto& m: sp.mesh())
                {
                    std::string fn = "asm-";
                    fn.append(std::to_string(m.tag()()));
                    fn.append("-");
                    fn.append(std::to_string(comm_->rank()));
                    fn.append("-");
                    fn.append(std::to_string(sp.tag()()));
                    fn.append(".vtk");
                    m.print_as_vtk_geometry(fn);
                }
            }

            for (const auto& sp: partition_->spc().sp())
            {
                for (const auto& m: sp.mesh())
                {
                    std::string fn = "asmwall-";
                    fn.append(std::to_string(m.tag()()));
                    fn.append("-");
                    fn.append(std::to_string(comm_->rank()));
                    fn.append("-");
                    fn.append(std::to_string(sp.tag()()));
                    fn.append(".vtk");
                    m.print_wall_as_vtk(fn);
                }
            }
        }
    }

    void Assembler::reconnectivity()
    {
        if (profiler_ != nullptr) {profiler_->start("asm-cae");}
        partition_->spc_->connect_after_exchange(nullptr, "asm");
        if (profiler_ != nullptr) {profiler_->stop("asm-cae");}
    }

    void Assembler::reset_oga_status()
    {
        if (profiler_ != nullptr) {profiler_->start("asm-reset-oga");}
        partition_->spc_->reset_all_oga_status();
        if (profiler_ != nullptr) {profiler_->stop("asm-reset-oga");}
    }

    void Assembler::exchange()
    {
        if (comm_->size() == 0) {
            return;
        }

        assert(nassemble_ != 0);

        auto& spc = partition_->spc_;

        CellExchanger cell_exc(spc, comm_, mergebins_);
        cell_exc.exchange(false, "asm-cell-exc", profiler_);

        mem_usage(comm_, "asm-exc");

        if (profiler_ != nullptr) {profiler_->start("asm-dis");}
        Disconnecter disconnecter(spc, &cell_exc, nullptr);
        disconnecter.disconnect(false, comm_->rank());
        if (profiler_ != nullptr) {profiler_->stop("asm-dis");}

        mem_usage(comm_, "asm-dis");

        if (profiler_ != nullptr) {profiler_->start("asm-arr");}
        ArrivalConnecter connecter(spc, &cell_exc.arrival(), nullptr);
        connecter.add_arrival_cells();
        if (profiler_ != nullptr) {profiler_->stop("asm-arr");}

        mem_usage(comm_, "asm-arr");
    }

    void Assembler::donor_search()
    {
        assert(!partition_->spc().sp().empty());
        auto& spc = partition_->spc_;

        DonorSearcher donor_searcher(comm_, spc, profiler_, donor_search_algo_, mergebins_, false, print_ds_info_);
        //if (profiler_ != nullptr) {profiler_->start("asm-ds");}
        donor_searcher.donor_search(nassemble_, pseudo3D_);
        //if (profiler_ != nullptr) {profiler_->stop("asm-ds");}

        if (profiler_ != nullptr) {profiler_->start("asm-con-undeftofield");}
        donor_searcher.convert_undefined_to_field();
        if (profiler_ != nullptr) {profiler_->stop("asm-con-undeftofield");}

        if (profiler_ != nullptr) {profiler_->start("asm-upd-prev-donor");}
        spc->update_prev_cand_donor();
        if (profiler_ != nullptr) {profiler_->stop("asm-upd-prev-donor");}

        CandReqExchanger cand_req_exc(&(partition_->spc()), comm_, mergebins_);
        cand_req_exc.exchange(false, "asm-candreq-exc", profiler_);
        //cand_req_exc.exchange();

        CandDonorExchanger cand_donor_exc(&(partition_->spc()), &cand_req_exc.arrival(), comm_);
        cand_donor_exc.exchange(false, "asm-canddonor-exc", profiler_);
        //cand_donor_exc.exchange();

        if (profiler_ != nullptr) {profiler_->start("asm-ds-rtof");}
        donor_searcher.receptor_to_field(cand_donor_exc.arrival());
        if (profiler_ != nullptr) {profiler_->stop("asm-ds-rtof");}
        //if (profiler_ != nullptr) {profiler_->start("asm-ds-handle");}
        donor_searcher.handle_donor_conflict(cand_donor_exc.arrival());
        //if (profiler_ != nullptr) {profiler_->stop("asm-ds-handle");}

        //donor_searcher.check_donor_validity();

        //if (profiler_ != nullptr) {profiler_->bstart("asm-canddonor-upd");}
        cand_donor_exc.update(profiler_, "asm-canddonor-upd");
        //if (profiler_ != nullptr) {profiler_->bstop("asm-canddonor-upd");}

        if (profiler_ != nullptr) {profiler_->start("asm-ds-best");}
        donor_searcher.determine_best_donor(cand_donor_exc.arrival());
        if (profiler_ != nullptr) {profiler_->stop("asm-ds-best");}

        //donor_searcher.check_donor_validity();

        if (profiler_ != nullptr) {profiler_->start("asm-ds-orphan");}
        donor_searcher.determine_orphan(cand_donor_exc.arrival());
        if (profiler_ != nullptr) {profiler_->stop("asm-ds-orphan");}

        if (profiler_ != nullptr) {profiler_->start("asm-con-receptohole");}
        donor_searcher.convert_receptor_to_hole();
        if (profiler_ != nullptr) {profiler_->stop("asm-con-receptohole");}

        //donor_searcher.check_donor_validity();

        //for (const auto& m: partition_->spc_->sp().front().mesh())
        //{
            //for (const auto& mc: m.cell())
            //{
                //if (mc.oga_cell_type() == OGA_cell_type_t::non_resident || mc.oga_cell_type() == OGA_cell_type_t::ghost)
                //{
                    //Tag spt;
                    //assert(!partition_->spc_->is_resident(mc.poly().centroid(), mc.tag()(), spt));
                //}
            //}
        //}
    }

    Partition* Assembler::partition()
    {
        return partition_;
    }

    const Partition* Assembler::partition() const
    {
        return partition_;
    }

}
