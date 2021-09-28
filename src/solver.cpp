#include "solver.h"

namespace Tailor
{
    Freestream Solver::fs() const
    {
        return fs_;
    }

    const boost::mpi::communicator *Solver::comm() const
    {
        return comm_;
    }

    void Solver::set_partition(Partition *partition)
    {
        assert(partition != nullptr);
        partition_ = partition;
    }

    void Solver::set_comm(boost::mpi::communicator *comm)
    {
        assert(comm != nullptr);
        assert(comm_ == nullptr);
        comm_ = comm;

        if (partition_ != nullptr)
        {
            partition_->set_comm(comm_);
        }
    }

    void Solver::set_profiler(Profiler *prof)
    {
        profiler_ = prof;

        if (partition_ != nullptr)
        {
            partition_->set_profiler(profiler_);
        }
    }

    bool calc_cell(const MeshCell &mc)
    {
        if (mc.oga_cell_type() != OGA_cell_type_t::field)
        {
            return false;
        }

        return true;
    }

    Solver::Solver() : increase_cfl_(true),
                       cfl_multiplier_(100.),
                       repart_ratio_(0),
                       initratio_(0.),
                       print_imbalance_(false),
                       print_repart_info_(false),
                       rebalance_thres_(0.),
                       can_rebalance_(false),
                       force_rebalance_(false),
                       profiler_(nullptr),
                       //save_solution_(false),
                       //restore_solution_(false),
                       half_cfl_(false),
                       print_map_(false),
                       pseudo3D_(false),
                       omega_(0.),
                       show_inner_res_(false),
                       show_inner_norm_(false),
                       print_residual_(false),
                       sorder_(0),
                       torder_(0),
                       maxtimestep_(0),
                       tol_(0.),
                       nsweep_(0.),
                       cfl_(0.),
                       delta_cfl_(0.),
                       cfl_increase_freq_(0),
                       ncfl_increase_(0),
                       progressive_cfl_(false),
                       dt_(0.),
                       steady_(false),
                       verbose_(false),
                       printfreq_(0),
                       finaltime_(0.),
                       temporal_discretization_(""),
                       nsolve_(0),
                       partition_(nullptr),
                       comm_(nullptr),
                       load_estim_type_(LoadEstim::solver),
                       make_load_balance_(false),
                       var_exc_(nullptr),
                       initial_global_residual_(0.),
                       last_global_residual_(0.),
                       donor_var_exc_(nullptr),
                       riemann_solver_type_(RiemannSolverType::roe),
                       dual_ts_(false),
                       global_nmesh_(0)
    {
    }

    Solver::Solver(boost::mpi::communicator *comm, const std::vector<std::string> &filename, Profiler *profiler, Partition *partition) : verbose_(true), maxtimestep_(10000), comm_(comm), var_exc_(nullptr), donor_var_exc_(nullptr), nsolve_(0), profiler_(profiler), partition_(partition),
    initial_global_residual_(0.),
    last_global_residual_(0.),
    increase_cfl_(true),
    cfl_multiplier_(100.),
    global_nmesh_(0)
    {
        read_settings();

        //if (steady_)
        if (use_local_time_step_)
        {
            dt_ = TAILOR_BIG_POS_NUM;
        } 

        if (!steady_)
        {
            assert(!use_local_time_step_);
        }

        ncfl_increase_ = 0;

        if (partition_ == nullptr)
        {
            partition_ = new Partition(comm_, load_estim_type_, pseudo3D_, "sol", profiler_);
            partition_->read_mesh(filename);

            partition_->make(true, make_load_balance_, nsolve_);

            partition_->clear_mesh();

            if (print_map_)
            {
                partition_->spc_->global_rm().print("sol");
            }
        }

        global_nmesh_ = partition_->spc_->mesh_system_size();

        if (comm_->size() != 1)
        {
            connect_partition_cells();
        }

        if (profiler_ != nullptr)
        {
            profiler_->start("sol-init");
        }
        partition_->spc_->init_flow();
        //partition_->spc_->init_sod();
        if (profiler_ != nullptr)
        {
            profiler_->stop("sol-init");
        }

        //double tl, fl;
        //imbalance(partition_->spc(), tl, fl);
        //double max, min;
        //double  = cell_imbalance(*partition_->spc_, max, min);
        initratio_ = imbalance(partition_->spc());

        if (print_imbalance_)
        {
            if (comm_->rank() == 0)
            {
                std::ofstream out;
                std::string s = "sol";
                s.append("-imbalance.dat");
                out.open(s, std::ios_base::app);
                out << nsolve_ << " " << initratio_ << std::endl;
                out.close();
            }
        }

        if (print_vtk_init_)
        {
            print_mesh_vtk("sol-init");
        }

        //for (const auto& sp: partition_->spc().sp())
        //{
        //    for (const auto& m: sp.mesh())
        //    {
        //        std::string fn = "soll-";
        //        fn.append(std::to_string(m.tag()()));
        //        fn.append("-");
        //        fn.append(std::to_string(comm_->rank()));
        //        fn.append("-");
        //        fn.append(std::to_string(sp.tag()()));
        //        fn.append(".vtk");
        //        m.print_wall_as_vtk(fn);
        //    }
        //}
    }

    void Solver::reconnectivity()
    {
        if (comm_->size() == 1)
        {
            return;
        }

        if (profiler_ != nullptr)
        {
            profiler_->start("sol-cae");
        }
        partition_->spc_->connect_after_exchange(nullptr, "sol");
        if (profiler_ != nullptr)
        {
            profiler_->stop("sol-cae");
        }
        //if (profiler_ != nullptr) {profiler_->start("sol-cpc");}
        connect_partition_cells();
        //if (profiler_ != nullptr) {profiler_->stop("sol-cpc");}
    }

    void Solver::reset_oga_status()
    {
        if (global_nmesh_ == 1)
        {
            return;
        }

        if (profiler_ != nullptr)
        {
            profiler_->start("sol-reset-oga");
        }
        partition_->spc_->reset_all_oga_status();
        if (profiler_ != nullptr)
        {
            profiler_->stop("sol-reset-oga");
        }
    }

    bool Solver::repartition()
    {
        bool repart = false;

        if (nsolve_ < 1)
        {
            return repart;
        }

        auto aa = [&]() {
            {
                std::string mems = "repart-a";
                mem_usage(comm_, mems);
            }

            if (profiler_ != nullptr)
            {
                profiler_->start("sol-rem-ghostnresi");
            }
            partition_->spc_->remove_ghosts();
            partition_->spc_->remove_nonresidents();
            if (profiler_ != nullptr)
            {
                profiler_->stop("sol-rem-ghostnresi");
            }

            {
                std::string mems = "repart-b";
                mem_usage(comm_, mems);
            }

            if (profiler_ != nullptr)
            {
                profiler_->start("sol-reset-status");
            }
            reset_oga_status();
            if (profiler_ != nullptr)
            {
                profiler_->stop("sol-reset-status");
            }

            {
                std::string mems = "repart-c";
                mem_usage(comm_, mems);
            }

            if (profiler_ != nullptr)
            {
                profiler_->start("sol-reset-connect");
            }
            partition_->spc_->reset_mesh_connectivity();
            if (profiler_ != nullptr)
            {
                profiler_->stop("sol-reset-connect");
            }

            {
                std::string mems = "repart-c";
                mem_usage(comm_, mems);
            }

            partition_->set_mesh(std::move(partition_->spc_->sp_.front().mesh_));
            {
                std::string mems = "repart-d";
                mem_usage(comm_, mems);
            }
            partition_->repartition(partition_->mesh(), make_load_balance_, RegType::aabb, true, false, load_estim_type_, nsolve_);
            {
                std::string mems = "repart-e";
                mem_usage(comm_, mems);
            }

            partition_->clear_mesh();
            {
                std::string mems = "repart-f";
                mem_usage(comm_, mems);
            }
            partition_->connect_cells();
            {
                std::string mems = "repart-g";
                mem_usage(comm_, mems);
            }

            if (comm_->size() != 1)
            {
                connect_partition_cells();
            }

            {
                std::string mems = "repart-h";
                mem_usage(comm_, mems);
            }

            repart = true;

            if (print_repart_info_)
            {
                if (comm_->rank() == 0)
                {
                    std::ofstream out;
                    std::string s = "sol";
                    s.append("-repart.dat");
                    out.open(s, std::ios_base::app);
                    out << nsolve_ << std::endl;
                    out.close();
                }
            }
        };

        //double tl, fl;
        //imbalance(partition_->spc(), tl, fl);
        double cratio = imbalance(partition_->spc());
        //tl = partition_->spc().tl_;
        //fl = partition_->spc().fl_;

        if (print_imbalance_)
        {
            if (comm_->rank() == 0)
            {
                std::ofstream out;
                std::string s = "sol";
                s.append("-imbalance.dat");
                out.open(s, std::ios_base::app);
                out << nsolve_ << " " << cratio << std::endl;
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
                //if (fl > rebalance_thres_) // or tl not sure yet
                if ((cratio / initratio_) > repart_ratio_) // or tl not sure yet
                {
                    //if (comm_->rank() == 0)
                    //{
                    //std::cout << "initratio_: " << initratio_ << std::endl;
                    //std::cout << "cratio: " << cratio << std::endl;
                    //}
                    aa();
                }
            }
        }

        return repart;
    }

    void Solver::exchange()
    {
        if (comm_->size() == 1)
        {
            return;
        }

        auto &spc = partition_->spc_;

        if (profiler_ != nullptr)
        {
            profiler_->start("sol-rem-ghost");
        }
        spc->remove_ghosts();
        if (profiler_ != nullptr)
        {
            profiler_->stop("sol-rem-ghost");
        }

        CellExchanger cell_exc(spc, comm_, true);
        cell_exc.exchange(false, "sol-cell-exc", profiler_, -1);
        //cell_exc.exchange();

        if (profiler_ != nullptr)
        {
            profiler_->start("sol-dis");
        }
        Disconnecter disconnecter(spc, &cell_exc, nullptr);
        disconnecter.disconnect(true, comm_->rank());
        if (profiler_ != nullptr)
        {
            profiler_->stop("sol-dis");
        }

        if (profiler_ != nullptr)
        {
            profiler_->start("sol-arr");
        }
        ArrivalConnecter connecter(spc, &cell_exc.arrival(), nullptr);
        connecter.add_arrival_cells();
        if (profiler_ != nullptr)
        {
            profiler_->stop("sol-arr");
        }
    }

    void Solver::rotate(const Tag &mesh, double ang, int axis, const Vector3 &pivot)
    {
        partition_->rotate(mesh, ang, axis, pivot);
    }

    void Solver::move(const Tag& mesh, const Vector3& v)
    {
        partition_->move(mesh, v);
    }

    void Solver::set_oga_cell_type_all_field()
    {
        for (Mesh &m : partition_->spc_->sp_.front().mesh_)
        {
            for (MeshCell &mc : m.cell_)
            {
                if (mc.oga_cell_type() != OGA_cell_type_t::non_resident && mc.oga_cell_type() != OGA_cell_type_t::ghost)
                {
                    mc.set_oga_cell_type(OGA_cell_type_t::field);
                }
            }
        }
        
        for (const MeshCell &mc : partition_->spc().sp().front().mesh().front().cell())
        {
            assert(mc.oga_cell_type() != OGA_cell_type_t::undefined);
        }
    }

    void Solver::transfer_oga_cell_type(const ArrCon<DonorInfo2> &di)
    {
        //assert(!di.empty());
        DonorRequestExchanger donor_req_exc(&partition_->spc(), &(partition_->spc_->sp_.front()), &di, comm_);
        donor_req_exc.exchange(false, "donor-req-exc", profiler_);
        //donor_req_exc.exchange();

        if (donor_var_exc_ != nullptr)
        {
            delete donor_var_exc_;
            donor_var_exc_ = nullptr;
        }

        donor_var_exc_ = new VarExchanger(&(partition_->spc().sp()), &donor_req_exc, comm_);
        donor_var_exc_->exchange(false, "donor-var-exc", profiler_);
        //donor_var_exc_->exchange();
    }

    void Solver::read_mesh(const std::vector<std::string> &file_name)
    {
        partition_->read_mesh(file_name);
    }

    void Solver::connect_partition_cells()
    {
        if (comm_->size() == 0)
        {
            return;
        }

        RequestExchanger req_exc(&partition_->spc(), comm_);
        //FaceRequestExchanger facereq_exc(&partition_->spc(), comm_);
        //facereq_exc.exchange(false, "req", profiler_);
        req_exc.exchange(false, "sol-cpc-req", profiler_);
        //req_exc.exchange();

        //NeiExchanger nei_exc(&(partition_->spc()), &facereq_exc.receiver(), comm_);
        NeiExchanger nei_exc(&(partition_->spc()), &req_exc.arrival(), comm_);
        nei_exc.exchange(false, "sol-cpc-nei", profiler_);
        //nei_exc.exchange();

        if (profiler_ != nullptr)
        {
            profiler_->start("sol-cpc-dup");
        }
        nei_exc.remove_dup();
        if (profiler_ != nullptr)
        {
            profiler_->stop("sol-cpc-dup");
        }

        if (profiler_ != nullptr)
        {
            profiler_->start("sol-cpc");
        }
        partition_->spc_->connect_partition_cells(nei_exc.arrival());
        if (profiler_ != nullptr)
        {
            profiler_->stop("sol-cpc");
        }
    }

    void Solver::read_settings()
    {
        namespace po = boost::program_options;

        po::options_description all_options;

        po::options_description desc{"Solver options"};
        desc.add_options()
            ("solver.increase_cfl", po::value<bool>()->default_value(false), "")
            ("solver.repart-ratio", po::value<int>()->default_value(100), "")
            ("solver.cfl_multiplier", po::value<double>()->default_value(10.), "")
            ("solver.cfl_ratio", po::value<double>()->default_value(10.), "")
            ("solver.rebalance-thres", po::value<double>()->default_value(40.), "")
            ("solver.show_inner_res", po::value<bool>()->default_value(true), "Show inner loop residual")
            ("solver.show_inner_norm", po::value<bool>()->default_value(true), "Show inner loop norm")
            ("solver.print_residual", po::value<bool>()->default_value(true), "Print outer norm")
            ("solver.steady", po::value<bool>()->default_value(false), "Steady state")
            ("solver.progressive_cfl", po::value<bool>()->default_value(false), "Progressive CFL increase")
            ("solver.temporal_discretization", po::value<std::string>()->default_value("forward_euler"), "Temporal discretization")
            ("solver.dt", po::value<double>()->default_value(0.001), "Real time step")
            ("solver.nsweep", po::value<int>()->default_value(1), "Number of sweeps in SOR")
            ("solver.omega", po::value<double>()->default_value(1), "Relaxation parameter in SOR")
            ("solver.tol", po::value<double>()->default_value(1e-12), "Error tolerance")
            ("solver.sorder", po::value<int>()->default_value(1), "Spatial order of accuracy")
            ("solver.torder", po::value<int>()->default_value(1), "Temporal order of accuracy")
            ("solver.printfreq", po::value<int>()->default_value(100), "Printting frequency")
            ("solver.cfl", po::value<double>()->default_value(0.5), "CFL number")
            ("solver.delta_cfl", po::value<double>()->default_value(0.), "Delta CFL number for progressive cfl increase")
            ("solver.cfl_increase_freq", po::value<int>()->default_value(5), "CFL increase frequency for progressive cfl increase")
            ("solver.finaltime", po::value<double>()->default_value(0.5), "Final time")
            ("solver.make-load-balance", po::value<bool>()->default_value(true), "Make load balance (for area or hybrid load estimation)")
            ("solver.load-estim", po::value<int>()->default_value(2), "Load estimation method")
            ("general.pseudo3D", po::value<bool>()->default_value(false), "2.5D simulation with 1 layer in depth")
            ("solver.print-map", po::value<bool>()->default_value(false), "Print map.")
            ("solver.max-time-step", po::value<int>()->default_value(10000), "")
            ("solver.half-cfl", po::value<bool>()->default_value(true), "")
            //("solver.restore-solution", po::value<bool>()->default_value(false), "")
            ("solver.can-rebalance", po::value<bool>()->default_value(true), "Make load rebalance if needed.")
            ("solver.force-rebalance", po::value<bool>()->default_value(false), "Force load rebalance even if relabalance is not needed.")
            ("solver.print-vtk-init", po::value<bool>()->default_value(false), "")
            ("solver.print-repart-info", po::value<bool>()->default_value(false), "")
            ("solver.print-imbalance", po::value<bool>()->default_value(false), "")
            ("solver.riemann-solver", po::value<int>()->default_value(0), "")
            ("solver.dual-ts", po::value<bool>()->default_value(false), "")
            ("solver.flow-init-type", po::value<int>()->default_value(0), "")
            ("solver.use-local-time-step", po::value<bool>()->default_value(false), "")
            ("solver.print-vtk-only-last-step", po::value<bool>()->default_value(true), "")
            ("solver.print-vtk-every-step", po::value<bool>()->default_value(false), "")
            ;

        po::options_description linear_solver{"Linear solver options"};
        linear_solver.add_options()
            ("linear-solver.max-iteration", po::value<int>()->default_value(0), "")
            ("linear-solver.max-restart", po::value<int>()->default_value(0), "")
            ("linear-solver.abs-error", po::value<double>()->default_value(0), "")
            ("linear-solver.rel-error", po::value<double>()->default_value(0), "")
            ("linear-solver.print-error", po::value<bool>()->default_value(false), "")
            ;

        all_options.add(desc);
        all_options.add(linear_solver);

        std::ifstream settings_file("settings.ini");

        boost::program_options::variables_map vm;
        po::store(po::parse_config_file(settings_file, all_options, true), vm);
        po::notify(vm);

        increase_cfl_ = vm["solver.increase_cfl"].as<bool>();
        print_residual_ = vm["solver.print_residual"].as<bool>();
        show_inner_res_ = vm["solver.show_inner_res"].as<bool>();
        show_inner_norm_ = vm["solver.show_inner_norm"].as<bool>();
        progressive_cfl_ = vm["solver.progressive_cfl"].as<bool>();
        steady_ = vm["solver.steady"].as<bool>();
        temporal_discretization_ = vm["solver.temporal_discretization"].as<std::string>();
        dt_ = vm["solver.dt"].as<double>();
        nsweep_ = vm["solver.nsweep"].as<int>();
        cfl_multiplier_ = vm["solver.cfl_multiplier"].as<double>();
        cfl_ratio_ = vm["solver.cfl_ratio"].as<double>();
        omega_ = vm["solver.omega"].as<double>();
        tol_ = vm["solver.tol"].as<double>();
        sorder_ = vm["solver.sorder"].as<int>();
        torder_ = vm["solver.torder"].as<int>();
        printfreq_ = vm["solver.printfreq"].as<int>();
        cfl_ = vm["solver.cfl"].as<double>();
        delta_cfl_ = vm["solver.delta_cfl"].as<double>();
        cfl_increase_freq_ = vm["solver.cfl_increase_freq"].as<int>();
        finaltime_ = vm["solver.finaltime"].as<double>();
        load_estim_type_ = static_cast<LoadEstim>(vm["solver.load-estim"].as<int>());
        make_load_balance_ = vm["solver.make-load-balance"].as<bool>();
        pseudo3D_ = vm["general.pseudo3D"].as<bool>();
        print_map_ = vm["solver.print-map"].as<bool>();
        maxtimestep_ = vm["solver.max-time-step"].as<int>();
        half_cfl_ = vm["solver.half-cfl"].as<bool>();
        //restore_solution_ = vm["solver.restore-solution"].as<bool>();
        can_rebalance_ = vm["solver.can-rebalance"].as<bool>();
        force_rebalance_ = vm["solver.force-rebalance"].as<bool>();
        rebalance_thres_ = vm["solver.rebalance-thres"].as<double>();
        print_vtk_init_ = vm["solver.print-vtk-init"].as<bool>();
        print_repart_info_ = vm["solver.print-repart-info"].as<bool>();
        print_imbalance_ = vm["solver.print-imbalance"].as<bool>();
        repart_ratio_ = vm["solver.repart-ratio"].as<int>();
        if (vm["solver.riemann-solver"].as<int>() == 0) {
            riemann_solver_type_ = RiemannSolverType::roe;
        }
        else if (vm["solver.riemann-solver"].as<int>() == 1) {
            riemann_solver_type_ = RiemannSolverType::hllc;
            assert(temporal_discretization_ != "backward_euler");
        }
        else {
            assert(false);
        }
        dual_ts_ = vm["solver.dual-ts"].as<bool>();
        use_local_time_step_ = vm["solver.use-local-time-step"].as<bool>();
        linear_solver_max_iteration_ = vm["linear-solver.max-iteration"].as<int>();
        linear_solver_max_restart_ = vm["linear-solver.max-restart"].as<int>();
        linear_solver_abs_error_ = vm["linear-solver.abs-error"].as<double>();
        linear_solver_rel_error_ = vm["linear-solver.rel-error"].as<double>();
        print_linear_solver_error_ = vm["linear-solver.print-error"].as<bool>();
        print_vtk_only_last_step_ = vm["solver.print-vtk-only-last-step"].as<bool>();
        print_vtk_every_step_ = vm["solver.print-vtk-every-step"].as<bool>();

        if (print_vtk_only_last_step_)
        {
            assert(!print_vtk_every_step_);
        }
    }

    double Solver::dt() const
    {
        return dt_;
    }

    std::tuple<Matrix5, Matrix5> get_rotation_matrix(const Vector3& normal)
    {
        Matrix5 rotation_matrix = make_rot_matrix(normal);
        Matrix5 inv_rotation_matrix = rotation_matrix.transpose(); // (inverve = transpose) since coordinate frames are orthogonal.

        return std::make_tuple(rotation_matrix, inv_rotation_matrix);
    }

    std::tuple<State, State> left_and_right_states(const Vector5& left_conservative_var, const Vector5& right_conservative_var, double gamma, const Matrix5& rotation_matrix, const Vector3& face_velocity)
    {
        State left(left_conservative_var , gamma, rotation_matrix, face_velocity);
        State right(right_conservative_var, gamma, rotation_matrix, face_velocity);

        return std::make_tuple(left, right);
    }

    std::tuple<Vector5, double, Matrix5> compute_flux(RiemannSolverType riemann_solver_type, double face_area, const Matrix5& inverse_rotation_matrix, const State& rotated_left_state, const State& rotated_right_state, double gamma, bool calculate_roe_jacobian)
    {
        Vector5 flux;
        double max_eigen;
        Matrix5 Aroe;

        RiemannSolver riemann_solver(riemann_solver_type, rotated_left_state, rotated_right_state, face_area, gamma, max_eigen, flux, Aroe, calculate_roe_jacobian, SpeedEstimateHLLC::roe);

        flux = inverse_rotation_matrix * flux;

        assert(!flux.isnan());

        return std::make_tuple(flux, max_eigen, Aroe);
    }

    void Solver::compute_sum_of_fluxes(Mesh &mesh)
    {
        double gamma = fs_.gamma_;

        mesh.reset_R_checked();
        mesh.reset_R();
        mesh.reset_D();

        for (MeshCell &mc : mesh.cell_)
        {
            assert(mc.oga_cell_type() != OGA_cell_type_t::undefined);
        }

        for (auto &mc : mesh.cell_)
        {
            mc.sumarea_ = 0.;
        }

        for (auto &mc : mesh.cell_)
        {
            if (!calc_cell(mc))
            {
                continue;
            }

            for (auto &mf : mc.face_p())
            {
                if (mf.R_checked())
                {
                    continue;
                }

                if (mf.btype() == BouType::empty)
                {
                    continue;
                }

                if (comm_->size() != 1)
                {
                    if (mf.parent_cell().size() == 1) // should be face of nonresi
                    {
                        assert(!mf.is_boundary());
                        continue;
                    }
                }

                Vector3 normal = mf.face().normal();
                Vector3 face_velocity = mf.vf();
                auto face_area = mf.face().signed_area();

                auto [left_cell, right_cell] = left_and_right_cells(mesh, mf, mc.tag());

                Vector5 left_cons;
                Vector5 right_cons;

                if (sorder_ == 1)
                {
                    left_cons = left_cell->cons_sp1();
                    right_cons = right_cell->cons_sp1();
                }
                else if (sorder_ == 2)
                {
                    //apply_limiter(mesh, *left_cell , mf);
                    //apply_limiter(mesh, *right_cell, mf);

                    left_cons = apply_limiter(mesh, *left_cell, mf);
                    right_cons = apply_limiter(mesh, *right_cell, mf);
                }

                MeshFace *commonface = nullptr;
                if (!mf.is_boundary())
                {
                    commonface = mesh.common_face(*right_cell, mf.tag());
                    assert(commonface != nullptr);
                }

                auto [rotation_matrix, inv_rotation_matrix] = get_rotation_matrix(normal);
                //auto [rotated_left_state, rotated_right_state] = left_and_right_states(left_cell->cons_sp1_, right_cell->cons_sp1_, gamma, rotation_matrix, face_velocity);
                auto [rotated_left_state, rotated_right_state] = left_and_right_states(left_cons, right_cons, gamma, rotation_matrix, face_velocity);
                bool calculate_roe_jacobian = false;
                if (temporal_discretization_ == "backward_euler")
                {
                    calculate_roe_jacobian = true;
                }
                auto [flux, max_eigen, Aroe] = compute_flux(riemann_solver_type_, face_area, inv_rotation_matrix, rotated_left_state, rotated_right_state, gamma, calculate_roe_jacobian);

                if (mf.btype() != BouType::empty)
                {
                    double t = std::abs(face_area * max_eigen);
                    left_cell->sumarea_ += t;
                    right_cell->sumarea_ += t;
                }

                left_cell->R_ -= flux;
                right_cell->R_ += flux;

                left_cell->R_mid_ -= flux;
                right_cell->R_mid_ += flux;

                if (temporal_discretization_ == "backward_euler")
                {
                    update_matrices(&mf, commonface, *left_cell, *right_cell, face_area, face_velocity, gamma, rotation_matrix, inv_rotation_matrix, Aroe, rotated_left_state, rotated_right_state);
                }

                if (commonface != nullptr)
                {
                    commonface->set_R_checked(true);
                }
            }
        }
    }

    void Solver::evolve_solution_in_time(Mesh& mesh)
    {
        for (MeshCell &mc : mesh.cell_)
        {
            if (!calc_cell(mc)) {
                assert(mc.dQ_ == 0.);
                continue;
            }

            if (dual_ts_)
            {
                mc.cons_sp1_ = mc.cons_s_ + mc.dQ_;
            }
            else
            {
                mc.cons_sp1_ = mc.cons_n_ + mc.dQ_;
            }

            mc.prim_ = cons_to_prim(mc.cons_sp1_, fs_.gamma_);

            if (mc.cons_sp1_(4) < 0.)
            {
                std::cout << "cons[0]: " << mc.cons_sp1_(0) << std::endl;
                std::cout << "cons[1]: " << mc.cons_sp1_(1) << std::endl;
                std::cout << "cons[2]: " << mc.cons_sp1_(2) << std::endl;
                std::cout << "cons[3]: " << mc.cons_sp1_(3) << std::endl;
                std::cout << "cons[4]: " << mc.cons_sp1_(4) << std::endl;

                std::cout << "oldcons[0]: " << mc.cons_s_(0) << std::endl;
                std::cout << "oldcons[1]: " << mc.cons_s_(1) << std::endl;
                std::cout << "oldcons[2]: " << mc.cons_s_(2) << std::endl;
                std::cout << "oldcons[3]: " << mc.cons_s_(3) << std::endl;
                std::cout << "oldcons[4]: " << mc.cons_s_(4) << std::endl;

                std::cout << "dQ[0]: " << mc.dQ_(0) << std::endl;
                std::cout << "dQ[1]: " << mc.dQ_(1) << std::endl;
                std::cout << "dQ[2]: " << mc.dQ_(2) << std::endl;
                std::cout << "dQ[3]: " << mc.dQ_(3) << std::endl;
                std::cout << "dQ[4]: " << mc.dQ_(4) << std::endl;
            }
        }   
    }

    void Solver::evolve_old_solution_in_time(Mesh& mesh, int runge_kutta_stage)
    {
        if (temporal_discretization_ == "runge_kutta_4")
        {
            if (runge_kutta_stage == 3)
            {
                for (MeshCell &mc : mesh.cell_)
                {
                    mc.cons_s_ = mc.cons_sp1_;
                }
            }

            return;
        }

        for (MeshCell &mc : mesh.cell_)
        {
            if (dual_ts_)
            {
                mc.cons_s_ = mc.cons_sp1_;
            }
            else if (steady_)
            {
                mc.cons_nm1_ = mc.cons_n_;
                mc.cons_n_ = mc.cons_sp1_;
            }
        }
    }

    void Solver::calc_change_in_conserved_var(Mesh& mesh, int runge_kutta_stage)
    {
        for (MeshCell &mc : mesh.cell_)
        {
            mc.dQ_ = 0.;
        }

        if (temporal_discretization_ == "backward_euler")
        {
            linear_solver(mesh);
        }
        else if (temporal_discretization_ == "forward_euler")
        {
            for (MeshCell &mc : mesh.cell_)
            {
                if (!calc_cell(mc)) {
                    continue;
                }

                double volume = mc.poly().volume();
                mc.dQ_ = mc.R_ / volume;
                        //std::cout << "R: " << mc.R_(0) << std::endl;
                        //std::cout << "R after: " << mc.R_(0) << std::endl;
                        //std::cout << "volume: " << volume << std::endl;
                        //std::cout << "cons_n(0)" << mc.cons_n(0) << std::endl;
                        //std::cout << "cons_nm1(0)" << mc.cons_nm1(0) << std::endl;
                        //assert(false);
                        //std::cout << "dQ_ before: " << mc.dQ_(0) << std::endl;

                //if (dual_ts_ || steady_) {
                if (dual_ts_ || use_local_time_step_) {
                    mc.dQ_ *= mc.dtao_;
                }
                else
                {
                    if (torder_ == 1)
                    {
                        mc.dQ_ *= dt_;
                        //std::cout << "dt: " << dt_ << std::endl;
                        //std::cout << "dQ_: " << mc.dQ_(0) << std::endl;
                        //for (int j = 0; j < NVAR; ++j)
                        //{
                            //assert(std::abs(mc.dQ_(j)) > 1e-6);
                        //}
                    }
                    else if (torder_ == 2)
                    {
                        mc.dQ_ *= (2. / 3.) * dt_;
                    }
                }
            }
        }
        else if (temporal_discretization_ == "runge_kutta_4")
        {
            for (MeshCell &mc : mesh.cell_)
            {
                if (!calc_cell(mc))
                {
                    continue;
                }

                double volume = mc.poly().volume();

                runge_kutta_coef_[runge_kutta_stage] = mc.dtao_ * mc.R_ / volume;

                if (runge_kutta_stage < 2)
                {
                    mc.dQ_ = 0.5 * runge_kutta_coef_[runge_kutta_stage];
                }
                else if (runge_kutta_stage == 2)
                {
                    mc.dQ_ = runge_kutta_coef_[runge_kutta_stage];
                }
                else
                {
                    mc.dQ_ = runge_kutta_coef_[0] / 6. + runge_kutta_coef_[1] / 3. + runge_kutta_coef_[2] / 3. + runge_kutta_coef_[3] / 6.;
                }
            }
        }
    }

    void Solver::init_partitioned_mesh_exchanger()
    {
        if (comm_->size() != 1)
        {
            NonResiExchanger nonresi_exc(&(partition_->spc()), comm_);
            nonresi_exc.exchange(false, "sol-nonresi", profiler_);

            if (var_exc_ != nullptr)
            {
                delete var_exc_;
                var_exc_ = nullptr;
            }

            var_exc_ = new VarExchanger(&(partition_->spc().sp()), &nonresi_exc, comm_);
            var_exc_->exchange(false, "sol-varexc", profiler_);
        }
    }

    void Solver::reset_overset_mesh_exchanger()
    {
        if (donor_var_exc_ != nullptr)
        {
            delete donor_var_exc_;
            donor_var_exc_ = nullptr;
        }
    }

    void Solver::reset_partitioned_mesh_exchanger()
    {
        if (var_exc_ != nullptr)
        {
            delete var_exc_;
            var_exc_ = nullptr;
        }
    }

    int Solver::nsolve() const
    {
        return nsolve_;
    }

    void Solver::solve(int max_time_step)
    {
        std::cout << "Solving" << std::endl;
        partition_->print_cell_dist(comm_, nsolve_);

        if (nsolve_ == 0)
        {
            fs_.read();
        }

        for (const MeshCell &mc : partition_->spc().sp().front().mesh().front().cell())
        {
            assert(mc.oga_cell_type() != OGA_cell_type_t::undefined);
        }

        init_partitioned_mesh_exchanger();
        calc_mesh_velocities();
        compute_gradient_coef();
        init_old_conservative_var();

        for (const MeshCell &mc : partition_->spc().sp().front().mesh().front().cell())
        {
            assert(mc.oga_cell_type() != OGA_cell_type_t::undefined);
        }

        auto residual = non_linear_iteration();
        print_residual(residual);
        if (print_vtk_only_last_step_)
        {
            if (nsolve_ == max_time_step)
            {
                print_mesh_vtk("sol");
            }
        }
        if (print_vtk_every_step_)
        {
            print_mesh_vtk("sol");
        }
        reset_overset_mesh_exchanger();
        reset_partitioned_mesh_exchanger();

        if (!steady_) {
            ++nsolve_;
        }
    }

    void Solver::print_mesh_vtk(std::string fn0)
    {
        for (const auto &sp : partition_->spc().sp())
        {
            for (const auto &m : sp.mesh())
            {
                std::string fn = fn0;
                fn.append("-");
                fn.append(std::to_string(m.tag()()));
                fn.append("-");
                fn.append(std::to_string(comm_->rank()));
                fn.append("-");
                fn.append(std::to_string(nsolve_));
                fn.append(".vtk");
                m.print_as_vtk(fn);
            }
        }
        for (const auto &sp : partition_->spc().sp())
        {
            for (const auto &m : sp.mesh())
            {
                std::string fn = fn0;
                fn.append("-");
                fn.append("wall-");
                fn.append(std::to_string(m.tag()()));
                fn.append("-");
                fn.append(std::to_string(comm_->rank()));
                fn.append("-");
                fn.append(std::to_string(nsolve_));
                fn.append(".vtk");
                m.print_wall_as_vtk(fn);
            }
        }
    }

    //void Solver::print_mesh_vtk()
    //{
    //    if (print_vtk_)
    //    {
    //        for (const auto &sp : partition_->spc().sp())
    //        {
    //            for (const auto &m : sp.mesh())
    //            {
    //                std::string fn = "sol-";
    //                fn.append(std::to_string(m.tag()()));
    //                fn.append("-");
    //                fn.append(std::to_string(comm_->rank()));
    //                fn.append("-");
    //                fn.append(std::to_string(sp.tag()()));
    //                fn.append(".vtk");
    //                m.print_as_vtk(fn);
    //            }
    //        }
    //        for (const auto &sp : partition_->spc().sp())
    //        {
    //            for (const auto &m : sp.mesh())
    //            {
    //                std::string fn = "solwall-";
    //                fn.append(std::to_string(m.tag()()));
    //                fn.append("-");
    //                fn.append(std::to_string(comm_->rank()));
    //                fn.append("-");
    //                fn.append(std::to_string(sp.tag()()));
    //                fn.append(".vtk");
    //                m.print_wall_as_vtk(fn);
    //            }
    //        }
    //    }
    //}

    //void Solver::restore_solution()
    //{
    //    // useless since boost serialization.
    //    for (auto &mesh : partition_->spc_->sp_.front().mesh_)
    //    {
    //        std::string fn = "save-";
    //        fn.append(std::to_string(comm_->rank()));
    //        fn.append("-");
    //        fn.append(std::to_string(mesh.tag()()));
    //        fn.append(".sav");

    //        std::ifstream in;
    //        in.open(fn);
    //        assert(in.is_open());

    //        for (auto &mc : mesh.cell_)
    //        {
    //            in >> mc.cons_sp1_(0);
    //            in >> mc.cons_sp1_(1);
    //            in >> mc.cons_sp1_(2);
    //            in >> mc.cons_sp1_(3);
    //            in >> mc.cons_sp1_(4);

    //            mc.prim_ = cons_to_prim(mc.cons_sp1(), fs_.gamma_);
    //        }

    //        in.close();
    //    }
    //}

    //void Solver::save_solution()
    //{
    //    // useless since boost serialization.
    //    assert(false);
    //    for (const auto &mesh : partition_->spc().sp().front().mesh())
    //    {
    //        std::string fn = "save-";
    //        fn.append(std::to_string(comm_->rank()));
    //        fn.append("-");
    //        fn.append(std::to_string(mesh.tag()()));
    //        fn.append(".sav");

    //        std::ofstream out;
    //        out.open(fn);

    //        for (const auto &mc : mesh.cell())
    //        {
    //            out << mc.cons_sp1(0) << " " << mc.cons_sp1(1) << " " << mc.cons_sp1(2) << " " << mc.cons_sp1(3) << " " << mc.cons_sp1(4) << "\n";
    //        }

    //        out.close();
    //    }
    //}

    void Solver::first_order_residual(Vector5& res, const MeshCell& mc)
    {
        for (int i = 0; i < mc.R_.nelm(); ++i)
        {
            res(i) = std::max(res(i), std::abs(mc.R_(i) - mc.poly().volume() * mc.dQ_(i) / dt_));
        }
    }

    Vector5 Solver::compute_residual(Mesh& mesh)
    {
        Vector5 res(TAILOR_BIG_NEG_NUM);

        for (MeshCell &mc : mesh.cell_)
        {
            if (!calc_cell(mc))
            {
                continue;
            }

            if (steady_ || dual_ts_)
            {
                for (int i = 0; i < mc.R_.nelm(); ++i)
                {
                    res(i) = std::max(res(i), std::abs(mc.R_(i)));
                }
            }
            else
            {
                if (torder_ == 1)
                {
                    first_order_residual(res, mc);
                }
                else if (torder_ == 2)
                {
                    //if (dual_ts_)
                    //{
                        //first_order_residual(res, mc);
                    //}
                    //else
                    {
                        if (nsolve_ > 1)
                        {
                            for (int i = 0; i < mc.R_.nelm(); ++i)
                            {
                                res(i) = std::max(res(i), std::abs(mc.R_(i) - 0.5 * mc.poly().volume() * (3. * mc.cons_sp1_(i) - 4. * mc.cons_n_(i) + mc.cons_nm1_(i)) / dt_));
                            }
                        }
                        else
                        {
                            first_order_residual(res, mc);
                        }
                    }
                }
            }
        }

        return res;
    }

    void Solver::oga_interpolate(Mesh& mesh)
    {
        Freestream fs;
        fs.read();

        auto &meshes = partition_->spc_->sp_.front().mesh_;

        assert(!meshes.empty());

        for (MeshCell &mc : mesh.cell_)
        {
            if (mc.oga_cell_type() == OGA_cell_type_t::receptor || mc.oga_cell_type() == OGA_cell_type_t::mandat_receptor)
            {
                assert(mc.donor().mesh_tag_.isvalid());
                assert(mc.donor().cell_tag_.isvalid());

                auto m = std::find_if(meshes.begin(), meshes.end(), [&](const auto &mm) { return mm.tag() == mc.donor().mesh_tag_; });
                assert(m != meshes.end());
                assert(m->tag() != mesh.tag());
                const auto &donor_cell = m->cell(mc.donor().cell_tag_);

                if (sorder_ == 2)
                {
                    //if (donor_cell.oga_cell_type() != OGA_cell_type_t::field || donor_cell.oga_cell_type() != OGA_cell_type_t::non_resident)
                    //{
                        //std::cout << "mc oga: " << static_cast<int>(donor_cell.oga_cell_type()) << std::endl;
                    //}

                    //assert(donor_cell.oga_cell_type() == OGA_cell_type_t::field || donor_cell.oga_cell_type() == OGA_cell_type_t::non_resident);
                    assert(donor_cell.oga_cell_type() == OGA_cell_type_t::field);

                    auto grad = gradient_.ls_grad(*m, donor_cell);

                    Limiter limiter(LimiterType::venka);
                    limiter.limit(mesh, mc, grad);

                    auto d = mc.poly().centroid() - donor_cell.poly().centroid();

                    for (int i = 0; i < NVAR; ++i)
                    {
                        mc.prim_(i) = donor_cell.prim(i) + limiter(i) * dot(grad[i], d);
                    }

                    //mc.prim_(1) += mc.vgn(0);
                    //mc.prim_(2) += mc.vgn(1);
                    //mc.prim_(3) += mc.vgn(2);

                    mc.cons_sp1_ = prim_to_cons(mc.prim_, fs.gamma_);
                }
                else
                {
                    mc.prim_ = donor_cell.prim();

                    //mc.prim_(1) += mc.vgn(0);
                    //mc.prim_(2) += mc.vgn(1);
                    //mc.prim_(3) += mc.vgn(2);

                    mc.dQ_ = donor_cell.dQ(); // ?
                    mc.cons_sp1_ = prim_to_cons(mc.prim_, fs.gamma_);
                    //for (int i = 0; i < NVAR; ++i)
                    //{
                        //mc.prim_[i] = donor_cell.prim(i);
                        //if (mesh.tag()() == 0)
                        //{
                        //std::cout << "u: " << donor_cell.prim(1) << std::endl;
                        //std::cout << "v: " << donor_cell.prim(2) << std::endl;
                        //std::cout << "w: " << donor_cell.prim(3) << std::endl;
                        //}
                    //}
                }

                //auto vg = donor_cell.vgn();
                //assert(vg(0) == 0.);
                //assert(vg(1) == 0.);
                //assert(vg(2) == 0.);
                //mc.prim_(1) -= vg(0);
                //mc.prim_(2) -= vg(1);
                //mc.prim_(3) -= vg(2);

                assert(!mc.prim_.isnan());
            }
        }
    }

    void Solver::print_settings() const
    {
        if (comm_->rank() != 0) {
            return;
        }

        std::ofstream out;
        out.open("solver_settings.log");

        out << "cfl_multiplier = " << cfl_multiplier_ << std::endl;
        out << "cfl_ratio = " << cfl_ratio_ << std::endl;
        out << "increase_cfl = " << increase_cfl_ << std::endl;
        out << "print_residual = " << print_residual_ << std::endl;
        out << "show_inner_res = " << show_inner_res_ << std::endl;
        out << "show_inner_norm = " << show_inner_norm_ << std::endl;
        out << "progressive_cfl = " << progressive_cfl_ << std::endl;
        out << "steady = " << steady_ << std::endl;
        out << "temporal_discretization = " << temporal_discretization_ << std::endl;
        out << "dt = " << dt_ << std::endl;
        out << "nsweep = " << nsweep_ << std::endl;
        out << "omega = " << omega_ << std::endl;
        out << "tol = " << tol_ << std::endl;
        out << "sorder = " << sorder_ << std::endl;
        out << "torder = " << torder_ << std::endl;
        out << "printfreq = " << printfreq_ << std::endl;
        out << "cfl = " << cfl_ << std::endl;
        out << "delta_cfl = " << delta_cfl_ << std::endl;
        out << "cfl_increase_freq = " << cfl_increase_freq_ << std::endl;
        out << "final_time = " << finaltime_ << std::endl;
        out << "load_estim_type = " << static_cast<int>(load_estim_type_) << std::endl;
        out << "make_load_balance = " << make_load_balance_ << std::endl;
        out << "pseudo3D = " << pseudo3D_ << std::endl;
        out << "print_map = " << print_map_ << std::endl;
        out << "maxtimestep = " << maxtimestep_ << std::endl;
        out << "half_cfl = " << half_cfl_ << std::endl;
        out << "can_rebalance = " << can_rebalance_ << std::endl;
        out << "force_rebalance = " << force_rebalance_ << std::endl;
        out << "rebalance_thres = " << rebalance_thres_ << std::endl;
        out << "print_repart_info = " << print_repart_info_ << std::endl;
        out << "print_imbalance = " << print_imbalance_ << std::endl;
        out << "repart_ratio = " << repart_ratio_ << std::endl;
        out << "riemann_solver = " << static_cast<int>(riemann_solver_type_) << std::endl;
        out << "dual_ts = " << dual_ts_ << std::endl;
        out << "linear_solver_max_restart_ = " << linear_solver_max_restart_ << std::endl;
        out << "linear_solver_max_iteration_ = " << linear_solver_max_iteration_ << std::endl;
        out << "linear_solver_abs_error_ = " << linear_solver_abs_error_ << std::endl;
        out << "linear_solver_rel_error_ = " << linear_solver_rel_error_ << std::endl;
        out << "use_local_time_step_ = " << use_local_time_step_ << std::endl;

        out.close();
    }

    void Solver::update_partitioned_mesh_exchanger()
    {
        //auto &sp = partition_->spc_->sp_.front();

        //auto global_nmesh = partition_->spc_->mesh_system_size();

        //for (int i = 0; i < global_nmesh; ++i)
        {
            //Mesh* mesh_ptr = nullptr;
            //auto meshp = std::find_if(sp.mesh_.begin(), sp.mesh_.end(), [i](Mesh& m){return m.tag()() == i;});
            //auto& mesh = *meshp;

            //if (meshp != sp.mesh_.end())
            //{
                //mesh_ptr = &mesh;
            //}

            if (var_exc_ != nullptr)
            {
                //var_exc_->update(mesh_ptr, profiler_, "sol-ghost-exc");
                var_exc_->update(profiler_, "sol-ghost-exc");
            }
        }
    }

    void Solver::update_ghosts()
    {
        auto &sp = partition_->spc_->sp_.front();

        //if (sp.mesh_.size() == 1) {
            //return;
        //}
        
        if (comm_->size() == 1) {
            return;
        }

        for (Mesh& mesh: sp.mesh_)
        {
            mesh.update_ghost_primitives(var_exc_->arrival(), comm_->rank(), fs_.gamma_);
        }
    }

    void Solver::update_overset_mesh_exchanger()
    {
        //auto &sp = partition_->spc_->sp_.front();

        //for (Mesh& mesh: sp.mesh_)
        {
            if (donor_var_exc_ != nullptr)
            {
                //donor_var_exc_->update(&mesh, profiler_, "sol-donor-exc");
                donor_var_exc_->update(profiler_, "sol-donor-exc");
            }
        }
    }

    void Solver::update_donors(Mesh& mesh)
    {
        auto &sp = partition_->spc_->sp_.front();

        if (sp.mesh_.size() == 1) {
            return;
        }

        //for (Mesh& mesh: sp.mesh_)
        {
            if (donor_var_exc_ == nullptr)
            {
                oga_interpolate(mesh);
            }
            else
            {
                mesh.oga_interpolate(donor_var_exc_->arrival(), comm_->rank());
            }
        }
    }

    void Solver::set_boundary_conditions(Mesh& mesh)
    {
        bc_.set_bc(mesh, profiler_);
    }

    void Solver::compute_sum_of_fluxes(Mesh& mesh, int ntimestep)
    {
        if (ntimestep == 0 || steady_ || dual_ts_)
        {
            compute_sum_of_fluxes(mesh);
        }
        else
        {
            mesh.reset_to_mid();
        }
    }

    Vector5 Solver::non_linear_iteration()
    {
        // loop is needed
        // for several communication across overset and partitioned meshes.
        // for steady state solution.
        // for dual time stepping.

        Vector5 global_residual(-1.);
        //auto global_nmesh = partition_->spc_->mesh_system_size();

        for (int ntimestep = 0; ntimestep < maxtimestep_; ++ntimestep)
        {
            Vector5 local_residual(TAILOR_BIG_NEG_NUM);
            auto &sp = partition_->spc_->sp_.front();

            int r_max = 1;
            if (temporal_discretization_ == "runge_kutta_4") {
                r_max = 4;
            }

            for (int runge_kutta_stage = 0; runge_kutta_stage < r_max; ++runge_kutta_stage)
            {
                update_ghosts();
                //for (int i = 0; i < global_nmesh_; ++i)
                //{
                //    Mesh* mesh_ptr = nullptr;
                //    auto meshp = std::find_if(sp.mesh_.begin(), sp.mesh_.end(), [i](Mesh& m){return m.tag()() == i;});
                //    auto& mesh = *meshp;

                //    if (meshp != sp.mesh_.end())
                //    {
                //        update_donors(mesh);
                //    }
                //}

                for (int i = 0; i < global_nmesh_; ++i)
                {
                    Mesh* mesh_ptr = nullptr;
                    auto meshp = std::find_if(sp.mesh_.begin(), sp.mesh_.end(), [i](Mesh& m){return m.tag()() == i;});
                    auto& mesh = *meshp;
        for (const MeshCell &mc : partition_->spc().sp().front().mesh().front().cell())
        {
            assert(mc.oga_cell_type() != OGA_cell_type_t::undefined);
        }

                    if (meshp != sp.mesh_.end())
                    {
                        update_donors(mesh);

                        set_boundary_conditions(mesh);
                        compute_sum_of_fluxes(mesh, ntimestep);
                        temporal_discretization(mesh);

                        local_residual = compute_residual(mesh);

                        calc_change_in_conserved_var(mesh, runge_kutta_stage);

                        evolve_solution_in_time(mesh);
                        evolve_old_solution_in_time(mesh, runge_kutta_stage);

                        mesh_ptr = &mesh;
                    }
                }

                update_partitioned_mesh_exchanger();
                update_overset_mesh_exchanger();
            }

            global_residual = get_global_residual(local_residual, ntimestep);
            increase_cfl(global_residual, ntimestep);
            print_sub_solver_residual(ntimestep, global_residual);

            if (ntimestep != 0)
            {
                bool all_good = true;
                for (int i = 0; i < global_residual.nelm(); ++i)
                {
                    if (global_residual(i) > tol_)
                    {
                        all_good = false;
                    }
                }

                if (all_good) {
                    return global_residual;
                }
            }
        }

        return global_residual;
    }

    void Solver::print_sub_solver_residual(int ntimestep, const Vector5& residual)
    {
        std::ofstream of;

        std::string fn = "residual-";
        fn.append(std::to_string(nsolve_));
        fn.append(".dat");
        of.open(fn, std::ios_base::app);

        if (comm_->rank() == 0)
        {
            of << ntimestep << " ";
            for (int i = 0; i < residual.nelm(); ++i)
            {
                of << std::scientific << residual(i) << std::fixed;
                of << " ";
            }
            of << std::endl;
        }

        of.close();
    }

    Vector5 Solver::get_global_residual(const Vector5& local_residual, int ntimestep)
    {
        Vector5 global_residual;

        for (int i = 0; i < local_residual.nelm(); ++i)
        //for (int i = 0; i < 1; ++i)
        {
            //double d;
            boost::mpi::all_reduce(*comm_, local_residual(i), global_residual(i), boost::mpi::maximum<double>());
            //boost::mpi::all_reduce(*comm_, local_residual(i), d, boost::mpi::maximum<double>());
            //global_residual(i) = d;
        }

        if (ntimestep ==  0)
        {
            initial_global_residual_ = global_residual;
            last_global_residual_ = global_residual;
        }

        return global_residual;
    }

    void Solver::increase_cfl(const Vector5& global_residual, int timestep)
    {
        if (increase_cfl_)
        {
            for (int i = 0; i < global_residual.nelm(); ++i)
            {
                if (last_global_residual_(i) / global_residual(i) < cfl_ratio_)
                {
                    return;
                }
            }

            if (comm_->rank() == 0)
            {
                std::ofstream out;
                out.open("cfl.dat", std::ios_base::app);

                out << timestep << " "; 
                out << cfl_ << " "; 
                out << cfl_ * cfl_multiplier_; 
                out << std::endl;

                out.close();
            }

            cfl_ *= cfl_multiplier_;
            last_global_residual_ = global_residual;
        }
    }

    void Solver::init_old_conservative_var()
    {
        auto &sp = partition_->spc_->sp_.front();

        for (Mesh &mesh : sp.mesh_)
        {
            if (nsolve_ == 0)
            {
                for (MeshCell &mc : mesh.cell_)
                {
                    mc.cons_s_ = mc.cons_sp1_;
                    mc.cons_nm1_ = mc.cons_sp1_;
                    mc.cons_n_ = mc.cons_sp1_;
                }
            }
            else
            {
                for (MeshCell &mc : mesh.cell_)
                {
                    mc.cons_s_ = mc.cons_sp1_;
                    mc.cons_nm1_ = mc.cons_n_;
                    mc.cons_n_ = mc.cons_sp1_;
                }
            }
        }
    }

    void Solver::compute_gradient_coef()
    {
        auto &sp = partition_->spc_->sp_.front();

        for (Mesh &mesh : sp.mesh_)
        {
            if (nsolve_ == 0)
            {
                if (sorder_ == 2)
                {
                    gradient_.calc_ls_coef(mesh);
                }
            }
        }
    }

    void Solver::calc_mesh_velocities()
    {
        auto &sp = partition_->spc_->sp_.front();

        for (Mesh &mesh : sp.mesh_)
        {
            mesh.calc_mesh_velocities(fs_, comm_->rank());
        }
    }

    void Solver::print_residual(const Vector5& residual)
    {
        if (print_residual_)
        {
            if (comm_->rank() == 0)
            {
                std::ofstream out;
                out.open("residual.dat", std::ios_base::app);

                out << nsolve_ << " "; 
                for (int i = 0; i < residual.nelm(); ++i)
                {
                    out << std::scientific << residual(i) << std::fixed;
                    out << " "; 
                }
                out << std::endl;

                out.close();
            }
        }
    }

    void Solver::update_fs(const Freestream &fs)
    {
        fs_ = fs;
        bc_.update_fs(fs_);
    }

    //Matrix5 Jacobian(const Vector5& prim, double gamma, const Vector3& vf)
    Matrix5 Jacobian(const State &state, double gamma)
    {
        double rho = state.rho;
        double u = state.u;
        double v = state.v;
        double w = state.w;
        double p = state.p;
        double e = state.e;
        double k = state.k;
        double E = state.E;
        double H = state.H;
        double a = state.a;

        //double rho = prim(0);
        //double u = prim(1) - vf(0);
        //double v = prim(2) - vf(1);
        //double w = prim(3) - vf(2);
        //double p = prim(4);

        //double e =  spec_inter_energy(rho, p, gamma);
        //double k = spec_kine_energy(u, v, w);
        //double E = total_energy(rho, k, e);
        //double H = total_enthalpy(rho, p, E);
        //double a = speed_of_sound(rho, p, gamma);

        double gs = gamma - 1.;

        Matrix5 A;

        A(0, 0) = 0.;
        A(0, 1) = 1.;
        A(0, 2) = 0.;
        A(0, 3) = 0.;
        A(0, 4) = 0.;

        A(1, 0) = gs * H - std::pow(u, 2.) - std::pow(a, 2.);
        A(1, 1) = (3. - gamma) * u;
        A(1, 2) = -gs * v;
        A(1, 3) = -gs * w;
        A(1, 4) = gs;

        A(2, 0) = -u * v;
        A(2, 1) = v;
        A(2, 2) = u;
        A(2, 3) = 0.;
        A(2, 4) = 0.;

        A(3, 0) = -u * w;
        A(3, 1) = w;
        A(3, 2) = 0.;
        A(3, 3) = u;
        A(3, 4) = 0.;

        A(4, 0) = 0.5 * u * ((gamma - 3.) * H - std::pow(a, 2.));
        A(4, 1) = H - gs * std::pow(u, 2.);
        A(4, 2) = -gs * u * v;
        A(4, 3) = -gs * u * w;
        A(4, 4) = gamma * u;

        return A;
    }

    void Solver::update_matrices(MeshFace *this_face, MeshFace *common_face, MeshCell& left_cell, MeshCell& right_cell, double facearea, const Vector3& face_velocity, double gamma, const Matrix5& rotation_matrix, const Matrix5& inv_rotation_matrix, const Matrix5& Aroe, const State& left_state_rotated, const State& right_state_rotated)
    {
        //auto [left_state, right_state] = left_and_right_states(left_cell.cons_sp1(), right_cell.cons_sp1(), gamma, unit_matrix<NVAR, NVAR, double>(), face_velocity);
        //auto [left_state, right_state] = left_and_right_states(left_cell.cons_sp1(), right_cell.cons_sp1(), gamma, rotation_matrix, face_velocity);

        Matrix5 JL = Jacobian(left_state_rotated , gamma);
        Matrix5 JR = Jacobian(right_state_rotated, gamma);

        this_face->M_ = inv_rotation_matrix * (JL + Aroe) * rotation_matrix * 0.5 * facearea;
        //this_face->M_ = inv_rotation_matrix * JL * 0.5 * facearea;
        //this_face->M_ = inv_rotation_matrix * (JL + Aroe) * 0.5 * facearea;
        //this_face->M_ = (JL + Aroe) * rotation_matrix * 0.5 * facearea;
        //this_face->M_ = (JL + Aroe) * 0.5 * facearea;
        if (common_face != nullptr)
        {
            //common_face->M_ = inv_rotation_matrix * (JR - Aroe) * 0.5 * facearea * -1;
            //common_face->M_ = (JR - Aroe) * rotation_matrix * 0.5 * facearea * -1;
            //common_face->M_ = (JR - Aroe) * 0.5 * facearea * -1;
            //common_face->M_ = inv_rotation_matrix * JR * 0.5 * facearea * -1;
            common_face->M_ = inv_rotation_matrix * (JR - Aroe) * rotation_matrix * 0.5 * facearea * -1;
        }

        left_cell.D_ += this_face->M_;
        left_cell.D_mid_ += this_face->M_;

        if (common_face != nullptr)
        {
            right_cell.D_ += common_face->M_;
            right_cell.D_mid_ += common_face->M_;
        }
    }

    Matrix5 Solver::slipwall_M(const Vector3 &n)
    {
        double nx = n(0);
        double ny = n(1);
        double nz = n(2);

        Matrix5 M;

        M(0, 0) = 1.;
        M(0, 1) = 0.;
        M(0, 2) = 0.;
        M(0, 3) = 0.;
        M(0, 4) = 0.;

        M(1, 0) = 0.;
        M(1, 1) = 1. - 2. * std::pow(nx, 2);
        M(1, 2) = -2. * nx * ny;
        M(1, 3) = -2. * nx * nz;
        M(1, 4) = 0.;

        M(2, 0) = 0.;
        M(2, 1) = -2. * nx * ny;
        M(2, 2) = 1. - 2. * std::pow(ny, 2);
        M(2, 3) = -2. * ny * nz;
        M(2, 4) = 0.;

        M(3, 0) = 0.;
        M(3, 1) = -2. * nx * nz;
        M(3, 2) = -2. * ny * nz;
        M(3, 3) = 1. - 2. * std::pow(nz, 2);
        M(3, 4) = 0.;

        M(4, 0) = 0.;
        M(4, 1) = 0.;
        M(4, 2) = 0.;
        M(4, 3) = 0.;
        M(4, 4) = 1.;

        return M;
    }

    Vector5 Solver::apply_limiter(const Mesh &mesh, const MeshCell &mc, const MeshFace &mf)
    {
        auto prim = mc.prim();
        auto cons = prim_to_cons(prim, fs_.gamma_);

        if (mc.oga_cell_type() == OGA_cell_type_t::non_resident || mc.oga_cell_type() == OGA_cell_type_t::ghost)
        {
            //prim = mc.prim();
            return cons;
        }

        assert(mc.oga_cell_type() != OGA_cell_type_t::undefined);

        //for (const auto &mff : mc.face())
        //{
        //    if (mff.btype() != BouType::interior) // TODO
        //    {
        //        //prim = mc.prim();
        //        return cons;
        //    }
        //    else
        //    {
        //        if (mff.parent_cell().size() != 2)
        //        {
        //            std::cout << "mc oga: " << static_cast<int>(mc.oga_cell_type()) << std::endl;
        //            std::cout << "mc btype: " << static_cast<int>(mc.btype()) << std::endl;
        //        }
        //        assert(mff.parent_cell().size() == 2);
        //    }
        //}

        for (const auto& f: mc.pnei())
        {
            const auto& nei = mesh.cell(f);
            if (mc.oga_cell_type() == OGA_cell_type_t::field)
            {
                assert(nei.oga_cell_type() != OGA_cell_type_t::hole);
            }
        }

        auto grad = gradient_.ls_grad(mesh, mc);

        Limiter limiter(LimiterType::venka);
        limiter.limit(mesh, mc, grad);

        auto d = mf.face().centroid() - mc.poly().centroid();

        assert(!mc.prim().isnan());

        for (int i = 0; i < NVAR; ++i)
        {
            assert(!std::isnan(dot(grad[i], d)));
            prim(i) = mc.prim(i) + limiter(i) * dot(grad[i], d);
            std::cout << "grad: " << limiter(i) << " " << dot(grad[i], d) << std::endl;
        }

        assert(!mc.prim().isnan());

        cons = prim_to_cons(prim, fs_.gamma_);

        return cons;
    }

    void Solver::sweep(Mesh &mesh, MeshCell &mc, Vector5 &r1, Vector5 &r2, Vector5 &r3, int &maxcell)
    {
        mc.dQ_ = (1. - omega_) * mc.old_dQ_ + omega_ * mc.R() / mc.D();

        //if (mc.tag()() == 1)
        //std::cout << "dq, old_dq, diff: " << mc.dQ_[0] << " " << mc.old_dQ_[0] << " " << std::abs(mc.dQ_[0] - mc.old_dQ_[0]) << std::endl;

        //if (mc.tag()() == 1369)
        //if (calc_r)
        {
            //if (std::abs(mc.dQ_[0] - mc.old_dQ_[0]) > r1[0])
            //maxcell = mc.tag()();
            r1(0) = std::max(r1(0), std::abs(mc.dQ_(0) - mc.old_dQ_(0)));
            r1(1) = std::max(r1(1), std::abs(mc.dQ_(1) - mc.old_dQ_(1)));
            r1(2) = std::max(r1(2), std::abs(mc.dQ_(2) - mc.old_dQ_(2)));
            r1(3) = std::max(r1(3), std::abs(mc.dQ_(3) - mc.old_dQ_(3)));
            r1(4) = std::max(r1(4), std::abs(mc.dQ_(4) - mc.old_dQ_(4)));
            r2(0) = std::pow(mc.dQ_(0) - mc.old_dQ_(0), 2.);
            r2(1) = std::pow(mc.dQ_(1) - mc.old_dQ_(1), 2.);
            r2(2) = std::pow(mc.dQ_(2) - mc.old_dQ_(2), 2.);
            r2(3) = std::pow(mc.dQ_(3) - mc.old_dQ_(3), 2.);
            r2(4) = std::pow(mc.dQ_(4) - mc.old_dQ_(4), 2.);
            //r3[0] = std::max(r3[0], std::abs(mc.dQ_[0]));
            //r3[1] = std::max(r3[1], std::abs(mc.dQ_[1]));
            //r3[2] = std::max(r3[2], std::abs(mc.dQ_[2]));
            //r3[3] = std::max(r3[3], std::abs(mc.dQ_[3]));
            //r3[4] = std::max(r3[4], std::abs(mc.dQ_[4]));
            r3(0) = std::max(r3(0), std::abs(mc.R(0)));
            r3(1) = std::max(r3(1), std::abs(mc.R(1)));
            r3(2) = std::max(r3(2), std::abs(mc.R(2)));
            r3(3) = std::max(r3(3), std::abs(mc.R(3)));
            r3(4) = std::max(r3(4), std::abs(mc.R(4)));
        }

        /*{
                if (std::abs(r1[3]) >= TAILOR_ZERO)
                {
                    std::cout << " r1[3]: " << r1[3] << std::endl;
                    std::cout << " mc.dQ_[3]: " << mc.dQ_[3] << std::endl;
                    std::cout << " mc.old_dQ[3]: " << mc.old_dQ_[3] << std::endl;
                    std::cout << " mc.R_[0]: " << mc.R_[0] << std::endl;
                    std::cout << " mc.R_[1]: " << mc.R_[1] << std::endl;
                    std::cout << " mc.R_[2]: " << mc.R_[2] << std::endl;
                    std::cout << " mc.R_[3]: " << mc.R_[3] << std::endl;
                    std::cout << " mc.R_[4]: " << mc.R_[4] << std::endl;
                }
                assert(std::abs(r1[3]) < TAILOR_ZERO);
            }*/

        auto ddQ = mc.dQ() - mc.old_dQ();
        mc.old_dQ_ = mc.dQ();

        //if (mc.tag()() == 1369)
        //std::cout << "new old dq: " << mc.old_dQ_[0] << std::endl;

        //for (const MeshFace& cmf: mc.face())
        for (MeshFace &mf : mc.face_)
        {
            //if (cmf.is_boundary())
            if (mf.is_boundary())
            {
                continue;
            }

            //MeshFace& mf = mesh.face_p(cmf.tag());

            //MeshCell *LC = nullptr;
            //MeshCell *RC = nullptr;

            auto [LC, RC] = left_and_right_cells(mesh, mf, mc.tag());
            //LC = &mc;
            //RC = &(mesh.cell_p(mf.left_cell()));

            //if (mc.tag() == LC->tag())
            //{
            //RC->R_ += mf.ML() * ddQ;
            //if (mc.tag()() == 1)
            //std::cout << "RC->R: " << RC->R_[0] << std::endl;
            //}
            //else if (mc.tag() == RC->tag())
            //{
            //LC->R_ -= mf.MR() * ddQ;
            //if (mc.tag()() == 1)
            //std::cout << "LC->R: " << LC->R_[0] << std::endl;
            //}
            //else
            //{
            //assert(false);
            //}
            assert(LC->tag() == mc.tag());
            RC->R_ += mf.M() * ddQ;
        }
        //if (mc.tag()() == 1)
        //assert(false);
    }

    void make_mat_st(const Mesh &mesh, int &n, int &nz_num, std::vector<int> &ia, std::vector<int> &ja, std::vector<double> &a, std::vector<double> &rhs, int rank)
    {
        std::map<int, int> cell_index;

        n = 0;
        for (const auto &mc : mesh.cell())
        {
            if (calc_cell(mc))
            {
                cell_index.insert(std::make_pair(mc.tag()(), n));
                ++n;
            }
        }

        n *= NVAR;
        rhs.resize(n);

        for (const auto &mc : mesh.cell())
        {
            if (calc_cell(mc))
            {
                auto block_row_index = cell_index[mc.tag()()];

                std::map<int, Matrix5> order;
                order.insert(std::make_pair(block_row_index, mc.D()));

                for (int j = 0; j < NVAR; ++j)
                {
                    rhs[block_row_index * NVAR + j] = mc.R(j);
                }

                for (const auto &mf : mc.face())
                {
                    if (mf.btype() == BouType::interior || mf.btype() == BouType::partition)
                    {
                        auto [left_cell, right_cell] = left_and_right_cells(mesh, mf, mc.tag());
                        assert(left_cell != nullptr);
                        assert(right_cell != nullptr);

                        auto common_face = mesh.common_face(*right_cell, mf.tag());
                        assert(common_face != nullptr);

                        if (calc_cell(*right_cell))
                        {
                            assert(mf.btype() == BouType::interior);
                            auto block_column_index = cell_index[right_cell->tag()()];
                            order.insert(std::make_pair(block_column_index, common_face->M() * -1));
                        }
                        else
                        {
                            assert(mf.btype() == BouType::interior || mf.btype() == BouType::partition);
                            Vector5 temp = common_face->M() * right_cell->dQ();

                            for (int j = 0; j < NVAR; ++j)
                            {
                                rhs[block_row_index * NVAR + j] += temp(j);
                            }
                        }
                    }
                }

                //std::sort(order.begin(), order.end(), [&](const auto& a, const auto& b){a.first < b.first;}); // map is instrisicly already ordered.

                for (int i = 0; i < NVAR; ++i)
                {
                    auto row_index = block_row_index * NVAR + i;
                    for (const auto &entry : order)
                    {
                        int block_column_index = entry.first;

                        for (int j = 0; j < NVAR; ++j)
                        {
                            auto column_index = block_column_index * NVAR + j;
                            double val = entry.second(i, j);

                            //if (std::abs(val) > TAILOR_ZERO)
                            {
                                ia.push_back(row_index);
                                ja.push_back(column_index);
                                a.push_back(val);
                            }
                        }
                    }
                }
            }
        }

        nz_num = a.size();
    }

    void make_mat_cr(const Mesh &mesh, int &n, int &nz_num, std::vector<int> &ia, std::vector<int> &ja, std::vector<double> &a, std::vector<double> &rhs, int rank)
    {
        make_mat_st(mesh, n, nz_num, ia, ja, a, rhs, rank);
        //if (rank == 0)
        //{
            //for (int i = 0; i < nz_num; ++i)
            //{
                //std::cout << ia[i] << " " << ja[i] << " " << a[i] << " " << rhs[i] << std::endl;
            //}
            //assert(false);
        //}

        std::vector<int> rp;

        for (int i = 0; i < nz_num; ++i)
        {
            if (i == 0)
            {
                rp.push_back(i);
            }
            else
            {
                if (ia[i] > ia[i - 1])
                {
                    rp.push_back(i);
                }
            }
        }

        rp.push_back(nz_num);

        if (rp.size() != (n + 1))
        {
            std::cout << "n: " << n << std::endl;
            std::cout << "rp size: " << rp.size() << std::endl;
            std::cout << "nz_num: " << nz_num << std::endl;
            std::cout << "ia.size(): " << ia.size() << std::endl;
        }
        assert(rp.size() == (n + 1));

        ia = rp;
        assert(ia.size() == rp.size());
        //if (rank == 0)
        //{
            //for (int i = 0; i < nz_num; ++i)
            //{
                //std::cout << rp[i] << " " << ja[i] << " " << a[i] << " " << rhs[i] << std::endl;
            //}
            //assert(false);
        //}
    }

    std::vector<double> Solver::amgcl(int n, int nz_num, const std::vector<int>& ia, const std::vector<int>& ja, const std::vector<double>& a, const std::vector<double>& rhs)
    {
        std::vector<double> x(n, 0.);

        //boost::property_tree::ptree prm;

        //prm.put("solver.type", "bicgstab");
        //prm.put("solver.type", "gmres");
        //prm.put("solver.M", 100);
        //if (n < 100)
        //{
        //    prm.put("solver.M", n-1);
        //}
        //prm.put("solver.tol", TAILOR_BIG_POS_NUM);
        //prm.put("solver.abstol", TAILOR_ZERO);
        //prm.put("solver.maxiter", 10000);
        //prm.put("precond.type", "dummy");
        //prm.put("precond.coarsening.type", "smoothed_aggregation");
        //prm.put("precond.relax.type", "spai0");
        //prm.put("precond.relax.type", "ilu0");

        typedef amgcl::make_solver<
            //amgcl::relaxation::as_preconditioner<Backend, amgcl::relaxation::ilu0>,
            //amgcl::relaxation::as_preconditioner<Backend, amgcl::relaxation::chebyshev>,
            //amgcl::relaxation::as_preconditioner<Backend, amgcl::relaxation::gauss_seidel>,
            //amgcl::relaxation::as_preconditioner<Backend, amgcl::relaxation::spai0>,
            //amgcl::preconditioner::dummy<Backend>,
            amgcl::amg<
            Backend,
            amgcl::coarsening::smoothed_aggregation,
            //amgcl::relaxation::chebyshev
            amgcl::relaxation::ilu0
            //amgcl::relaxation::spai0
                >,
            amgcl::solver::gmres<Backend>
            //amgcl::solver::preonly<Backend>
            //amgcl::solver::bicgstab<Backend>
                > AMGCLSolver;

        AMGCLSolver::params prm;
        if (linear_solver_max_restart_ != 0)
        {
            prm.solver.M = linear_solver_max_restart_;
        }
        if (linear_solver_max_iteration_ != 0)
        {
            prm.solver.maxiter = linear_solver_max_iteration_;
        }
        if (linear_solver_abs_error_ != 0)
        {
            prm.solver.abstol = linear_solver_abs_error_;
        }
        if (linear_solver_rel_error_ != 0)
        {
            prm.solver.tol = linear_solver_rel_error_;
        }

        AMGCLSolver amgcl_solver(std::tie(n, ia, ja, a), prm);
        //amgcl::make_solver<
        //    amgcl::relaxation::as_preconditioner<Backend, amgcl::relaxation::ilu0>,
        //    amgcl::solver::gmres<Backend>
        //        > amgcl_solver(
        //        amgcl::adapter::zero_copy(
        //            n,
        //            ia, 
        //            ja, 
        //            a
        //            ), 
        //        prm);
        auto [iters, error] = amgcl_solver(rhs, x);

        if (print_linear_solver_error_)
        {
            std::ofstream out;
            std::string s = "linear-error-";
            s.append(std::to_string(comm_->rank()));
            s.append(".dat");
            out.open(s, std::ios_base::app);
            out << iters << " " << error << std::endl;
            out.close();
        }

        //assert(error < 1e-6);

        return x;
    }

    void Solver::linear_solver(Mesh &mesh)
    {
        int n, nz_num;
        std::vector<int> ia, ja;
        std::vector<double> a;
        std::vector<double> rhs;

        make_mat_cr(mesh, n, nz_num, ia, ja, a, rhs, comm_->rank());

        if (n == 0)
        {
            return;
        }

        auto x = amgcl(n, nz_num, ia, ja, a, rhs);
        //auto x = gmres(n, nz_num, ia, ja, a, rhs);

        int i = 0;
        for (auto &mc : mesh.cell_)
        {
            if (calc_cell(mc))
            {
                for (int j = 0; j < NVAR; ++j)
                {
                    mc.dQ_(j) = x[i * NVAR + j];
                    //if (std::abs(mc.dQ_(j) < 1e-6))
                    //{
                        //std::cout << "R: " << mc.R(j) << std::endl;
                        //std::cout << "rhs: " << rhs[i * NVAR + j] << std::endl;
                        //assert(false);
                    //}
                }
                ++i;
            }
        }
    }

    std::vector<double> Solver::gmres(int n, int nz_num, std::vector<int>& ia, std::vector<int>& ja, std::vector<double>& a, std::vector<double>& rhs)
    {
        std::vector<double> x(n, 0.);

        int itr_max = 100;
        int mr = 100;

        if (n < mr)
        {
            mr = n - 1;
        }
        double tol_abs = TAILOR_ZERO;
        double tol_rel = TAILOR_BIG_POS_NUM;

        bool verbose = true;

        //mgmres_st(n, nz_num, ia.data(), ja.data(), a.data(), x.data(), rhs.data(), itr_max, mr, tol_abs, tol_rel, verbose);
        pmgmres_ilu_cr(n, nz_num, ia.data(), ja.data(), a.data(), x.data(), rhs.data(), itr_max, mr, tol_abs, tol_rel, verbose);

        return x;
    }

    bool Solver::sor(Mesh &mesh, int ntimestep)
    {
        for (auto mc = mesh.cell_.begin(); mc != mesh.cell_.end(); ++mc)
        {
            if (calc_cell(*mc))
            //if (mc->oga_cell_type() == OGA_cell_type_t::field)
            {
                mc->old_dQ_ = 0.;
            }
        }

        int maxcell;
        Vector5 r1;
        Vector5 r3;
        for (int i = 0; i < nsweep_; ++i)
        {
            Vector5 r2;
            r2 = 0.;
            r1 = TAILOR_BIG_NEG_NUM;
            r3 = TAILOR_BIG_NEG_NUM;
            for (auto mc = mesh.cell_.begin(); mc != std::prev(mesh.cell_.end()); ++mc)
            {
                if (calc_cell(*mc))
                //if (mc->oga_cell_type() == OGA_cell_type_t::field)
                {
                    sweep(mesh, *mc, r1, r2, r3, maxcell);
                }
            }

            r1 = TAILOR_BIG_NEG_NUM;

            for (auto mc = mesh.cell_.rbegin(); mc != mesh.cell_.rend(); ++mc)
            {
                if (calc_cell(*mc))
                //if (mc->oga_cell_type() == OGA_cell_type_t::field)
                {
                    sweep(mesh, *mc, r1, r2, r3, maxcell);
                }
            }

            if (r1(0) == TAILOR_BIG_NEG_NUM)
            {
                // no field in mesh.
                return true;
            }
            if (show_inner_norm_)
            {
                std::string fn = "inner-norm-max-";
                fn.append(std::to_string(comm_->rank()));
                fn.append(".dat");
                std::ofstream outmax;
                outmax.open(fn, std::ios_base::app);

                //std::string fn = "inner-norm-";
                //fn.append(std::to_string(comm_->rank()));
                //fn.append(".dat");
                //std::ofstream out;
                //out.open(fn, std::ios_base::app);

                //std::cout << "maxcell: " << maxcell << std::scientific;
                outmax << "ntimestep: " << ntimestep << " ";
                outmax << "cfl: " << cfl_ << " ";
                outmax << "nsweep: " << i << " ";
                outmax << " | ";
                outmax << std::scientific;
                //std::cout << std::sqrt(r1[0]);
                //std::cout << " ";
                //std::cout << std::sqrt(r1[1]);
                //std::cout << " ";
                //std::cout << std::sqrt(r1[2]);
                //std::cout << " ";
                //std::cout << std::sqrt(r1[3]);
                //std::cout << " ";
                //std::cout << std::sqrt(r1[4]);
                outmax << (r1(0));
                outmax << " ";
                outmax << (r1(1));
                outmax << " ";
                outmax << (r1(2));
                outmax << " ";
                outmax << (r1(3));
                outmax << " ";
                outmax << (r1(4));
                outmax << std::endl;
                //out << "r3: " << std::scientific;
                //std::cout << std::sqrt(r1[0]);
                //std::cout << " ";
                //std::cout << std::sqrt(r1[1]);
                //std::cout << " ";
                //std::cout << std::sqrt(r1[2]);
                //std::cout << " ";
                //std::cout << std::sqrt(r1[3]);
                //std::cout << " ";
                //std::cout << std::sqrt(r1[4]);
                //out << (r3[0]);
                //out << " ";
                //out << (r3[1]);
                //out << " ";
                //out << (r3[2]);
                //out << " ";
                //out << (r3[3]);
                //out << " ";
                //out << (r3[4]);
                //out << std::endl;
                //out.close();
                //std::cout << "r2: " << std::scientific;
                //std::cout << std::sqrt(r2[0]);
                //std::cout << " ";
                //std::cout << std::sqrt(r2[1]);
                //std::cout << " ";
                //std::cout << std::sqrt(r2[2]);
                //std::cout << " ";
                //std::cout << std::sqrt(r2[3]);
                //std::cout << " ";
                //std::cout << std::sqrt(r2[4]);
                //std::cout << std::endl;;
                outmax.close();
            }

            if (std::abs(r1(0)) < TAILOR_ZERO && std::abs(r1(1)) < TAILOR_ZERO && std::abs(r1(2)) < TAILOR_ZERO && std::abs(r1(3)) < TAILOR_ZERO && std::abs(r1(4)) < TAILOR_ZERO)
            {
                return true;
            }
        }

        return false;
    }

    double calc_local_time_step(double sumarea, double volume, double cfl)
    {
        double charlen = volume / sumarea;
        double dtau = cfl * charlen;

        return dtau;
    }

    void Solver::temporal_discretization(Mesh& mesh)
    {
        for (MeshCell &mc : mesh.cell_)
        {
            if (!calc_cell(mc))
            {
                continue;
            }

            double volume = mc.poly().volume();
            mc.dtao_ = calc_local_time_step(mc.sumarea_, volume, cfl_);

            if (torder_ == 1)
            {
                //if (dual_ts_ || steady_)
                if (dual_ts_)
                {
                    assert(false);
                    double t = volume / mc.dtao_;
                    mc.D_.add_diag(t);

                    mc.R_ -= volume * (mc.cons_sp1() - mc.cons_n()) / dt_;
                }
                else if (use_local_time_step_)
                {
                    double t = volume / mc.dtao_;
                    mc.D_.add_diag(t);
                }
                else
                {
                    double t = volume / dt_;
                    mc.D_.add_diag(t);
                }
            }
            else if (torder_ == 2)
            {
                //double t = volume * (1. / mc.dtao_ + 1.5 / dt_);
                //mc.D_.add_diag(t);

                if (dual_ts_)
                {
                    assert(false);
                    //double t = volume / mc.dtao_;
                    double t = volume / mc.dtao_;
                    mc.D_.add_diag(t);

                    if (nsolve_ > 1)
                    {
                        mc.R_ -= 0.5 * volume * (3. * mc.cons_sp1() - 4. * mc.cons_n() + mc.cons_nm1()) / dt_;
                    }
                    else
                    {
                        mc.R_ -= volume * (mc.cons_sp1() - mc.cons_n()) / dt_;
                    }
                }
                else if (use_local_time_step_)
                {
                    double t = 1.5 * volume / mc.dtao_;
                    mc.D_.add_diag(t);
                    mc.R_ += 0.5 * volume * (mc.cons_n() - mc.cons_nm1()) / mc.dtao_;
                }
                else
                {
                    double t = 1.5 * volume / dt_;
                    mc.D_.add_diag(t);

                    if (nsolve_ > 1)
                    {
                        mc.R_ += 0.5 * volume * (mc.cons_n() - mc.cons_nm1()) / dt_;
                    }
                }
            }
        }
    }

    Partition *Solver::partition()
    {
        return partition_;
    }
    const Partition *Solver::partition() const
    {
        return partition_;
    }
}
