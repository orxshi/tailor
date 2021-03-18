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

    Solver::Solver() : repart_ratio_(0),
                       initratio_(0.),
                       print_imbalance_(false),
                       print_vtk_(false),
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
                       print_outer_norm_(false),
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
                       init_max_res_(0.),
                       last_max_res_(0.),
                       donor_var_exc_(nullptr),
                       riemann_solver_type_(RiemannSolverType::roe),
                       dual_ts_(false)
    {
    }

    Solver::Solver(boost::mpi::communicator *comm, const std::vector<std::string> &filename, Profiler *profiler, Partition *partition) : verbose_(true), maxtimestep_(10000), comm_(comm), var_exc_(nullptr), donor_var_exc_(nullptr), nsolve_(0), profiler_(profiler), partition_(partition), 
    init_max_res_(0.),
    last_max_res_(0.)
    {
        if (steady_)
        {
            dt_ = TAILOR_BIG_POS_NUM;
            maxtimestep_ = int(1e6);
        }

        read_settings();

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

        if (comm_->size() != 1)
        {
            connect_partition_cells();
        }

        if (profiler_ != nullptr)
        {
            profiler_->start("sol-init");
        }
        partition_->spc_->init(); // TODO
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
            ("solver.repart-ratio", po::value<int>()->default_value(100), "")
            ("solver.rebalance-thres", po::value<double>()->default_value(40.), "")
            ("solver.show_inner_res", po::value<bool>()->default_value(true), "Show inner loop residual")
            ("solver.show_inner_norm", po::value<bool>()->default_value(true), "Show inner loop norm")
            ("solver.print_outer_norm", po::value<bool>()->default_value(true), "Print outer norm")
            ("solver.steady", po::value<bool>()->default_value(false), "Steady state")
            ("solver.progressive_cfl", po::value<bool>()->default_value(false), "Progressive CFL increase")
            ("solver.temporal-discretization", po::value<std::string>()->default_value("forward_euler"), "Temporal discretization")
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
            ("solver.print-vtk", po::value<bool>()->default_value(false), "")
            ("solver.print-repart-info", po::value<bool>()->default_value(false), "")
            ("solver.print-imbalance", po::value<bool>()->default_value(false), "")
            ("solver.riemann-solver", po::value<int>()->default_value(0), "")
            ("solver.dual-ts", po::value<bool>()->default_value(false), "");

        all_options.add(desc);

        std::ifstream settings_file("settings.ini");

        boost::program_options::variables_map vm;
        po::store(po::parse_config_file(settings_file, all_options, true), vm);
        po::notify(vm);

        print_outer_norm_ = vm["solver.print_outer_norm"].as<bool>();
        show_inner_res_ = vm["solver.show_inner_res"].as<bool>();
        show_inner_norm_ = vm["solver.show_inner_norm"].as<bool>();
        progressive_cfl_ = vm["solver.progressive_cfl"].as<bool>();
        steady_ = vm["solver.steady"].as<bool>();
        temporal_discretization_ = vm["solver.temporal-discretization"].as<std::string>();
        dt_ = vm["solver.dt"].as<double>();
        nsweep_ = vm["solver.nsweep"].as<int>();
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
        print_vtk_ = vm["solver.print-vtk"].as<bool>();
        print_repart_info_ = vm["solver.print-repart-info"].as<bool>();
        print_imbalance_ = vm["solver.print-imbalance"].as<bool>();
        repart_ratio_ = vm["solver.repart-ratio"].as<int>();
        if (vm["solver.riemann-solver"].as<int>() == 0) {
            riemann_solver_type_ = RiemannSolverType::roe;
        }
        else if (vm["solver.riemann-solver"].as<int>() == 1) {
            riemann_solver_type_ = RiemannSolverType::hllc;
        }
        else {
            assert(false);
        }
        dual_ts_ = vm["solver.dual-ts"].as<bool>();
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

    std::tuple<Vector5, double> compute_flux(RiemannSolverType riemann_solver_type, double face_area, const Matrix5& inverse_rotation_matrix, const State& rotated_left_state, const State& rotated_right_state, double gamma)
    {
        Vector5 flux;
        double max_eigen;

        RiemannSolver riemann_solver(riemann_solver_type, rotated_left_state, rotated_right_state, face_area, gamma, max_eigen, flux);

        flux = inverse_rotation_matrix * flux;

        assert(!flux.isnan());
        assert(max_eigen > 0.);

        return std::make_tuple(flux, max_eigen);
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

                if (sorder_ == 2)
                {
                    apply_limiter(mesh, *left_cell , mf);
                    apply_limiter(mesh, *right_cell, mf);
                }

                MeshFace *commonface = nullptr;
                if (!mf.is_boundary())
                {
                    commonface = mesh.common_face(*right_cell, mf.tag());
                    assert(commonface != nullptr);
                }

                auto [rotation_matrix, inv_rotation_matrix] = get_rotation_matrix(normal);
                auto [rotated_left_state, rotated_right_state] = left_and_right_states(left_cell->cons_sp1_, right_cell->cons_sp1_, gamma, rotation_matrix, face_velocity);
                auto [flux, max_eigen] = compute_flux(riemann_solver_type_, face_area, inv_rotation_matrix, rotated_left_state, rotated_right_state, gamma);

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
                    update_matrices(&mf, commonface, *left_cell, *right_cell, face_area, face_velocity, gamma);
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
            mc.cons_sp1_ = mc.cons_s_ + mc.dQ_;
            mc.prim_ = cons_to_prim(mc.cons_sp1_, fs_.gamma_);
        }   
    }

    void Solver::evolve_old_solution_in_time(Mesh& mesh)
    {
        for (MeshCell &mc : mesh.cell_)
        {
            mc.cons_s_ = mc.cons_sp1_;
            if (steady_)
            {
                mc.cons_nm1_ = mc.cons_n_;
                mc.cons_n_ = mc.cons_sp1_;
            }
        }
    }

    void Solver::calc_change_in_conserved_var(Mesh &mesh)
    {
        if (temporal_discretization_ == "backward_euler")
        {
            gmres(mesh);
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

                if (dual_ts_ || steady_) {
                    mc.dQ_ *= mc.dtao_;
                }
                else
                {
                    if (torder_ == 1)
                    {
                        mc.dQ_ *= dt_;
                    }
                    else if (torder_ == 2)
                    {
                        mc.dQ_ *= (2. / 3.) * dt_;
                    }
                }
            }
        }
    }

    //void Solver::update_cons_implicitly(Mesh &mesh, int ntimestep)
    //{
    //    //tempo_discre(mesh);
    //    gmres(mesh);
    //    //sor(mesh, ntimestep);
    //    /*bool converged = false;
    //    double orig_cfl = cfl_;
    //    while(!converged)
    //    {
    //        //std::string fn = "inner-norm-";
    //        //fn.append(std::to_string(comm_->rank()));
    //        //fn.append(".dat");
    //        //std::ofstream out;
    //        //out.open(fn, std::ios_base::app);
    //        //out << "cfl = " << cfl_ << std::endl;
    //        //out.close();
    //        tempo_discre(mesh);
    //        converged = sor(mesh, ntimestep);
    //        if (!half_cfl_) {
    //            break;
    //        }
    //        if (!converged)
    //        {
    //            // actually whole calc_R isn't needed. just need to revert R and D after calc_R, before update_matrices state.
    //            //calc_R(mesh);
    //            mesh.equate_RD_to_RDmid();
    //            cfl_ = cfl_ / 2;
    //        }
    //    }
    //    cfl_ = orig_cfl;*/

    //    for (MeshCell &mc : mesh.cell_)
    //    {
    //        if (!calc_cell(mc))
    //        {
    //            continue;
    //        }

    //        update_solution(mc, gamma);
    //    }
    //}

    void Solver::solve()
    {
        partition_->print_cell_dist(comm_, nsolve_);

        if (nsolve_ == 0)
        {
            fs_.read();
        }

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

        //if (nsolve_ == 0)
        //{
        //    if (restore_solution_)
        //    {
        //        restore_solution();
        //    }
        //}

        solve_(comm_->rank());
        if (!steady_) {
            ++nsolve_;
        }

        //if (save_solution_)
        //{
            //save_solution();
        //}

        if (print_vtk_)
        {
            for (const auto &sp : partition_->spc().sp())
            {
                for (const auto &m : sp.mesh())
                {
                    std::string fn = "sol-";
                    fn.append(std::to_string(m.tag()()));
                    fn.append("-");
                    fn.append(std::to_string(comm_->rank()));
                    fn.append("-");
                    fn.append(std::to_string(sp.tag()()));
                    fn.append(".vtk");
                    m.print_as_vtk(fn);
                }
            }
            for (const auto &sp : partition_->spc().sp())
            {
                for (const auto &m : sp.mesh())
                {
                    std::string fn = "solwall-";
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

        if (donor_var_exc_ != nullptr)
        {
            delete donor_var_exc_;
            donor_var_exc_ = nullptr;
        }

        if (var_exc_ != nullptr)
        {
            delete var_exc_;
            var_exc_ = nullptr;
        }
    }

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

    Vector5 Solver::resi(Mesh &mesh)
    {
        // TODO refactor as L1 norm, abs, etc.
        Vector5 res;
        res = 0.;
        //res = TAILOR_BIG_NEG_NUM;
        for (MeshCell &mc : mesh.cell_)
        {
            if (!calc_cell(mc))
            {
                continue;
            }

            res += abs(mc.R_);

            //for (int i = 0; i < NVAR; ++i)
            //{
                //res[i] = std::max(res[i], std::abs(mc.cons_sp1_[i] - mc.cons_s_[i]) * mc.poly().volume() / mc.dtao_);
                //res[i] = std::max(res[i], std::abs(mc.cons_sp1_[i] - mc.cons_s_[i]) * (1./mc.dtao_ + 1./dt_));
                //res[i] = std::max(res[i], std::abs(-mc.R_[i] - mc.poly().volume() * (mc.cons_s_[i] - mc.cons_n_[i]) / dt_));
                //assert(std::abs(mc.R_[i]) < 1e3);
                //res[i] += std::abs(mc.poly().volume() * (mc.cons_sp1_[i] - mc.cons_s_[i]) / dt_ + mc.R_[i]);
                //res[i] += std::abs(mc.cons_sp1(i) - mc.cons_s(i));
                //std::cout << "resi: " << mc.tag()() << " " << mc.R_[i] << std::endl;
                //assert(std::abs(res[i]) < 1e3);
                //res[i] += std::pow((mc.cons_sp1_(i) - mc.cons_s_(i)), 2);
                //res[i] += std::pow((mc.R_[i]), 2.);
            //}
        }
        //for (int i=0; i<NVAR; ++i)
        //{
        //res[i] = std::sqrt(res[i]);
        //}
        //averes_ = res[0];
        //double averes = TAILOR_BIG_NEG_NUM;
        res /= mesh.cell().size();
        //for (int i = 0; i < NVAR; ++i)
        //{
            //averes = std::max(averes, res[i]);
            //res[i] = res[i] / mesh.cell().size();
            //res[i] = std::sqrt(res[i]);
        //}

        return res;

        //return averes;
    }

    void Solver::oga_interpolate(Mesh &mesh)
    {
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
                    auto grad = gradient_.ls_grad(*m, donor_cell);

                    Limiter limiter(LimiterType::venka);
                    limiter.limit(mesh, mc, grad);

                    auto d = mc.poly().centroid() - donor_cell.poly().centroid();

                    for (int i = 0; i < NVAR; ++i)
                    {
                        mc.prim_(i) = donor_cell.prim(i) + limiter(i) * dot(grad[i], d);
                    }
                }
                else
                {
                    mc.prim_ = donor_cell.prim();
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

                auto vg = donor_cell.vgn();
                assert(vg(0) == 0.);
                assert(vg(1) == 0.);
                assert(vg(2) == 0.);
                mc.prim_(1) -= vg(0);
                mc.prim_(2) -= vg(1);
                mc.prim_(3) -= vg(2);

                //mc.prim_ = donor_cell.prim_;

                assert(!mc.prim_.isnan());
            }
        }
    }

    void Solver::calc_steady()
    {
        auto &sp = partition_->spc_->sp_.front();

        std::vector<double> maxres_init(sp.mesh_.size(), 0.);

        for (int ntimestep = 0; ntimestep < maxtimestep_; ++ntimestep)
        {

            if (var_exc_ != nullptr)
            {
                var_exc_->update(profiler_, "sol-ghost-exc");
            }

            if (donor_var_exc_ != nullptr)
            {
                donor_var_exc_->update(profiler_, "sol-donor-exc");
            }

            if (profiler_ != nullptr)
            {
                profiler_->start("sol-1");
            }

            double local_res = TAILOR_BIG_NEG_NUM;

            for (int i = 0; i < sp.mesh_.size(); ++i)
            {
                Mesh &mesh = sp.mesh_[i];


                if (sp.mesh_.size() != 1)
                {
                    if (donor_var_exc_ == nullptr)
                    {
                        //mesh.oga_interpolate(sp.mesh(), comm_->rank(), gradient_);
                        oga_interpolate(mesh);
                    }
                    else
                    {
                        mesh.oga_interpolate(donor_var_exc_->arrival(), comm_->rank());
                    }
                }

                bc_.set_bc(mesh, profiler_);

                if (comm_->size() != 1)
                {
                    mesh.update_ghost_primitives(var_exc_->arrival(), comm_->rank(), fs_.gamma_);
                }

                if (ntimestep == 0 || steady_)
                {
                    compute_sum_of_fluxes(mesh);
                }
                else
                {
                    mesh.reset_to_mid();
                }

                if (progressive_cfl_)
                {
                    if (ncfl_increase_ == cfl_increase_freq_)
                    {
                        cfl_ += delta_cfl_;
                        ncfl_increase_ = 0.;
                    }
                    ++ncfl_increase_;
                }

                temporal_discretization(mesh);

                auto res = resi(mesh);
                local_res = std::max(local_res, max(res));

                calc_change_in_conserved_var(mesh);
                evolve_solution_in_time(mesh);
                evolve_old_solution_in_time(mesh);
            }

            if (profiler_ != nullptr)
            {
                profiler_->stop("sol-1");
            }

            if (profiler_ != nullptr)
            {
                profiler_->bstart("sol-reduce-resi");
            }
            double maxres;
            //boost::mpi::all_reduce(*comm_, local_res, maxres, std::plus<double>());
            boost::mpi::all_reduce(*comm_, local_res, maxres, boost::mpi::maximum<double>());

            if (nsolve_ ==  0) {
                init_max_res_ = maxres;
                last_max_res_ = maxres;
            }

            if (maxres / last_max_res_ > 100.)
            {
                cfl_ *= 10.;
                last_max_res_ = maxres;
            }

            if (profiler_ != nullptr)
            {
                profiler_->bstop("sol-reduce-resi");
            }

            if (profiler_ != nullptr)
            {
                profiler_->start("sol-print-norm");
            }
            
            
            if (print_outer_norm_)
            {
                if (comm_->rank() == 0)
                {
                    std::ofstream outer_norm;
                    std::string fn = "norm-";
                    //fn.append(std::to_string(mesh.tag()()));
                    //fn.append("-");
                    fn.append(std::to_string(nsolve_));
                    fn.append(".dat");
                    outer_norm.open(fn, std::ios_base::app);
                    //outer_norm << ntimestep << std::scientific << " " << maxres << " " << res[0] << " " << res[1] << " " << res[2] << " " << res[3] << " " << res[4] << std::fixed << std::endl;
                    outer_norm << ntimestep << std::scientific << " " << maxres << std::fixed << std::endl;
                    outer_norm.close();
                }
            }
            if (profiler_ != nullptr)
            {
                profiler_->stop("sol-print-norm");
            }

            if (ntimestep != 0)
            {
                if (maxres < tol_)
                {
                    //all_meshes_solved = false;
                    break;
                }
            }

            //if (all_meshes_solved) {
            //break;
            //}
        }

        //assert(all_meshes_solved);
    }

    void Solver::solve_(int rank)
    {
        assert(partition_->spc_->sp().size() == 1);
        auto &sp = partition_->spc_->sp_.front();

        for (Mesh &mesh : sp.mesh_)
        {
            for (MeshCell &mc : mesh.cell_)
            {
                assert(mc.oga_cell_type() != OGA_cell_type_t::undefined);
            }
        }

        if (profiler_ != nullptr)
        {
            profiler_->start("sol-facevel");
        }
        for (Mesh &mesh : sp.mesh_)
        {
            mesh.calc_mesh_velocities(fs_, comm_->rank());

            if (nsolve_ == 0)
            {
                if (sorder_ == 2)
                {
                    gradient_.calc_ls_coef(mesh);
                }
            }

            for (MeshCell &mc : mesh.cell_)
            {
                assert(mc.oga_cell_type() != OGA_cell_type_t::undefined);
            }
            
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
        if (profiler_ != nullptr)
        {
            profiler_->stop("sol-facevel");
        }

        //for (double currenttime = 0.; currenttime < finaltime_; currenttime += dt_)
        {
            //std::cout << "time: " << currenttime + dt_ << std::endl;

            calc_steady();

            //if (steady_)
            //{
            //    for (Mesh &mesh : sp.mesh_)
            //    {
            //        for (MeshCell &mc : mesh.cell_)
            //        {
            //            mc.cons_nm1_ = mc.cons_n_;
            //            mc.cons_n_ = mc.cons_sp1_;
            //        }
            //    }
            //}

            //if (!steady_)
            //{
                //return;
            //}
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

    /*Matrix5 Jacobian(const Vector5& q, const Vector3& n, double vgn, double gamma)
    {
        // this is not Roe jacobian.

        Matrix5 A;

        double gs = gamma - 1.;

        double nx = n(0);
        double ny = n(1);
        double nz = n(2);

        double kq = spec_kine_energy(q[1], q[2], q[3]);
        double qn = q[1]*nx + q[2]*ny + q[3]*nz;

        A(0,0) = -vgn;
        A(0,1) = nx;
        A(0,2) = ny;
        A(0,3) = nz;
        A(0,4) = 0.;

        A(1,0) = ( -1. / std::pow(q[0],2) ) * ( q[1]*qn - gs*kq*nx );
        A(1,1) = ( 1. / q[0] ) * (q[1]*nx*(3.-gamma) + q[2]*ny + q[3]*nz ) - vgn;
        A(1,2) = ( 1. / q[0] ) * (q[1]*ny - gs*q[2]*nx );
        A(1,3) = ( 1. / q[0] ) * (q[1]*nz - gs*q[3]*nx );
        A(1,4) = gs * nx;

        A(2,0) = ( -1. / std::pow(q[0],2) ) * ( q[2]*qn - gs*kq*ny );
        A(2,1) = ( 1. / q[0] ) * (q[2]*nx - gs*q[1]*ny );
        A(2,2) = ( 1. / q[0] ) * (q[1] * nx + q[2]*ny*(3.-gamma) + q[3]*nz ) - vgn;
        A(2,3) = ( 1. / q[0] ) * (q[2]*nz - gs*q[3]*ny );
        A(2,4) = gs * ny;

        A(3,0) = ( -1. / std::pow(q[0],2) ) * ( q[3]*qn - gs*kq*nz );
        A(3,1) = ( 1. / q[0] ) * (q[3]*nx - gs*q[1]*nz );
        A(3,2) = ( 1. / q[0] ) * (q[3]*ny - gs*q[2]*nz );
        A(3,3) = ( 1. / q[0] ) * (q[1]*nx + q[2]*ny + q[3]*nz*(3.-gamma) ) - vgn;
        A(3,4) = gs * nz;

        A(4,0) = -( q[4]*qn / std::pow(q[0],2) ) * gamma + 2.*qn*gs*kq/std::pow(q[0],3);
        A(4,1) = ( q[4]/q[0] ) * gamma * nx - gs*q[1]/q[0];
        A(4,2) = ( q[4]/q[0] ) * gamma * ny - gs*q[2]/q[0];
        A(4,3) = ( q[4]/q[0] ) * gamma * nz - gs*q[3]/q[0];
        A(4,4) = ( qn / q[0] ) * gamma - vgn;

        return A;
    }*/

    /*Matrix5 Jacobian(const Vector5& q, const Vector3& n, double vgn)
    {
        // this is not Roe jacobian.

        Matrix5 A;

        double gs = GAMMA - 1.;

        double nx = n(0);
        double ny = n(1);
        double nz = n(2);

        double kq = spec_kine_energy(q[1], q[2], q[3]);
        double qn = q[1]*nx + q[2]*ny + q[3]*nz;

        double q0_r = 1. / q[0];
        double q0_rs  = 1. / std::pow(q[0],2);
        double g40 = ( q[4]/q[0] ) * GAMMA;
        double g0 = gs / q[0];

        A(0,0) = -vgn;
        A(0,1) = nx;
        A(0,2) = ny;
        A(0,3) = nz;
        A(0,4) = 0.;

        A(1,0) = -q0_rs * ( q[1]*qn - gs*kq*nx );
        A(1,1) = q0_r * (q[1]*nx*(3.-GAMMA) + q[2]*ny + q[3]*nz ) - vgn;
        A(1,2) = q0_r * (q[1]*ny - gs*q[2]*nx );
        A(1,3) = q0_r * (q[1]*nz - gs*q[3]*nx );
        A(1,4) = gs * nx;

        A(2,0) = -q0_rs * ( q[2]*qn - gs*kq*ny );
        A(2,1) = q0_r * (q[2]*nx - gs*q[1]*ny );
        A(2,2) = q0_r * (q[1] * nx + q[2]*ny*(3.-GAMMA) + q[3]*nz ) - vgn;
        A(2,3) = q0_r * (q[2]*nz - gs*q[3]*ny );
        A(2,4) = gs * ny;

        A(3,0) = -q0_rs * ( q[3]*qn - gs*kq*nz );
        A(3,1) = q0_r * (q[3]*nx - gs*q[1]*nz );
        A(3,2) = q0_r * (q[3]*ny - gs*q[2]*nz );
        A(3,3) = q0_r * (q[1]*nx + q[2]*ny + q[3]*nz*(3.-GAMMA) ) - vgn;
        A(3,4) = gs * nz;

        A(4,0) = -( q[4]*qn / std::pow(q[0],2) ) * GAMMA + 2.*qn*gs*kq/std::pow(q[0],3);
        A(4,1) = g40 * nx - g0*q[1];
        A(4,2) = g40 * ny - g0*q[2];
        A(4,3) = g40 * nz - g0*q[3];
        A(4,4) = ( qn / q[0] ) * GAMMA - vgn;

        return A;
    }*/

    /*void Solver::rhhl_update_matrices(const Matrix<NVAR, NVAR>& Aroe, MeshFace* myface, MeshFace* commonface, MeshCell& LC, MeshCell& RC, const Vector3& n, double vgn, double facearea, double SLm, double SRp, double vfn)
    {
        auto JL = Jacobian(LC.prim(), fs_.gamma_, vfn);
        auto JR = Jacobian(RC.prim(), fs_.gamma_, vfn);

        double tempL = SRp / (SRp - SLm);
        double tempR = -SLm / (SRp - SLm);
        myface->M_ = (JL * tempL + Aroe * 0.5) * facearea;
        if (commonface != nullptr)
        {
            commonface->M_ = (JR * tempR - Aroe * 0.5) * facearea * -1;
        }

        LC.D_ += myface->M_;
        LC.D_mid_ += myface->M_;
        if (commonface != nullptr)
        {
            RC.D_ += commonface->M_;
            RC.D_mid_ += commonface->M_;
        }
    }*/

    /*void Solver::hhl_update_matrices(MeshFace* myface, MeshFace* commonface, MeshCell& LC, MeshCell& RC, double facearea, double SL, double SR, const Matrix5& T, double vfn)
    {
        auto JL = Jacobian(T * LC.prim(), fs_.gamma_, vfn);
        auto JR = Jacobian(T * RC.prim(), fs_.gamma_, vfn);

        Matrix5 I;
        I = 0.;
        I(0,0) = 1.;
        I(1,1) = 1.;
        I(2,2) = 1.;
        I(3,3) = 1.;
        I(4,4) = 1.;

        //myface->M_ = (JL - I * SLm) * (SRp * facearea / (SRp - SLm));
        //myface->M_ = JL * (facearea * (SR - SL * SR) / (SR - SL));
        myface->M_ = (JL - I * SL) * (facearea * SR / (SR - SL));
        if (commonface != nullptr)
        {
            //commonface->M_ = (JR * -1 + I * SRp) * (SLm * facearea * -1 / (SRp - SLm));
            //commonface->M_ = JR * (facearea * (SL * SR - SL) / (SR - SL));
            commonface->M_ = (JR - I * SR) * (facearea * SL / (SR - SL));
        }

        LC.D_ += myface->M_;
        LC.D_mid_ += myface->M_;
        if (commonface != nullptr)
        {
            RC.D_ += commonface->M_;
            RC.D_mid_ += commonface->M_;
        }
    }*/

    void Solver::update_matrices(MeshFace *this_face, MeshFace *common_face, MeshCell& left_cell, MeshCell& right_cell, double facearea, const Vector3& face_velocity, double gamma)
    {
        auto [left_state, right_state] = left_and_right_states(left_cell.cons_sp1(), right_cell.cons_sp1(), gamma, unit_matrix<NVAR, NVAR, double>(), face_velocity);

        Matrix5 JL = Jacobian(left_state , gamma);
        Matrix5 JR = Jacobian(right_state, gamma);

        this_face->M_ = JL * 0.5 * facearea;
        if (common_face != nullptr)
        {
            common_face->M_ = JR * 0.5 * facearea * -1;
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

    void Solver::apply_limiter(Mesh &mesh, MeshCell &mc, const MeshFace &mf)
    {
        if (mc.oga_cell_type() == OGA_cell_type_t::non_resident || mc.oga_cell_type() == OGA_cell_type_t::ghost)
        {
            //prim = mc.prim();
            return;
        }

        for (const auto &mff : mc.face())
        {
            if (mff.btype() != BouType::interior)
            {
                //prim = mc.prim();
                return;
            }
            else
            {
                if (mff.parent_cell().size() != 2)
                {
                    std::cout << "mc oga: " << static_cast<int>(mc.oga_cell_type()) << std::endl;
                    std::cout << "mc btype: " << static_cast<int>(mc.btype()) << std::endl;
                }
                assert(mff.parent_cell().size() == 2);
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
            mc.prim_(i) = mc.prim(i) + limiter(i) * dot(grad[i], d);
        }

        assert(!mc.prim().isnan());

        mc.cons_sp1_ = prim_to_cons(mc.prim(), fs_.gamma_);
    }

    void Solver::RK4(Mesh &mesh)
    {
        Vector5 k1, k2, k3, k4;

        for (MeshCell &mc : mesh.cell_)
        {
            if (!calc_cell(mc))
            {
                continue;
            }

            double vol = mc.poly().volume();

            k1 = mc.dtao_ * mc.R_ / vol;
            mc.cons_sp1_ = mc.cons_s_ + 0.5 * k1;
            mc.prim_ = cons_to_prim(mc.cons_sp1_, fs_.gamma_);
        }

        /*bc_.set_bc(mesh, profiler_);
        if (comm_->size() != 1)
        {
            mesh.update_ghost_primitives(var_exc_->receiver().front().arrival(), comm_->rank(), fs_.gamma_);
        }
        calc_R(mesh);
        tempo_discre(mesh, false);

        for (MeshCell& mc: mesh.cell_)
        {
            if (mc.oga_cell_type() != OGA_cell_type_t::field) {
                continue;
            }

            double vol = mc.poly().volume();

            k2 = mc.dtao_ * mc.R_ / vol;
            mc.cons_sp1_ = mc.cons_s_ + 0.5 * k2;
            mc.prim_ = cons_to_prim(mc.cons_sp1_, fs_.gamma_);
        }

        bc_.set_bc(mesh, profiler_);
        if (comm_->size() != 1)
        {
            mesh.update_ghost_primitives(var_exc_->receiver().front().arrival(), comm_->rank(), fs_.gamma_);
        }
        calc_R(mesh);
        tempo_discre(mesh, false);

        for (MeshCell& mc: mesh.cell_)
        {
            if (mc.oga_cell_type() != OGA_cell_type_t::field) {
                continue;
            }

            double vol = mc.poly().volume();

            k3 = mc.dtao_ * mc.R_ / vol;
            mc.cons_sp1_ = mc.cons_s_ + k3;
            mc.prim_ = cons_to_prim(mc.cons_sp1_, fs_.gamma_);
        }

        bc_.set_bc(mesh, profiler_);
        if (comm_->size() != 1)
        {
            mesh.update_ghost_primitives(var_exc_->receiver().front().arrival(), comm_->rank(), fs_.gamma_);
        }
        calc_R(mesh);
        tempo_discre(mesh, false);*/

        for (MeshCell &mc : mesh.cell_)
        {
            if (!calc_cell(mc))
            {
                continue;
            }

            double vol = mc.poly().volume();

            k4 = mc.dtao_ * mc.R_ / vol;
            //mc.cons_sp1_ = mc.cons_s_ + (k1 + 2. * k2 + 2. * k3 + k4) / 6.;
            mc.cons_sp1_ = mc.cons_s_ + k1;
            mc.prim_ = cons_to_prim(mc.cons_sp1_, fs_.gamma_);
        }
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
            //if (mc.oga_cell_type() == OGA_cell_type_t::field)
            {
                cell_index.insert(std::make_pair(mc.tag()(), n));
                ++n;
            }
        }

        n = n * NVAR;

        nz_num = 0;
        for (const auto &mc : mesh.cell())
        {
            if (calc_cell(mc))
            //if (mc.oga_cell_type() == OGA_cell_type_t::field)
            {
                int rindex = cell_index[mc.tag()()];

                std::map<int, Matrix5> order;
                order.insert(std::make_pair(rindex, mc.D()));

                for (const auto &mf : mc.face())
                {
                    if (mf.btype() == BouType::interior)
                    {
                        auto nct = std::find_if(mf.parent_cell().begin(), mf.parent_cell().end(), [&](const auto &pc) { return pc != mc.tag(); });
                        assert(nct != mf.parent_cell().end());
                        auto &nc = mesh.cell(*nct);

                        if (calc_cell(nc))
                        //if (nc.oga_cell_type() == OGA_cell_type_t::field)
                        {
                            order.insert(std::make_pair(cell_index[nc.tag()()], mf.M() * -1));
                        }
                    }
                }

                //std::sort(order.begin(), order.end(), [&](const auto& a, const auto& b){a.first < b.first;}); // map is instrisicly already ordered.

                for (int i = 0; i < NVAR; ++i)
                {
                    for (int j = 0; j < NVAR; ++j)
                    {
                        for (const auto &entry : order)
                        {
                            int cindex = entry.first;
                            double val = entry.second(i, j);
                            if (std::abs(val) > TAILOR_ZERO)
                            {
                                ia.push_back(rindex * NVAR + i);
                                ja.push_back(cindex * NVAR + j);
                                a.push_back(val);
                                ++nz_num;
                            }
                        }
                    }
                }
            }
        }

        rhs.resize(n);
        int i = 0;
        for (const auto &mc : mesh.cell())
        {
            if (calc_cell(mc))
            //if (mc.oga_cell_type() == OGA_cell_type_t::field)
            {
                for (int j = 0; j < NVAR; ++j)
                {
                    rhs[i * NVAR + j] = mc.R(j);
                }
                ++i;
            }
        }

        assert(ia.size() == nz_num);
        assert(ja.size() == nz_num);
        assert(a.size() == nz_num);

        //auto max_ia = std::max_element(ia.begin(), ia.end());
        //if (*max_ia > n)
        //{
        //    std::cout << "max_ia: " << *max_ia << std::endl;
        //    std::cout << "n: " << n << std::endl;
        //}
        //assert(*max_ia <= n);

        //if (rank == 1)
        //{
        //    auto temp = ia;
        //    auto it = std::unique(temp.begin(), temp.end());
        //    temp.erase(it, temp.end());
        //    if (temp.size() != n)
        //    {
        //        std::cout << "temp size: " << temp.size() << std::endl;
        //        std::cout << "ia size: " << ia.size() << std::endl;
        //        std::cout << "n: " << n << std::endl;

        //        for (int i: temp)
        //        {
        //            std::cout << "temp: " << i << std::endl;
        //        }
        //    }
        //    assert(temp.size() == n);
        //}
    }

    void make_mat_cr(const Mesh &mesh, int &n, int &nz_num, std::vector<int> &ia, std::vector<int> &ja, std::vector<double> &a, std::vector<double> &rhs, int rank)
    {
        make_mat_st(mesh, n, nz_num, ia, ja, a, rhs, rank);

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
        }
        assert(rp.size() == (n + 1));

        ia = rp;
        assert(ia.size() == rp.size());
    }

    void Solver::gmres(Mesh &mesh)
    {
        int n, nz_num;
        std::vector<int> ia, ja;
        std::vector<double> a;
        std::vector<double> rhs;

        //make_mat_st(mesh, n, nz_num, ia, ja, a, rhs);
        make_mat_cr(mesh, n, nz_num, ia, ja, a, rhs, comm_->rank());

        boost::property_tree::ptree prm;

        prm.put("solver.type", "bicgstab");
        prm.put("solver.tol", 1e-3);
        prm.put("solver.maxiter", 10);
        prm.put("precond.coarsening.type", "smoothed_aggregation");
        prm.put("precond.relax.type", "spai0");

        AMGCLSolver solve(std::tie(n, ia, ja, a), prm);

        if (n == 0)
        {
            int nfield = 0;
            for (const auto &mc : mesh.cell())
            {
                if (calc_cell(mc))
                //if (mc.oga_cell_type() == OGA_cell_type_t::field)
                {
                    ++nfield;
                }
            }
            assert(nfield == 0);
        }

        if (n == 0)
        {
            return;
        }

        std::vector<double> x(n, 0.);
        //int itr_max = 1;
        //int mr = 100;
        //if (n < mr)
        //{
        //    mr = n - 1;
        //}
        //double tol_abs = TAILOR_ZERO;
        ////double tol_rel = 1e6;
        //double tol_rel = TAILOR_BIG_POS_NUM;

        //bool verbose = true;
        ////if (comm_->rank() == 1) {
        ////verbose = true;
        ////}

        ////mgmres_st(n, nz_num, ia.data(), ja.data(), a.data(), x.data(), rhs.data(), itr_max, mr, tol_abs, tol_rel, verbose);
        //pmgmres_ilu_cr(n, nz_num, ia.data(), ja.data(), a.data(), x.data(), rhs.data(), itr_max, mr, tol_abs, tol_rel, verbose);

        //int i = 0;
        //for (auto &mc : mesh.cell_)
        //{
        //    if (calc_cell(mc))
        //    //if (mc.oga_cell_type() == OGA_cell_type_t::field)
        //    {
        //        for (int j = 0; j < NVAR; ++j)
        //        {
        //            mc.dQ_(j) = x[i * NVAR + j];
        //        }
        //        ++i;
        //    }
        //}

        auto [iters, error] = solve(rhs, x);
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

    void Solver::temporal_discretization(Mesh &mesh)
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
                if (dual_ts_ || steady_)
                {
                    double t = volume / mc.dtao_;
                    mc.D_.add_diag(t);

                    mc.R_ -= volume * (mc.cons_sp1() - mc.cons_n()) / dt_;
                }
                else
                {
                    double t = volume / dt_;
                    mc.D_.add_diag(t);
                }
            }
            else if (torder_ == 2)
            {
                double t = volume * (1. / mc.dtao_ + 1.5 / dt_);
                mc.D_.add_diag(t);

                if (dual_ts_)
                {
                    if (nsolve_ > 1)
                    {
                        mc.R_ -= 0.5 * volume * (3. * mc.cons_sp1() - 4. * mc.cons_n() + mc.cons_nm1()) / dt_;
                    }
                    else
                    {
                        mc.R_ -= volume * (mc.cons_sp1() - mc.cons_n()) / dt_;
                    }
                }
                else
                {
                    if (nsolve_ > 1)
                    {
                        mc.R_ -= 0.5 * volume * (-4. * mc.cons_n() + mc.cons_nm1()) / dt_;
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
