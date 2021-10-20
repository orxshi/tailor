#include "tailor.h"

namespace Tailor
{
    const boost::mpi::communicator& Tailor::comm() const
    {
        return comm_;
    }

    std::string make_serialization_folder(int time_step, std::string save_folder)
    {
        std::string sf = save_folder;
        sf.append(std::to_string(time_step));
        std::filesystem::create_directory(sf);
        return sf;
    }

    void serialize(const Assembler& assembler, int rank, std::string file_name)
    {
        file_name.append("/assembler-");
        file_name.append(std::to_string(rank));
        file_name.append(".ser");

        std::ofstream ofs(file_name);
        //boost::archive::text_oarchive oa(ofs);
        boost::archive::binary_oarchive oa(ofs);
        oa << assembler;
    }

    void deserialize(Assembler& assembler, boost::mpi::communicator* comm, Profiler* profiler, std::string file_name)
    {
        assert(comm != nullptr);

        file_name.append("/assembler-");
        file_name.append(std::to_string(comm->rank()));
        file_name.append(".ser");

        std::ifstream ifs(file_name);
        //boost::archive::text_iarchive ia(ifs);
        boost::archive::binary_iarchive ia(ifs);
        ia >> assembler;

        assembler.set_comm(comm);

        if (profiler != nullptr) {
            assembler.set_profiler(profiler);
        }

        assembler.read_settings();
        assembler.partition()->spc().global_rm().update_address();
    }

    void serialize(const Solver& solver, int rank, std::string file_name)
    {
        file_name.append("/solver-");
        file_name.append(std::to_string(rank));
        file_name.append(".ser");

        std::ofstream ofs(file_name);
        //boost::archive::text_oarchive oa(ofs);
        boost::archive::binary_oarchive oa(ofs);
        oa << solver;
    }

    void deserialize(Assembler* assembler, Solver& solver, boost::mpi::communicator* comm, Profiler* profiler, bool use_shared_partition, std::string file_name)
    {
        assert(comm != nullptr);

        file_name.append("/solver-");
        file_name.append(std::to_string(comm->rank()));
        file_name.append(".ser");

        std::ifstream ifs(file_name);
        //boost::archive::text_iarchive ia(ifs);
        boost::archive::binary_iarchive ia(ifs);
        ia >> solver;

        if (use_shared_partition)
        {
            assert(assembler != nullptr);
            solver.set_partition(assembler->partition());
        }

        solver.set_comm(comm);

        if (profiler != nullptr) {
            solver.set_profiler(profiler);
        }

        solver.read_settings();

        if (!use_shared_partition)
        {
            assert(solver.partition() != nullptr);
            solver.partition()->spc().global_rm().update_address();
        }

        solver.read_settings();
    }

    void Tailor::rotate(const Tag& mesh, double ang, int axis, const Vector3& pivot)
    {
        if (assembler_on_)
        {
            if (!use_shared_partition_ || !solver_on_)
            {
                assembler_->rotate(mesh, ang, axis, pivot);
            }
        }

        if (solver_on_)
        {
            solver_->rotate(mesh, ang, axis, pivot);
        }
    }

    void Tailor::move(const Tag& mesh, const Vector3& v)
    {
        if (assembler_on_)
        {
            if (!use_shared_partition_ || !solver_on_)
            {
                assembler_->move(mesh, v);
            }
        }

        if (solver_on_)
        {
            solver_->move(mesh, v);
        }
    }

    Tailor::Tailor(): assembler_on_(true), solver_on_(true)
    {
        read_settings();

        assert(!mesh_folder_.empty());

        if (mesh_folder_.size() == 1)
        {
            assembler_on_ = false;
            use_shared_partition_ = false;
        }

        if (profiler_on_)
        {
            profiler_ = std::make_unique<Profiler>(&comm_, false);
        }

        if (load_)
        {
            if (assembler_on_)
            {
                assembler_ = std::make_unique<Assembler>();
                deserialize(*assembler_, &comm_, nullptr, load_folder_);
            }
            if (solver_on_)
            {
                solver_ = std::make_unique<Solver>();
                deserialize(assembler_.get(), *solver_, &comm_, nullptr, use_shared_partition_, load_folder_);
            }
        }
        else
        {
            if (assembler_on_)
            {
                profiler_start("create_assembler");
                assembler_ = std::make_unique<Assembler>(&comm_, nullptr, mesh_folder_);
                profiler_stop("create_assembler");
            }

            if (solver_on_)
            {
                profiler_start("create_solver");
                if (use_shared_partition_)
                {
                    solver_ = std::make_unique<Solver>(&comm_, mesh_folder_, nullptr, assembler_->partition());
                }
                else
                {
                    solver_ = std::make_unique<Solver>(&comm_, mesh_folder_, nullptr);
                }
                std::cout << "constructured solver" << std::endl;
                profiler_stop("create_solver");
            }
        }

        if (!assembler_on_ && solver_on_)
        {
            solver_->set_oga_cell_type_all_field();
        }
    }

    void Tailor::profiler_start(std::string s)
    {
        if (profiler_on_) {
            assert(profiler_ != nullptr);
            profiler_->start(s);
        }
    }

    void Tailor::profiler_stop(std::string s)
    {
        if (profiler_on_) {
            assert(profiler_ != nullptr);
            profiler_->stop(s);
        }
    }

    void Tailor::compute_aerodyn_coef(std::vector<AeroCoefPara> (*compute_para)())
    {
        if (solver_ == nullptr) {
            return;
        }

        if (compute_para == nullptr) {
            return;
        }

        auto aero_para = compute_para();

        solver_->partition()->spc().get_coef(aero_para, solver_->nsolve(), solver_->dt());
    }

    void Tailor::pre(int time_step, std::vector<AeroCoefPara> (*compute_para)())
    {
        if (assembler_on_)
        {
            assembler_->assemble();
            mem_usage(&comm_, "assemble");

            if (!use_shared_partition_ && solver_on_ && assembler_on_)
            {
                DIExchanger di_exc(&comm_, *assembler_, *solver_);
                di_exc.exchange(false, "asm-di", nullptr);

                solver_->transfer_oga_cell_type(di_exc.arrival());
            }
        }

        if (solver_on_)
        {
            if (solver_->nsolve() > 0)
            {
                if (assembler_on_)
                {
                    solver_->update_flow_field_after_mesh_motion();
                }
            }

            profiler_start("solve");
            solver_->print_settings(); 
            solver_->solve(max_time_step_);
            compute_aerodyn_coef(compute_para);
            mem_usage(&comm_,  "solve");
            profiler_stop("solve");
        }
    }

    void Tailor::post()
    {
        if (solver_on_)
        {
            profiler_start("exchange");
            solver_->exchange();
            profiler_stop("exchange");
            mem_usage(&comm_,  "exchange");

            solver_->reset_oga_status();

            profiler_start("reconnect");
            solver_->reconnectivity();
            profiler_stop("reconnect");
            mem_usage(&comm_,  "reconnect");

            if (!assembler_on_)
            {
                solver_->set_oga_cell_type_all_field();
            }
        }

        if (assembler_on_)
        {
            if (!use_shared_partition_ || !solver_on_)
            {
                assembler_->exchange();
                assembler_->reset_oga_status();
                assembler_->reconnectivity();
            }
        }

        if (solver_on_)
        {
            profiler_start("repartition");
            solver_->repartition();
            profiler_stop("repartition");
        }
        mem_usage(&comm_,  "repartition");

        if (assembler_on_)
        {
            if (!use_shared_partition_ || !solver_on_)
            {
                assembler_->repartition();
            }
        }
    }

    void Tailor::save(int time_step, int& save_counter)
    {
        if (save_)
        {
            ++save_counter;

            if (save_counter == save_interval_)
            {
                std::string save_folder;

                if (assembler_)
                {
                    save_folder = make_serialization_folder(assembler_->nassemble(), save_folder_);
                }
                else if (solver_)
                {
                    save_folder = make_serialization_folder(solver_->nsolve(), save_folder_);
                }

                if (assembler_)
                {
                    serialize(*assembler_, comm_.rank(), save_folder);
                }
                if (solver_)
                {
                    serialize(*solver_, comm_.rank(), save_folder);
                }
            }

            if (save_counter > save_interval_) {
                save_counter = 0;
            }
        }
    }

    void Tailor::read_settings()
    {
        namespace po = boost::program_options;

        po::options_description all_options;

        po::options_description desc{"Tailor options"};
        desc.add_options()
            ("tailor.load", po::value<bool>()->default_value(false), "")
            ("tailor.save", po::value<bool>()->default_value(false), "")
            ("tailor.use_shared_partition", po::value<bool>()->default_value(true), "")
            ("tailor.save_interval", po::value<int>()->default_value(int(1e15)), "")
            ("tailor.max_time_step", po::value<int>()->default_value(1), "")
            ("tailor.load_folder", po::value<std::string>()->default_value(""), "")
            ("tailor.save_folder", po::value<std::string>()->default_value("sv-"), "")
            ("tailor.mesh_folder", po::value<std::vector<std::string>>()->multitoken(), "")
            ("tailor.profiler", po::value<bool>()->default_value(false), "")
            ("tailor.solver", po::value<bool>()->default_value(true), "")
            ("tailor.only-motion", po::value<bool>()->default_value(false), "")
            ;

        all_options.add(desc);

        std::ifstream settings_file("settings.ini");

        boost::program_options::variables_map vm;
        po::store(po::parse_config_file(settings_file, all_options, true), vm);
        po::notify(vm);

        save_interval_ = vm["tailor.save_interval"].as<int>();
        max_time_step_ = vm["tailor.max_time_step"].as<int>();
        mesh_folder_ = vm["tailor.mesh_folder"].as<std::vector<std::string>>();
        load_ = vm["tailor.load"].as<bool>();
        save_ = vm["tailor.save"].as<bool>();
        use_shared_partition_ = vm["tailor.use_shared_partition"].as<bool>();
        load_folder_ = vm["tailor.load_folder"].as<std::string>();
        save_folder_ = vm["tailor.save_folder"].as<std::string>();
        profiler_on_ = vm["tailor.profiler"].as<bool>();
        solver_on_ = vm["tailor.solver"].as<bool>();
        only_motion_ = vm["tailor.only-motion"].as<bool>();
    }

    void Tailor::make(void (*callback)(Tailor&), std::vector<AeroCoefPara> (*compute_para)())
    {
        int save_counter = 0;

        assert(max_time_step_ > 0);

        for (int time_step = 0; time_step < max_time_step_; ++time_step)
        {
            if (!only_motion_)
            {
                pre(time_step, compute_para);
            }
            if (max_time_step_ > 1)
            {
                if (callback != nullptr) {
                    callback(*this);
                }
            }
            if (!only_motion_)
            {
                post();
            }
            save(time_step, save_counter);

            if (profiler_on_)
            {
                profiler_->print_iter(time_step);
                profiler_->clear_iter();
            }
        }

        if (solver_on_)
        {
            solver_->print_settings();
        }
    }

    const Solver* Tailor::solver() const
    {
        assert(solver_);
        return solver_.get();
    }

    const Assembler* Tailor::assembler() const
    {
        assert(assembler_);
        return assembler_.get();
    }
}
