#include "assembler.h" 
#include "solver.h" 
#include "profiler.h" 
#include "di_exchanger.h" 
#include "coef.h" 
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <filesystem>

const int SAVE_THRES = 0; // 74

namespace bmpi = boost::mpi;

Tailor::Solver* solver; // made global since signal handler cannot take arguments.
Tailor::Assembler* assembler; // made global since signal handler cannot take arguments.
bmpi::communicator* comm;

void save(const Tailor::Assembler& assembler, const Tailor::Solver& solver, boost::mpi::communicator* comm, std::string savefolder)
{
    {
        std::string sername = savefolder;
        sername.append("/assembler-");
        sername.append(std::to_string(comm->rank()));
        sername.append(".ser");

        std::ofstream ofs(sername);
        boost::archive::text_oarchive oa(ofs);
        oa << assembler;
    }

    {
        std::string sername = savefolder;
        sername.append("/solver-");
        sername.append(std::to_string(comm->rank()));
        sername.append(".ser");

        std::ofstream ofs(sername);
        boost::archive::text_oarchive oa(ofs);
        oa << solver;
    }
}

void load(Tailor::Assembler& assembler, Tailor::Solver& solver, boost::mpi::communicator* comm, Tailor::Profiler* profiler, bool use_shared_partition, std::string savefolder)
{
    {
        std::string sername = savefolder;
        sername.append("/assembler-");
        sername.append(std::to_string(comm->rank()));
        sername.append(".ser");

        std::ifstream ifs(sername);
        boost::archive::text_iarchive ia(ifs);
        ia >> assembler;

        assembler.set_comm(comm);
        if (profiler != nullptr) {
            assembler.set_profiler(profiler);
        }
        assembler.read_settings();
        assembler.partition()->spc().global_rm().update_address();
    }

    {
        std::string sername = savefolder;
        sername.append("/solver-");
        sername.append(std::to_string(comm->rank()));
        sername.append(".ser");

        std::ifstream ifs(sername);
        boost::archive::text_iarchive ia(ifs);
        ia >> solver;

        if (use_shared_partition)
        {
            solver.set_partition(assembler.partition());
        }

        solver.set_comm(comm);
        if (profiler != nullptr) {
            solver.set_profiler(profiler);
        }

        solver.read_settings();

        if (!use_shared_partition)
        {
            solver.partition()->spc().global_rm().update_address();
        }
    }
}

std::string make_save_folder(int i)
{
    std::string save_folder = "sv-";
    save_folder.append(std::to_string(i));
    std::filesystem::create_directory(save_folder);
    return save_folder;
}

//void signal_handler(int signum)
//{
//    assert(assembler != nullptr);
//    assert(solver != nullptr);
//    std::cout << "Interrupt signal (" << signum << ") received. Saving the latest solution." << std::endl;
//    save(*assembler, *solver, comm);
//    exit(signum);
//}

int main()
{
    //signal(SIGTERM, signal_handler);
    bmpi::environment env;
    comm = new bmpi::communicator;

    {
        std::string mems = "main";
        Tailor::mem_usage(comm, mems);
    }

    Tailor::Profiler profiler(comm, false);

    std::vector<std::string> fn;
    fn.push_back("msh/64/fuspyl");
    fn.push_back("msh/64/wing0");
    fn.push_back("msh/64/wing1");
    fn.push_back("msh/64/wing2");
    fn.push_back("msh/64/wing3");
    fn.push_back("msh/64/hubshaft");

    Tailor::Freestream fs;
    fs.read();

    Tailor::Component compo;
    compo.read(Tailor::Tag(1));

    bool use_shared_partition = true;

    //assembler = new Tailor::Assembler();
    //solver = new Tailor::Solver();
    //load(*assembler, *solver, comm, nullptr, use_shared_partition, "sv-2249");

    int counter = 0;

    for (int i=0; i<1; ++i)
    {
        std::string a = "iter-";
        a.append(std::to_string(i));
        profiler.start(a);

        if (i == 0)
        {
            //profiler.start("asm-ctr");
            //assembler = new Tailor::Assembler(comm, &profiler, fn);
            assembler = new Tailor::Assembler(comm, nullptr, fn);
            //assembler = new Tailor::Assembler();
            //profiler.stop("asm-ctr");

            //profiler.start("sol-ctr");
            if (use_shared_partition)
            {
                //solver = new Tailor::Solver(comm, fn, &profiler, assembler->partition());
                solver = new Tailor::Solver(comm, fn, nullptr, assembler->partition());
                //solver = new Tailor::Solver();
            }
            else
            {
                //solver = new Tailor::Solver(comm, fn, &profiler);
                solver = new Tailor::Solver(comm, fn, nullptr);
            }
            //profiler.stop("sol-ctr");
        }

        {
            std::string mems = "start-iter-";
            mems.append(std::to_string(i));
            Tailor::mem_usage(comm, mems);
        }

        //profiler.start("asm-asm");
        assembler->assemble();
        //profiler.stop("asm-asm");

        if (!use_shared_partition)
        {
            //profiler.start("asm-di");
            Tailor::DIExchanger di_exc(comm, *assembler, *solver);
            di_exc.exchange(false, "asm-di", &profiler);
            //profiler.stop("asm-di");

            //profiler.start("sol-di");
            solver->transfer_oga_cell_type(di_exc.arrival());
            //profiler.stop("sol-di");
        }

        {
            std::string mems = "di-iter-";
            mems.append(std::to_string(i));
            Tailor::mem_usage(comm, mems);
        }

        //profiler.start("sol-sol");
        solver->solve();
        solver->partition()->spc().get_coef(fs, i, solver->dt());
        //profiler.stop("sol-sol");
        //save(*assembler, *solver, comm);

        double rpm = compo.rpm_;
        double om = rpm * 2. * Tailor::PI / 60.; // rad/s
        double aoa = om * solver->dt(); // rad/s * time step
        int axis = compo.rotaxis_;
        Tailor::vec3<double> pivot(compo.pivot_(0), compo.pivot_(1), compo.pivot_(2));

        if (!use_shared_partition)
        {
            //profiler.start("asm-rot");
            assembler->rotate(Tailor::Tag(1), aoa, axis, pivot);
            assembler->rotate(Tailor::Tag(2), aoa, axis, pivot);
            assembler->rotate(Tailor::Tag(3), aoa, axis, pivot);
            assembler->rotate(Tailor::Tag(4), aoa, axis, pivot);
            assembler->rotate(Tailor::Tag(5), aoa, axis, pivot);
            //profiler.stop("asm-rot");
        }

        //profiler.start("sol-rot");
        solver->rotate(Tailor::Tag(1), aoa, axis, pivot);
        solver->rotate(Tailor::Tag(2), aoa, axis, pivot);
        solver->rotate(Tailor::Tag(3), aoa, axis, pivot);
        solver->rotate(Tailor::Tag(4), aoa, axis, pivot);
        solver->rotate(Tailor::Tag(5), aoa, axis, pivot);
        //profiler.stop("sol-rot");

        {
            std::string mems = "rotate-iter-";
            mems.append(std::to_string(i));
            Tailor::mem_usage(comm, mems);
        }

        //profiler.start("sol-exc");
        solver->exchange();
        //profiler.stop("sol-exc");

        {
            std::string mems = "exchange-";
            mems.append(std::to_string(i));
            Tailor::mem_usage(comm, mems);
        }

        //profiler.start("sol-reset");
        solver->reset_oga_status();
        //profiler.stop("sol-reset");

        {
            std::string mems = "reset-";
            mems.append(std::to_string(i));
            Tailor::mem_usage(comm, mems);
        }

        //profiler.start("sol-recon");
        solver->reconnectivity();
        //profiler.stop("sol-recon");

        {
            std::string mems = "recon-";
            mems.append(std::to_string(i));
            Tailor::mem_usage(comm, mems);
        }

        if (!use_shared_partition)
        {
            //profiler.start("asm-exc");
            assembler->exchange();
            //profiler.stop("asm-exc");

            //profiler.start("asm-reset");
            assembler->reset_oga_status();
            //profiler.stop("asm-reset");

            //profiler.start("asm-recon");
            assembler->reconnectivity();
            //profiler.stop("asm-recon");
        }

        bool sol_repart = false;
        bool asm_repart = false;

        //profiler.start("sol-repart");
        sol_repart = solver->repartition();
        //profiler.stop("sol-repart");

        {
            std::string mems = "repart-";
            mems.append(std::to_string(i));
            Tailor::mem_usage(comm, mems);
        }

        if (!use_shared_partition)
        {
            //profiler.start("asm-repart");
            asm_repart = assembler->repartition();
            //profiler.stop("asm-repart");
        }

        {
            std::string mems = "end-iter-";
            mems.append(std::to_string(i));
            Tailor::mem_usage(comm, mems);
        }

        profiler.stop(a);

        profiler.print_iter(i);
        profiler.clear_iter();

        if (counter == SAVE_THRES) {
            auto save_folder = make_save_folder(i);
            save(*assembler, *solver, comm, save_folder);
        }

        ++counter;

        if (counter > SAVE_THRES) {
            counter = 0;
        }
    }

    return 0;
}
