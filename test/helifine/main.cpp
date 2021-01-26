#include "assembler.h" 
#include "solver.h" 
#include "profiler.h" 
#include "di_exchanger.h" 
#include "coef.h" 
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

namespace bmpi = boost::mpi;

Tailor::Solver* solver; // made global since signal handler cannot take arguments.

void save(const Tailor::Assembler& assembler, const Tailor::Solver& solver, boost::mpi::communicator* comm)
{
    {
        std::string sername = "assembler-";
        sername.append(std::to_string(comm->rank()));
        sername.append(".ser");

        std::ofstream ofs(sername);
        boost::archive::text_oarchive oa(ofs);
        oa << assembler;
    }

    {
        std::string sername = "solver-";
        sername.append(std::to_string(comm->rank()));
        sername.append(".ser");

        std::ofstream ofs(sername);
        boost::archive::text_oarchive oa(ofs);
        oa << solver;
    }
}

void load(Tailor::Assembler& assembler, Tailor::Solver& solver, boost::mpi::communicator* comm, Tailor::Profiler* profiler, bool use_shared_partition)
{
    {
        std::string sername = "assembler-";
        sername.append(std::to_string(comm->rank()));
        sername.append(".ser");

        std::ifstream ifs(sername);
        boost::archive::text_iarchive ia(ifs);
        ia >> assembler;

        assembler.set_comm(comm);
        assembler.set_profiler(profiler);
        assembler.partition()->spc().global_rm().update_address();
    }

    {
        std::string sername = "solver-";
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
        solver.set_profiler(profiler);

        if (!use_shared_partition)
        {
            solver.partition()->spc().global_rm().update_address();
        }
    }
}

int main()
{
    //std::cout << "inmain" << "\n"; 

    bmpi::environment env;
    bmpi::communicator comm;

    {
        std::string mems = "main";
        Tailor::mem_usage(&comm, mems);
    }

    Tailor::Profiler profiler(&comm, false);

    //std::cout << "profctr" << "\n"; 

    //std::cout << "profctr" << "\n"; 

    //std::vector<std::string> fn;
    //fn.push_back("msh/96/fuspyl");
    //fn.push_back("msh/96/wing0");
    //fn.push_back("msh/96/wing1");
    //fn.push_back("msh/96/wing2");
    //fn.push_back("msh/96/wing3");
    //fn.push_back("msh/96/hubshaft");

    Tailor::Freestream fs;
    fs.read();

    Tailor::Component compo;
    compo.read(Tailor::Tag(1));

    bool use_shared_partition = true;

    //profiler.start("asm-ctr");
    ////Tailor::Assembler assembler(&comm, &profiler, fn);
    //Tailor::Assembler assembler(&comm, nullptr, fn);
    //profiler.stop("asm-ctr");
    ////std::cout << "asmctr" << "\n"; 

    //profiler.start("sol-ctr");
    //if (use_shared_partition)
    //{
    //    //solver = new Tailor::Solver(&comm, fn, &profiler, assembler.partition());
    //    solver = new Tailor::Solver(&comm, fn, nullptr, assembler.partition());
    //}
    //else
    //{
    //    //solver = new Tailor::Solver(&comm, fn, &profiler);
    //    solver = new Tailor::Solver(&comm, fn, nullptr);
    //}
    ////std::cout << "solctr" << "\n"; 
    //profiler.stop("sol-ctr");

    Tailor::Assembler assembler;
    load(assembler, *solver, &comm, &profiler, use_shared_partition);

    for (int i=2; i<3; ++i)
    {
        //std::cout << "start-" << i << "\n"; 
        //profiler.start("asm-ctr");

        {
            std::string mems = "start-iter-";
            mems.append(std::to_string(i));
            Tailor::mem_usage(&comm, mems);
        }

        profiler.start("asm-asm");
        assembler.assemble();
        //std::cout << "assembled-" << i << "\n"; 
        profiler.stop("asm-asm");

        if (!use_shared_partition)
        {
            //profiler.start("asm-di");
            Tailor::DIExchanger di_exc(&comm, assembler, *solver);
            di_exc.exchange(false, "asm-di", &profiler);
            //profiler.stop("asm-di");

            //profiler.start("sol-di");
            solver->transfer_oga_cell_type(di_exc.arrival());
            //profiler.stop("sol-di");
        }

        {
            std::string mems = "di-iter-";
            mems.append(std::to_string(i));
            Tailor::mem_usage(&comm, mems);
        }

        profiler.start("sol-sol");
        solver->solve();
        //std::cout << "solved-" << i << "\n"; 
        solver->partition()->spc().get_coef(fs, i, solver->dt());
        //std::cout << "coef-" << i << "\n"; 
        profiler.stop("sol-sol");

        double rpm = compo.rpm_;
        double om = rpm * 2. * Tailor::PI / 60.; // rad/s
        double aoa = om * solver->dt(); // rad/s * time step
        int axis = compo.rotaxis_;
        Tailor::vec3<double> pivot(compo.pivot_(0), compo.pivot_(1), compo.pivot_(2));

        if (!use_shared_partition)
        {
            profiler.start("asm-rot");
            assembler.rotate(Tailor::Tag(1), aoa, axis, pivot);
            assembler.rotate(Tailor::Tag(2), aoa, axis, pivot);
            assembler.rotate(Tailor::Tag(3), aoa, axis, pivot);
            assembler.rotate(Tailor::Tag(4), aoa, axis, pivot);
            assembler.rotate(Tailor::Tag(5), aoa, axis, pivot);
            profiler.stop("asm-rot");
        }

        profiler.start("sol-rot");
        solver->rotate(Tailor::Tag(1), aoa, axis, pivot);
        solver->rotate(Tailor::Tag(2), aoa, axis, pivot);
        solver->rotate(Tailor::Tag(3), aoa, axis, pivot);
        solver->rotate(Tailor::Tag(4), aoa, axis, pivot);
        solver->rotate(Tailor::Tag(5), aoa, axis, pivot);
        profiler.stop("sol-rot");

        //std::cout << "rotated-" << i << "\n"; 

        {
            std::string mems = "rotate-iter-";
            mems.append(std::to_string(i));
            Tailor::mem_usage(&comm, mems);
        }

        profiler.start("sol-exc");
        solver->exchange();
        //std::cout << "exchanged-" << i << "\n"; 
        profiler.stop("sol-exc");

        {
            std::string mems = "exchange-";
            mems.append(std::to_string(i));
            Tailor::mem_usage(&comm, mems);
        }

        //profiler.start("sol-reset");
        solver->reset_oga_status();
        //std::cout << "reseted-" << i << "\n"; 
        //profiler.stop("sol-reset");

        {
            std::string mems = "reset-";
            mems.append(std::to_string(i));
            Tailor::mem_usage(&comm, mems);
        }

        profiler.start("sol-recon");
        solver->reconnectivity();
        //std::cout << "reconnected-" << i << "\n"; 
        profiler.stop("sol-recon");

        {
            std::string mems = "recon-";
            mems.append(std::to_string(i));
            Tailor::mem_usage(&comm, mems);
        }

        if (!use_shared_partition)
        {
            //profiler.start("asm-exc");
            assembler.exchange();
            //profiler.stop("asm-exc");

            //profiler.start("asm-reset");
            assembler.reset_oga_status();
            //profiler.stop("asm-reset");

            //profiler.start("asm-recon");
            assembler.reconnectivity();
            //profiler.stop("asm-recon");
        }

        bool sol_repart = false;
        bool asm_repart = false;

        profiler.start("sol-repart");
        sol_repart = solver->repartition();
        //std::cout << "repartitioned-" << i << "\n"; 
        profiler.stop("sol-repart");

        {
            std::string mems = "repart-";
            mems.append(std::to_string(i));
            Tailor::mem_usage(&comm, mems);
        }

        if (!use_shared_partition)
        {
            //profiler.start("asm-repart");
            asm_repart = assembler.repartition();
            //profiler.stop("asm-repart");
        }

        profiler.print_iter(i);
        profiler.clear_iter();
        //std::cout << "profiler-" << i << "\n"; 

        {
            std::string mems = "end-iter-";
            mems.append(std::to_string(i));
            Tailor::mem_usage(&comm, mems);
        }
    }

    //save(assembler, *solver, &comm);

    return 0;
}
