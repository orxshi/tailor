#ifndef TAILOR_TAILOR_H
#define TAILOR_TAILOR_H

#include "assembler.h"
#include "solver.h"
#include "di_exchanger.h"
#include <filesystem>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

namespace Tailor
{
    class Tailor
    {
        public:

            Tailor();

            void make(void (*callback)(Tailor&) = nullptr, std::vector<AeroCoefPara> (*compute_para)() = nullptr);
            void rotate(const Tag& mesh, double ang, const Vector3& axis, const Vector3& pivot);
            void move(const Tag& mesh, const Vector3& v);
            const Solver* solver() const;
            const Assembler* assembler() const;
            const boost::mpi::communicator& comm() const;

        private:

            bool only_motion_;
            int max_time_step_;
            int save_interval_;
            bool assembler_on_;
            bool solver_on_;
            bool profiler_on_;
            bool load_;
            bool save_;
            bool use_shared_partition_;
            std::vector<std::string> mesh_folder_;
            std::string load_folder_;
            std::string save_folder_;
            std::unique_ptr<Assembler> assembler_;
            std::unique_ptr<Solver> solver_;
            boost::mpi::communicator comm_;
            std::unique_ptr<Profiler> profiler_;
            boost::mpi::environment env_;

            void pre(int time_step, std::vector<AeroCoefPara> (*compute_para)());
            void post();
            void save(int time_step, int& save_counter);
            void read_settings();
            void profiler_start(std::string s);
            void profiler_stop(std::string s);
            void compute_aerodyn_coef(std::vector<AeroCoefPara> (*compute_para)() = nullptr);
    };

    std::string make_serialization_folder(int time_step, std::string save_folder);
    void serialize(const Assembler& assembler, int rank, std::string file_name);
    void deserialize(Assembler& assembler, boost::mpi::communicator* comm, Profiler* profiler, std::string file_name);
    void serialize(const Solver& solver, int rank, std::string file_name);
    void deserialize(Assembler* assembler, Solver& solver, boost::mpi::communicator* comm, Profiler* profiler, bool use_shared_partition, std::string file_name);
}

#endif
