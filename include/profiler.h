#ifndef PROFILER_H
#define	PROFILER_H

#include "timer.h"
//#include "settings.h"
#include <boost/mpi.hpp>
#include <fstream>
#include <iomanip>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/trivial.hpp>
#include <boost/filesystem.hpp>
#include <sys/time.h>
#include <sys/resource.h>

namespace Tailor
{
    struct Function
    {
        std::string name;
        int elapsed_index_;
        int tier_index_;
    };

    class Tier
    {
        std::string tag_;
        std::vector<std::shared_ptr<Tier>> child_;

        public:

        Tier(std::string tier_string);
        Tier(std::vector<int>::const_iterator begin, std::vector<int>::const_iterator end);

        void add(std::vector<int>::const_iterator begin, std::vector<int>::const_iterator end);
        void add(std::string tier_string);
        std::string tag() const;
        static std::vector<int> to_int(std::string s);
        static std::string to_string(std::vector<int>::const_iterator begin, std::vector<int>::const_iterator end);
    };

    class Profiler
    {
        //std::shared_ptr<Settings> settings_;
        bool with_barrier_;
        std::vector<std::string> excluded_;
        std::vector<std::string> included_;
        //std::vector<Tier> tier_;
        std::vector<Tier> tier_;
        size_t nprobe_;
        std::vector<double> total_elapsed_;
        std::vector<double> total_elapsed_iter_;
        const boost::mpi::communicator* world_;
        //Timer timer;
        std::vector<double> elapsed_batch_;
        std::vector<std::vector<double>> elapsed_;
        std::vector<std::vector<double>> elapsed_main_;
        std::vector<std::vector<double>> elapsed_main_iter_;
        //std::vector<time_point> start_;
        std::vector<double> start_;
        //std::vector<time_point> stop_;
        std::vector<double> stop_;
        std::vector<std::pair<std::string, int>> func_;
        std::map<int, std::string> func_rev_;
        //std::tuple<std::string, int, int> func_; // name, index in elapsed_main_iter_, index in tier_
        //std::vector<Function> func_;
        std::map<std::string, std::shared_ptr<Tier>> tier_tag_address_;
        std::map<std::string, std::string> tier_func_;
        void add_tier(std::string tier_string);
        bool exclude(std::string func_name) const;
        bool include(std::string func_name) const;
        double skewness(const std::pair<std::string, int>& func);
        void gather(int f);
        void update();
        void gather_batch();
        void update_batch();

        public:

        //Profiler(const MPI_Comm& comm, bool with_barrier, const Settings& settings);
        //Profiler(const MPI_Comm& comm, bool with_barrier);
        Profiler(const boost::mpi::communicator* comm, bool with_barrier);
        Profiler() = default;

        void clear_iter();
        void print_iter(int iteration);
        const size_t nprobe() const;
        void print() const;
        void start(std::string func_name);
        void stop(std::string func_name);
        void bstart(std::string func_name);
        void bstop(std::string func_name);
        const std::vector<std::pair<std::string, int>>& func() const;
        //const Function& func() const;

        /*template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            ar & with_barrier_;
            ar & nprobe_;
            ar & total_elapsed_;
            ar & total_elapsed_iter_;
            ar & elapsed_batch_;
            ar & elapsed_;
            ar & elapsed_main_;
            ar & elapsed_main_iter_;
            ar & start_;
            ar & stop_;
            ar & func_;
            ar & func_rev_;
        }*/
    };

    void mem_usage(const boost::mpi::communicator* comm, std::string s);
}

#endif
