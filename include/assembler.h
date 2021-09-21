#ifndef TAILOR_ASSEMBLER_H
#define TAILOR_ASSEMBLER_H

#include "partition.h"
#include "donor_search.h"
#include "cell_exchanger.h"
#include "disconnecter.h"
#include "arrival_connecter.h"
//#include "donor_info.h"
#include "cand_request_exchanger.h"
#include "cand_donor_exchanger.h"

namespace Tailor
{
    class Assembler
    {
        public:

            Assembler(boost::mpi::communicator* comm, Profiler*, const std::vector<std::string>&);
            Assembler();
            ~Assembler();
            void assemble();
            const Partition* partition() const;
            Partition* partition();
            void rotate(const Tag& mesh, double ang, int axis, const Vector3& pivot);
            void move(const Tag& mesh, const Vector3& v);
            void exchange();
            void reset_oga_status();
            void reconnectivity();
            bool repartition();
            void set_comm(boost::mpi::communicator* comm);
            void set_profiler(Profiler* prof);
            void read_settings();
            int nassemble() const;

            template<class Archive> void serialize(Archive & ar, const unsigned int version)
            {
                ar & rebalance_thres_;
                ar & print_ds_info_;
                ar & print_vtk_;
                ar & print_pre_vtk_;
                ar & print_repart_info_;
                ar & print_imbalance_;
                ar & mergebins_;
                ar & donor_search_algo_;
                ar & nassemble_;
                ar & load_estim_type_;
                ar & can_rebalance_;
                ar & force_rebalance_;
                ar & make_load_balance_;
                ar & pseudo3D_;
                ar & print_map_;
                ar & partition_;
                ar & nlayer_of_overlap_;
                ar & overlap_minimization_;
            }

        private:

            //DonorInfo* donor_info_;
            //std::deque<Mesh> wall_;
            //std::deque<Mesh> dirichlet_;
            //std::deque<Mesh> farfield_;
            //std::deque<Mesh> empty_;
            //std::deque<Mesh> mesh_;
            double rebalance_thres_;
            bool print_ds_info_;
            bool print_vtk_;
            bool print_pre_vtk_;
            bool print_repart_info_;
            bool print_imbalance_;
            bool mergebins_;
            DonorSearchAlgo donor_search_algo_;
            Profiler* profiler_;
            Partition* partition_;
            int nassemble_;
            boost::mpi::communicator* comm_;
            LoadEstim load_estim_type_;
            bool can_rebalance_;
            bool force_rebalance_;
            bool make_load_balance_;
            bool pseudo3D_;
            bool print_map_;
            int nlayer_of_overlap_;
            bool overlap_minimization_;

            void donor_search();
            void print_mesh_vtk(std::string fn0);
    };
}

#endif
