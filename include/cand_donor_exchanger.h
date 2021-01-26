#ifndef CAND_DONOR_EXCHANGER_H
#define CAND_DONOR_EXCHANGER_H

#include "exchanger.h"
#include "spc.h"
#include "donor_info2.h"
#include "request.h"

namespace Tailor
{
    class CandDonorExchanger: public Exchanger<DonorInfo2>
    {
        public:

            CandDonorExchanger(const SpatialPartitionContainer* spc, const ArrCon<Request>*, boost::mpi::communicator* comm);
            void update(Profiler* profiler=nullptr, std::string s="");

        private:

            boost::mpi::communicator* comm_;
            const SpatialPartitionContainer* spc_;
            const ArrCon<Request>* arrival_;

            void prepare_storage();
    };
}

#endif
