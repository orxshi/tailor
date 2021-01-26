#ifndef DONOR_REQUEST_EXCHANGER_H
#define DONOR_REQUEST_EXCHANGER_H

#include "exchanger.h"
#include "meshcell.h"
#include "spc.h"
#include "request.h"
#include "donor_info2.h"

namespace Tailor
{
    class DonorRequestExchanger: public Exchanger<Request>
    {
        public:

            DonorRequestExchanger(const SpatialPartitionContainer* spc, SpatialPartition* sp, const ArrCon<DonorInfo2>* di, boost::mpi::communicator* comm);

        private:

            boost::mpi::communicator* comm_;
            SpatialPartition* sp_;
            const SpatialPartitionContainer* spc_;
            const ArrCon<DonorInfo2>* di_;

            void prepare_storage();
    };
}

#endif
