#ifndef TAILOR_TCR_EXCHANGER_H
#define TAILOR_TCR_EXCHANGER_H

#include "exchanger.h"
#include "spc.h"
#include "request.h"
#include "donor_info2.h"

namespace Tailor
{
    class TCRExchanger: public Exchanger<Request>
    {
        public:

            TCRExchanger(SpatialPartitionContainer* spc, boost::mpi::communicator* comm, const ArrCon<DonorInfo2>* arrival);

        private:

            boost::mpi::communicator* comm_;
            SpatialPartitionContainer* spc_;
            const ArrCon<DonorInfo2>* arrival_;

            void prepare_storage();
    };
}

#endif
