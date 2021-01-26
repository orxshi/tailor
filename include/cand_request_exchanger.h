#ifndef TAILOR_CAND_REQUEST_EXCHANGER_H
#define TAILOR_CAND_REQUEST_EXCHANGER_H

#include "exchanger.h"
#include "meshcell.h"
#include "spc.h"
#include "request.h"

namespace Tailor
{
    class CandReqExchanger: public Exchanger<Request>
    {
        public:

            CandReqExchanger(const SpatialPartitionContainer* spc, boost::mpi::communicator*, bool mergebins);

        private:

            bool mergebins_;
            boost::mpi::communicator* comm_;
            const SpatialPartitionContainer* spc_;

            void prepare_storage();
    };
}

#endif
