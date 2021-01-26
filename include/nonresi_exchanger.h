#ifndef NON_RESI_EXCHANGER_H
#define NON_RESI_EXCHANGER_H

#include "exchanger.h"
#include "spc.h"
#include "request.h"

namespace Tailor
{
    class NonResiExchanger: public Exchanger<Request>
    {
        public:

            NonResiExchanger(const SpatialPartitionContainer* spc, boost::mpi::communicator* comm);

        private:

            boost::mpi::communicator* comm_;
            const SpatialPartitionContainer* spc_;

            void prepare_storage();
    };
}

#endif
