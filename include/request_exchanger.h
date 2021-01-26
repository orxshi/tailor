#ifndef REQUEST_EXCHANGER_H
#define REQUEST_EXCHANGER_H

#include "exchanger.h"
#include "meshcell.h"
#include "spc.h"
#include "request.h"

namespace Tailor
{
    class RequestExchanger: public Exchanger<Request>
    {
        public:

            RequestExchanger(const SpatialPartitionContainer* spc, boost::mpi::communicator*, Profiler* profiler=nullptr);

        private:

            Profiler* profiler_;
            boost::mpi::communicator* comm_;
            const SpatialPartitionContainer* spc_;

            void prepare_storage();
    };
}

#endif
