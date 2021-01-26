#ifndef TAILOR_SP_EXCHANGER_H
#define TAILOR_SP_EXCHANGER_H

#include "exchanger.h"
#include "sp.h"
#include "loadmap.h"

namespace Tailor
{
    class SpExchanger: public Exchanger<SpatialPartition>
    {
        public:

            SpExchanger(const Loadmap* lm, boost::mpi::communicator* comm);

        private:

            boost::mpi::communicator* comm_;
            const Loadmap* lm_;

            void prepare_storage();
    };
}

#endif
