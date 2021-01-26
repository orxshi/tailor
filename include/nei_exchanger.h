#ifndef NEI_EXCHANGER_H
#define NEI_EXCHANGER_H

#include "exchanger.h"
//#include "facerequest.h"
#include "request.h"
#include "spc.h"
#include "nei.h"

namespace Tailor
{
    class NeiExchanger: public Exchanger<Nei>
    {
        public:

            NeiExchanger(const SpatialPartitionContainer* spc, ArrCon<Request>*, boost::mpi::communicator*);

            void remove_dup();

        private:

            boost::mpi::communicator* comm_;
            const SpatialPartitionContainer* spc_;
            ArrCon<Request>* other_arrival_;

            void prepare_storage();
    };
}

#endif
