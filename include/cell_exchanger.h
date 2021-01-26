#ifndef TAILOR_CELL_EXCHANGER_H
#define TAILOR_CELL_EXCHANGER_H

#include "exchanger.h"
#include "spc.h"
#include "cell.h"

namespace Tailor
{
    class CellExchanger: public Exchanger<Cell>
    {
        public:

            CellExchanger(SpatialPartitionContainer* spc, boost::mpi::communicator* comm, bool merge_bins);

        private:

            bool merge_bins_;
            const boost::mpi::communicator* comm_;
            SpatialPartitionContainer* spc_;

            void prepare_storage();
    };
}

#endif
