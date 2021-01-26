#ifndef DI_EXCHANGER_H
#define DI_EXCHANGER_H

#include "exchanger.h"
#include "assembler.h"
#include "solver.h"
#include "donor_info2.h"

namespace Tailor
{
    class DIExchanger: public Exchanger<DonorInfo2>
    {
        public:

            DIExchanger(boost::mpi::communicator* comm, const Assembler& assembler, const Solver& solver);

        private:

            boost::mpi::communicator* comm_;
            const Assembler& assembler_;
            const Solver& solver_;

            void prepare_storage();
    };
}

#endif
