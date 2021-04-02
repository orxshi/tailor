#ifndef VAR_EXCHANGER_H
#define VAR_EXCHANGER_H

#include "exchanger.h"
#include "request.h"
#include "sp.h"
#include "var.h"

namespace Tailor
{
    class VarExchanger: public Exchanger<Var>
    {
        public:

            //VarExchanger(const std::vector<SpatialPartition>* sp, const RecvCon<Request>*, const boost::mpi::communicator*);
            VarExchanger(const std::vector<SpatialPartition>* sp, const Exchanger<Request>*, const boost::mpi::communicator*);
            void update(Mesh& mesh, Profiler* profiler=nullptr, std::string fname="");

        private:

            const boost::mpi::communicator* comm_;
            const std::vector<SpatialPartition>* sp_;
            //const RecvCon<Request>* other_receiver_;
            const Exchanger<Request>* other_exc_;

            void prepare_storage();
    };
}

#endif
