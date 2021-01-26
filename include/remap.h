#ifndef TAILOR_REMAP_H
#define TAILOR_REMAP_H

#include "loadmap.h"
#include "spc.h"
#include "sp_exchanger.h"

namespace Tailor
{
    class Remapper
    {
        public:

            Remapper(SpatialPartitionContainer* spc, boost::mpi::communicator* comm, std::string name_);
            void remap(Loadmap& lm, bool mergebins, int iter, std::string, Profiler*);

        private:

            std::string name_;
            boost::mpi::communicator* comm_;
            SpatialPartitionContainer* spc_;
    };
}

#endif
