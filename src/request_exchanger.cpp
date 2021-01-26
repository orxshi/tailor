#include "request_exchanger.h"

namespace Tailor
{
    Request::Request(int mesh, int cell, int source_rank, int source_tag): source_rank_(source_rank), source_tag_(source_tag)
    {
        //mesh_cell_.push_back(std::make_pair(mesh, cell));
        mesh_cell_ = std::make_pair(mesh, cell);
    }

    RequestExchanger::RequestExchanger(const SpatialPartitionContainer* spc, boost::mpi::communicator* comm, Profiler* profiler): comm_(comm), Exchanger(comm), spc_(spc), profiler_(profiler)
    {
    }

    void RequestExchanger::prepare_storage()
    {
        assert(!spc_->sp().empty());
        assert(spc_->sp().size() == 1);

        const SpatialPartition& sp = spc_->sp().front();
        assert(sp.tag()() == comm_->rank());

        //sender_.push_back(Sender<Request>(sp.tag()()));
        //receiver_.push_back(Receiver<Request>(sp.tag()()));

        for (const Mesh& m: sp.mesh())
        {
            for (const MeshCell& mc: m.cell())
            {
                if (mc.oga_cell_type() == OGA_cell_type_t::non_resident || mc.btype() != BouType::partition) {
                    continue;
                }

                std::pair<Tag, int> resbin;
                auto map = spc_->search_rm(mc.poly(), mc.poly().centroid(), resbin, false);

                for (const auto& pair: map)
                {
                    int dest_rank = pair.second;
                    if (dest_rank == comm_->rank())
                    {
                        continue;
                    }

                    //sender_.front().add(dest_rank, dest_rank, Request(m.tag()(), mc.tag()(), comm_->rank(), sp.tag()()));
                    add(dest_rank, dest_rank, Request(m.tag()(), mc.tag()(), comm_->rank(), sp.tag()()));
                }
            }
        }
    }
}
