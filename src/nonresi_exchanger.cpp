#include "nonresi_exchanger.h"

namespace Tailor
{
    NonResiExchanger::NonResiExchanger(const SpatialPartitionContainer* spc, boost::mpi::communicator* comm): comm_(comm), Exchanger(comm), spc_(spc)
    {
    }

    void NonResiExchanger::prepare_storage()
    {
        assert(spc_->sp().size() == 1);

        const SpatialPartition& sp = spc_->sp().front();
        //storage_.push_back(Storage<Request>(sp.tag()()));
        //sender_.push_back(Sender<Request>(sp.tag()()));
        //receiver_.push_back(Receiver<Request>(sp.tag()()));

        for (const Mesh& m: sp.mesh())
        {
            for (const MeshCell& mc: m.cell())
            {
                //if (spc_->is_resident(mc.poly().centroid(), mc.tag()(), spt)) {
                //if (mc.oga_cell_type() != OGA_cell_type_t::non_resident) {
                if (mc.oga_cell_type() != OGA_cell_type_t::non_resident && mc.oga_cell_type() != OGA_cell_type_t::ghost) {
                //if (sp.is_resident(mc.poly().centroid())) {
                    continue;
                }

                //Tag spt;
                //if (spc_->is_resident(mc.poly().centroid(), mc.tag()(), spt))
                //{
                    //std::cout << "mc oga: " << static_cast<int>(mc.oga_cell_type()) << std::endl;
                //}
                //assert(!spc_->is_resident(mc.poly().centroid(), mc.tag()(), spt));

                //int temp;
                Tag tag;
                //spc_->global_rm().get_bintag_adaptive_unique(mc.poly().centroid(), tag, comm_->rank(), mc.tag()());
                int dest_rank;
                bool resu = spc_->search_rm(mc.poly().centroid(), tag, dest_rank);
                assert(resu);
                assert(tag.isvalid());
                //assert(tag.isvalid());

                //auto it = spc_->bintag_proc_map().find(tag);
                //assert(it != spc_->bintag_proc_map().end());
                //int dest_rank = it->second; // dest_rank = bintag

                //const Bin& _bin = spc_->global_rm().bin(tag);
                if (sp.tag()() == dest_rank) {
                    continue;
                }

                const Bin& _bin = spc_->global_rm().bin(tag);
                assert(_bin.aabb().do_intersect(AABB(mc.poly())));
                //if (_bin.aabb().do_intersect(AABB(mc.poly())))
                {
                    //storage_.front().add(dest_rank, dest_rank, Request(m.tag()(), mc.tag()(), comm_->rank(), sp.tag()()));
                    //sender_.front().add(dest_rank, dest_rank, Request(m.tag()(), mc.tag()(), comm_->rank(), sp.tag()()));
                    add(dest_rank, dest_rank, Request(m.tag()(), mc.tag()(), comm_->rank(), sp.tag()()));
                }
            }
        }
    }
}
