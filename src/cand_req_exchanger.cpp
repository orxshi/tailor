#include "cand_request_exchanger.h"

namespace Tailor
{
    CandReqExchanger::CandReqExchanger(const SpatialPartitionContainer* spc, boost::mpi::communicator* comm, bool mergebins): comm_(comm), Exchanger(comm), spc_(spc), mergebins_(mergebins)
    {
    }

    void CandReqExchanger::prepare_storage()
    {
        // if a cell is receptor, if one of its cand donor is non-resident, transfer it to get its type.
        // first we need to exchange requests
        // then exchange cand donors based on request.

        if (!mergebins_)
        {
            //sender_.reserve(spc_->global_rm().size());
            //receiver_.reserve(spc_->global_rm().size());
        }

        if (mergebins_)
        {
            assert(spc_->sp().size() == 1);
        }

        //for (const auto& sp: spc_->sp())
        //{
            //sender_.push_back(Sender<Request>(sp.tag()()));
            //receiver_.push_back(Receiver<Request>(sp.tag()()));
        //}

        for (int i=0; i<spc_->sp().size(); ++i)
        {
            const auto& sp = spc_->sp()[i];

            for (const Mesh& m: sp.mesh())
            {
                for (const MeshCell& mc: m.cell())
                {
                    if (mc.oga_cell_type() != OGA_cell_type_t::receptor && mc.oga_cell_type() != OGA_cell_type_t::mandat_receptor) {
                        continue;
                    }

                    //for (const Donor& cand: mc.cand_donor())
                    for (auto cand = mc.cand_donor().begin(); cand != mc.cand_donor().end(); ++cand)
                    {
                        const MeshCell* addr = cand->addr_;
                        assert(addr != nullptr);

                        //if (sp.is_resident(addr->poly().centroid())) {
                            //continue;
                        //}

                        //BinRMTag tag;
                        //int tag = -1;
                        //spc_->global_rm().get_bintag_adaptive_unique(addr->poly().centroid(), tag, comm_->rank(), -1);
                        //assert(tag.isvalid());
//
                        //auto it = spc_->bintag_proc_map().find(tag);
                        //assert(it != spc_->bintag_proc_map().end());
                        //int dest_rank = it->second;

                        int dest_rank;
                        Tag tag;
                        bool resu = spc_->search_rm(addr->poly().centroid(), tag, dest_rank); 
                        assert(resu);

                        // if resident continue.
                        if (mergebins_)
                        {
                            if (dest_rank == comm_->rank()) {
                                continue;
                            }
                        }
                        else
                        {
                            if (tag == sp.tag()) {
                                continue;
                            }
                        }

                        //sender_[i].add(dest_rank, tag.bintag()(), Request(addr->parent_mesh()(), addr->tag()(), comm_->rank(), sp.tag()()));
                        int desttag;
                        if (mergebins_) {
                            desttag = dest_rank;
                        }
                        else {
                            desttag = tag();
                        }
                        //sender_[i].add(dest_rank, desttag, Request(addr->parent_mesh()(), addr->tag()(), comm_->rank(), sp.tag()()));
                        add(dest_rank, desttag, Request(addr->parent_mesh()(), addr->tag()(), comm_->rank(), sp.tag()()));
                    }
                }
            }
        }
    }
}
