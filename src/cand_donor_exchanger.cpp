#include "cand_donor_exchanger.h"

namespace Tailor
{
    CandDonorExchanger::CandDonorExchanger(const SpatialPartitionContainer* spc, const ArrCon<Request>* arrival, boost::mpi::communicator* comm): comm_(comm), Exchanger(comm), spc_(spc), arrival_(arrival)
    {
    }

    void CandDonorExchanger::prepare_storage()
    {
        // if a cell is receptor, if one of its cand donor is non-resident, transfer it.
        // first we need to exchange requests
        // then exchange cand donors based on request.

        //sender_.reserve(spc_->global_rm().size());
        //receiver_.reserve(spc_->global_rm().size());

        //for (const auto& sp: spc_->sp())
        //{
            //sender_.push_back(Sender<DonorInfo2>(sp.tag()()));
            //receiver_.push_back(Receiver<DonorInfo2>(sp.tag()()));
        //}

        for (int i=0; i<spc_->sp().size(); ++i)
        //for (const auto& sp: spc_->sp())
        {
            const auto& sp = spc_->sp()[i];

            //for (const auto& stor: *other_receiver_)
            {
                //if (stor.tag() != sp.tag()()) {
                    //continue;
                //}

                //for (const Request& r: stor.arrival())
                for (const Request& r: *arrival_)
                {
                    auto m = std::find_if(sp.mesh().begin(), sp.mesh().end(), [&](const Mesh& mm){return mm.tag()() == r.mesh_cell_.first;});
                    if (m == sp.mesh().end())
                    {
                        std::cout << "sp mesh size:" << sp.mesh().size() << std::endl;
                        std::cout << "mesh tag:" << r.mesh_cell_.first << std::endl;
                    }
                    assert(m != sp.mesh().end());
                    auto mc = std::find_if(m->cell().begin(), m->cell().end(), [&](const MeshCell& cc){return cc.tag()() == r.mesh_cell_.second;});
                    assert(mc != m->cell().end());

                    assert(m->tag().isvalid());
                    assert(mc->tag().isvalid());
                    //sender_[i].add(r.source_rank_, r.source_tag_, DonorInfo2(m->tag()(), mc->tag()(), mc->oga_cell_type(), mc->donor().mesh_tag_(), mc->donor().cell_tag_(), comm_->rank(), comm_->rank(), sp.tag()()));
                    add(r.source_rank_, r.source_tag_, DonorInfo2(m->tag()(), mc->tag()(), mc->oga_cell_type(), mc->donor().mesh_tag_(), mc->donor().cell_tag_(), comm_->rank(), comm_->rank(), sp.tag()()));
                }
            }
        }
    }

    void CandDonorExchanger::update(Profiler* profiler, std::string s)
    {
        if (profiler != nullptr) {profiler->start(s);}
        //assert(sender_.size() == spc_->sp().size());
        for (int i=0; i<spc_->sp().size(); ++i)
        {
            const auto& sp = spc_->sp()[i];
            //auto& s = sender_[i];

            //for (auto& gr: s.group())
            for (auto& gr: group())
            {
                for (auto& v: gr.data_)
                {
                    auto m = std::find_if(sp.mesh().begin(), sp.mesh().end(), [&](const Mesh& mm){return mm.tag()() == v.meshtag_;});
                    assert(m != sp.mesh().end());
                    const auto& mc = m->cell(Tag(v.celltag_));
                    v.oga_cell_type_ = mc.oga_cell_type();
                    v.donor_mesh_ = mc.donor().mesh_tag_();
                    v.donor_cell_ = mc.donor().cell_tag_();
                }
            }
        }
        if (profiler != nullptr) {profiler->stop(s);}

        send_recv(false, -1, s, profiler);
    }
}
