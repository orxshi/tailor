#include "nei_exchanger.h"

namespace Tailor
{
    Nei::Nei(int source_rank, int source_tag, const MeshCell& cell): source_rank_(source_rank), source_tag_(source_tag), cell_(cell)
    {
        //cell_.push_back(cell);
    }

    //const std::vector<MeshCell>& Nei::cell() const
    //{
    //return cell_;
    //}

    const MeshCell& Nei::cell() const
    {
        return cell_;
    }

    //void Nei::remove_dup()
    //{
    //std::sort(cell_.begin(), cell_.end(), [&](const MeshCell& a, const MeshCell& b){return a.tag() < b.tag();});
    //cell_.erase(std::unique(cell_.begin(), cell_.end()), cell_.end());
    //}

    bool Nei::operator==(const Nei& other) const
    {
        if (cell_.parent_mesh() == other.cell().parent_mesh())
        {
            if (cell_.tag() == other.cell().tag())
            {
                return true;
            }
        }

        return false;
    }

    bool Nei::operator<(const Nei& other) const
    {
        if (cell_.parent_mesh() < other.cell().parent_mesh())
        {
            return true;
        }
        else if (cell_.parent_mesh() == other.cell().parent_mesh())
        {
            if (cell_ < other.cell())
            {
                return true;
            }
        }

        return false;
    }

    void NeiExchanger::remove_dup()
    {
        //for (auto& s: storage_)
        //for (auto& s: receiver_)
        {
            //{
            //for (const auto& cell: s.arrival())
            //{
            //    const auto& mc = cell.cell_;
            //    {
            //        if (mc.parent_mesh()() == 0)
            //        {
            //            if (mc.tag()() == 4929)
            //            {
            //                {
            //                    for (const auto& ff: mc.face())
            //                    {
            //                        if (ff.tag() == FaceTag(std::vector<int>{521, 1319, 1566}))
            //                        {
            //                            std::cout << "mmmmvvvvvvvvvvvvvvvvv btype: " << static_cast<int>(ff.btype()) << std::endl;
            //                            std::cout << "mmmmrrrrrrrrrrrrrrrank: " << comm_->rank() << std::endl;
            //                            //assert(ff.btype() == BouType::interior);
            //                        }
            //                    }
            //                }
            //            }
            //        }
            //    }
            //}
            //}
            //for (auto& nei: s.arrival_p())
            //{
            // TODO Don't delete mesh cells which have same tag but belong to different parent meshes!
            //std::sort(s.arrival_p().begin(), s.arrival_p().end());
            std::sort(arrival_.begin(), arrival_.end());
            //s.arrival_p().erase(std::unique(s.arrival_p().begin(), s.arrival_p().end()), s.arrival_p().end());
            arrival_.erase(std::unique(arrival_.begin(), arrival_.end()), arrival_.end());

            //nei.remove_dup();
            //}
        }
    }

    NeiExchanger::NeiExchanger(const SpatialPartitionContainer* spc, ArrCon<Request>* other_arrival, boost::mpi::communicator* comm): comm_(comm), Exchanger(comm), spc_(spc), other_arrival_(other_arrival)
    {
    }

    void NeiExchanger::prepare_storage()
    {
        const auto* sp_ = &(spc_->sp());

        assert(sp_->size() == 1);
        //assert(other_receiver_->size() == sp_->size());

        const SpatialPartition& sp = sp_->front();
        //sender_.push_back(Sender<Nei>(sp.tag()()));
        //receiver_.push_back(Receiver<Nei>(sp.tag()()));
        //const Receiver<Request>& os = other_receiver_->front();

        //assert(os.tag() == (sp.tag()()));

        for (const auto& r: *other_arrival_)
        {
            auto m = std::find_if(sp.mesh().begin(), sp.mesh().end(), [&](const Mesh& mm){return mm.tag()() == r.mesh_cell_.first;});
            if (m == sp.mesh().end()) {
                continue;
            }
            auto mc = std::find_if(m->cell().begin(), m->cell().end(), [&](const MeshCell& cc){return cc.tag()() == r.mesh_cell_.second;});
            assert(mc != m->cell().end());

            for (const Tag& inei: mc->pnei())
            {
                const MeshCell& nc = m->cell(inei);
                //sender_.front().add(r.source_rank_, r.source_tag_, Nei(comm_->rank(), sp.tag()(), nc));
                add(r.source_rank_, r.source_tag_, Nei(comm_->rank(), sp.tag()(), nc));
            }
        }
    }
}
