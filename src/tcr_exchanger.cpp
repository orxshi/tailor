#include "tcr_exchanger.h"

namespace Tailor
{
    TCRExchanger::TCRExchanger(SpatialPartitionContainer* spc, boost::mpi::communicator* comm, const ArrCon<DonorInfo2>* arrival): comm_(comm), Exchanger(comm), spc_(spc), arrival_(arrival)
    {
    }

    bool no_valid_cand(const MeshCell& mc, const SpatialPartition& sp, const ArrCon<DonorInfo2>& arrival)
    {
        for (auto cand_donor = mc.cand_donor().begin(); cand_donor != mc.cand_donor().end(); ++cand_donor)
        {
            const auto& cdm = cand_donor->mesh_tag_;
            const auto& cdc = cand_donor->cell_tag_;

            auto mesh_it = std::find_if(sp.mesh().begin(), sp.mesh().end(), [&](const Mesh& m){return m.tag() == cdm;});
            assert(mesh_it != sp.mesh().end());
            //assert(mesh_it->tag() != m.tag());

            //auto cit = mesh_it->query(Tag(cdc));
            assert(mesh_it->query(Tag(cdc)) != nullptr);
            const auto& cit = mesh_it->cell(Tag(cdc));

            OGA_cell_type_t type = OGA_cell_type_t::undefined;
            //assert(cit != nullptr);

            if (cit.oga_cell_type() == OGA_cell_type_t::non_resident)
            {
                //assert(!sp.is_resident(cit.poly().centroid()));
                //for (const auto& receiver: nonresi)
                {
                    //if (receiver.tag() == sp.tag()())
                    {
                        //auto nr = std::find_if(receiver.arrival().begin(), receiver.arrival().end(), [&](const auto& n){return (n.meshtag_ == cdm() && n.celltag_ == cdc());});
                        auto nr = std::find_if(arrival.begin(), arrival.end(), [&](const auto& n){return (n.meshtag_ == cdm() && n.celltag_ == cdc());});
                        //assert(nr != receiver.arrival().end());
                        assert(nr != arrival.end());
                        type = nr->oga_cell_type_;
                        //break;
                    }
                }
            }
            else
            {
                type = cit.oga_cell_type();
            }

            assert(type != OGA_cell_type_t::undefined);

            if (type != OGA_cell_type_t::receptor && type != OGA_cell_type_t::mandat_receptor && type != OGA_cell_type_t::orphan)
            {
                return false;
            }
        }

        return true;
    }

    void TCRExchanger::prepare_storage()
    {
        assert(!spc_->sp().empty());

        //sender_.reserve(spc_->sp().size());
        //receiver_.reserve(spc_->sp().size());

        for (int i=0; i<spc_->sp().size(); ++i)
        {
            SpatialPartition& sp = spc_->sp_[i];

            //sender_.push_back(Sender<Request>(sp.tag()()));
            //receiver_.push_back(Receiver<Request>(sp.tag()()));

            for (Mesh& m: sp.mesh_)
            {
                for (MeshCell& mc: m.cell_)
                {
                    //if (mc.tag()() == 25688 && m.tag()() == 1 && comm_->rank() == 1)
                    //{
                        //std::cout << "rankkkkkkkkkkkkkk: " << comm_->rank() << std::endl;
                        //std::cout << "typeeeeeeeeeeeeee: " << static_cast<int>(mc.oga_cell_type()) << std::endl;
                        //assert(false);
                    //}
                    auto type = mc.oga_cell_type();
                    if (type == OGA_cell_type_t::mandat_receptor)
                    {
                        assert(!mc.cand_donor().empty());

                        if (no_valid_cand(mc, sp, *arrival_))
                        {
                            for (auto cand_donor = mc.cand_donor().begin(); cand_donor != mc.cand_donor().end();)
                            {
                                const auto& cdm = cand_donor->mesh_tag_;
                                const auto& cdc = cand_donor->cell_tag_;

                                auto mesh_it = std::find_if(sp.mesh_.begin(), sp.mesh_.end(), [&](const Mesh& m){return m.tag() == cdm;});
                                assert(mesh_it != sp.mesh_.end());
                                assert(mesh_it->tag() != m.tag());

                                assert(mesh_it->query(Tag(cdc)) != nullptr);
                                auto& cit = mesh_it->cell_p(Tag(cdc));

                                OGA_cell_type_t type = OGA_cell_type_t::undefined;
                                //assert(cit != nullptr);

                                if (cit.oga_cell_type() == OGA_cell_type_t::non_resident)
                                {
                                    //assert(!sp.is_resident(cit.poly().centroid()));
                                    //for (const auto& receiver: *nonresi_)
                                    {
                                        //if (receiver.tag() == sp.tag()())
                                        {
                                            //auto nr = std::find_if(receiver.arrival().begin(), receiver.arrival().end(), [&](const auto& n){return (n.meshtag_ == cdm() && n.celltag_ == cdc());});
                                            auto nr = std::find_if(arrival_->begin(), arrival_->end(), [&](const auto& n){return (n.meshtag_ == cdm() && n.celltag_ == cdc());});
                                            //sender_[i].add(nr->source_rank_, nr->source_tag_, Request(nr->meshtag_, nr->celltag_, comm_->rank(), sp.tag()()));
                                            add(nr->source_rank_, nr->source_tag_, Request(nr->meshtag_, nr->celltag_, comm_->rank(), sp.tag()()));
                                            //sender_[i].add(nr->source_rank_, nr->source_tag_, Request(nr->donor_mesh_, nr->donor_cell_, comm_->rank(), sp.tag()()));
                                            //break;
                                        }
                                    }
                                    break;
                                }
                                else
                                {
                                    cit.set_oga_cell_type(OGA_cell_type_t::field);
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
