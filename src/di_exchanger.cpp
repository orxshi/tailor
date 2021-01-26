#include "di_exchanger.h"

namespace Tailor
{
    DonorInfo2::DonorInfo2(int meshtag, int celltag, OGA_cell_type_t type, int donor_mesh, int donor_cell, int donor_rank, int source_rank, int source_tag): meshtag_(meshtag), celltag_(celltag), oga_cell_type_(type), donor_mesh_(donor_mesh), donor_cell_(donor_cell), donor_rank_(donor_rank), source_rank_(source_rank), source_tag_(source_tag)
    {
    }

    DIExchanger::DIExchanger(boost::mpi::communicator* comm, const Assembler& assembler, const Solver& solver):
        comm_(comm),
        Exchanger(comm),
        assembler_(assembler),
        solver_(solver)
    {
    }

    void DIExchanger::prepare_storage()
    {
        const auto& a_partition = assembler_.partition();
        const auto& s_partition = solver_.partition();

        //sender_.reserve(a_partition->spc().global_rm().size());
        //receiver_.reserve(s_partition->spc().global_rm().size());

        //for (const auto& sp: a_partition->spc().sp())
        //{
        //    sender_.push_back(Sender<DonorInfo2>(sp.tag()() + comm_->size()));
        //}

        //for (const auto& sp: s_partition->spc().sp())
        //{
        //    receiver_.push_back(Receiver<DonorInfo2>(sp.tag()()));
        //}

        //for (const auto& se: sender_)
        //{
        //    for (const auto& re: receiver_)
        //    {
        //        assert(se.tag() != re.tag());
        //    }
        //}

        for (int i=0; i<a_partition->spc().sp().size(); ++i)
        {
            const auto& sp = a_partition->spc().sp()[i];

            for (const Mesh& mesh: sp.mesh())
            {
                for (const MeshCell& mc: mesh.cell())
                {
                    const auto& cnt = mc.poly().centroid();
                    //if (sp.is_resident(cnt))
                    int celltag = -1;
                    Tag brmt;
                    if (a_partition->spc().is_resident(cnt, celltag, brmt))
                    {
                        std::pair<Tag, int> resbin;
                        auto map = s_partition->spc().search_rm(mc.poly(), mc.poly().centroid(), resbin, false);

                        assert(!map.empty());

                        Tag donor;
                        int donor_rank = -1;
                        if (mc.oga_cell_type() == OGA_cell_type_t::receptor || mc.oga_cell_type() == OGA_cell_type_t::mandat_receptor)
                        {
                            assert(mc.donor().addr_ != nullptr);
                            bool resu = s_partition->spc().search_rm(mc.donor().addr_->poly().centroid(), donor, donor_rank);
                            assert(resu);
                            assert(donor.isvalid());
                        }

                        for (const auto& pair: map)
                        {
                            int dest_rank = pair.second;
                            assert(mc.oga_cell_type() != OGA_cell_type_t::undefined);
                            if (mc.oga_cell_type() == OGA_cell_type_t::receptor || mc.oga_cell_type() == OGA_cell_type_t::mandat_receptor)
                            {
                                assert(mc.donor().mesh_tag_.isvalid());
                            }
                            add(dest_rank, dest_rank, DonorInfo2(mesh.tag()(), mc.tag()(), mc.oga_cell_type(), mc.donor().mesh_tag_(), mc.donor().cell_tag_(), donor_rank, comm_->rank(), sp.tag()()));
                        }
                    }
                }
            }
        }
    }
}
