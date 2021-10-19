#include "donor_request_exchanger.h"

namespace Tailor
{
    DonorRequestExchanger::DonorRequestExchanger(const SpatialPartitionContainer* spc, SpatialPartition* sp, const ArrCon<DonorInfo2>* di, boost::mpi::communicator* comm): comm_(comm), Exchanger(comm), sp_(sp), di_(di), spc_(spc)
    {
    }

    void DonorRequestExchanger::prepare_storage()
    {
        auto sptag = sp_->tag()();
        //sender_.push_back(Sender<Request>(sptag));
        //receiver_.push_back(Receiver<Request>(sptag));

        //auto rece = std::find_if(di_->begin(), di_->end(), [&](const auto& re){return re.tag() == sptag;});
        //assert(rece != di_->end());
        //const auto& arr = rece->arrival();

        for (Mesh& mesh: sp_->mesh_)
        {
            //for (const MeshFace& mf: mesh.face())
            //{
            //    if (mf.parent_cell().size() == 1)
            //    {
            //        const auto& mc = mesh.cell(Tag(mf.parent_cell()[0]));
            //        if (mc.oga_cell_type() != OGA_cell_type_t::non_resident && mc.oga_cell_type() != OGA_cell_type_t::ghost)
            //        {
            //            std::cout << comm_->rank() << " mesh: " << mesh.tag()() << std::endl;
            //            std::cout << comm_->rank() << " mc: " << mc.tag()() << std::endl;
            //            std::cout << comm_->rank() << " oga type: " << static_cast<int>(mc.oga_cell_type()) << std::endl;
            //            std::cout << comm_->rank() << " mf btype: " << static_cast<int>(mf.btype()) << std::endl;
            //            std::cout << comm_->rank() << " mc btype: " << static_cast<int>(mc.btype()) << std::endl;
            //            std::cout << comm_->rank() << " mf: " << mf.tag() << std::endl;
            //        }
            //        //assert(mesh.cell(Tag(mf.parent_cell()[0])).oga_cell_type() == OGA_cell_type_t::non_resident || mesh.cell(Tag(mf.parent_cell()[0])).oga_cell_type() == OGA_cell_type_t::ghost);
            //        assert(mc.oga_cell_type() == OGA_cell_type_t::non_resident || mc.oga_cell_type() == OGA_cell_type_t::ghost);
            //        assert(!mf.is_boundary());
            //    }
            //}
            int meshtag = mesh.tag()();

            for (MeshCell& mc: mesh.cell_)
            {
                auto mctag = mc.tag()();
                auto ogatype = mc.oga_cell_type();

                if (ogatype == OGA_cell_type_t::non_resident || ogatype == OGA_cell_type_t::ghost)
                {
                    continue;
                }

                //if (rece != di_->end())
                {
                    assert(ogatype == OGA_cell_type_t::undefined);
                    //if (ogatype != OGA_cell_type_t::undefined)
                    //{
                    //    if (ogatype != OGA_cell_type_t::non_resident || ogatype != OGA_cell_type_t::ghost)
                    //    {
                    //        std::cout << "ogatype: " << static_cast<int>(ogatype) << std::endl;
                    //    }
                    //    assert(ogatype == OGA_cell_type_t::non_resident || ogatype == OGA_cell_type_t::ghost);
                    //    continue;
                    //}
                }

                //if (rece == di_->end())
                //{
                //    if (ogatype != OGA_cell_type_t::receptor && ogatype != OGA_cell_type_t::mandat_receptor) {
                //        continue;
                //    }

                //    int donor_rank = -1;
                //    Tag donor;
                //    assert(mc.donor().addr_ != nullptr);
                //    bool resu = spc_->search_rm(mc.donor().addr_->poly().centroid(), donor, donor_rank);
                //    assert(resu);
                //    assert(donor.isvalid());

                //    //sender_.front().add(donor_rank, donor_rank, Request(mc.donor().mesh_tag_(), mc.donor().cell_tag_(), comm_->rank(), sptag));
                //    add(donor_rank, donor_rank, Request(mc.donor().mesh_tag_(), mc.donor().cell_tag_(), comm_->rank(), sptag));
                //}
                //else
                {
                    //const auto& arr = rece->arrival();

                    auto mc_iter = std::find_if(di_->begin(), di_->end(), [&](const DonorInfo2& d){return ((d.meshtag_ == meshtag) && (d.celltag_ == mctag));});
                    if (mc_iter == di_->end())
                    {
                        std::cout << "comm_->rank(): " << comm_->rank() << std::endl;
                        std::cout << "meshtag: " << meshtag << std::endl;
                        std::cout << "mctag: " << mctag << std::endl;
                        //std::cout << "arr.size(): " << arr.size() << std::endl;
                        std::cout << "mesh.cell().size(): " << mesh.cell().size() << std::endl;
                        std::cout << "ogatype: " << static_cast<int>(ogatype) << std::endl;
                        mesh.print_as_vtk("weird.vtk");
                    }
                    assert(mc_iter != di_->end());

                    assert(mc_iter->oga_cell_type_ != OGA_cell_type_t::undefined);
                    mc.set_oga_cell_type(mc_iter->oga_cell_type_);
                    ogatype = mc.oga_cell_type();

                    if (ogatype == OGA_cell_type_t::receptor || ogatype == OGA_cell_type_t::mandat_receptor)
                    {
                        mc.set_donor(Tag(mc_iter->donor_mesh_), Tag(mc_iter->donor_cell_), nullptr);
                        // TODO set receptor.
                    }

                    if (ogatype != OGA_cell_type_t::receptor && ogatype != OGA_cell_type_t::mandat_receptor) {
                        continue;
                    }
                    
                    //sender_.front().add(mc_iter->donor_rank_, mc_iter->donor_rank_, Request(mc_iter->donor_mesh_, mc_iter->donor_cell_, comm_->rank(), sptag));
                    add(mc_iter->donor_rank_, mc_iter->donor_rank_, Request(mc_iter->donor_mesh_, mc_iter->donor_cell_, comm_->rank(), sptag));
                }
            }
        }
    }
}
