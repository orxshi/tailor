#include "cell_exchanger.h"

namespace Tailor
{
    CellExchanger::CellExchanger(SpatialPartitionContainer* spc, boost::mpi::communicator* comm, bool merge_bins): comm_(comm), Exchanger(comm), spc_(spc), merge_bins_(merge_bins)
    {
    }

    void CellExchanger::prepare_storage()
    {
        {
            std::string mems = "bef-cell-exc-prepare";
            mem_usage(comm_, mems);
        }
        //std::ofstream out;
        //std::string fn = "ce-";
        //fn.append(std::to_string(comm_->rank()));
        //fn.append(".dat");
        //out.open(fn);
        //int na = 0;
        //int nb = 0;
        //int nsearch = 0;
        //int nconstruct = 0;

        for (int i=0; i<spc_->sp_.size(); ++i)
        {
            auto& sp = spc_->sp_[i];
            sp.reset_erase_marks();
        }

        //sender_.reserve(spc_->sp().size());
        //receiver_.reserve(spc_->sp().size());

        for (int i=0; i<spc_->sp().size(); ++i)
        {
            SpatialPartition& sp = spc_->sp_[i];

            //sender_.push_back(Sender<Cell>(sp.tag()()));
            //receiver_.push_back(Receiver<Cell>(sp.tag()()));

            for (Mesh& m: sp.mesh_p())
            {
                for (MeshCell& mc: m.cell_p())
                {
                    //++na;

                    if (mc.oga_cell_type() == OGA_cell_type_t::ghost) {
                        continue;
                    }

                    AABB aabb = AABB(mc.poly());

                    //++nb;

                    //if (!merge_bins_)
                    //{
                    //    if (sp.aabb().do_intersect(mc.poly().centroid()))
                    //    {
                    //        if (mc.oga_cell_type() == OGA_cell_type_t::non_resident) {
                    //            mc.set_oga_cell_type(OGA_cell_type_t::undefined);
                    //        }
                    //    }
                    //    else
                    //    {
                    //        mc.set_oga_cell_type(OGA_cell_type_t::non_resident);
                    //    }
                    //}

                    if (!merge_bins_)
                    {
                        if (sp.aabb().do_contain(aabb, true)) {
                            continue;
                        }
                    }

                    //std::vector<BinRMTag> tag;
                    //int temp;
                    //spc_->global_rm_.get_bintag_adaptive(AABB(mc.poly()), tag, temp, mc.tag()());
                    bool overlap = false;
                    //auto map = spc_->search_rm(aabb); 
                    std::pair<Tag, int> resbin;
                    auto map = spc_->search_rm(aabb, mc.poly().centroid(), resbin, false);
                    //++nsearch;

                    //if (merge_bins_)
                    //{
                    //    if (resbin.second != comm_->rank())
                    //    {
                    //        mc.set_oga_cell_type(OGA_cell_type_t::non_resident);
                    //    }
                    //    else
                    //    {
                    //        mc.set_oga_cell_type(OGA_cell_type_t::undefined);
                    //    }
                    //}


                    //for (const auto& brmt: tag)
                    std::unique_ptr<Cell> cll;
                    std::list<MeshCell*> wallp, dirichletp, farfieldp, emptyp, interogp, symmetryp;

                    for (auto& pair: map)
                    {
                        const auto& brmt = pair.first;
                        int dest_rank = pair.second;
                        //const Bin& bin = spc_->global_rm().bin(brmt);

                        //auto it = spc_->bintag_proc_map().find(brmt);
                        //int dest_rank = it->second;

                        //assert(bin.aabb().do_intersect(aabb));
                        //if (bin.aabb().do_intersect(aabb))
                        int desttag;
                        if (merge_bins_)
                        {
                            assert(sp.tag()() == comm_->rank());
                            if (dest_rank == comm_->rank())
                            {
                                overlap = true;

                                if (map.size() == 1) {
                                    break;
                                }
                            }
                            desttag = dest_rank;
                        }
                        else
                        {
                            if (sp.tag() == brmt)
                            {
                                overlap = true;

                                if (map.size() == 1) {
                                    break;
                                }
                            }
                            desttag = brmt();
                        }

                        if (cll == nullptr)
                        {
                            std::list<MeshCell> wall, dirichlet, farfield, empty, interog, symmetry;

                            //for (const auto& bc: mc.wall_boundary())
                            for (auto bc = mc.wall_boundary().begin(); bc != mc.wall_boundary().end(); ++bc)
                            {
                                auto bb = m.wall_boundary_p(*bc);
                                //if (!overlap) {
                                //bb->mark_to_be_erased();
                                //}
                                assert(bb != nullptr);
                                wall.push_back(*bb);
                                wallp.push_back(bb);
                            }
                            //for (const auto& bc: mc.symmetry_boundary())
                            for (auto bc = mc.symmetry_boundary().begin(); bc != mc.symmetry_boundary().end(); ++bc)
                            {
                                auto bb = m.symmetry_boundary_p(*bc);
                                //if (!overlap) {
                                //bb->mark_to_be_erased();
                                //}
                                assert(bb != nullptr);
                                symmetry.push_back(*bb);
                                symmetryp.push_back(bb);
                            }
                            //for (const auto& bc: mc.dirichlet_boundary())
                            for (auto bc = mc.dirichlet_boundary().begin(); bc != mc.dirichlet_boundary().end(); ++bc)
                            {
                                auto bb = m.dirichlet_boundary_p(*bc);
                                //if (!overlap) {
                                //bb->mark_to_be_erased();
                                //}
                                assert(bb != nullptr);
                                dirichlet.push_back(*bb);
                                dirichletp.push_back(bb);
                            }
                            //for (const auto& bc: mc.farfield_boundary()) {
                            for (auto bc = mc.farfield_boundary().begin(); bc != mc.farfield_boundary().end(); ++bc)
                            {
                                auto bb = m.farfield_boundary_p(*bc);
                                //if (!overlap) {
                                //bb->mark_to_be_erased();
                                //}
                                assert(bb != nullptr);
                                farfield.push_back(*bb);
                                farfieldp.push_back(bb);
                            }
                            //for (const auto& bc: mc.empty_boundary()) {
                            for (auto bc = mc.empty_boundary().begin(); bc != mc.empty_boundary().end(); ++bc)
                            {
                                auto bb = m.empty_boundary_p(*bc);
                                //if (!overlap) {
                                //bb->mark_to_be_erased();
                                //}
                                assert(bb != nullptr);
                                empty.push_back(*bb);
                                emptyp.push_back(bb);
                            }
                            //for (const Tag& bc: mc.interog_boundary()) {
                            for (auto bc = mc.interog_boundary().begin(); bc != mc.interog_boundary().end(); ++bc)
                            {
                                auto bb = m.interog_boundary_p(*bc);
                                //if (!overlap) {
                                //bb->mark_to_be_erased();
                                //}
                                if (bb == nullptr)
                                {
                                    std::cout << "m.tag: " << m.tag()() << std::endl;
                                    std::cout << "mc.tag: " << mc.tag()() << std::endl;
                                    std::cout << "m interog size: " << m.interog_boundaries().size() << std::endl;
                                    std::cout << "mc interog size: " << mc.interog_boundary().size() << std::endl;
                                    std::cout << "oga: " << static_cast<int>(mc.oga_cell_type()) << std::endl;
                                    std::cout << "btype: " << static_cast<int>(mc.btype()) << std::endl;
                                }
                                assert(bb != nullptr);
                                interog.push_back(*bb);
                                interogp.push_back(bb);
                            }

                            //cll = new Cell(comm_->rank(), sp.tag()(), mc, std::move(wall), std::move(dirichlet), std::move(farfield), std::move(empty), std::move(interog));
                            cll = std::unique_ptr<Cell>(new Cell(comm_->rank(), sp.tag()(), mc, std::move(wall), std::move(dirichlet), std::move(farfield), std::move(empty), std::move(interog), std::move(symmetry)));
                            //if (!overlap) {
                            //mc.mark_to_be_erased();
                            //}
                        }


                        if (merge_bins_)
                        {
                            if (dest_rank == comm_->rank())
                            {
                                continue;
                            }
                        }
                        else
                        {
                            if (sp.tag() == brmt)
                            {
                                continue;
                            }
                        }

                        //int desttag;
                        //if (merge_bins_) {
                        //desttag = dest_rank;
                        //}
                        //else {
                        ////desttag = bin.tag()();
                        //desttag = brmt();
                        //}

                        //sender_[i].add(dest_rank, desttag, Cell(comm_->rank(), sp.tag()(), mc, std::move(wall), std::move(dirichlet), std::move(farfield), std::move(empty), std::move(interog)));
                        assert(cll != nullptr);
                        //sender_[i].add(dest_rank, desttag, *cll);
                        //std::cout << "merge_bins_: " << merge_bins_ << std::endl;
                        //std::cout << "dest_rank: " << dest_rank << std::endl;
                        //std::cout << "comm rank: " << comm_->rank() << std::endl;
                        add(dest_rank, desttag, *cll);
                        //++nconstruct;

                    }


                    if (!overlap)
                    {
                        //sp.mark_to_be_erased(m.tag()(), mc.tag()());
                        mc.mark_to_be_erased();
                        for (auto& ptr: wallp) {
                            ptr->mark_to_be_erased();
                        }
                        for (auto& ptr: symmetryp) {
                            ptr->mark_to_be_erased();
                        }
                        for (auto& ptr: dirichletp) {
                            ptr->mark_to_be_erased();
                        }
                        for (auto& ptr: farfieldp) {
                            ptr->mark_to_be_erased();
                        }
                        for (auto& ptr: emptyp) {
                            ptr->mark_to_be_erased();
                        }
                        for (auto& ptr: interogp) {
                            ptr->mark_to_be_erased();
                        }
                    }

                    //delete cll;
                    //delete dirichletp;
                    //delete farfieldp;
                    //delete emptyp;
                    //delete interogp;
                }
            }

            //out << "nsearch: " << nsearch << std::endl;
            //out << "nconstruct: " << nconstruct << std::endl;
            //out << "na: " << na << std::endl;
            //out << "nb: " << nb << std::endl;
            //out.close();
        }

        {
            std::string mems = "aft-cell-exc-prepare";
            mem_usage(comm_, mems);
        }
    }

    Cell::Cell(int source_rank, int source_tag, const MeshCell& cell, std::list<MeshCell>&& wall, std::list<MeshCell>&& dirichlet, std::list<MeshCell>&& farfield, std::list<MeshCell>&& empty, std::list<MeshCell>&& interog, std::list<MeshCell>&& symmetry): source_rank_(source_rank), source_tag_(source_tag), cell_(cell), wall_(wall), dirichlet_(dirichlet), farfield_(farfield), empty_(empty), interog_(interog), symmetry_(symmetry)
    {
        if (!cell_.farfield_boundary().empty())
        {
            assert(!farfield_.empty());
        }
    }

    const MeshCell& Cell::cell() const
    {
        return cell_;
    }
}
