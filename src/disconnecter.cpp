#include "disconnecter.h"

namespace Tailor
{
    Disconnecter::Disconnecter(SpatialPartitionContainer* spc, const Exchanger<Cell>* exc, Profiler* profiler): profiler_(profiler), exc_(exc), spc_(spc)
    {
        assert(spc_ != nullptr);
    }

    int Disconnecter::cellsize(int fn) const
    {
        int size = 0;
        for (int i=0; i<spc_->sp_.size(); ++i)
        {
            const auto& sp = spc_->sp_[i];
            for (const Mesh& mesh: sp.mesh())
            {
                if (mesh.tag()() == 0)
                {
                    size += mesh.cell().size();
                }
            }
        }

        return size;
    }

    void Disconnecter::remove_ghosts()
    {
        for (int i=0; i<spc_->sp_.size(); ++i)
        {
            auto& sp = spc_->sp_[i];
            sp.remove_ghosts();
        }
    }

    void Disconnecter::disconnect(bool mergebins, int rank)
    {
        if (profiler_ != nullptr) {profiler_->start("asm-dis-nover");}
        disconnect_nonoverlapping_cells(mergebins, rank);
        if (profiler_ != nullptr) {profiler_->stop("asm-dis-nover");}
        if (profiler_ != nullptr) {profiler_->start("asm-dis-empty");}
        disconnect_empty_meshes();
        if (profiler_ != nullptr) {profiler_->stop("asm-dis-empty");}
        if (profiler_ != nullptr) {profiler_->start("asm-dis-iso");}
        disconnect_isolated_points();
        if (profiler_ != nullptr) {profiler_->stop("asm-dis-iso");}
        if (profiler_ != nullptr) {profiler_->start("asm-dis-orphan");}
        disconnect_orphan_faces();
        if (profiler_ != nullptr) {profiler_->stop("asm-dis-orphan");}
    }

    void Disconnecter::disconnect_orphan_faces()
    {
        //for (SpatialPartition& sp: spc_->sp_)
        //{
        //    sp.disconnect_orphan_faces();
        //}
    }

    void Disconnecter::disconnect_neiless_cells(int rank)
    {
        for (int i=0; i<spc_->sp_.size(); ++i)
        {
            auto& sp = spc_->sp_[i];
            sp.reset_erase_marks();
        }

        //for (int i=0; i<spc_->sp_.size(); ++i)
        //{
        //    auto& sp = spc_->sp_[i];

        //    for (const auto& m: sp.mesh())
        //    {
        //        for (const auto& mc: m.cell())
        //        {
        //            for (const auto& mf: mc.face())
        //            {
        //                assert(mf.parent_cell().size() == m.face(mf.tag()).parent_cell().size());
        //            }
        //        }
        //    }
        //}

        for (int i=0; i<spc_->sp_.size(); ++i)
        {
            auto& sp = spc_->sp_[i];
            //auto& stor = (*storage_)[i];

            for (const auto& gr: exc_->group())
            {
                for (const auto& c: gr.data_)
                {
                    const MeshCell& mc = c.cell();

                    Tag ic = mc.tag();
                    Tag im = mc.parent_mesh();

                    const MeshCell* _cell = sp.mesh(im).cell_ptr(ic);
                    if (_cell == nullptr) {
                        continue;
                    }

                    if (_cell->pnei().empty())
                    {
                        std::cout << "icc: " << ic() << std::endl;
                        sp.mesh(im).print_as_vtk("empt.vtk");
                        assert(false);
                        sp.mark_to_be_erased(im, ic);
                    }
                }
            }
        }

        for (int i=0; i<spc_->sp_.size(); ++i)
        {
            auto& sp = spc_->sp_[i];
            sp.erase_marked_cells();
        }
    }

    void Disconnecter::disconnect_empty_meshes()
    {
        for (SpatialPartition& s: spc_->sp_)
        {
            for (int i=0; i<s.mesh_.size();)
            {
                if (s.mesh_[i].cell().empty())
                {
                    s.remove_mesh(s.mesh_[i].tag());
                }
                else {
                    ++i;
                }
            }
        }
    }

    void Disconnecter::disconnect_isolated_points()
    {
        for (SpatialPartition& s: spc_->sp_)
        {
            for (Mesh& m: s.mesh_) {
                m.remove_isolated_points();
            }
        }
    }

    void Disconnecter::disconnect_nonoverlapping_cells(bool mergebins, int rank)
    {
        /*
        for (int i=0; i<spc_->sp_.size(); ++i)
        {
            auto& sp = spc_->sp_[i];
            sp.reset_erase_marks();
        }

        int oversize = 0;
        for (int i=0; i<spc_->sp_.size(); ++i)
        {
            auto& sp = spc_->sp_[i];
            auto& stor = (*storage_)[i];

            for (const auto& gr: stor.group())
            {
                for (const auto& c: gr.data_)
                {
                    const MeshCell& mc = c.cell();

                    Tag ic = mc.tag();
                    Tag im = mc.parent_mesh();

                    const MeshCell* _cell = sp.mesh(im).cell_ptr(ic);
                    if (_cell == nullptr) {
                        continue;
                    }

                    for (const MeshPoint& p: _cell->point())
                    {
                        assert(p.parent_cell().size() > 0);
                    }

                    bool overlap = false;
                    if (mergebins)
                    {
                        //std::vector<BinRMTag> tag;
                        //int temp;
                        //spc_->global_rm().get_bintag_adaptive(AABB(mc.poly()), tag, temp, -1);
                        auto map = spc_->search_rm(mc.poly());

                        for (const auto& pair: map)
                        //for (const auto& brmt: tag)
                        {
                            //auto it = spc_->bintag_proc_map().find(brmt);
                            //int dest_rank = it->second;
                            auto brmt = pair.first;
                            int dest_rank = pair.second;

                            if (dest_rank != rank) {
                                continue;
                            }

                            const Bin& bin = spc_->global_rm().bin(brmt);
                            if (bin.aabb().do_intersect(AABB(mc.poly())))
                            {
                                overlap = true;
                                break;
                            }
                        }
                    }
                    else
                    {
                        overlap = sp.aabb().do_intersect(AABB(_cell->poly()));
                    }

                    if (!overlap)
                    {
                        ++oversize;
                        sp.mark_to_be_erased(im, ic);
                    }
                }
            }
        }
    */

        // but remove_if won't take advantage of sorted vector.

        for (int i=0; i<spc_->sp_.size(); ++i)
        {
            auto& sp = spc_->sp_[i];
            sp.erase_marked_cells();
        }
    }
}
