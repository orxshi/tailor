#include "remap.h"

namespace Tailor
{
    Remapper::Remapper(SpatialPartitionContainer* spc, boost::mpi::communicator* comm, std::string name): spc_(spc), comm_(comm), name_(name)
    {
    }

    void Remapper::remap(Loadmap& lm, bool mergebins, int iter, std::string s, Profiler* profiler_)
    {
        auto s1 = s;
        s1.append("-sp-exc");
        //if (profiler_ != nullptr) {profiler_->bstart(s1);}
        SpExchanger sp_exc(&lm, comm_);
            {
                std::string mems = "bef-sp-ex";
                mem_usage(comm_, mems);
            }
        sp_exc.exchange(false, s1, profiler_, -1);
        //sp_exc.exchange();
        //if (profiler_ != nullptr) {profiler_->stop(s1);}

            {
                std::string mems = "aft-sp-ex";
                mem_usage(comm_, mems);
            }

        auto s2 = s;
        s2.append("-reduce");
        if (profiler_ != nullptr) {profiler_->start(s2);}
        if (mergebins) {
            assert(spc_->sp().empty());
            spc_->reduce_sp(lm.aug_bin_tag(), sp_exc.arrival(), s);
        }
        else {
            spc_->add_sp(lm.aug_bin_tag(), sp_exc.arrival());
        }
        if (profiler_ != nullptr) {profiler_->stop(s2);}

            {
                std::string mems = "aft-reduce";
                mem_usage(comm_, mems);
            }

        auto s3 = s;
        s3.append("-sort-sp-meshes");
        if (profiler_ != nullptr) {profiler_->start(s3);}
        spc_->sort_sp_meshes();
        if (profiler_ != nullptr) {profiler_->stop(s3);}

            {
                std::string mems = "aft-sort";
                mem_usage(comm_, mems);
            }

        for (const auto& sp: spc_->sp())
        {
            assert(!sp.aabb().faces().empty());
        }

        auto s4 = s;
        s4.append("-sp-misc");
        if (profiler_ != nullptr) {profiler_->start(s4);}
        lm.remove_cells_from_rm();
        spc_->set_global_rm(lm.rm());
        spc_->make_mesh_aabbs();
        assert(spc_->global_rm().tag()() == 0);
        spc_->make_global_adt();
        if (profiler_ != nullptr) {profiler_->stop(s4);}
    }
}
