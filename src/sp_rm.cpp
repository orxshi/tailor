#include "sp.h"

namespace Tailor
{
    void SpatialPartition::make_regular_map_for_mesh_serial(const Tag& mt)
    {
        AABB mesh_aabb(mesh(mt).rawpoint());
        AABB sp_aabb = aabb_;
        sp_aabb.extend(mesh_aabb);

        rm_p(mt).set_props(mesh(mt), sp_aabb);
        //rm_p(mt).set_props(mesh(mt), AABB(mesh(mt).point()));
        rm_p(mt).insert_bins(0, 0, -1);
        assert(!mesh(mt).cell().empty());
        //assert(!rm_p(mt).bintag_index_map().left.empty());
        //std::cout << "registering mesh to map" << std::endl;
        rm_p(mt).register_mesh(mesh(mt), comm_->rank(), false);
    }
    
    void SpatialPartition::make_regular_map_for_mesh_uniproc(const Tag& mt, Profiler* profiler_)
    {
        //AABB mesh_aabb(mesh(mt).rawpoint());

        profiler_->start("mrm - aabb");
        vec3<double> minn;
        vec3<double> maxx;
        mesh(mt).bbox(minn, maxx);
        AABB mesh_aabb(minn, maxx);
        profiler_->stop("mrm - aabb");

        profiler_->start("mrm - props");
        rm_p(mt).set_props(mesh(mt), mesh_aabb);
        profiler_->stop("mrm - props");

        profiler_->start("mrm - insert");
        rm_p(mt).insert_bins(0, 0, -1);
        profiler_->stop("mrm - insert");

        for (const Bin& b: rm(mt).bin())
        {
            assert(b.rm() == nullptr);
        }

        assert(!mesh(mt).cell().empty());
        //assert(!rm_p(mt).bintag_index_map().left.empty());
        //profiler_->start("mrm - register");
        rm_p(mt).register_mesh(mesh(mt), comm_->rank(), false);
        //profiler_->stop("mrm - register");

        assert(!mesh(mt).cell().empty());

        profiler_->start("mrm - assert");
        bool all_empty = true;
        for (const Bin& b0: rm_p(mt).bin())
        {
            if (!b0.cell().empty())
            {
                all_empty = false;
                break;
            }
        }
        assert(!all_empty);
        profiler_->stop("mrm - assert");
    }

    void SpatialPartition::make_regular_map_for_mesh(const Tag& mt)
    {
        const Mesh& msh = mesh(mt);
        RegularMesh& rm = rm_p(mt);

        AABB mesh_aabb(msh.rawpoint());
        if (mesh_aabb.min(0) >= mesh_aabb.max(0))
        {
            std::cout << "meshaabb min 0: " << mesh_aabb.min(0)<< std::endl;
            std::cout << "meshaabb max 0: " << mesh_aabb.max(0)<< std::endl;
        }
        assert(mesh_aabb.min(0) < mesh_aabb.max(0));
        //AABB sp_aabb = aabb_;
        //sp_aabb.extend(mesh_aabb);

        rm.set_props(msh, mesh_aabb);
        rm.insert_bins(0, 0, -1);
        assert(!rm.bin().empty());
        assert(!msh.cell().empty());
        //assert(!rm.bintag_index_map().left.empty());
        rm.register_mesh(msh, comm_->rank(), false);

        assert(!msh.cell().empty());

        bool all_empty = true;
        for (const Bin& b0: rm.bin())
        {
            if (!b0.cell().empty())
            {
                all_empty = false;
                break;
            }
        }
        assert(!all_empty);

        assert(rm.aabb().min(0) < rm.aabb().max(0));
    }

    void SpatialPartition::make_regular_maps_for_mesh_serial()
    {
        {
            rm_.clear();
            rm_.resize(mesh_.size());
            for (int i=0; i<rm_.size(); ++i)
            {
                //rm_[i].set_tag(mesh_[i].tag());
                //rm_[i].insert_to_rmtag_address_map(rm_[i].tag(), &rm_[i]);
            }

            for (Mesh& m: mesh_)
            {
                make_regular_map_for_mesh_serial(m.tag());
            }
        }
    }

    void SpatialPartition::make_regular_maps_for_mesh_uniproc(Profiler* profiler_)
    {
        profiler_->start("mrm - push");
        rm_.clear();
        //rm_.resize(mesh_.size());
        rm_.reserve(mesh_.size());

        for (int i=0; i<mesh_.size(); ++i)
        {
            rm_.push_back(RegularMesh(mesh_[i].tag()));
        }
        profiler_->stop("mrm - push");

        for (Mesh& m: mesh_)
        {
            make_regular_map_for_mesh_uniproc(m.tag(), profiler_);
        }
    }

    void SpatialPartition::make_regular_maps_for_mesh()
    {
        rm_.clear();
        //rm_.resize(mesh_.size());
        rm_.reserve(mesh_.size());
        //for (int i=0; i<rm_.size(); ++i)
        for (int i=0; i<mesh_.size(); ++i)
        {
            //rm_[i] = RegularMesh(mesh_[i].tag());
            rm_.push_back(RegularMesh(mesh_[i].tag()));
            //rm_[i].set_tag(mesh_[i].tag());
            //rm_[i].insert_to_rmtag_address_map(rm_[i].tag(), &rm_[i]);
        }

        for (Mesh& m: mesh_) {
            make_regular_map_for_mesh(m.tag());
        }
        assert(mesh_.size() == rm_.size());
    }
}
