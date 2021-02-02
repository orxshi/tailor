#include "sp.h"

namespace Tailor
{
    const std::vector<RegularMesh>& SpatialPartition::rm() const
    {
        return rm_;
    }

    const Tag& SpatialPartition::tag() const
    {
        return tag_;
    }

    void SpatialPartition::set_tag(const Tag& t)
    {
        tag_ = t;
    }

    /*const std::vector<MeshCell> SpatialPartition::arrival_cell() const
    {
        return arrival_cell_;
    }*/
    const AABB& SpatialPartition::aabb() const
    {
        return aabb_;
    }

    const RegularMesh& SpatialPartition::rm(const Tag& mt) const
    {
        assert(mt.isvalid());

        auto it = std::find_if(rm_.begin(), rm_.end(), [&](const RegularMesh& r){return r.tag() == mt;});
        assert(it != rm_.end());
        return *it;

        //auto it = mesh_tag_index_map.find(mt());
        //assert(it != mesh_tag_index_map.end());
        //return rm_[it->second];
    }

    //const ADT& SpatialPartition::adt(const Tag& mt) const
    //{
        //assert(mt.isvalid());

        //auto it = mesh_tag_index_map.left.find(mt());
        //assert(it != mesh_tag_index_map.left.end());
        //if(it->second >= adt_.size())
        //{
            //std::cout << "it->second = " << it->second << std::endl;
            //std::cout << "adt_.size() = " << adt_.size() << std::endl;
        //}
        //assert(it->second < adt_.size());
        //return adt_[it->second];
    //}

    //ADT& SpatialPartition::adt_p(const Tag& mt)
    //{
        //assert(mt.isvalid());

        //auto it = mesh_tag_index_map.left.find(mt());
        //assert(it != mesh_tag_index_map.left.end());
        //assert(it->second < adt_.size());
        //return adt_[it->second];
    //}

    RegularMesh& SpatialPartition::rm_p(const Tag& mt)
    {
        assert(mt.isvalid());

        auto it = std::find_if(rm_.begin(), rm_.end(), [&](const RegularMesh& r){return r.tag() == mt;});
        if (it == rm_.end())
        {
            std::cout << "mt = " << mt() << std::endl;
            std::cout << "rm size = " << rm_.size() << std::endl;
        }
        assert(it != rm_.end());
        return *it;

        //auto it = mesh_tag_index_map.find(mt());
        //assert(it != mesh_tag_index_map.end());
        //if (it->second >= rm_.size())
        //{
            //std::cout << "index = " << it->second << std::endl;
            //std::cout << "rm size = " << rm_.size() << std::endl;
            //std::cout << "mesh size = " << mesh_.size() << std::endl;
        //}
        //assert(it->second < rm_.size());
        //return rm_[it->second];
    }

    bool SpatialPartition::fully_resident(const MeshCell& mc) const
    {
        for (const MeshPoint& mp: mc.point())
        {
            if (!aabb_.do_intersect(mp.p().r()))
            {
                return false;
            }
        }

        return true;
    }


    //bool SpatialPartition::is_resident(const vec3<double>& cnt) const
    //{
    //    assert(!aabb_.faces().empty());
    //    return aabb_.do_intersect(cnt);
    //}

    //bool SpatialPartition::is_resident(const vec3<double>& cnt, const Outline& outline) const
    //{
    //    assert(false);
    //    //return outline.do_contain(cnt);
    //}

    Tag SpatialPartition::generate_meshtag() const
    {
        int i = TAILOR_BIG_NEG_NUM;
        Tag t;

        if (mesh_.empty())
        {
            t.set(0);
            return t;
        }
        for (const Mesh& m: mesh_)
        {
            i = std::max(i, m.tag()()); 
        }

        t.set(i+1);

        return t;
    }

    const std::deque<Mesh>& SpatialPartition::mesh() const
    {
        return mesh_;
    }

    std::deque<Mesh>& SpatialPartition::mesh_p()
    {
        return mesh_;
    }

    const Mesh& SpatialPartition::mesh(int i) const
    {
        assert(i >= 0);
        assert(i < mesh_.size());

        return mesh_[i];
    }

    const Mesh& SpatialPartition::mesh(const Tag& t) const
    {
        assert(t.isvalid());

        //auto it = mesh_tag_index_map.find(t());
        auto it = std::find_if(mesh_.begin(), mesh_.end(), [&](const Mesh& m){return m.tag() == t;});

        //if(it == mesh_tag_index_map.end())
        //{
            //std::cout << "mesh tag = " << t() << std::endl;
            //std::cout << "mesh size = " << mesh_.size() << std::endl;
        //}
        //assert(it != mesh_tag_index_map.end());
        if (it == mesh_.end())
        {
            std::cout << "mesh size: " << mesh_.size() << std::endl;
            std::cout << "t: " << t() << std::endl;
        }
        assert(it != mesh_.end());

        //return mesh_[it->second];
        return *it;
    }

    Mesh& SpatialPartition::mesh_p(const Tag& t)
    {
        assert(t.isvalid());

        //auto it = mesh_tag_index_map.find(t());
        auto it = std::find_if(mesh_.begin(), mesh_.end(), [&](const Mesh& m){return m.tag() == t;});

        //assert(it != mesh_tag_index_map.end());
        assert(it != mesh_.end());

        //return mesh_[it->second];
        return *it;
    }
}
