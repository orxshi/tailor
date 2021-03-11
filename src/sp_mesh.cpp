#include "sp.h"

namespace Tailor
{
    void SpatialPartition::rotate_mesh(const Tag& _parent_mesh, double angle, int axis, const Vector3& rot_point)
    {
        //assert(mesh_.size() == adt_.size());

            assert(!std::isnan(rot_point(0)));
            assert(!std::isnan(rot_point(1)));
            assert(!std::isnan(rot_point(2)));

        for (int i=0; i<mesh_.size(); ++i)
        {
            //assert(mesh_[i].parent_mesh().isvalid());

            //if (mesh_[i].parent_mesh() == _parent_mesh)
            if (mesh_[i].tag() == _parent_mesh)
            {
                mesh_[i].rotate(angle, axis, rot_point);
                //rm_[i].rotate(angle, axis, rot_point);
                //adt_[i].rotate(ang, rot_axis);
            }
        }
    }

    void SpatialPartition::move_mesh(const Tag& _parent_mesh, const Vector3& v)
    {
        //assert(mesh_.size() == adt_.size());

        for (int i=0; i<mesh_.size(); ++i)
        {
            //assert(mesh_[i].parent_mesh().isvalid());

            //if (mesh_[i].parent_mesh() == _parent_mesh)
            if (mesh_[i].tag() == _parent_mesh)
            {
                mesh_[i].move(v);
                //rm_[i].move(v);
                //adt_[i].move(v);
            }
        }
    }

    void SpatialPartition::add_mesh(Mesh&& m)
    {
        add_mesh(m);
    }

    void SpatialPartition::add_mesh(const Mesh& m)
    {
        assert(m.tag().isvalid());

        //std::cout << "mesh size: " << m.mem() / 1e6 << std::endl;
        mesh_.push_back(m);
        //std::cout << "pushed mesh " << std::endl;
        //assert(mesh_tag_index_map.count(m.tag()()) == 0);
        //mesh_tag_index_map.insert(std::pair<int, int>(m.tag()(), mesh_.size() - 1));
    }

    void SpatialPartition::remove_mesh(Tag mt)
    {
        assert(mt.isvalid());
        //auto it = mesh_tag_index_map.find(mt());
        //assert(it != mesh_tag_index_map.end());
        //int thres = it->second;
        //mesh_.erase(mesh_.begin() + it->second);
        //
        auto it = std::find_if(mesh_.begin(), mesh_.end(), [&](const Mesh& m){return m.tag() == mt;});
        assert(it != mesh_.end());
        mesh_.erase(it);

        /*mesh_tag_index_map.erase(mt());
        assert(mesh_tag_index_map.count(thres) == 0);
        for (int j=thres+1; j<=mesh_.size(); ++j)
        {
            auto itt = mesh_tag_index_map.right.find(j);
            assert(itt != mesh_tag_index_map.right.end());
            int rep = itt->first - 1;
            bool successful_replace = mesh_tag_index_map.right.replace_key(itt, rep);
            assert(successful_replace);
            assert(mesh_[rep].tag()() == itt->second);
        }
        assert(mesh_tag_index_map.left.count(mt()) == 0);*/
    }
}
