#include "loadmap.h"

namespace Tailor
{
    LoadEstim Loadmap::load_estim_type() const
    {
        return load_estim_type_;
    }

    double Loadmap::refine_tol() const
    {
        return refine_tol_;
    }

    /*double Loadmap::final_dev() const
    {
        assert(master());
        return final_dev_;
    }*/

    const std::vector<std::vector<BinRMTag>>& Loadmap::aug_bin_tag() const
    {
        return aug_bin_tag_;
    }

    //const std::map<BinRMTag, int> Loadmap::bintag_proc_map() const
    const std::map<int, int> Loadmap::bintag_proc_map() const
    {
        return bintag_proc_map_;
    }

    /*const size_t Loadmap::nmesh() const
    {
        return mesh_->size();
    }*/

    const RegularMesh& Loadmap::rm() const
    {
        assert(rm_ != nullptr);
        return *rm_;
    }

    /*const Mesh& Loadmap::mesh(const Tag& tag) const
    {
        auto it = mesh_tag_index_map.left.find(tag());

        assert(it != mesh_tag_index_map.left.end());

        return mesh_[it->second];
    }*/

    const std::deque<Mesh>* Loadmap::mesh() const
    {
        return mesh_;
    }
}
