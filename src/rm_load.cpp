#include "regular_mesh.h"

namespace Tailor
{
    /*double Bin::load_solver(const std::deque<Mesh>& mesh)
    {
        double load_oga = 0;

        //if (mesh_load_.size() < 2) {
            //load_ = load_oga;
            //return load_oga;
        //}

        //for (double ml: mesh_load_) {
            //std::cout << "bin: " << tag_() << " ml: " << ml << std::endl;
            //load_oga += ml;
        //}

        for (auto z=mesh_tag_index_map_.left.begin(); z!=mesh_tag_index_map_.left.end(); ++z)
        {
            int mt = z->first;
            int ii = z->second;

            auto mm = std::find_if(mesh.begin(), mesh.end(), [&](const Mesh& m){return m.tag()() == mt;});
            assert(mm != mesh.end());

            //int nouter = mm->outer_boundaries().size();
            //int nwall = mm->wall_boundaries().size();
            //int npart = mm->npartition();
            //double nbou = nouter + nwall + npart;
            //double nbou = nouter + nwall;
            // cannot use npart without connectivity information.

            //load_oga += mesh_load_[ii] * nbou;
            load_oga += mesh_load_[ii];
        }

        load_ = load_oga;
        return load_oga;
    }

    double Bin::load_with_area(const std::deque<Mesh>& mesh)
    {
        assert(!mesh.empty());

        double load_oga = 0;

        if (cell_.empty()) {
            load_ = load_oga;
            return load_oga;
        }

        std::deque<AABB> aabbs = mesh_aabb(mesh);

        for (auto z=mesh_tag_index_map_.left.begin(); z!=mesh_tag_index_map_.left.end(); ++z)
        {
            int mt = z->first;
            int ii = z->second;

            auto mm = std::find_if(mesh.begin(), mesh.end(), [&](const Mesh& m){return m.tag()() == mt;});
            assert(mm != mesh.end());

            // how to incorporate nouter and nwall to load estimation?
            //int nouter = mm->outer_boundaries().size();
            //int nwall = mm->wall_boundaries().size();
            //int npart = mm->npartition();
            //double nbou = nouter + nwall + npart;
            //double nbou = nouter + nwall;

            int di = std::distance(mesh.begin(), mm);

            if (std::next(z) == mesh_tag_index_map_.left.end()) break;
            
            for (auto zz=std::next(z); zz!=mesh_tag_index_map_.left.end(); ++zz)
            {
                int mtt = zz->first;
                int jj = zz->second;

                auto mmm = std::find_if(mesh.begin(), mesh.end(), [&](const Mesh& m){return m.tag()() == mtt;});
                assert(mmm != mesh.end());

                int dj = std::distance(mesh.begin(), mmm);

                double oa = area_of_overlap(aabbs[di], aabbs[dj]);

                assert(oa >= 0.);
                assert((oa / aabbs[di].area()) <= 1.);
                assert((oa / aabbs[di].area()) >= 0.);
                assert((oa / aabbs[dj].area()) <= 1.);
                assert((oa / aabbs[dj].area()) >= 0.);

                double ncommon_i = (oa * mesh_load_[ii] / aabbs[di].area());
                double ncommon_j = (oa * mesh_load_[jj] / aabbs[dj].area());

                //load_oga += (ncommon_i + ncommon_j) * nbou;
                load_oga += (ncommon_i + ncommon_j);
                // so each potentially overlapped cell is subject to nbou boundary checks.
            }
        }

        load_ = load_oga;
        return load_oga;
    }
    
    double Bin::load(const std::deque<Mesh>& mesh, LoadEstimType load_estim_type)
    {
        if (load_estim_type == LoadEstimType::solver)
        {
            return load_solver(mesh);
        }
        else if (load_estim_type == LoadEstimType::area)
        {
            return load_with_area(mesh);
        }
        else
        {
            assert(false);
        }
    }*/
}
