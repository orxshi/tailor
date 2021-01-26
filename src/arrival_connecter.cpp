#include "arrival_connecter.h"

namespace Tailor
{
    ArrivalConnecter::ArrivalConnecter(SpatialPartitionContainer* spc, ArrCon<Cell>* arrival, Profiler* profiler): profiler_(profiler), arrival_(arrival), spc_(spc)
    {
    }

    void ArrivalConnecter::add_arrival_cells()
    {
        for (int i=0; i<spc_->sp_.size(); ++i)
        {
            if (profiler_ != nullptr) {profiler_->start("asm-arr-pre");}
            auto& sp = spc_->sp_[i];
            //auto& stor = (*storage_)[i];
            
            erase_interior_.clear();
            //erase_interior_.resize(stor.arrival().size(), false);
            erase_interior_.resize(arrival_->size(), false);

            //auto& arrival = stor.arrival_p();
            for (Cell& c: *arrival_)
            {
                MeshCell& ac = c.cell_;
                ac.reset_oga_status();
            }
            if (profiler_ != nullptr) {profiler_->stop("asm-arr-pre");}

            if (profiler_ != nullptr) {profiler_->start("asm-arr-det-cell");}
            determine_added_cells(sp, *arrival_);
            if (profiler_ != nullptr) {profiler_->stop("asm-arr-det-cell");}

            if (profiler_ != nullptr) {profiler_->start("asm-arr-add-cell");}
            add_arrival_cells(sp);
            if (profiler_ != nullptr) {profiler_->stop("asm-arr-add-cell");}

            if (profiler_ != nullptr) {profiler_->start("asm-arr-det-wall");}
            determine_added_bou(sp, *arrival_, BouType::wall);
            if (profiler_ != nullptr) {profiler_->stop("asm-arr-det-wall");}
            
            if (profiler_ != nullptr) {profiler_->start("asm-arr-add-wall");}
            add_arrival_bou(sp, BouType::wall);
            if (profiler_ != nullptr) {profiler_->stop("asm-arr-add-wall");}

            if (profiler_ != nullptr) {profiler_->start("asm-arr-det-dir");}
            determine_added_bou(sp, *arrival_, BouType::dirichlet);
            if (profiler_ != nullptr) {profiler_->stop("asm-arr-det-dir");}

            if (profiler_ != nullptr) {profiler_->start("asm-arr-add-dir");}
            add_arrival_bou(sp, BouType::dirichlet);
            if (profiler_ != nullptr) {profiler_->stop("asm-arr-add-dir");}

            if (profiler_ != nullptr) {profiler_->start("asm-arr-det-far");}
            determine_added_bou(sp, *arrival_, BouType::farfield);
            if (profiler_ != nullptr) {profiler_->stop("asm-arr-det-far");}

            if (profiler_ != nullptr) {profiler_->start("asm-arr-add-far");}
            add_arrival_bou(sp, BouType::farfield);
            if (profiler_ != nullptr) {profiler_->stop("asm-arr-add-far");}

            if (profiler_ != nullptr) {profiler_->start("asm-arr-det-empty");}
            determine_added_bou(sp, *arrival_, BouType::empty);
            if (profiler_ != nullptr) {profiler_->stop("asm-arr-det-empty");}

            if (profiler_ != nullptr) {profiler_->start("asm-arr-add-empty");}
            add_arrival_bou(sp, BouType::empty);
            if (profiler_ != nullptr) {profiler_->stop("asm-arr-add-empty");}

            if (profiler_ != nullptr) {profiler_->start("asm-arr-det-interog");}
            determine_added_bou(sp, *arrival_, BouType::interog);
            if (profiler_ != nullptr) {profiler_->stop("asm-arr-det-interog");}

            if (profiler_ != nullptr) {profiler_->start("asm-arr-add-interog");}
            add_arrival_bou(sp, BouType::interog);
            if (profiler_ != nullptr) {profiler_->stop("asm-arr-add-interog");}
        }
    }

    void ArrivalConnecter::determine_added_cells(SpatialPartition& sp, ArrCon<Cell>& arrival)
    {
        //auto& arrival = storage.arrival_p();

        //for (Cell& c: arrival)
        //{
        //    MeshCell& ac = c.cell_;

        //    for (auto& mf: ac.face_p())
        //    {
        //        mf.set_faceaddr(nullptr);
        //    }
        //}

        for (const Cell& c: arrival)
        {
            const MeshCell& ac = c.cell();
            //assert(ac.oga_cell_type() != OGA_cell_type_t::ghost);
            //if (ac.oga_cell_type() == OGA_cell_type_t::non_resident) {
                //continue;
            //}

            auto m = std::find_if(sp.mesh_.begin(), sp.mesh_.end(), [&](const Mesh& _m){return _m.tag() == ac.parent_mesh();});
            if (m != sp.mesh_.end())
            {
                int i = &c - &arrival.front();
                check_cell_exist2(*m, ac, i);
            }
        }

        auto it = std::remove_if(arrival.begin(), arrival.end(), [&](Cell& c){return erase_interior_[&c - &arrival.front()] == true;});
        arrival.erase(it, arrival.end());

        /*std::vector<int> meshtags;
        for (const auto& ar: arrival)
        {
            assert(ar.cell().parent_mesh()() == 0 || ar.cell().parent_mesh()() == 1);
            int c = std::count(meshtags.begin(), meshtags.end(), ar.cell().parent_mesh()());
            assert(c == 0 || c == 1);
            if (c == 0)
            {
                meshtags.push_back(ar.cell().parent_mesh()());
                if (meshtags.size() > 2)
                {
                    for (int ii: meshtags)
                    {
                        std::cout << "meshtags: " << ii << std::endl;
                    }
                }
                assert(meshtags.size() <= 2);
            }
        }*/

        arr_interior_.clear();
        arr_interior_.resize(sp.mesh_.size());
        //arr_interior_.resize(meshtags.size());

        for (const Cell& c: arrival)
        {
            const MeshCell& mc = c.cell();

            //if (mc.oga_cell_type() == OGA_cell_type_t::non_resident) {
                //continue;
            //}

            auto m = std::find_if(sp.mesh_.begin(), sp.mesh_.end(), [&](const Mesh& _m){return _m.tag() == mc.parent_mesh();});
            if (m == sp.mesh_.end())
            {
                Mesh amesh;
                amesh.set_tag(mc.parent_mesh());
                //amesh.set_parent_mesh(mc.parent_mesh());
                sp.add_mesh(amesh);
                auto temp = {std::move(mc)};
                arr_interior_.push_back(temp);
            }
            else {
                assert(std::distance(sp.mesh_.begin(), m) < arr_interior_.size());
                //arr_interior_[std::distance(sp.mesh_.begin(), m)].push_back(std::move(mc));
                arr_interior_[std::distance(sp.mesh_.begin(), m)].push_back(std::move(mc));
            }
        }
    }

    void ArrivalConnecter::determine_added_bou(SpatialPartition& sp, ArrCon<Cell>& arrival, BouType bou)
    {
        //auto& arrival = storage.arrival_p();

        std::vector<mcc>* container;
        if (bou == BouType::wall) {
            container = &arr_wall_;
        }
        else if (bou == BouType::dirichlet) {
            container = &arr_dirichlet_;
        }
        else if (bou == BouType::farfield) {
            container = &arr_farfield_;
        }
        else if (bou == BouType::empty) {
            container = &arr_empty_;
        }
        else if (bou == BouType::interog) {
            container = &arr_interog_;
        }
        else
        {
            assert(false);
        }

        container->clear();
        container->resize(sp.mesh_.size());

        if (bou == BouType::wall)
        {
            for (Cell& c: arrival)
            {
                for (MeshCell& mc: c.wall_)
                {
                    auto m = std::find_if(sp.mesh_.begin(), sp.mesh_.end(), [&](const Mesh& _m){return _m.tag() == mc.parent_mesh();});
                    assert(m != sp.mesh_.end());
                    if (m->query_bou(mc.tag(), bou) == nullptr) {
                        //for (auto& mf: mc.face_p())
                        //{
                        //    mf.set_faceaddr(nullptr);
                        //}
                        //(*container)[std::distance(sp.mesh_.begin(), m)].push_back(mc);
                        (*container)[std::distance(sp.mesh_.begin(), m)].push_back(mc);
                    }
                }
            }
        }
        else if (bou == BouType::dirichlet)
        {
            for (Cell& c: arrival)
            {
                for (MeshCell& mc: c.dirichlet_)
                {
                    //auto m = std::find_if(sp.mesh_.begin(), sp.mesh_.end(), [&](const Mesh& _m){return _m.parent_mesh() == mc.parent_mesh();});
                    auto m = std::find_if(sp.mesh_.begin(), sp.mesh_.end(), [&](const Mesh& _m){return _m.tag() == mc.parent_mesh();});
                    assert(m != sp.mesh_.end());
                    if (m->query_bou(mc.tag(), bou) == nullptr) {
                        //for (auto& mf: mc.face_p())
                        //{
                        //    mf.set_faceaddr(nullptr);
                        //}
                        //(*container)[std::distance(sp.mesh_.begin(), m)].push_back(mc);
                        (*container)[std::distance(sp.mesh_.begin(), m)].push_back(mc);
                    }
                }
            }
        }
        else if (bou == BouType::farfield)
        {
            for (Cell& c: arrival)
            {
                for (MeshCell& mc: c.farfield_)
                {
                    //auto m = std::find_if(sp.mesh_.begin(), sp.mesh_.end(), [&](const Mesh& _m){return _m.parent_mesh() == mc.parent_mesh();});
                    auto m = std::find_if(sp.mesh_.begin(), sp.mesh_.end(), [&](const Mesh& _m){return _m.tag() == mc.parent_mesh();});
                    assert(m != sp.mesh_.end());
                    if (m->query_bou(mc.tag(), bou) == nullptr) {
                        //for (auto& mf: mc.face_p())
                        //{
                        //    mf.set_faceaddr(nullptr);
                        //}
                        //(*container)[std::distance(sp.mesh_.begin(), m)].push_back(mc);
                        (*container)[std::distance(sp.mesh_.begin(), m)].push_back(mc);
                    }
                }
            }
        }
        else if (bou == BouType::empty)
        {
            for (Cell& c: arrival)
            {
                for (MeshCell& mc: c.empty_)
                {
                    //auto m = std::find_if(sp.mesh_.begin(), sp.mesh_.end(), [&](const Mesh& _m){return _m.parent_mesh() == mc.parent_mesh();});
                    auto m = std::find_if(sp.mesh_.begin(), sp.mesh_.end(), [&](const Mesh& _m){return _m.tag() == mc.parent_mesh();});
                    assert(m != sp.mesh_.end());
                    if (m->query_bou(mc.tag(), bou) == nullptr) {
                        //for (auto& mf: mc.face_p())
                        //{
                        //    mf.set_faceaddr(nullptr);
                        //}
                        //(*container)[std::distance(sp.mesh_.begin(), m)].push_back(mc);
                        (*container)[std::distance(sp.mesh_.begin(), m)].push_back(mc);
                    }
                }
            }
        }
        else if (bou == BouType::interog)
        {
            for (Cell& c: arrival)
            {
                for (MeshCell& mc: c.interog_)
                {
                    //auto m = std::find_if(sp.mesh_.begin(), sp.mesh_.end(), [&](const Mesh& _m){return _m.parent_mesh() == mc.parent_mesh();});
                    auto m = std::find_if(sp.mesh_.begin(), sp.mesh_.end(), [&](const Mesh& _m){return _m.tag() == mc.parent_mesh();});
                    assert(m != sp.mesh_.end());
                    if (m->query_bou(mc.tag(), bou) == nullptr) {
                        //for (auto& mf: mc.face_p())
                        //{
                        //    mf.set_faceaddr(nullptr);
                        //}
                        //(*container)[std::distance(sp.mesh_.begin(), m)].push_back(mc);
                        (*container)[std::distance(sp.mesh_.begin(), m)].push_back(mc);
                    }
                }
            }
        }
    }

    void ArrivalConnecter::add_arrival_cells(SpatialPartition& sp)
    {
        //if (profiler_ != nullptr) {profiler_->start("asm-arr-add-cell-pre");}
        for (auto& msh: arr_interior_)
        {
            for (MeshCell& ac: msh)
            {
                ac.remove_all_neighbors(); // we do this because arrival cells came from mover not after when modified in disconnecter where 
                //remove_cell_from_cellhood(ic) is called.

                ac.deparent_neis_from_faces();

                //for (MeshFace& mf: ac.face_p())
                //{
                    //assert(mf.parent_cell().size() < 2);
                    //mf.remove_parent_cells();
                //}
            }
        }
        //if (profiler_ != nullptr) {profiler_->stop("asm-arr-add-cell-pre");}

        for (int i=0; i<arr_interior_.size(); ++i)
        {
            mcc& vec = arr_interior_[i];

            if (vec.empty())
            {
                continue;
            }

            //if (profiler_ != nullptr) {profiler_->start("asm-arr-add-int");}
            sp.mesh_[i].add_interiors_nonsorted_nopoint(vec);
            //if (profiler_ != nullptr) {profiler_->stop("asm-arr-add-int");}

            //if (profiler_ != nullptr) {profiler_->start("asm-arr-rem-dup");}
            //sp.mesh_[i].remove_duplicate_cells();
            sp.mesh_[i].remove_duplicate_cells(profiler_, "arr");
            //if (profiler_ != nullptr) {profiler_->stop("asm-arr-rem-dup");}

            //if (profiler_ != nullptr) {profiler_->start("asm-arr-add-pts");}
            sp.mesh_[i].add_points(vec, spc_->comm()->rank());
            //if (profiler_ != nullptr) {profiler_->stop("asm-arr-add-pts");}

            //if (profiler_ != nullptr) {profiler_->start("asm-arr-rpc");}
            sp.mesh_[i].remove_parent_cells_of_vertices_of_all_cells();
            //if (profiler_ != nullptr) {profiler_->stop("asm-arr-rpc");}

            //if (profiler_ != nullptr) {profiler_->start("asm-arr-spc");}
            sp.mesh_[i].set_parent_cell_of_vertices_of_all_cells();
            //if (profiler_ != nullptr) {profiler_->stop("asm-arr-spc");}

            //if (profiler_ != nullptr) {profiler_->start("asm-rem-pts");}
            sp.mesh_[i].remove_parent_cells_of_all_points();
            //if (profiler_ != nullptr) {profiler_->stop("asm-rem-pts");}

            //if (profiler_ != nullptr) {profiler_->start("asm-upd");}
            sp.mesh_[i].update_points_from_cell_vertices(spc_->comm()->rank());
            //if (profiler_ != nullptr) {profiler_->stop("asm-upd");}
        }

        //if (profiler_ != nullptr) {profiler_->start("asm-sort");}
        //for (Mesh& m: sp.mesh_)
        //{
            //m.sort_cells();
        //}
        //if (profiler_ != nullptr) {profiler_->stop("asm-sort");}
    }

    void ArrivalConnecter::add_arrival_bou(SpatialPartition& sp, BouType bou)
    {
        if (bou == BouType::wall)
        {
            for (int i=0; i<arr_wall_.size(); ++i)
            {
                sp.mesh_[i].add_bous(arr_wall_[i], bou);
            }
        }
        else if (bou == BouType::dirichlet)
        {
            for (int i=0; i<arr_dirichlet_.size(); ++i)
            {
                sp.mesh_[i].add_bous(arr_dirichlet_[i], bou);
            }
        }
        else if (bou == BouType::farfield)
        {
            for (int i=0; i<arr_farfield_.size(); ++i)
            {
                sp.mesh_[i].add_bous(arr_farfield_[i], bou);
            }
        }
        else if (bou == BouType::empty)
        {
            for (int i=0; i<arr_empty_.size(); ++i)
            {
                sp.mesh_[i].add_bous(arr_empty_[i], bou);
            }
        }
        else if (bou == BouType::interog)
        {
            for (int i=0; i<arr_interog_.size(); ++i)
            {
                sp.mesh_[i].add_bous(arr_interog_[i], bou);
            }
        }
    }

    void ArrivalConnecter::check_cell_exist2(const Mesh& m, const MeshCell& ac, int i)
    {
        if (m.query_sorted(ac.tag()) != nullptr)
        {
            assert(erase_interior_.size() > i);
            erase_interior_[i] = true;
        }
    }
}
