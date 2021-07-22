#include "regular_mesh.h"
#include <iomanip>

namespace Tailor
{
    void RegularMesh::register_bincells(const std::deque<BinCell>& bincell, const std::deque<Mesh>& meshes, bool pseudo3D, int rank, RegType regtype)
    {
        size_t quarter;
        if (pseudo3D) {
            quarter = bincell.size() / 4;
        }
        else {
            quarter = bincell.size() / 8;
        }

        for (Bin& b: bin_)
        {
            b.reserve(quarter);
        }

        //int count_1 = 0;
        //int count_2 = 0;
        //int count_3 = 0;
        //int count_4 = 0;
        //int count_5 = 0;
        //int count_6 = 0;
        for (const BinCell& bc: bincell)
        {
            //++count_1;
            Tag mt = bc.mesh();
            Tag ct = bc.cell();


            AABB _aabb;
            //const Mesh* _mesh = nullptr;
            const MeshCell* _cell = nullptr;
            for (auto m=meshes.begin(); m!=meshes.end(); ++m)
            {
                if (mt == m->tag())
                {
                    assert(m->query(ct));
                    //_mesh = &(*m);
                    _cell = &m->cell(ct);
                    _aabb = AABB(_cell->poly());
                    break;
                }
            }
            
            //++count_2;

            assert(_cell != NULL);

            for (Bin& b: bin_)
            {
                //assert(f->cross_len() != 0.);
                assert(!b.aabb().degenerate());
                assert(b.aabb().face_cross_len() == false);
                //++count_3;
                if (regtype == RegType::aabb)
                {
                    if (b.aabb().do_intersect(_aabb))
                    {
                        b.copy_bincell(bc);
                        if (b.aabb().do_intersect(_cell->poly().centroid()))
                        {
                            b.increment_mesh_load(mt);
                        }
                        else
                        {
                            b.insert_to_mesh_tag_index_map_all(mt);
                        }
                    }
                }
                else if (regtype == RegType::centroid)
                {
                    if (b.aabb().do_intersect(_cell->poly().centroid()))
                    {
                        b.copy_bincell(bc);
                        b.increment_mesh_load(mt);
                        //b.insert_to_mesh_tag_index_map_all(mt);
                    }
                }
                else
                {
                    assert(false);
                }
                /*else
                {
                    if (rank == 3)
                    {
                    std::cout << "min(0): " << b.aabb().min(0) << std::endl;
                    std::cout << "min(1): " << b.aabb().min(1) << std::endl;
                    std::cout << "min(2): " << b.aabb().min(2) << std::endl;
                    std::cout << "max(0): " << b.aabb().max(0) << std::endl;
                    std::cout << "max(1): " << b.aabb().max(1) << std::endl;
                    std::cout << "max(2): " << b.aabb().max(2) << std::endl;
                    std::cout << "cnt(0): " << _cell->poly().centroid()(0) << std::endl;
                    std::cout << "cnt(1): " << _cell->poly().centroid()(1) << std::endl;
                    std::cout << "cnt(2): " << _cell->poly().centroid()(2) << std::endl;
                    _mesh->print_as_vtk("aaa.vtk");
                    assert(false);
                    }
                }*/
            }
        }

        for (Bin& b: bin_)
        {
            b.shrink();
        }

        for (const Bin& b: bin_)
        {
            if (!b.mesh_tag_index_map_res().left.empty())
            {
                if (b.cell().empty())
                {
                    std::cout << "bintag: " << b.tag()() << std::endl;
                    std::cout << "rmtag: " << b.parent_rm()() << std::endl;
                }
                assert(!b.cell().empty());
            }
        }

        //std::cout << rank << " count: " << count_1 << " " << count_2 << " " << count_3 << " " << count_4  << " " << count_5 << " " << count_6 << std::endl;

        /*for (const Bin& b: bin_)
        {
            std::cout << "rmtag: " << b.parent_rm()() << std::endl;
            std::cout << "bintag: " << b.tag()() << std::endl;
            std::cout << "ml size: " << b.mesh_load().size() << std::endl;
            for (const auto& ml: b.mesh_load())
            {
                std::cout << "ml: " << ml << std::endl;
            }
        }*/

            /*for (const auto& b: bin_)
            {
                if (b.tag()() == 17)
                {
                    for (const auto& bc: b.cell())
                    {
                        if (bc.mesh()() == 1 && bc.cell()() == 274)
                        {
                            std::cout << "oooooooooooooooo: " << b.parent_rm()() << " " << b.tag()() << " " << rank << std::endl;
                            assert(false);
                        }
                    }
                }
            }*/
    }

    void RegularMesh::register_resident_mesh(const Mesh& mesh, int rank)
    {
        //assert(rank != 0);
        // ghosts must be connected to interiors at this stage.
        assert(!mesh.cell().empty());
        for (const MeshCell& c: mesh.cell())
        {
            //if (mesh.tag()() == 0 && c.tag()() == 18003)
            //{
                //assert(false);
            //}
            assert(c.tag().isvalid());
            assert(mesh.tag().isvalid());

            {
                register_resident_cell(c, c.tag(), mesh.tag(), rank);
            }
        }
        {
            bool all_empty = true;
            for (const Bin& b0: bin())
            {
                if (!b0.cell().empty())
                {
                    all_empty = false;
                    break;
                }
            }
            assert(!all_empty);
        }
    }

    void RegularMesh::register_overlapping_mesh(const Mesh& mesh, int rank, bool adaptive)
    {
        assert(!mesh.cell().empty());
        for (const MeshCell& c: mesh.cell())
        {
            if (aabb_.do_intersect(AABB(c.poly())) == false) {
                continue;
            }
            assert(c.tag().isvalid());
            assert(mesh.tag().isvalid());

            assert(aabb_.do_intersect(AABB(c.poly())));

            for (const Bin& b: bin_)
            {
                assert(b.mesh_tag_index_map_res().size() == b.mesh_load().size());
            }

            register_cell(c, c.tag(), mesh.tag(), rank, adaptive);

            for (const Bin& b: bin_)
            {
                assert(b.mesh_tag_index_map_res().size() == b.mesh_load().size());
            }
        }
        {
            bool all_empty = true;
            for (const Bin& b0: bin())
            {
                if (!b0.cell().empty())
                {
                    all_empty = false;
                    break;
                }
            }
            assert(!all_empty);
        }
    }

    void RegularMesh::register_mesh(const Mesh& mesh, int rank, bool adaptive)
    {
        //assert(!mesh.cell().empty());

        //if (rank == 0)
        //{
        //    mesh.print_as_vtk_geometry("zimesh.vtk");
        //    print("zioctree.vtk");
        //    assert(false);
        //}
        //std::cout << "rankkkk: " << rank << std::endl;

        for (const MeshCell& c: mesh.cell())
        {
            assert(c.tag().isvalid());
            assert(mesh.tag().isvalid());

            assert(aabb_.do_intersect(AABB(c.poly())));

            for (const Bin& b: bin_)
            {
                assert(b.mesh_tag_index_map_res().size() == b.mesh_load().size());
            }

            assert(aabb_.do_intersect(AABB(c.poly())));
            
            register_cell(c, c.tag(), mesh.tag(), rank, adaptive);

            for (const Bin& b: bin_)
            {
                assert(b.mesh_tag_index_map_res().size() == b.mesh_load().size());
            }
        }
    }

    void RegularMesh::register_resident_cell(const MeshCell& cell, const Tag& cell_tag, const Tag& mesh_tag, int rank)
    {
        //bool overlap = aabb_.do_intersect(AABB(cell.poly()));
        bool overlap = aabb_.do_intersect(cell.poly().centroid());
        //if (!overlap)
        //{
        //    auto aabb = AABB(cell.poly());
        //    std::cout << "cell min0: " << aabb.min(0) << std::endl;
        //    std::cout << "cell min1: " << aabb.min(1) << std::endl;
        //    std::cout << "cell max0: " << aabb.max(0) << std::endl;
        //    std::cout << "cell max1: " << aabb.max(1) << std::endl;

        //    std::cout << "rm min0: " << aabb_.min(0) << std::endl;
        //    std::cout << "rm min1: " << aabb_.min(1) << std::endl;
        //    std::cout << "rm max0: " << aabb_.max(0) << std::endl;
        //    std::cout << "rm max1: " << aabb_.max(1) << std::endl;
        //}
        //assert(overlap);
        if (!overlap)
        {
            return;
        }

        BinRMTag tag;
        //get_bintag_adaptive(cell.poly().centroid(), tag, false); // centroid might be inbetween bins. so tags.size coudl be > 1.
        get_bintag_adaptive_unique(cell.poly().centroid(), tag, rank, cell.tag()());
        assert(tag.isvalid());

        //for (const BinRMTag& brmt: tag) 
        {
            //Bin& _bin = bin_p(tag.front()); // chose only one bin even if multiple tags exist.
            //Bin& _bin = bin_p(brmt);
            //Bin& _bin = bin_p(tag);
            Bin& _bin = bin_p(tag.bintag());

            assert(cell_tag.isvalid());
            assert(mesh_tag.isvalid());

            bool point_in_bin = _bin.is_resident(cell.poly().centroid());

            /*if (cell.tag()() == 18003 && mesh_tag() == 0)
            {
                std::cout << "bin min 0: " << _bin.aabb().min(0) << std::endl;
                std::cout << "bin min 1: " << _bin.aabb().min(1) << std::endl;
                std::cout << "bin min 2: " << _bin.aabb().min(2) << std::endl;

                std::cout << "bin max 0: " << _bin.aabb().max(0) << std::endl;
                std::cout << "bin max 1: " << _bin.aabb().max(1) << std::endl;
                std::cout << "bin max 2: " << _bin.aabb().max(2) << std::endl;

                std::cout << "cell cnt 0: " << cell.poly().centroid()(0) << std::endl;
                std::cout << "cell cnt 1: " << cell.poly().centroid()(1) << std::endl;
                std::cout << "cell cnt 2: " << cell.poly().centroid()(2) << std::endl;

                for (const auto& pp: cell.poly().vertices())
                {
                    std::cout << "cell vertex 0: " << pp.r(0) << std::endl;
                    std::cout << "cell vertex 1: " << pp.r(1) << std::endl;
                    std::cout << "cell vertex 2: " << pp.r(2) << std::endl;
                }

                assert(false);
            }*/

            if (point_in_bin)
            {
                _bin.add_bincell(cell_tag, mesh_tag, rank);
                assert(!_bin.cell().empty());
                //if (calcload)
                {
                    _bin.increment_mesh_load(mesh_tag);
                }

                //break;
            }
            //else
            //{
                //_bin.insert_to_mesh_tag_index_map_all(mesh_tag);
            //}
        }
    }

    void RegularMesh::register_cell(const MeshCell& cell, const Tag& cell_tag, const Tag& mesh_tag, int rank, bool adaptive)
    {
        bool overlap = aabb_.do_intersect(AABB(cell.poly()));
        assert(overlap);
        if (!overlap)
        {
            return;
        }

        std::vector<BinRMTag> tag;
        tag.reserve(64);
        int dummy;
        for (const Bin& b: bin_)
        {
            assert(b.rm() == nullptr);
        }

        get_bintag_adaptive(AABB(cell.poly()), tag, dummy, -1);
        assert(!tag.empty());

        for (const BinRMTag& _tag: tag)
        {
            Bin& _bin = bin_p(_tag.bintag());

            if (_bin.aabb().do_intersect(AABB(cell.poly())))
            {
                assert(cell_tag.isvalid());
                assert(mesh_tag.isvalid());
                _bin.add_bincell(cell_tag, mesh_tag, rank);
                assert(!_bin.cell().empty());

                bool point_in_bin = _bin.is_resident(cell.poly().centroid());
                if (point_in_bin)
                {
                    //if (rank == 2)
                    //{
                    //    if (mesh_tag() == 1)
                    //    {
                    //        if (cell_tag() == 26447)
                    //        {
                    //            auto aabb = AABB(cell.poly());
                    //            std::cout << "rank: " << rank << std::endl;
                    //            std::cout << "min(0): " << aabb.min(0) << std::endl;
                    //            std::cout << "min(1): " << aabb.min(1) << std::endl;
                    //            std::cout << "min(2): " << aabb.min(2) << std::endl;
                    //            std::cout << "max(0): " << aabb.max(0) << std::endl;
                    //            std::cout << "max(1): " << aabb.max(1) << std::endl;
                    //            std::cout << "max(2): " << aabb.max(2) << std::endl;
                    //            std::cout << "bin: " << _bin.tag()() << std::endl;
                    //        }
                    //    }
                    //}
                    _bin.increment_mesh_load(mesh_tag);
                    assert(_bin.mesh_tag_index_map_res().size() == _bin.mesh_load().size());
                    assert(!_bin.mesh_tag_index_map_res().empty());
                }
                else
                {
                    _bin.insert_to_mesh_tag_index_map_all(mesh_tag);
                }
            }
        }
    }
}
