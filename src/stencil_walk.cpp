#include "stencil_walk.h"

namespace Tailor
{
    std::ostream& operator<<(std::ostream& os, const StencilWalkResult& res)
    {
        if (res == StencilWalkResult::inside_cell)
        {
            os << "inside cell";
        }
        else if (res == StencilWalkResult::inside_hole)
        {
            os << "inside hole";
        }
        else if (res == StencilWalkResult::outside_mesh)
        {
            os << "outside mesh";
        }
        else
        {
            assert(false);
        }
        return os;
    }

    // TODO: instead of checking all partition cells, checking only nearby partition cells is better by using ADT. 
    void check_partition(const Mesh& donor_mesh, const MeshCell*& current_cell, const Point& target, const MeshCell* closest_cell, bool& found, Tag& prev_nei)
    {
        Segment s(current_cell->poly().centroid(), target);

        for (const MeshCell& mc: donor_mesh.cell())
        {
            if (mc.tag() == current_cell->tag()) {
                continue;
            }
            //if (!mc.is_partition()) {
            if (mc.btype() != BouType::partition) {
                continue;
            }

            for (const auto& f: mc.face())
            {
                bool just_on_face;
                if (f.face().do_intersect(s, just_on_face))
                {
                    const Tag& temp_prev_nei = current_cell->tag();
                    prev_nei = temp_prev_nei;
                    current_cell = &mc;
                    closest_cell = current_cell;
                    found = true;
                    if (current_cell->oga_cell_type() == OGA_cell_type_t::ghost)
                    {
                        found = false;
                        return;
                    }
                    break;
                }
            }

            if (found) {break;}
        }
    }

    void check_boundary(const mcc& container, const Mesh& donor_mesh, const MeshCell*& current_cell, const Point& target, const MeshCell* closest_cell, bool& found, Tag& prev_nei/*, const ADT& outer_adt*/, const MeshCell& starting_cell)
    {
        if (container.empty()) return;

        Segment s(current_cell->poly().centroid(), target);
        for (const MeshCell& mc: container)
        {
            bool just_on_face;
            if (mc.face()[0].face().do_intersect(s, just_on_face))
            {
                //std::vector<int> ftt = {1, 2005, 4163, 4428};
                //if (mc.face()[0].tag() == FaceTag(ftt))
                //{
                    //if (current_cell->tag()() == 472)
                    //{
                        //std::cout << "target(0): " << target.r(0) << std::endl;
                        //std::cout << "target(1): " << target.r(1) << std::endl;
                        //assert(mc.face()[0].tag() != FaceTag(ftt));
                    //}
                //}
                const Tag& temp_prev_nei = current_cell->tag();
                const MeshCell& gc = mc;
                assert(gc.interior_boundary().isvalid());
                const MeshCell* temp_cell = &donor_mesh.cell(gc.interior_boundary());
                //const MeshCell* temp_cell = gc.interior_boundary().const_addr();
                if (temp_cell->tag() == temp_prev_nei)
                {
                    if (current_cell->tag() != starting_cell.tag()) {
                        continue;
                    }
                }
                prev_nei = temp_prev_nei;
                current_cell = temp_cell;
                closest_cell = current_cell;
                found = true;
                break;
            }
        }
    }

    BouType find_intersected_face(const Mesh& mesh, const MeshCell*& current_cell, Tag& prev_nei, Tag& prev_prev_nei, const MeshCell*& closest_cell, const MeshFace*& inter_mf, const MeshCell*& nei, const Segment& start_to_target, bool verbose)
    {
        assert(current_cell->oga_cell_type() != OGA_cell_type_t::ghost);

        inter_mf = NULL;
        nei = NULL;
        assert(prev_nei.isvalid());

        if (verbose)
        {
            std::cout << "start to target 0: " << start_to_target.vertex(0).r(0) << " " << start_to_target.vertex(0).r(1) << " " << start_to_target.vertex(0).r(2) << std::endl;
            std::cout << "start to target 1: " << start_to_target.vertex(1).r(0) << " " << start_to_target.vertex(1).r(1) << " " << start_to_target.vertex(1).r(2) << std::endl;
        }

        //for (const MeshFace& imf: current_cell->face())
        for (const MeshFace& mf: current_cell->face())
        {
            //auto fiter = std::find_if(mesh.face().begin(), mesh.face().end(), [&](const auto& ff){return ff.tag() == imf.tag();});
            //if (fiter == mesh.face().end())
            //{
            //    std::cout << "current cell: " << current_cell->tag()() << std::endl;
            //    std::cout << "imf: " << imf.tag() << std::endl;
            //    std::cout << "face type: " << static_cast<int>(imf.btype()) << std::endl;
            //    std::cout << "mesh face size: " << mesh.face().size() << std::endl;
            //    mesh.print_as_vtk("probmesh.vtk");
            //}
            //assert(fiter != mesh.face().end());

            //const MeshFace& mf = mesh.face(imf.tag());

                /*if (verbose) {
                    std::cout << "face tag: " << mf.tag()() << std::endl;
                    if (verbose)
                    {
                        std::cout << "face vertex: " << mf.face().vertex(0).r(0) << " " << mf.face().vertex(0).r(1) << " " << mf.face().vertex(0).r(2) << std::endl;
                        std::cout << "face vertex: " << mf.face().vertex(1).r(0) << " " << mf.face().vertex(1).r(1) << " " << mf.face().vertex(1).r(2) << std::endl;
                        std::cout << "face vertex: " << mf.face().vertex(2).r(0) << " " << mf.face().vertex(2).r(1) << " " << mf.face().vertex(2).r(2) << std::endl;
                        std::cout << "face vertex: " << mf.face().vertex(3).r(0) << " " << mf.face().vertex(3).r(1) << " " << mf.face().vertex(3).r(2) << std::endl;
                    }
                }*/

            const Polygon& s = mf.face();
            //assert(current_cell->btype() == BouType::interior || current_cell->btype() == BouType::partition);
            assert(!mf.parent_cell().empty());
            Tag inei;
            for (const auto& pc: mf.parent_cell())
            {
                //if (!mf.is_boundary())
                //{
                    //assert(mesh.query(pc) != nullptr);
                    //assert(pc.const_addr() != nullptr);
                    //assert(mesh.cell(pc).btype() == BouType::interior || mesh.cell(pc).btype() == BouType::partition);
                    //assert(pc.const_addr()->btype() == BouType::interior || pc.const_addr()->btype() == BouType::partition);
                //}
                if (pc == current_cell->tag()) continue;
                inei = pc;
                assert(inei.isvalid());
            }

            //if (mf.is_boundary())
            //{
                //assert(mesh.query(mf.parent_cell()[1]) != nullptr);
                //assert(mf.const_parent_cell()[1].const_addr() != nullptr);
            //}

            bool just_on_face;
            bool inter = s.do_intersect(start_to_target, just_on_face);
            if (inter)
            {
                if (verbose) {
                    std::cout << "mf inter" << std::endl;
                }
                inter_mf = &mf;
                assert(inter_mf != NULL);
                if (mf.is_boundary() || mf.btype() == BouType::partition) {
                    if (verbose) {
                        std::cout << "mf is boundary or partition" << std::endl;
                    }
                    return mf.btype();
                }
                if (!inei.isvalid())
                {
                    std::cout << "current cell btype: " << static_cast<int>(current_cell->btype()) << std::endl;
                    std::cout << "mf btype: " << static_cast<int>(mf.btype()) << std::endl;
                    std::cout << "current cell: " << current_cell->tag()() << std::endl;
                    std::cout << "mf.parent_cell().size(): " << mf.parent_cell().size() << std::endl;
                    for (const auto& pcc: mf.parent_cell())
                    {
                        std::cout << "pc: " << pcc() << std::endl;
                    }

                }
                assert(inei.isvalid());
                if (mf.btype() != BouType::interior)
                {
                    std::cout << "btype: " << static_cast<int>(mf.btype()) << std::endl;
                }
                assert(mf.btype() == BouType::interior);
                //if (!inei.isvalid()) // possible if current cell is partition or near bou cell.
                    //continue;

                //assert(mesh.query(inei) != nullptr);
                nei = &mesh.cell(inei);
                if (nei->oga_cell_type() == OGA_cell_type_t::ghost)
                {
                    return BouType::partition;
                }
                //if (nei->is_ghost())
                    //return;
                prev_prev_nei = prev_nei;
                prev_nei = current_cell->tag();
                current_cell = nei;
                assert(current_cell->tag() != prev_nei);
                //closest_cell = current_cell->tag();
                closest_cell = current_cell;
                return BouType::interior;
                //break;
            }
            else
            {
                if (verbose) {
                    std::cout << "mf no inter" << std::endl;
                }
            }
        }

        return BouType::undefined;
    }

    StencilWalkResult stencil_walk(const Mesh& donor_mesh, const Point& target, const MeshCell& starting_cell, const MeshCell*& closest_cell/*, const ADT& wall_adt, const ADT& outer_adt*/, int dummyrank, int& iter, bool verbose)
    {
        assert(starting_cell.oga_cell_type() != OGA_cell_type_t::ghost);
        assert(!donor_mesh.cell().empty());
        const MeshCell* current_cell = &starting_cell;
        //assert(current_cell->is_interior());
        //assert(current_cell->btype() == BouType::interior || current_cell->btype() == BouType::partition);
        const MeshFace* inter_mf = NULL;
        const MeshCell* nei = NULL;
        closest_cell = current_cell;
        Tag prev_nei = current_cell->tag();
        Tag prev_prev_nei = current_cell->tag();
        std::vector<Tag> prev_partition;
        iter = 0;

        std::vector<int> current_cell_list;

        while (true)
        {
                if (verbose) {
                    std::cout << "current cell: " << current_cell->tag()() << std::endl;
                }
            if (iter >= 1000)
            {
                if (dummyrank == 6)
                {
                    std::cout << "starting cell: " << starting_cell.tag()() << std::endl;
                    std::cout << "current cell: " << current_cell->tag()() << std::endl;
                    std::cout << "target(0): " << target.r(0) << std::endl;
                    std::cout << "target(1): " << target.r(1) << std::endl;
                    std::cout << "target(2): " << target.r(2) << std::endl;
                    donor_mesh.print_as_vtk("donormesh.vtk");
                    for (int ii: current_cell_list)
                    {
                        std::cout << "cell list: " << ii << std::endl;
                    }
                }
            }
            assert(iter < 1000);
            current_cell_list.push_back(current_cell->tag()());
            Segment start_to_target(current_cell->poly().centroid(), target);
            if (iter != 0)
            {
                if (current_cell->tag() == prev_prev_nei)
                {
                    closest_cell = current_cell;
                    assert(current_cell->oga_cell_type() != OGA_cell_type_t::ghost);
                    return StencilWalkResult::inside_cell;
                }
            }
            assert(iter == 0 || current_cell->tag() != prev_nei);
            BouType btype = find_intersected_face(donor_mesh, current_cell, prev_nei, prev_prev_nei, closest_cell, inter_mf, nei, start_to_target, verbose);
            if (inter_mf == NULL)
            {
                if (verbose) {
                    std::cout << "inter_mf is NULL" << std::endl;
                }
                closest_cell = current_cell;
                //if (!current_cell->poly().do_intersect(target.r()))
                //{
                //    std::cout << "target.r(0): " << target.r(0) << std::endl;
                //    std::cout << "target.r(1): " << target.r(1) << std::endl;
                //    std::cout << "target.r(2): " << target.r(2) << std::endl;
                //    for (const auto& v: current_cell->poly().vertices())
                //    {
                //        std::cout << "v(0): " << v.r(0) << " " << v.r(1) << " " << v.r(2) << std::endl;
                //    }
                //    current_cell->poly().do_intersect(target.r(), true);
                //}
                //assert(current_cell->poly().do_intersect(target.r()));
                return StencilWalkResult::inside_cell;
            }
            else if (btype == BouType::wall)
            {
                bool found = false;
                check_boundary(donor_mesh.wall_boundaries(), donor_mesh, current_cell, target, closest_cell, found, prev_nei/*, wall_adt*/, starting_cell);

                if (verbose) {
                    std::cout << "btype is wall. found=" << found << std::endl;
                }

                if (!found)
                {
                    //if (current_cell->poly().do_intersect(target.r()))
                    //{
                    //    donor_mesh.print_as_vtk("domesh.vtk");
                    //    std::cout << "target(0):" << target.r(0) << std::endl;
                    //    std::cout << "target(1):" << target.r(1) << std::endl;
                    //    std::cout << "target(2):" << target.r(2) << std::endl;
                    //    std::cout << "current cell" << current_cell->tag()() << std::endl;
                    //    std::cout << "intermf:" << inter_mf->tag() << std::endl;
                    //    std::cout << "starting cell:" << starting_cell.tag()() << std::endl;
                    //    for (auto iii: current_cell_list)
                    //    {
                    //        std::cout << "prev cell:" << iii << std::endl;
                    //    }
                    //}

                    //assert(!current_cell->poly().do_intersect(target.r()));
                    return StencilWalkResult::inside_hole;
                }
            }
            else if (btype == BouType::dirichlet)
            {
                bool found = false;
                check_boundary(donor_mesh.dirichlet_boundaries(), donor_mesh, current_cell, target, closest_cell, found, prev_nei/*, outer_adt*/, starting_cell);

                if (verbose) {
                    std::cout << "btype is dirichlet. found=" << found << std::endl;
                }

                if (!found)
                {
                    //assert(!current_cell->poly().do_intersect(target.r()));
                    return StencilWalkResult::outside_mesh;
                }
            }
            else if (btype == BouType::farfield)
            {
                bool found = false;
                check_boundary(donor_mesh.farfield_boundaries(), donor_mesh, current_cell, target, closest_cell, found, prev_nei/*, outer_adt*/, starting_cell);

                if (verbose) {
                    std::cout << "btype is farfield. found=" << found << std::endl;
                }

                if (!found)
                {
                    //assert(!current_cell->poly().do_intersect(target.r()));
                    return StencilWalkResult::outside_mesh;
                }
            }
            else if (btype == BouType::interog)
            {
                bool found = false;
                check_boundary(donor_mesh.interog_boundaries(), donor_mesh, current_cell, target, closest_cell, found, prev_nei/*, outer_adt*/, starting_cell);

                if (verbose) {
                    std::cout << "btype is interog. found=" << found << std::endl;
                }

                if (!found)
                {
                    //assert(!current_cell->poly().do_intersect(target.r()));
                    return StencilWalkResult::outside_mesh;
                }
            }
            else if (btype == BouType::empty)
            {
                bool found = false;
                check_boundary(donor_mesh.empty_boundaries(), donor_mesh, current_cell, target, closest_cell, found, prev_nei/*, outer_adt*/, starting_cell);

                if (verbose) {
                    std::cout << "btype is empty. found=" << found << std::endl;
                }

                if (!found)
                {
                    //assert(!current_cell->poly().do_intersect(target.r()));
                    return StencilWalkResult::outside_mesh;
                }
            }
            else if (btype == BouType::partition)
            {
                bool found = false;
                check_partition(donor_mesh, current_cell, target, closest_cell, found, prev_nei);
                if (verbose) {
                    std::cout << "btype is partition. found=" << found << std::endl;
                }
                if (!found) {
                    //assert(!current_cell->poly().do_intersect(target.r()));
                    return StencilWalkResult::outside_mesh;
                }
            }
            else if (btype == BouType::interior)
            {
                if (verbose) {
                    std::cout << "btype is interior" << std::endl;
                }
                // just continue. nothing to do.
            }
            else
            {
                std::cout << "btypee: " << static_cast<int>(btype) << std::endl;
                assert(false);
            }

            ++iter;
        }
    }
}
