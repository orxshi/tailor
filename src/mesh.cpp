#include "mesh.h"

namespace Tailor
{
    void Mesh::increase_overlap_thickness_(MeshCell& mc, int& count, int nlayer, const ADT& passive_cell_adt, Mesh& passive_mesh)
    {
        assert(mc.oga_cell_type() == OGA_cell_type_t::mandat_receptor || mc.oga_cell_type() == OGA_cell_type_t::undefined);

        if (count >= nlayer)
        {
            return;
        }

        ++count;

        for (const auto& inei: mc.pnei())
        {
            auto& nei = cell_p(inei);

            if (nei.oga_cell_type() == OGA_cell_type_t::field || nei.oga_cell_type() == OGA_cell_type_t::undefined)
            {
                const Vector3& target = nei.poly().centroid();
                ADTPoint targetadt(target, nei.tag()());
                std::vector<int> res = passive_cell_adt.search(targetadt);

                assert(!res.empty());

                for (int r: res)
                {
                    MeshCell& passive_cell = passive_mesh.cell_p(Tag(r));

                    if (passive_cell.poly().do_intersect(target))
                    {
                        nei.set_oga_cell_type(OGA_cell_type_t::undefined); // works if convert_undefined_to_field() is called before increase_overlap_thickness.
                        assert(passive_cell.oga_cell_type() != OGA_cell_type_t::hole);
                        passive_cell.set_oga_cell_type(OGA_cell_type_t::field);
                        nei.add_cand_donor(passive_cell.parent_mesh(), passive_cell.tag(), &passive_cell);
                    }
                }

            }
        }

        for (const auto& inei: mc.pnei())
        {
            if (mc.oga_cell_type() == OGA_cell_type_t::mandat_receptor)
            {
                count = 0;
            }

            auto& nei = cell_p(inei);

            if (nei.oga_cell_type() == OGA_cell_type_t::undefined)
            {
                increase_overlap_thickness_(nei, count, nlayer, passive_cell_adt, passive_mesh);
            }
        }
    }

    void Mesh::increase_overlap_thickness(int nlayer, const ADT& passive_cell_adt, Mesh& passive_mesh)
    {
        for (auto& mc: cell_)
        {
            if (mc.oga_cell_type() != OGA_cell_type_t::mandat_receptor) {
                continue;
            }

            //if (!mc.near_interog()) {
                //continue;
            //}

            int count = 0;

            increase_overlap_thickness_(mc, count, nlayer, passive_cell_adt, passive_mesh);
        }

        for (auto& mc: cell_)
        {
            if (mc.oga_cell_type() == OGA_cell_type_t::undefined)
            {
                mc.set_oga_cell_type(OGA_cell_type_t::mandat_receptor);
            }
        }
    }

    void Mesh::reset_to_mid()
    {
        for (auto& mc: cell_)
        {
            mc.R_ = mc.R_mid_;
            mc.D_ = mc.D_mid_;
        }        
    }

    Mesh::Mesh(): priority_(-1)
    {
    }
    /*void Mesh::update_cell_pnei_addresses()
    {
        for (auto& mc: cell_)
        {
            if (mc.oga_cell_type() == OGA_cell_type_t::ghost) {
                continue;
            }

            for (auto& pnei: mc.pnei_)
            {
                if (pnei.tag() == mc.tag())
                {
                    pnei.set_addr(&mc);
                }
                else
                {
                    if (mc.btype() == BouType::interior || mc.btype() == BouType::partition)
                    {
                        pnei.set_addr(&cell_p(pnei.tag()));
                    }
                    else
                    {
                        pnei.set_addr(&boundary(mc.btype(), pnei.tag()));
                    }
                }
            }
        }
        
        auto fn = [&](auto& container)
        {
            for (auto& mc: container)
            {
                mc.interior_boundary_.set_addr(&cell_p(mc.interior_boundary().tag()));
            }
        };

        fn(wall_boundaries_);
        fn(dirichlet_boundaries_);
        fn(farfield_boundaries_);
        fn(empty_boundaries_);
        fn(interog_boundaries_);
    }

    void Mesh::update_face_parent_addresses()
    {
        for (auto& mf: face_)
        {
            if (mf.is_boundary())
            {
                auto& pc = mf.parent_cell(); 
                pc[0].set_addr(&boundary(mf.btype(), pc[0].tag()));
                pc[1].set_addr(&cell_p(pc[1].tag()));
            }
            else
            {
                for (auto& pc: mf.parent_cell())
                {
                    pc.set_addr(&cell_p(pc.tag()));
                }
            }
        }

        for (auto& mc: cell_)
        {
            if (mc.oga_cell_type() == OGA_cell_type_t::ghost) {
                continue;
            }

            for (auto& mf: mc.face_)
            {
                mf.parent_cell() = face_p(mf.tag()).parent_cell();
            }
        }
       
        auto fn = [&](auto& container)
        {
            for (auto& mc: container)
            {
                for (auto& mf: mc.face_)
                {
                    mf.parent_cell() = face_p(mf.tag()).parent_cell();
                }
            }
        };

        fn(wall_boundaries_);
        fn(dirichlet_boundaries_);
        fn(farfield_boundaries_);
        fn(empty_boundaries_);
        fn(interog_boundaries_);
    }

    void Mesh::update_cell_vertex_addresses()
    {
        for (auto& mc: cell_)
        {
            for (auto& mp: mc.point_)
            {
                for (auto& pc: mp.parent_cell())
                {
                    if (pc.tag() == mc.tag())
                    {
                        pc.set_addr(&mc);
                        assert(pc.addr() == &cell(pc.tag()));
                    }
                }
            }
        }

        auto fn = [&](auto& container, BouType btype)
        {
            for (auto& mc: container)
            {
                for (auto& mp: mc.point_)
                {
                    for (auto& pc: mp.parent_cell())
                    {
                        if (pc.tag() == mc.tag())
                        {
                            pc.set_addr(&mc);
                            assert(pc.addr() == &boundary(btype, pc.tag()));
                        }
                    }
                }
            }
        };

        fn(wall_boundaries_, BouType::wall);
        fn(dirichlet_boundaries_, BouType::dirichlet);
        fn(farfield_boundaries_, BouType::farfield);
        fn(empty_boundaries_, BouType::empty);
        fn(interog_boundaries_, BouType::interog);
    }*/

    /*void Mesh::set_oga_cell_type(const std::vector<Receiver<DonorInfo2>>& di)
    {
        for (MeshCell& mc: cell_)
        {
            auto iter = std::find_if(di.begin(), di.end(), [&](const auto& d){return (di.meshtag_ == tag_()() && di.celltag_ == mc.tag()());});
            assert(iter != di.end());
            mc.set_oga_cell_type(di.oga_cell_type_);
        }
    }*/

    void Mesh::determine_non_residents(const AABB& aabb)
    {
        for (MeshCell& mc: cell_)
        {
            if (!aabb.do_intersect(mc.poly().centroid()))
            {
                mc.set_oga_cell_type(OGA_cell_type_t::non_resident);
            }
        }
    }

    void Mesh::oga_interpolate(const std::deque<Mesh>& mesh, int rank)
    {
        assert(!mesh.empty());

        for (MeshCell& mc: cell_)
        {
            if (mc.oga_cell_type() == OGA_cell_type_t::receptor || mc.oga_cell_type() == OGA_cell_type_t::mandat_receptor)
            {
                assert(mc.donor().mesh_tag_.isvalid());
                assert(mc.donor().cell_tag_.isvalid());

                auto m = std::find_if(mesh.begin(), mesh.end(), [&](const auto& mm){return mm.tag() == mc.donor().mesh_tag_;});
                assert(m != mesh.end());
                assert(m->tag() != tag_);
                const auto& donor_cell = m->cell(mc.donor().cell_tag_);

                mc.prim_ = donor_cell.prim_;
                mc.prim_(1) -= donor_cell.vgn()(0);
                mc.prim_(2) -= donor_cell.vgn()(1);
                mc.prim_(3) -= donor_cell.vgn()(2);

                assert(!mc.prim_.isnan());
            }
        }
    }

    void Mesh::oga_interpolate(const ArrCon<Var>& arrival, int rank)
    {
        if (arrival.empty())
        {
            for (MeshCell& mc: cell_)
            {
                assert(mc.oga_cell_type() != OGA_cell_type_t::receptor && mc.oga_cell_type() != OGA_cell_type_t::mandat_receptor);
                assert(mc.oga_cell_type() != OGA_cell_type_t::undefined);
            }
        }

        //assert(!arrival.empty());

        for (MeshCell& mc: cell_)
        {
            if (mc.oga_cell_type() == OGA_cell_type_t::receptor || mc.oga_cell_type() == OGA_cell_type_t::mandat_receptor)
            {
                assert(mc.donor().mesh_tag_.isvalid());
                assert(mc.donor().cell_tag_.isvalid());

                auto iter = std::find_if(arrival.begin(), arrival.end(), [&](const auto& arr){return (arr.mesh_cell_.first == mc.donor().mesh_tag_() && arr.mesh_cell_.second == mc.donor().cell_tag_());});
                if (iter == arrival.end())
                {
                    std::cout << "rank: " << rank << std::endl;
                    std::cout << "parent mesh: " << mc.parent_mesh()() << std::endl;
                    std::cout << "mc: " << mc.tag()() << std::endl;
                    std::cout << "donor mesh: " << mc.donor().mesh_tag_() << std::endl;
                    std::cout << "donor cell: " << mc.donor().cell_tag_() << std::endl;
                    std::cout << "arrival.size(): " << arrival.size() << std::endl;
                    for (const auto& arr: arrival)
                    {
                        std::cout << "inside: " << arr.mesh_cell_.first << " " << arr.mesh_cell_.second << std::endl;
                    }
                }
                assert(iter != arrival.end());

                assert(!iter->var_.isnan());

                mc.prim_ = iter->var_;

                //assert(mc.prim(0) > 0.);
            }
        }
    }

    void Mesh::update_ghost_primitives(const ArrCon<Var>& arrival, int rank, double gamma)
    {
        if (arrival.empty())
        {
            std::cout << "arrival empty rank: " << rank << std::endl;
            print_as_vtk("empty.vtk");
        }
        assert(!arrival.empty());

        for (const Var& v: arrival)
        {
            if (v.mesh_cell_.first != tag_()) {
                continue;
            }

            assert(query(Tag(v.mesh_cell_.second)) != nullptr);
            MeshCell& mc = cell_p(Tag(v.mesh_cell_.second));

            assert(mc.oga_cell_type() == OGA_cell_type_t::non_resident || mc.oga_cell_type() == OGA_cell_type_t::ghost);
            //if (mc.oga_cell_type() != OGA_cell_type_t::ghost)
            //{
                //std::cout << "type: " << static_cast<int>(mc.oga_cell_type()) << std::endl;
            //}
            //assert(mc.oga_cell_type() == OGA_cell_type_t::ghost);

            assert(!v.var_.isnan());

            mc.prim_ = v.var_;
            mc.cons_sp1_ = prim_to_cons(mc.prim_, gamma);
            mc.dQ_ = v.dQ_;

            assert(mc.prim(0) > 0.);
        }
    }

    std::vector<Tag> ordered_bou_points(std::vector<MeshCell> input)
    {
        std::vector<Tag> output;

        for (const MeshPoint& mp: input.front().point())
        {
            if (mp.p().r(2) == 0.)
            {
                output.push_back(mp.tag());
                //std::cout << "output.back: " << output.back()() << std::endl;
            }
        }
        input.erase(input.begin());

        while(!input.empty())
        {
            bool loop = true;
            bool found = false;
            for (auto iter = input.begin(); iter != input.end(); ++iter)
            {
                bool consequtive = false;
                for (const MeshPoint& mp: iter->point())
                {
                    if (mp.p().r(2) != 0.) {
                        continue;
                    }
                    if (mp.tag() == output.back())
                    {
                        consequtive = true;
                        break;
                    }
                }
                if (!consequtive) {
                    continue;
                }
                for (const MeshPoint& mp: iter->point())
                {
                    if (mp.p().r(2) != 0.) {
                        continue;
                    }
                    if (mp.tag() == output.back())
                    {
                        continue;
                    }
                    {
                        output.push_back(mp.tag());
                        //std::cout << "output.back: " << output.back()() << std::endl;
                        input.erase(iter);
                        found = true;
                        loop = false;
                        break;
                    }
                }
                if (found) {
                    break;
                }
            }

            if (loop)
            {
                assert(false);
            }
        }

        return output;
    }

    //void convert_mesh_to_fortran(const Mesh& mesh)
    //{
    //    std::ofstream out;
    //    out.open("fortran.grid");

    //    out << mesh.point().size() << " 0 " << mesh.cell().size() << "\n";

    //    for (const MeshPoint& mp: mesh.point())
    //    {
    //        out << mp.p().r(0) << " " << mp.p().r(1) << "\n";
    //    }

    //    for (const MeshCell& mc: mesh.cell())
    //    {
    //        for (const MeshFace& mf: mc.face())
    //        {
    //            const MeshFace& rmf = mesh.face(mf.tag());

    //            bool allzero = true;
    //            for (const Point& p: rmf.face().vertex())
    //            {
    //                if (p.r(2) != 0.) {
    //                    allzero = false;
    //                    break;
    //                }
    //            }
    //            if (!allzero) {
    //                continue;
    //            }

    //            for (const Tag& t: rmf.mesh_point())
    //            {
    //                out << t() << " ";
    //            }
    //            out << "\n";
    //        }
    //    }

    //    out << "2\n"; // number of boundaries
    //    out << mesh.wall_boundaries().size() + 1 << "\n";
    //    out << mesh.dirichlet_boundaries().size() + 1 << "\n";

    //    std::cout << "wall size: " << mesh.wall_boundaries().size() << std::endl;
    //    std::cout << "dirichlet size: " << mesh.dirichlet_boundaries().size() << std::endl;

    //    auto container0 = ordered_bou_points(mesh.wall_boundaries());
    //    auto container1 = ordered_bou_points(mesh.dirichlet_boundaries());

    //    std::cout << "container0 size: " << container0.size() << std::endl;
    //    assert(container0.size() == mesh.wall_boundaries().size() + 1);
    //    std::cout << "container1 size: " << container1.size() << std::endl;
    //    assert(container1.size() == mesh.dirichlet_boundaries().size() + 1);

    //    for (auto iter = container0.begin(); iter != container0.end(); ++iter)
    //    {
    //        out << (*iter)() << "\n";
    //    }

    //    for (auto iter = container1.begin(); iter != container1.end(); ++iter)
    //    {
    //        out << (*iter)() << "\n";
    //    }

    //    out.close();
    //}
    const MeshCell* Mesh::boundary(BouType btype, const Tag& t) const
    {
        if (btype == BouType::wall)
        {
            return wall_boundary(t);
        }
        else if (btype == BouType::symmetry)
        {
            return symmetry_boundary(t);
        }
        else if (btype == BouType::dirichlet)
        {
            return dirichlet_boundary(t);
        }
        else if (btype == BouType::farfield)
        {
            return farfield_boundary(t);
        }
        else if (btype == BouType::empty)
        {
            return empty_boundary(t);
        }
        else if (btype == BouType::interog)
        {
            return interog_boundary(t);
        }
        else
        {
            assert(false);
        }
    }

    MeshCell* Mesh::boundary(BouType btype, const Tag& t)
    {
        if (btype == BouType::wall)
        {
            return wall_boundary_p(t);
        }
        else if (btype == BouType::symmetry)
        {
            return symmetry_boundary_p(t);
        }
        else if (btype == BouType::dirichlet)
        {
            return dirichlet_boundary_p(t);
        }
        else if (btype == BouType::farfield)
        {
            return farfield_boundary_p(t);
        }
        else if (btype == BouType::empty)
        {
            return empty_boundary_p(t);
        }
        else if (btype == BouType::interog)
        {
            return interog_boundary_p(t);
        }
        else
        {
            std::cout << "btype: " << static_cast<int>(btype) << std::endl;
            assert(false);
        }
    }
    void Mesh::set_all_cells_as_interior()
    {       
        for (MeshCell& mc: cell_)
        {   
            mc.set_btype(BouType::interior);
        }   
    }

    std::tuple<const MeshCell*, const MeshCell*> left_and_right_cells(const Mesh& mesh, const MeshFace& mf, const Tag& mctag)
    {
        const MeshCell *LC = nullptr;
        const MeshCell *RC = nullptr;

        if (mf.is_boundary())
        {
            assert(mf.parent_cell().size() == 2);
            //if (mesh.query(mf.left_cell()) == nullptr)
            //if (mf.left_cell().addr() == nullptr)
            //{
            //    std::cout << "btype" << static_cast<int>(mf.btype()) << std::endl;
            //    std::cout << "mf.left_cell(): " << mf.left_cell().tag()() << std::endl;
            //    std::cout << "mf.right_cell(): " << mf.right_cell().tag()() << std::endl;
            //    mesh.print_as_vtk("ard.vtk");
            //}
            //assert(mesh.query(mf.left_cell()) != nullptr);
            //assert(mf.left_cell().addr() != nullptr);
            LC = &mesh.cell(mf.left_cell());
            //LC = mf.left_cell().addr();
            RC = mesh.boundary(mf.btype(), mf.right_cell());
            assert(RC != nullptr);
            //RC = mf.right_cell().addr();
        }
        else if (mf.btype() == BouType::partition)
        {
            assert(mf.parent_cell().size() == 2);
            //if (mesh.query(mf.left_cell()) == nullptr || mesh.query(mf.right_cell()) == nullptr)
            //if (mf.left_cell().addr() == nullptr || mf.right_cell().addr() == nullptr)
            //{
            //    std::cout << "btype" << static_cast<int>(mf.btype()) << std::endl;
            //    std::cout << "mf.left_cell(): " << mf.left_cell().tag()() << std::endl;
            //    std::cout << "mf.right_cell(): " << mf.right_cell().tag()() << std::endl;
            //    mesh.print_as_vtk("ard.vtk");
            //}
            //assert(mesh.query(mf.left_cell()) != nullptr);
            //assert(mf.left_cell().addr() != nullptr);
            //assert(mesh.query(mf.right_cell()) != nullptr);
            //assert(mf.right_cell().addr() != nullptr);
            //LC = &mesh.cell_p(mf.left_cell());
            LC = &mesh.cell(mf.right_cell());
            //LC = mf.left_cell().addr();
            //RC = &mesh.cell_p(mf.right_cell());
            RC = &mesh.cell(mf.left_cell());
            //RC = mf.right_cell().addr();
            if (RC->oga_cell_type() == OGA_cell_type_t::non_resident)
            {
                std::cout << "LC oga: " << static_cast<int>(LC->oga_cell_type()) << std::endl;
                std::cout << "RC oga: " << static_cast<int>(RC->oga_cell_type()) << std::endl;
            }
            assert(RC->oga_cell_type() != OGA_cell_type_t::non_resident);
        }
        else
        {
            assert(mf.parent_cell().size() == 2);
            //if (mesh.query(mf.left_cell()) == nullptr || mesh.query(mf.right_cell()) == nullptr)
            //if (mf.left_cell().addr() == nullptr || mf.right_cell().addr() == nullptr)
            //{
            //    std::cout << "btype" << static_cast<int>(mf.btype()) << std::endl;
            //    std::cout << "mf.left_cell(): " << mf.left_cell().tag()() << std::endl;
            //    std::cout << "mf.right_cell(): " << mf.right_cell().tag()() << std::endl;
            //    mesh.print_as_vtk("ard.vtk");
            //}

            //LC = &mesh.cell_p(mf.left_cell());
            LC = &mesh.cell(mf.right_cell());
            //LC = mf.left_cell().addr();
            //RC = &mesh.cell_p(mf.right_cell());
            RC = &mesh.cell(mf.left_cell());
            //RC = mf.right_cell().addr();

        }

        assert(LC->tag() == mctag);
        assert(RC->tag() != mctag);

        return std::make_tuple(LC, RC);
    }

    std::tuple<MeshCell*, MeshCell*> left_and_right_cells(Mesh& mesh, MeshFace& mf, const Tag& mctag)
    {
        MeshCell *LC = nullptr;
        MeshCell *RC = nullptr;

        if (mf.is_boundary())
        {
            assert(mf.parent_cell().size() == 2);
            //if (mesh.query(mf.left_cell()) == nullptr)
            //if (mf.left_cell().addr() == nullptr)
            //{
            //    std::cout << "btype" << static_cast<int>(mf.btype()) << std::endl;
            //    std::cout << "mf.left_cell(): " << mf.left_cell().tag()() << std::endl;
            //    std::cout << "mf.right_cell(): " << mf.right_cell().tag()() << std::endl;
            //    mesh.print_as_vtk("ard.vtk");
            //}
            //assert(mesh.query(mf.left_cell()) != nullptr);
            //assert(mf.left_cell().addr() != nullptr);
            LC = &mesh.cell_p(mf.left_cell());
            //LC = mf.left_cell().addr();
            RC = mesh.boundary(mf.btype(), mf.right_cell());
            assert(RC != nullptr);
            //RC = mf.right_cell().addr();
        }
        else if (mf.btype() == BouType::partition)
        {
            assert(mf.parent_cell().size() == 2);
            //if (mesh.query(mf.left_cell()) == nullptr || mesh.query(mf.right_cell()) == nullptr)
            //if (mf.left_cell().addr() == nullptr || mf.right_cell().addr() == nullptr)
            //{
            //    std::cout << "btype" << static_cast<int>(mf.btype()) << std::endl;
            //    std::cout << "mf.left_cell(): " << mf.left_cell().tag()() << std::endl;
            //    std::cout << "mf.right_cell(): " << mf.right_cell().tag()() << std::endl;
            //    mesh.print_as_vtk("ard.vtk");
            //}
            //assert(mesh.query(mf.left_cell()) != nullptr);
            //assert(mf.left_cell().addr() != nullptr);
            //assert(mesh.query(mf.right_cell()) != nullptr);
            //assert(mf.right_cell().addr() != nullptr);
            //LC = &mesh.cell_p(mf.left_cell());
            LC = &mesh.cell_p(mf.right_cell());
            //LC = mf.left_cell().addr();
            //RC = &mesh.cell_p(mf.right_cell());
            RC = &mesh.cell_p(mf.left_cell());
            //RC = mf.right_cell().addr();
            if (RC->oga_cell_type() == OGA_cell_type_t::non_resident)
            {
                std::cout << "LC oga: " << static_cast<int>(LC->oga_cell_type()) << std::endl;
                std::cout << "RC oga: " << static_cast<int>(RC->oga_cell_type()) << std::endl;
            }
            assert(RC->oga_cell_type() != OGA_cell_type_t::non_resident);
        }
        else
        {
            assert(mf.parent_cell().size() == 2);
            //if (mesh.query(mf.left_cell()) == nullptr || mesh.query(mf.right_cell()) == nullptr)
            //if (mf.left_cell().addr() == nullptr || mf.right_cell().addr() == nullptr)
            //{
            //    std::cout << "btype" << static_cast<int>(mf.btype()) << std::endl;
            //    std::cout << "mf.left_cell(): " << mf.left_cell().tag()() << std::endl;
            //    std::cout << "mf.right_cell(): " << mf.right_cell().tag()() << std::endl;
            //    mesh.print_as_vtk("ard.vtk");
            //}

            //LC = &mesh.cell_p(mf.left_cell());
            LC = &mesh.cell_p(mf.right_cell());
            //LC = mf.left_cell().addr();
            //RC = &mesh.cell_p(mf.right_cell());
            RC = &mesh.cell_p(mf.left_cell());
            //RC = mf.right_cell().addr();

        }

        assert(LC->tag() == mctag);
        assert(RC->tag() != mctag);

        return std::make_tuple(LC, RC);
    }

    void left_right_cells(const Mesh& mesh, const MeshFace& mf, const MeshCell*& LC, const MeshCell*& RC)
    {
        if (mf.is_boundary())
        {
            LC = &mesh.cell(mf.left_cell());
            //LC = mf.left_cell().const_addr();
            RC = mesh.boundary(mf.btype(), mf.right_cell());
            //RC = mf.right_cell().const_addr();
        }
        else if (mf.btype() == BouType::partition)
        {
            //RC = &mesh.cell(mf.right_cell());
            LC = &mesh.cell(mf.right_cell());
            RC = &mesh.cell(mf.left_cell());
            //RC = mf.right_cell().const_addr();
        }
        else
        {
            //LC = &mesh.cell(mf.left_cell());
            LC = &mesh.cell(mf.right_cell());
            //LC = mf.left_cell().const_addr();
            //RC = &mesh.cell(mf.right_cell());
            RC = &mesh.cell(mf.left_cell());
            //RC = mf.right_cell().const_addr();
        }
    }

    const MeshCell* opposing_nei(const Mesh& mesh, const MeshFace& mf, const Tag& me)
    {
        if (mf.is_boundary())
        {
            const MeshCell* LC = &mesh.cell(mf.left_cell());
            //const MeshCell* LC = mf.left_cell().const_addr();
            const MeshCell* RC = mesh.boundary(mf.btype(), mf.right_cell());
            //const MeshCell* RC = mf.right_cell().const_addr();

            return RC;

        }
        else
        {
            const MeshCell* LC = &mesh.cell(mf.left_cell());
            //const MeshCell* LC = mf.left_cell().const_addr();
            const MeshCell* RC = &mesh.cell(mf.right_cell());
            //const MeshCell* RC = mf.right_cell().const_addr();

            if (me == LC->tag())
            {
                return RC;
            }
            else if (me == RC->tag())
            {
                return LC;
            }
            else
            {
                assert(false);
            }
        }
    }

    bool Mesh::is_gcl_satisfied(int rank) const
    {
        for (const MeshCell& mc: cell())
        {
            if (mc.oga_cell_type() == OGA_cell_type_t::non_resident || mc.oga_cell_type() == OGA_cell_type_t::ghost) {
                continue;
            }
            double sum = 0.;

            //for (const MeshFace& cmf: mc.face())
            for (const MeshFace& mf: mc.face())
            {
                //if (!cmf.tag().isvalid())
                //{
                    //std::cout << "mc: " << mc.tag()() << std::endl;
                    //std::cout << "oga type: " << static_cast<int>(mc.oga_cell_type()) << std::endl;
                    //std::cout << "btype: " << static_cast<int>(mc.btype()) << std::endl;
                    //std::cout << "face size: " << mc.face().size() << std::endl;
                    //std::cout << "pnei size: " << mc.pnei().size() << std::endl;
                    //std::cout << std::abs(sum) << std::endl;
                    //print_as_vtk("nnn.vtk");
                //}
                //assert(cmf.tag().isvalid());
                //const MeshFace* mf;
                //if (cmf.btype() == BouType::partition)
                //{
                    //mf = &cmf;
                //}
                //else
                //{
                    //auto it = std::find_if(face_.begin(), face_.end(), [&](const auto& f){return f.tag() == cmf.tag();});
                    //if (it == face_.end())
                    //{
                        //std::cout << "mc: " << mc.tag()() << std::endl;
                        //std::cout << "oga type: " << static_cast<int>(mc.oga_cell_type()) << std::endl;
                        //std::cout << "mc btype: " << static_cast<int>(mc.btype()) << std::endl;
                        //std::cout << "cmf btype: " << static_cast<int>(cmf.btype()) << std::endl;
                    //}
                    //assert(it != face_.end());
                    //mf = &face(cmf.tag());
                //}

                //double dp = mf->vgn() * mf->face().signed_area();
                //double dp = mf.vgn() * mf.face().signed_area();
                //double dp = mf.vgn().len() * mf.face().signed_area();
                //double dp = dotp(mf.vgn(), mf.face().normal()) * mf.face().signed_area();
                double dp = mf.vgn() * mf.face().signed_area();

                //if (dotp(mf.face().normal(), cmf.face().normal()) >= 0.)
                //{
                    sum += dp;
                //}
                //else
                //{
                    //sum -= dp;
                //}
            }

            //if (rank == 0)
            {
                //if (std::abs(sum) > TAILOR_ZERO)
                if (std::abs(sum) > 1e-5)
                {
                    std::cout << "GCL was not satisfied" << std::endl;
                    std::cout << "mc: " << mc.tag()() << std::endl;
                    std::cout << "oga type: " << static_cast<int>(mc.oga_cell_type()) << std::endl;
                    std::cout << "btype: " << static_cast<int>(mc.btype()) << std::endl;
                    std::cout << "face size: " << mc.face().size() << std::endl;
                    std::cout << "pnei size: " << mc.pnei().size() << std::endl;
                    std::cout << "sum: " << std::abs(sum) << std::endl;
                    print_as_vtk("nnn.vtk");

                    for (const MeshPoint& mp: mc.point())
                    {
                        std::cout << "mp r0: " << mp.p().r(0) << std::endl;
                        std::cout << "mp r1: " << mp.p().r(1) << std::endl;
                        std::cout << "mp r2: " << mp.p().r(2) << std::endl;
                    }

                    //for (const MeshFace& cmf: mc.face())
                    for (const MeshFace& mf: mc.face())
                    {
                        //const MeshFace& mf = face(cmf.tag());
                        std::cout << "normal: " << mf.face().normal()(0) << " " << mf.face().normal()(1) << " " << mf.face().normal()(2) << std::endl;
                        std::cout << "area: " << mf.face().signed_area() << std::endl;
                        //std::cout << "dotp: " << dotp(mf.vgn(), mf.face().normal()) << std::endl;
                        std::cout << "dotp: " << mf.vgn() << std::endl;
                        //std::cout << "full: " << dotp(mf.vgn(), mf.face().normal()) * mf.face().signed_area() << std::endl;
                        std::cout << "full: " << mf.vgn() * mf.face().signed_area() << std::endl;
                    }

                    //assert(false); // TODO should be made per cell basis.
                    //return false;
                    return true;
                }
            }
        }

        return true;
    }

    void Mesh::init_flow()
    {
        // Read init file for this mesh.
        // init accordingly.
        
        FlowInit finit;
        finit.read(tag_);

        Freestream fs;
        fs.read();

        if (finit.type == "uniform")
        {
             init_uniform(finit, fs.gamma_);
        }
        else if (finit.type == "gaussian")
        {
            init_gaussian(finit, fs.gamma_);
        }
        else if (finit.type == "xsplit")
        {
            init_xsplit(finit, fs.gamma_);
        }
        else
        {
            assert(false);
        }
    }

    void Mesh::init_xsplit(const FlowInit& finit, double gamma)
    {
        Vector5 priml;
        priml(0) = finit.rhol;
        priml(1) = finit.ul;
        priml(2) = 0.;
        priml(3) = 0.;
        priml(4) = finit.pl;

        Vector5 primr;
        primr(0) = finit.rhor;
        primr(1) = finit.ur;
        primr(2) = 0.;
        primr(3) = 0.;
        primr(4) = finit.pr;

        double x = finit.x;

        for (auto& mc: cell_)
        {
            if (mc.poly().centroid()(0) <= x)
            {
                mc.prim_ = priml;
            }
            else
            {
                mc.prim_ = primr;
            }

            mc.cons_sp1_ = prim_to_cons(mc.prim_, gamma);
        }

        for (auto& mc: dirichlet_boundaries_)
        {
            if (mc.poly().centroid()(0) <= x)
            {
                mc.prim_ = priml;
            }
            else
            {
                mc.prim_ = primr;
            }

            mc.cons_sp1_ = prim_to_cons(mc.prim_, gamma);
        }
    }

    void Mesh::init_gaussian(const FlowInit& finit, double gamma)
    {
        //GaussianInit ginit;
        //ginit.read();

        auto uni_prim = uniform_prim(finit);

        double b = finit.strength;
        Vector3 vcnt;
        vcnt(0) = finit.cnt_x;
        vcnt(1) = finit.cnt_y;
        vcnt(2) = finit.cnt_z;

        double rhoinf = uni_prim(0);
        double uinf = uni_prim(1);
        double vinf = uni_prim(2);
        double winf = uni_prim(3);
        double pinf = uni_prim(4);

        for (auto& mc: cell_)
        {
            const auto& cnt = mc.poly().centroid();
            double rsq = std::pow(cnt(0) - vcnt(0), 2.) + std::pow(cnt(1) - vcnt(1), 2.);

            double rho = std::pow(rhoinf - (gamma - 1.) * b * b * std::exp(1 - rsq) / (8. * gamma * PI * PI), (1./(gamma - 1.)));
            double p = std::pow(rho, gamma);

            double a1 = b * std::exp(0.5 * (1. - rsq)) / (2. * PI);
            double u = uinf - a1 * (cnt(1) - vcnt(1));
            double v = vinf + a1 * (cnt(0) - vcnt(0));

            Vector5 prim;
            prim(0) = rho;
            prim(1) = u;
            prim(2) = v;
            prim(3) = 0.;
            prim(4) = p;

            mc.set_prim_cons(prim, gamma);
        }

        for (auto& mc: dirichlet_boundaries_)
        {
            mc.set_prim_cons(uni_prim, gamma);
        }
    }

    //Vector5 Mesh::uniform_prim(const Freestream& fs)
    //{
    //    double cinf = std::sqrt(fs.gamma_ * fs.pinf_ / fs.rhoinf_);
    //    Vector3 vinf_air;

    //    if (fs.velair_ != 0.)
    //    {
    //        vinf_air = Vector3(
    //                fs.velair_ * std::cos(deg_to_rad(fs.aoa_air_x_)),
    //                fs.velair_ * std::cos(deg_to_rad(90. - fs.aoa_air_x_)),
    //                fs.velair_ * std::cos(deg_to_rad(fs.aoa_air_z_)));
    //    }
    //    else
    //    {
    //        vinf_air = Vector3(
    //                fs.machair_ * cinf * std::cos(deg_to_rad(fs.aoa_air_x_)),
    //                fs.machair_ * cinf * std::cos(deg_to_rad(90. - fs.aoa_air_x_)),
    //                fs.machair_ * cinf * std::cos(deg_to_rad(fs.aoa_air_z_)));
    //    }

    //    Vector5 prim;
    //    prim(0) = fs.rhoinf_;
    //    prim(1) = vinf_air(0);
    //    prim(2) = vinf_air(1);
    //    prim(3) = vinf_air(2);
    //    prim(4) = fs.pinf_;

    //    return prim;
    //}

    Vector5 Mesh::uniform_prim(const FlowInit& finit)
    {
        Vector5 prim;
        prim(0) = finit.rho;
        prim(1) = finit.u;
        prim(2) = finit.v;
        prim(3) = finit.w;
        prim(4) = finit.p;

        return prim;
    }

    void Mesh::init_uniform(const FlowInit& finit, double gamma)
    {
        auto prim = uniform_prim(finit);

        for (auto& mc: cell_)
        {
            //mc.init(vinf_air, fs, compo);
            mc.set_prim_cons(prim, gamma);
        }

        for (auto& mc: dirichlet_boundaries_)
        {
            //mc.init(vinf_air, fs, compo);
            mc.set_prim_cons(prim, gamma);
        }
    }

    //void Mesh::calc_face_velocities(double dt, const Freestream& fs, int rank)
    void Mesh::calc_mesh_velocities(const Freestream& fs, int rank, double real_time)
    {
        Component compo;
        compo.read(tag_);

        for (auto& mc: cell_)
        {
            //mc.mesh_velocity(fs, compo);

            for (auto& mf: mc.face_p())
            {
                mf.face_velocity(fs, compo, real_time);
            }
        }

        assert(is_gcl_satisfied(rank));
    }

    void Mesh::set_prim_cons(BouType btype, const Vector5& prim, double gamma)
    {
        if (btype == BouType::interior)
        {
            for (MeshCell& mc: cell_)
            {
                mc.prim_ = prim;
                mc.cons_sp1_ = prim_to_cons(prim, gamma);
            }
        }
        else if (btype == BouType::dirichlet)
        {
            for (MeshCell& mc: dirichlet_boundaries_)
            {
                mc.prim_ = prim;
                mc.cons_sp1_ = prim_to_cons(prim, gamma);
            }
        }
        else if (btype == BouType::farfield)
        {
            for (MeshCell& mc: dirichlet_boundaries_)
            {
                mc.prim_ = prim;
                mc.cons_sp1_ = prim_to_cons(prim, gamma);
            }
        }
        else
        {
            assert(false);
        }
    }

    void Mesh::reset_R_checked()
    {
        for (MeshCell& mc: cell_)
        {
            for (auto& mf: mc.face_)
            {
                mf.set_R_checked(false);
            }
        }
    }

    void Mesh::reset_R()
    {
        for (MeshCell& mc: cell_)
        {
            mc.R_ = 0.;
            mc.R_mid_ = 0.;
        }
    }

    void Mesh::equate_RD_to_RDmid()
    {
        for (MeshCell& mc: cell_)
        {
            mc.R_ = mc.R_mid_;
            mc.D_ = mc.D_mid_;
        }
    }

    void Mesh::reset_D()
    {
        for (MeshCell& mc: cell_)
        {
            mc.D_ = 0.;
            mc.D_mid_ = 0.;
        }
    }

    bool are_common_faces(const MeshFace& a, const MeshFace& b)
    {
        if (a.mesh_point().size() != b.mesh_point().size()) {
            return false;
        }

        int count = 0;
        for (const Tag& p: a.mesh_point())
        {
            for (const Tag& op: b.mesh_point())
            {
                if (p == op) {
                    ++count;
                    break;
                }
            }
        }

        assert(a.mesh_point().size() == b.mesh_point().size());

        if (count == a.mesh_point().size())
        {
            return true;
        }

        return false;
    }

    MeshFace* Mesh::common_face(MeshCell& mc, const FaceTag& ft)
    {
        auto it = std::find_if(mc.face_.begin(), mc.face_.end(), [&](const auto& f){return f.tag() == ft;});
        if (it != mc.face_.end())
        {
            return &(*it);
        }

        return nullptr;
    }

    const MeshFace* Mesh::common_face(const MeshCell& mc, const FaceTag& ft) const
    {
        auto it = std::find_if(mc.face_.begin(), mc.face_.end(), [&](const auto& f){return f.tag() == ft;});
        if (it != mc.face_.end())
        {
            return &(*it);
        }

        return nullptr;
    }

    MeshFace* Mesh::common_face(MeshCell& mc, const Tag& celltag)
    {
        for (MeshFace& mf: mc.face_)
        {
            if (is_potential_parent_cell(mf, celltag))
            {
                assert(&mf != nullptr);
                return &mf;
            }
        }


        return nullptr;
    }

    bool Mesh::is_potential_parent_cell(const MeshFace& mf, const Tag& celltag)
    {
        // checks if celltag could be a parent cell of mf.

        int counter = 0;
        for (int i=0; i<mf.mesh_point().size(); ++i)
        {
            auto pp = point(mf.mesh_point(i));
            assert(pp != nullptr);
            //const MeshPoint& p = point(mf.mesh_point(i));
            const MeshPoint& p = *pp;
            for (const auto& _t: p.parent_cell())
            {
                if (_t == celltag)
                {
                    ++counter;
                    break;
                }
            }
        }

        if (counter == mf.mesh_point().size())
        {
            return true;
        }

        return false;
    }

    void Mesh::update_prev_cand_donor()
    {
        for (MeshCell& mc: cell_)
        {
            mc.update_prev_cand_donor();
        }
    }

    void Mesh::determine_unreached_hole_cells(const AABB& cutter_hole_aabb)
    {
        // it is possible that when determining hole cells from seeds to interior hole, some cells left in undefined oga_cell_type since all neighbords of undefined cell all tagged as hole cells. use this function to determine the rest of the hole cells.

        for (MeshCell& mc: cell_)
        {
            if (mc.oga_cell_type() != OGA_cell_type_t::hole) {
                continue;
            }

            std::ofstream out;
            flood_fill(mc.tag(), cutter_hole_aabb, 0, false, out);
        }
    }
    /*void Mesh::determine_unreached_hole_cells()
    {
        // it is possible that when determining hole cells from seeds to interior hole, some cells left in undefined oga_cell_type since all neighbords of undefined cell all tagged as hole cells. use this function to determine the rest of the hole cells.

        for (MeshCell& mc: cell_)
        {
            if (mc.oga_cell_type() == OGA_cell_type_t::hole) {
                continue;
            }

            assert(mc.oga_cell_type() == OGA_cell_type_t::undefined);

            bool neis_all_hole = true;
            for (const Tag& inei: mc.pnei())
            {
                const MeshCell& nei = cell(inei);

                if (nei.oga_cell_type() != OGA_cell_type_t::hole) {
                    neis_all_hole = false;
                    break;
                }
            }

            if (!neis_all_hole) {
                continue;
            }

            mc.set_oga_cell_type(OGA_cell_type_t::hole);
        }
    }

    void Mesh::determine_unreached_hole_cells_via_aabb(const AABB& cutter_aabb)
    {
        // it is possible that when determining hole cells from seeds to interior hole, some cells left in undefined oga_cell_type since all neighbords of undefined cell all tagged as hole cells. use this function to determine the rest of the hole cells.

        for (MeshCell& mc: cell_)
        {
            if (mc.oga_cell_type() == OGA_cell_type_t::hole) {
                continue;
            }

            assert(mc.oga_cell_type() == OGA_cell_type_t::undefined);

            if (cutter_aabb.do_contain(AABB(mc.poly())))
            {
                mc.set_oga_cell_type(OGA_cell_type_t::hole);
            }
        }
    }*/

    void Mesh::print_donor_info() const
    {
        std::ofstream out;
        out.open("donor-info.dat");

        for (const MeshCell& mc: cell_)
        {
            out << tag_();
            out << " ";
            out << mc.tag()();
            out << " ";
            out << mc.donor().mesh_tag_();
            out << " ";
            out << mc.donor().cell_tag_();
        }

        out.close();
    }

    //const mfc& Mesh::face() const
    //{
        //return face_;
    //}

    void Mesh::bbox(Vector3& min_, Vector3& max_) const
    {
        min_ = TAILOR_BIG_POS_NUM;
        max_ = TAILOR_BIG_NEG_NUM;

        for (const MeshPoint& mp: point_)
        {
            min_(0) = std::min(min_(0), mp.p().r(0));
            min_(1) = std::min(min_(1), mp.p().r(1));
            min_(2) = std::min(min_(2), mp.p().r(2));

            max_(0) = std::max(max_(0), mp.p().r(0));
            max_(1) = std::max(max_(1), mp.p().r(1));
            max_(2) = std::max(max_(2), mp.p().r(2));
        }
    }

    /*size_t Mesh::mem() const
    {
        size_t size = 0;

        size += sizeof(wall_boundaries_);
        for (const auto& mc: wall_boundaries_)
        {
            size += mc.mem();
        }
        size += sizeof(dirichlet_boundaries_);
        for (const auto& mc: dirichlet_boundaries_)
        {
            size += mc.mem();
        }

        size += sizeof(cell_);
        size_t size_interior = 0;
        for (const auto& mc: cell_)
        {
            size_interior += mc.mem();
        }

        std::cout << "mem interior (MB): " << size_interior/1e6 << std::endl;
        std::cout << "ave cell mem (MB): " << size_interior/cell_.size()/1e6 << std::endl;
        std::cout << "ncell: " << cell_.size() << std::endl;

        size += size_interior;
        size += sizeof(point_);
        for (const auto& mp: point_)
        {
            size += mp.mem();
        }
        size += sizeof(int) * 2;
        size += hole_aabb_.mem();


        return size;

        //bimap_int wall_tag_index_map_; // ignored
        //bimap_int dirichlet_tag_index_map_; // ignored
    }*/

    bool Mesh::operator<(const Mesh& other) const
    {
        return tag_ < other.tag();
    }

    /*void Mesh::add_points_from_cells()
      {
      for (const MeshPoint& p: mc->point())
      {
      auto it = std::lower_bound(point_.begin(), point_.end(), p, [](const MeshPoint& left, const MeshPoint& right){return left.tag()() < right.tag()();});
      if (it == point_.end())
      { 
      point_.push_back(p);
      }
      else if (it->tag() == p.tag())
      {
      it->add_parent_cell(mc->tag());
      }
      else
      {
      point_.insert(it, p);
      }
      }
      }*/


    size_t Mesh::npartition() const
    {
        // Mesh cells should with set_partition() in connect_cells() before this function.

        size_t count = 0;
        for (const MeshCell& mc: cell_)
        {
            if (mc.btype() == BouType::partition)
            {
                ++count;
            }
        }

        return count;
    }

    /*void Mesh::batch_merge(Mesh& m)
      {
      {
      int size = cell_.size();
      for (auto mc = m.cell().rbegin(); mc != m.cell().rend(); ++mc)
      {
      auto tit = m.cell_tag_index_map.left.find(mc->tag()());
      assert(tit != m.cell_tag_index_map.left.end());
      int rep = tit->second + size;
      bool successful_replace = m.cell_tag_index_map.left.replace_data(tit, rep);
      assert(successful_replace);
      }

      cell_tag_index_map.insert(m.cell_tag_index_map.begin(), m.cell_tag_index_map.end());
      }
      {
      int size = wall_boundaries_.size();
      for (auto mc = m.wall_boundaries().rbegin(); mc != m.wall_boundaries().rend(); ++mc)
      {
      auto tit = m.wall_tag_index_map_.left.find(mc->tag()());
      assert(tit != m.wall_tag_index_map_.left.end());
      int rep = tit->second + size;
      bool successful_replace = m.wall_tag_index_map_.left.replace_data(tit, rep);
      assert(successful_replace);
      }

      wall_tag_index_map_.insert(m.wall_tag_index_map_.begin(), m.wall_tag_index_map_.end());
      }
      {
      int size = dirichlet_boundaries_.size();
      for (auto mc = m.dirichlet_boundaries().rbegin(); mc != m.dirichlet_boundaries().rend(); ++mc)
      {
      auto tit = m.dirichlet_tag_index_map_.left.find(mc->tag()());
      assert(tit != m.dirichlet_tag_index_map_.left.end());
      int rep = tit->second + size;
      bool successful_replace = m.dirichlet_tag_index_map_.left.replace_data(tit, rep);
      assert(successful_replace);
      }

      dirichlet_tag_index_map_.insert(m.dirichlet_tag_index_map_.begin(), m.dirichlet_tag_index_map_.end());
      }

      cell_.insert(cell_.end(), m.cell().begin(), m.cell().end());
      wall_boundaries_.insert(wall_boundaries_.end(), m.wall_boundaries().begin(), m.wall_boundaries().end());
      dirichlet_boundaries_.insert(dirichlet_boundaries_.end(), m.dirichlet_boundaries().begin(), m.dirichlet_boundaries().end());

      int psize = point_.size();
      point_.resize(psize + m.point().size());
      auto it = point_.begin() + psize;
      std::copy_if(m.point().begin(), m.point().end(), it, [&](const MeshPoint& mp){return query_point(mp.tag()) == nullptr;});

      for (auto it = point_.begin() + psize; it != point_.end(); ++it)
      {
      point_tag_index_map.insert(boost::bimap<int, int>::value_type(it->tag()(), std::distance(point_.begin(), it)));
      }

    // now there should not be any duplicates.
    //std::sort(cell_.begin(), cell_.end(), [](const MeshCell& left, const MeshCell& right){return left.tag() < right.tag();});
    //auto itt = std::adjacent_find(cell_.begin(), cell_.end(), [](const MeshCell& left, const MeshCell& right){return left.tag() == right.tag();});
    //assert(itt == cell_.end());
    }*/

    //const bimap_int& Mesh::dirichlet_tag_index_map() const
    //{
        //return dirichlet_tag_index_map_;
    //}

    void Mesh::destroy_cell_hood()
    {
        for (MeshCell& mc: cell_)
        {
            mc.remove_all_neighbors();
            assert(mc.pnei().empty());
        }
    }

    const mcc& Mesh::wall_boundaries() const
    {
        return wall_boundaries_;
    }

    const mcc& Mesh::symmetry_boundaries() const
    {
        return symmetry_boundaries_;
    }

    const mcc& Mesh::dirichlet_boundaries() const
    {
        return dirichlet_boundaries_;
    }

    const mcc& Mesh::farfield_boundaries() const
    {
        return farfield_boundaries_;
    }

    const mcc& Mesh::empty_boundaries() const
    {
        return empty_boundaries_;
    }

    const mcc& Mesh::interog_boundaries() const
    {
        return interog_boundaries_;
    }

    const MeshCell* Mesh::wall_boundary(const Tag& t) const
    {
        assert(t.isvalid());
        auto it = std::find_if(wall_boundaries_.begin(), wall_boundaries_.end(), [&](const MeshCell& tmc){return tmc.tag() == t;});
        //if(it == wall_boundaries_.end())
        //{
            //std::cout << "t: " << t() << std::endl;
            //std::cout << "size: " << wall_boundaries_.size() << std::endl;
            //std::cout << "mesh tag: " << tag_() << std::endl;
        //}
        if (it != wall_boundaries_.end())
        {
            return &(*it);
        }

        return nullptr;
    }

    const MeshCell* Mesh::symmetry_boundary(const Tag& t) const
    {
        assert(t.isvalid());
        auto it = std::find_if(symmetry_boundaries_.begin(), symmetry_boundaries_.end(), [&](const MeshCell& tmc){return tmc.tag() == t;});
        if (it != symmetry_boundaries_.end())
        {
            return &(*it);
        }

        return nullptr;
    }

    const MeshCell* Mesh::dirichlet_boundary(const Tag& t) const
    {
        auto it = std::find_if(dirichlet_boundaries_.begin(), dirichlet_boundaries_.end(), [&](const MeshCell& tmc){return tmc.tag() == t;});
        if (it != dirichlet_boundaries_.end())
        {
            return &(*it);
        }

        return nullptr;
    }

    const MeshCell* Mesh::farfield_boundary(const Tag& t) const
    {
        auto it = std::find_if(farfield_boundaries_.begin(), farfield_boundaries_.end(), [&](const MeshCell& tmc){return tmc.tag() == t;});
        if (it != farfield_boundaries_.end())
        {
            return &(*it);
        }

        return nullptr;
    }

    const MeshCell* Mesh::interog_boundary(const Tag& t) const
    {
        auto it = std::find_if(interog_boundaries_.begin(), interog_boundaries_.end(), [&](const MeshCell& tmc){return tmc.tag() == t;});
        if (it != interog_boundaries_.end())
        {
            return &(*it);
        }

        return nullptr;
    }

    const MeshCell* Mesh::empty_boundary(const Tag& t) const
    {
        auto it = std::find_if(empty_boundaries_.begin(), empty_boundaries_.end(), [&](const MeshCell& tmc){return tmc.tag() == t;});
        if (it != empty_boundaries_.end())
        {
            return &(*it);
        }

        return nullptr;
    }

    MeshCell* Mesh::wall_boundary_p(const Tag& t)
    {
        auto it = std::find_if(wall_boundaries_.begin(), wall_boundaries_.end(), [&](const MeshCell& tmc){return tmc.tag() == t;});
        if (it != wall_boundaries_.end())
        {
            return &(*it);
        }

        return nullptr;
    }

    MeshCell* Mesh::symmetry_boundary_p(const Tag& t)
    {
        auto it = std::find_if(symmetry_boundaries_.begin(), symmetry_boundaries_.end(), [&](const MeshCell& tmc){return tmc.tag() == t;});
        if (it != symmetry_boundaries_.end())
        {
            return &(*it);
        }

        return nullptr;
    }

    MeshCell* Mesh::dirichlet_boundary_p(const Tag& t)
    {
        auto it = std::find_if(dirichlet_boundaries_.begin(), dirichlet_boundaries_.end(), [&](const MeshCell& tmc){return tmc.tag() == t;});
        if (it != dirichlet_boundaries_.end())
        {
            return &(*it);
        }

        return nullptr;
    }

    MeshCell* Mesh::farfield_boundary_p(const Tag& t)
    {
        auto it = std::find_if(farfield_boundaries_.begin(), farfield_boundaries_.end(), [&](const MeshCell& tmc){return tmc.tag() == t;});
        if (it != farfield_boundaries_.end())
        {
            return &(*it);
        }

        return nullptr;
    }

    MeshCell* Mesh::empty_boundary_p(const Tag& t)
    {
        auto it = std::find_if(empty_boundaries_.begin(), empty_boundaries_.end(), [&](const MeshCell& tmc){return tmc.tag() == t;});
        if (it == empty_boundaries_.end())
        {
            std::cout << "tag: " << t() << std::endl;
            std::cout << "empty bou size: " << empty_boundaries_.size() << std::endl;
        }
        if (it != empty_boundaries_.end())
        {
            return &(*it);
        }

        return nullptr;
    }

    MeshCell* Mesh::interog_boundary_p(const Tag& t)
    {
        auto it = std::find_if(interog_boundaries_.begin(), interog_boundaries_.end(), [&](const MeshCell& tmc){return tmc.tag() == t;});
        if (it != interog_boundaries_.end())
        {
            return &(*it);
        }

        return nullptr;
    }

    /*void Mesh::connect_add_wall_to_interior(const Mesh& wm)
      {
      for (MeshCell mc: wm.cell())
      {
    //auto pit = std::find_if(point_.begin(), point_.end(), [&](const MeshPoint& mp){return mp.tag() == mc.point()[0].tag();});
    //auto pit = std::find_if(point_.begin(), point_.end(), [&](const MeshPoint& mp){return mp.p() == mc.point()[0].p();});
    auto pit = std::find_if(point_.begin(), point_.end(), [&](const MeshPoint& mp){return (std::abs(mp.p().r(0) - mc.point(0).p().r(0)) <= TAILOR_ZERO && std::abs(mp.p().r(1) - mc.point(0).p().r(1)) <= TAILOR_ZERO);});

    if (pit == point_.end()) continue;
    mc.set_point_tag(0, pit->tag());

    for (const Tag& pc: pit->parent_cell())
    {
    //auto pit2 = std::find_if(cell(pc).point().begin(), cell(pc).point().end(), [&](const MeshPoint& mp){return mp.tag() == mc.point()[1].tag();});
    //auto pit2 = std::find_if(cell(pc).point().begin(), cell(pc).point().end(), [&](const MeshPoint& mp){return mp.p() == mc.point(1).p();});
    auto pit2 = std::find_if(cell(pc).point().begin(), cell(pc).point().end(), [&](const MeshPoint& mp){return (std::abs(mp.p().r(0) - mc.point(1).p().r(0)) <= TAILOR_ZERO && std::abs(mp.p().r(1) - mc.point(1).p().r(1)) <= TAILOR_ZERO);});
    if (pit2 != cell(pc).point().end())
    {
    mc.set_point_tag(1, pit2->tag());
    cell_p(pc).add_wall_boundary(mc.tag());
    add_wall_boundary(std::move(mc));
    wall_boundaries_.back().set_interior_boundary(pc);
    break;
    }
    }
    }

    if (wall_boundaries_.size() == 0)
    {
    for (const MeshCell& mc: cell_)
    {
    assert(mc.wall_boundary().empty());
    }
    }
    }*/

    void Mesh::connect_add_bou_to_interior(BouType boutype, int rank)
    {
        mcc* container = nullptr;

        if (boutype == BouType::wall)
        {
            container = &wall_boundaries_;
        }
        else if (boutype == BouType::symmetry)
        {
            container = &symmetry_boundaries_;
        }
        else if (boutype == BouType::dirichlet)
        {
            container = &dirichlet_boundaries_;
        }
        else if (boutype == BouType::farfield)
        {
            container = &farfield_boundaries_;
        }
        else if (boutype == BouType::empty)
        {
            container = &empty_boundaries_;
        }
        else if (boutype == BouType::interog)
        {
            container = &interog_boundaries_;
        }
        else
        {
            assert(false);
        }

        for (MeshCell& cll: *container)
        {
            assert(cll.btype() == boutype);

            auto pit1 = std::find_if(point_.begin(), point_.end(), [&](const MeshPoint& mp){return mp.tag() == cll.point(0).tag();});
            //auto pit1 = std::find_if(point_.begin(), point_.end(), [&](const MeshPoint& mp){return (std::abs(mp.p().r(0) - cll.point(0).p().r(0)) <= TAILOR_ZERO && std::abs(mp.p().r(1) - cll.point(0).p().r(1)) <= TAILOR_ZERO && std::abs(mp.p().r(2) - cll.point(0).p().r(2)) <= TAILOR_ZERO);});
            if (pit1 == point_.end()) continue;

            for (auto& pc: pit1->parent_cell())
            {
                auto& cellpc = cell_p(pc);
                                //if (tag_() == 0 && cellpc.tag()() == 33711)
                                //{
                                    //std::cout << "cellpc face size: " << cellpc.face_.size() << std::endl;
                                    //for (MeshFace& mf: cellpc.face_)
                                    //{
                                    //    for (const Tag& mp: mf.mesh_point())
                                    //    {
                                    //        std::cout << "mpp: " << mp() << std::endl;
                                    //    }
                                    //}
                                    //assert(false);
                                //}
                auto pit2 = std::find_if(cellpc.point().begin(), cellpc.point().end(), [&](const MeshPoint& mp){return mp.tag() == cll.point(1).tag();});
                //auto pit2 = std::find_if(cellpc.point().begin(), cellpc.point().end(), [&](const MeshPoint& mp){return mp.p() == cll.point(1).p();});
                if (pit2 != cellpc.point().end())
                {
                    auto pit3 = std::find_if(cellpc.point().begin(), cellpc.point().end(), [&](const MeshPoint& mp){return mp.tag() == cll.point(2).tag();});
                    //auto pit3 = std::find_if(cellpc.point().begin(), cellpc.point().end(), [&](const MeshPoint& mp){return mp.p() == cll.point(2).p();});
                    if (pit3 != cellpc.point().end())
                    {
                        cll.set_point_tag(0, pit1->tag());
                        cll.set_point_tag(1, pit2->tag());
                        cll.set_point_tag(2, pit3->tag());

                        if (cll.point().size() > 3)
                        {
                            for (int i=3; i<cll.point().size(); ++i)
                            {
                                auto pit = std::find_if(cellpc.point().begin(), cellpc.point().end(), [&](const MeshPoint& mp){return mp.p() == cll.point(i).p();});
                                cll.set_point_tag(i, pit->tag());
                            }
                        }
                        cellpc.add_boundary(cll.tag(), boutype);

                        assert(cll.face_.size() == 1);
                        auto& facezero = cll.face_[0];
                        assert(facezero.parent_cell().size() == 0);
                        facezero.add_parent_cell(cll.tag());
                        assert(facezero.parent_cell().size() <= 2);
                        facezero.set_btype(boutype);

                        int fcount = 0;
                        for (MeshFace& mf: cellpc.face_)
                        {
                            int countt = 0;
                            for (const Tag& mp: mf.mesh_point())
                            {
                                if (mp == pit1->tag() || mp == pit2->tag() || mp == pit3->tag())
                                {
                                    ++countt;
                                }
                            }

                            if (countt >= 3)
                            {
                                auto ft = gen_face_tag(mf);
                                mf.set_tag(ft);
                                assert(mf.tag().isvalid());
                                facezero.set_tag(ft); // temporary face tag to connect face and wall face.
                                mf.add_parent_cell(cll.tag());
                                assert(mf.parent_cell().size() <= 2);
                                mf.add_parent_cell(pc);
                                assert(pc.isvalid());
                                mf.set_btype(boutype);
                                facezero.add_parent_cell(pc);
                                cll.set_interior_boundary(cellpc.tag());
                                assert(mf.parent_cell().size() <= 2);

                                assert(mf.tag().isvalid());
                                break;
                            }
                            ++fcount;
                        }

                        break;
                    }
                }
            }
        }

        for (const auto& mc: cell_)
        {
            for (const auto& mf: mc.face())
            {
                if (mf.is_boundary())
                {
                    if (mf.parent_cell().size() != 2)
                    {
                        //std::cout << "rank: " << rank << std::endl;
                        std::cout << "mf: " << mf.tag() << std::endl;
                        std::cout << "btype: " << static_cast<int>(mf.btype()) << std::endl;
                        std::cout << "mf pc size: " << mf.parent_cell().size() << std::endl;
                    }
                    assert(mf.parent_cell().size() == 2);
                }
            }
        }
    }

    void Mesh::connect_add_bou_to_interior(Mesh& boumesh, BouType boutype, int rank)
    {
        mcc* container = nullptr;

        if (boutype == BouType::wall)
        {
            container = &wall_boundaries_;
        }
        else if (boutype == BouType::symmetry)
        {
            container = &symmetry_boundaries_;
        }
        else if (boutype == BouType::dirichlet)
        {
            container = &dirichlet_boundaries_;
        }
        else if (boutype == BouType::farfield)
        {
            container = &farfield_boundaries_;
        }
        else if (boutype == BouType::empty)
        {
            container = &empty_boundaries_;
        }
        else if (boutype == BouType::interog)
        {
            container = &interog_boundaries_;
        }
        else
        {
            assert(false);
        }

        for (MeshCell mc: boumesh.cell())
        {
        //{
        //    auto& ff = mc.face_[0];
        //    auto ft = gen_face_tag(ff);
        //    {
        //        if (tag_() == 0)
        //        {
        //            {
        //                {
        //                    if (ft == FaceTag(std::vector<int>{10043,10044,10179}))
        //                    {
        //                        std::cout << "mccccccccccccccccc: " << mc.tag()() << std::endl;
        //                        assert(false);
        //                    }
        //                }
        //            }
        //        }
        //    }
        //}
            assert(mc.btype() == boutype);

            auto pit1 = std::find_if(point_.begin(), point_.end(), [&](const MeshPoint& mp){return (std::abs(mp.p().r(0) - mc.point(0).p().r(0)) <= TAILOR_ZERO && std::abs(mp.p().r(1) - mc.point(0).p().r(1)) <= TAILOR_ZERO && std::abs(mp.p().r(2) - mc.point(0).p().r(2)) <= TAILOR_ZERO);});
            if (pit1 == point_.end()) continue;

            for (auto& pc: pit1->parent_cell())
            {
                auto& cellpc = cell_p(pc);
                auto pit2 = std::find_if(cellpc.point().begin(), cellpc.point().end(), [&](const MeshPoint& mp){return mp.p() == mc.point(1).p();});
                if (pit2 != cellpc.point().end())
                {
                    auto pit3 = std::find_if(cellpc.point().begin(), cellpc.point().end(), [&](const MeshPoint& mp){return mp.p() == mc.point(2).p();});
                    if (pit3 != cellpc.point().end())
                    {
                        mc.set_point_tag(0, pit1->tag());
                        mc.set_point_tag(1, pit2->tag());
                        mc.set_point_tag(2, pit3->tag());

                        if (mc.point().size() > 3)
                        {
                            for (int i=3; i<mc.point().size(); ++i)
                            {
                                auto pit = std::find_if(cellpc.point().begin(), cellpc.point().end(), [&](const MeshPoint& mp){return mp.p() == mc.point(i).p();});
                                mc.set_point_tag(i, pit->tag());
                            }
                        }
                        cellpc.add_boundary(mc.tag(), boutype);
                        add_element(mc);
                        assert(!container->empty());

                        MeshCell& cll = container->back();
                        assert(cll.face_.size() == 1);
                        auto& facezero = cll.face_[0];
                        assert(facezero.parent_cell().size() == 0);
                        facezero.add_parent_cell(cll.tag());
                        assert(facezero.parent_cell().size() <= 2);
                        facezero.set_btype(boutype);

                        int fcount = 0;
                        for (MeshFace& mf: cellpc.face_)
                        {
                            int countt = 0;
                            for (const Tag& mp: mf.mesh_point())
                            {
                                if (mp == pit1->tag() || mp == pit2->tag() || mp == pit3->tag())
                                {
                                    ++countt;
                                }
                            }

                            if (countt >= 3)
                            {
                                auto ft = gen_face_tag(mf);
                                mf.set_tag(ft);
                                assert(mf.tag().isvalid());
                                facezero.set_tag(ft); // temporary face tag to connect face and wall face.
                                mf.add_parent_cell(cll.tag());
                                assert(mf.parent_cell().size() <= 2);
                                mf.add_parent_cell(pc);
                                assert(pc.isvalid());
                                mf.set_btype(boutype);
                                facezero.add_parent_cell(pc);
                                cll.set_interior_boundary(cellpc.tag());
                                assert(mf.parent_cell().size() <= 2);
                                break;
                            }
                            ++fcount;
                        }

                        break;
                    }
                }
            }
        }

        for (const auto& mc: cell_)
        {
            for (const auto& mf: mc.face())
            {
                if (mf.is_boundary())
                {
                    if (mf.parent_cell().size() != 2)
                    {
                        std::cout << "rank: " << rank << std::endl;
                        std::cout << "mf: " << mf.tag() << std::endl;
                        std::cout << "btype: " << static_cast<int>(mf.btype()) << std::endl;
                        std::cout << "mf pc size: " << mf.parent_cell().size() << std::endl;
                    }
                    assert(mf.parent_cell().size() == 2);
                }
            }
        }
    }

    /*void Mesh::connect_add_wall_to_interior(const Mesh& wm)
      {
      for (MeshCell mc: wm.cell())
      {
      auto pit1 = std::find_if(point_.begin(), point_.end(), [&](const MeshPoint& mp){return (std::abs(mp.p().r(0) - mc.point(0).p().r(0)) <= TAILOR_ZERO && std::abs(mp.p().r(1) - mc.point(0).p().r(1)) <= TAILOR_ZERO);});
      if (pit1 == point_.end()) continue;

      for (const Tag& pc: pit1->parent_cell())
      {
      auto pit2 = std::find_if(cell(pc).point().begin(), cell(pc).point().end(), [&](const MeshPoint& mp){return mp.p() == mc.point(1).p();});
      if (pit2 != cell(pc).point().end())
      {
      auto pit3 = std::find_if(cell(pc).point().begin(), cell(pc).point().end(), [&](const MeshPoint& mp){return mp.p() == mc.point(2).p();});
      if (pit3 != cell(pc).point().end())
      {
      auto pit4 = std::find_if(cell(pc).point().begin(), cell(pc).point().end(), [&](const MeshPoint& mp){return mp.p() == mc.point(3).p();});

      mc.set_point_tag(0, pit1->tag());
      mc.set_point_tag(1, pit2->tag());
      mc.set_point_tag(2, pit3->tag());
      mc.set_point_tag(3, pit4->tag());

      assert(cell_p(pc).wall_boundary().size() < 6);
      cell_p(pc).add_wall_boundary(mc.tag());
      add_wall_boundary(std::move(mc));
      wall_boundaries_.back().set_interior_boundary(pc);

      int countt = 0;
      for (MeshFace& mf: cell_p(pc).face_)
      {
      for (const Tag& mp: mf.mesh_point())
      {
      if (mp == pit1->tag() || mp == pit2->tag() || mp == pit3->tag())
      {
      ++countt;
      }
      }

      if (countt >= 3)
      {
      mf.add_parent_cell(mc.tag());
      mf.set_as_boundary();
      mf.set_btype(boundary_t::wall);
      }
      }

      break;
      }
      }
      }
      }
      }*/

    /*void Mesh::connect_add_dirichlet_to_interior(const Mesh& wm, std::string dummyfilename)
      {
      for (MeshCell mc: wm.cell())
      {
      auto pit1 = std::find_if(point_.begin(), point_.end(), [&](const MeshPoint& mp){return (std::abs(mp.p().r(0) - mc.point(0).p().r(0)) <= TAILOR_ZERO && std::abs(mp.p().r(1) - mc.point(0).p().r(1)) <= TAILOR_ZERO);});
      if (pit1 == point_.end()) continue;

      for (const Tag& pc: pit1->parent_cell())
      {
      auto pit2 = std::find_if(cell(pc).point().begin(), cell(pc).point().end(), [&](const MeshPoint& mp){return mp.p() == mc.point(1).p();});
      if (pit2 != cell(pc).point().end())
      {
      auto pit3 = std::find_if(cell(pc).point().begin(), cell(pc).point().end(), [&](const MeshPoint& mp){return mp.p() == mc.point(2).p();});
      if (pit3 != cell(pc).point().end())
      {
      auto pit4 = std::find_if(cell(pc).point().begin(), cell(pc).point().end(), [&](const MeshPoint& mp){return mp.p() == mc.point(3).p();});

      mc.set_point_tag(0, pit1->tag());
      mc.set_point_tag(1, pit2->tag());
      mc.set_point_tag(2, pit3->tag());
      mc.set_point_tag(3, pit4->tag());

      assert(cell_p(pc).dirichlet_boundary().size() < 6);
      cell_p(pc).add_dirichlet_boundary(mc.tag());
      add_dirichlet_boundary(std::move(mc));
      dirichlet_boundaries_.back().set_interior_boundary(pc);

      int countt = 0;
      for (MeshFace& mf: cell_p(pc).face_)
      {
      for (const Tag& mp: mf.mesh_point())
      {
      if (mp == pit1->tag() || mp == pit2->tag() || mp == pit3->tag())
      {
      ++countt;
      }
      }

      if (countt >= 3)
      {
      mf.add_parent_cell(mc.tag());
      mf.set_as_boundary();
      mf.set_btype(boundary_t::dirichlet);
      }
      }

      break;
      }
      }
      }
      }
      }*/

    /*void Mesh::connect_add_dirichlet_to_interior(const Mesh& wm)
      {
      for (MeshCell mc: wm.cell())
      {
    //auto pit = std::find_if(point_.begin(), point_.end(), [&](const MeshPoint& mp){return mp.tag() == mc.point()[0].tag();});
    //auto pit = std::find_if(point_.begin(), point_.end(), [&](const MeshPoint& mp){return mp.p() == mc.point()[0].p();});
    auto pit = std::find_if(point_.begin(), point_.end(), [&](const MeshPoint& mp){return (std::abs(mp.p().r(0) - mc.point(0).p().r(0)) <= TAILOR_ZERO && std::abs(mp.p().r(1) - mc.point(0).p().r(1)) <= TAILOR_ZERO);});
    if (pit == point_.end()) continue;
    mc.set_point_tag(0, pit->tag());

    for (const Tag& pc: pit->parent_cell())
    {
    //int count = std::count_if(cell(pc).point().begin(), cell(pc).point().end(), [&](const MeshPoint& mp){return mp.tag() == mc.point()[1].tag();});
    auto pit2 = std::find_if(cell(pc).point().begin(), cell(pc).point().end(), [&](const MeshPoint& mp){return mp.p() == mc.point(1).p();});
    if (pit2 != cell(pc).point().end())
    {
    mc.set_point_tag(1, pit2->tag());
    if (cell_p(pc).dirichlet_boundary().size() >= 6)
    {
    std::cout << "dirichlet boundary size: " << cell_p(pc).dirichlet_boundary().size() << std::endl;
    }
    assert(cell_p(pc).dirichlet_boundary().size() < 6);
    cell_p(pc).add_dirichlet_boundary(mc.tag());
    add_dirichlet_boundary(std::move(mc));
    dirichlet_boundaries_.back().set_interior_boundary(pc);
    break;
    }
    }
    }
    }*/

    void Mesh::add_element(MeshCell&& mc)
    {
        if (mc.btype() == BouType::wall)
        {
            add_wall_boundary(std::move(mc));
        }
        else if (mc.btype() == BouType::symmetry)
        {
            add_symmetry_boundary(std::move(mc));
        }
        else if (mc.btype() == BouType::dirichlet)
        {
            add_dirichlet_boundary(std::move(mc));
        }
        else if (mc.btype() == BouType::farfield)
        {
            add_farfield_boundary(std::move(mc));
        }
        else if (mc.btype() == BouType::empty)
        {
            add_empty_boundary(std::move(mc));
        }
        else if (mc.btype() == BouType::interior)
        {
            add_interior_cell(std::move(mc));
        }
        else if (mc.btype() == BouType::interog)
        {
            add_interog_boundary(std::move(mc));
        }
        else
        {
            assert(false);
        }
    }

    void Mesh::add_element(const MeshCell& mc)
    {
        if (mc.btype() == BouType::wall)
        {
            add_wall_boundary(std::move(mc));
        }
        else if (mc.btype() == BouType::symmetry)
        {
            add_symmetry_boundary(std::move(mc));
        }
        else if (mc.btype() == BouType::dirichlet)
        {
            add_dirichlet_boundary(std::move(mc));
        }
        else if (mc.btype() == BouType::farfield)
        {
            add_farfield_boundary(std::move(mc));
        }
        else if (mc.btype() == BouType::empty)
        {
            add_empty_boundary(std::move(mc));
        }
        else if (mc.btype() == BouType::interog)
        {
            add_interog_boundary(std::move(mc));
        }
        else if (mc.btype() == BouType::interior)
        {
            add_interior_cell(std::move(mc));
        }
        else
        {
            assert(false);
        }
    }

    void Mesh::add_element_nonsorted(const MeshCell& mc)
    {
        if (mc.btype() == BouType::wall)
        {
            add_wall_boundary(mc);
        }
        else if (mc.btype() == BouType::symmetry)
        {
            add_symmetry_boundary(mc);
        }
        else if (mc.btype() == BouType::dirichlet)
        {
            add_dirichlet_boundary(mc);
        }
        else if (mc.btype() == BouType::farfield)
        {
            add_farfield_boundary(mc);
        }
        else if (mc.btype() == BouType::empty)
        {
            add_empty_boundary(mc);
        }
        else if (mc.btype() == BouType::interog)
        {
            add_interog_boundary(mc);
        }
        else if (mc.btype() == BouType::interior)
        {
            add_interior_cell_nonsorted(mc);
        }
        else
        {
            assert(false);
        }
    }

    void Mesh::remove_merge_duplicate_points()
    {
        remove_merge_dup(point_);
    }

    void Mesh::add_element_no_check(const MeshCell& mc)
    {
        if (mc.btype() == BouType::wall)
        {
            add_wall_boundary(mc);
        }
        else if (mc.btype() == BouType::symmetry)
        {
            add_symmetry_boundary(mc);
        }
        else if (mc.btype() == BouType::dirichlet)
        {
            add_dirichlet_boundary(mc);
        }
        else if (mc.btype() == BouType::farfield)
        {
            add_farfield_boundary(mc);
        }
        else if (mc.btype() == BouType::empty)
        {
            add_empty_boundary(mc);
        }
        else if (mc.btype() == BouType::interog)
        {
            add_interog_boundary(mc);
        }
        else if (mc.btype() == BouType::interior)
        {
            add_interior_cell_no_check(mc);
        }
        else
        {
            assert(false);
        }
    }

    void Mesh::add_bous(mcc& cells, BouType type)
    {
        mcc* container;
        bimap_int* bimap;
        if (type == BouType::wall) {
            container = &wall_boundaries_;
            //bimap = &wall_tag_index_map_;
        }
        else if (type == BouType::symmetry) {
            container = &symmetry_boundaries_;
        }
        else if (type == BouType::dirichlet) {
            container = &dirichlet_boundaries_;
            //bimap = &dirichlet_tag_index_map_;
        }
        else if (type == BouType::farfield) {
            container = &farfield_boundaries_;
            //bimap = &farfield_tag_index_map_;
        }
        else if (type == BouType::empty) {
            container = &empty_boundaries_;
            //bimap = &empty_tag_index_map_;
        }
        else if (type == BouType::interog) {
            container = &interog_boundaries_;
            //bimap = &interog_tag_index_map_;
        }
        else
        {
            assert(false);
        }

        //assert(container->size() == bimap->left.size());
        //for (int i=0; i<cells.size(); ++i)
        //{
        //    assert(std::count_if(container->begin(), container->end(), [&](const MeshCell& mc){return mc.tag() == cells[i].tag();}) == 0);
        //    assert(bimap->left.count(cells[i].tag()()) == 0);
        //}
        for (auto mc = cells.begin(); mc != cells.end(); ++mc)
        {
            //if (bimap->left.count(cells[i].tag()()) != 0) {
                //continue;
            //}
            if (std::count_if(container->begin(), container->end(), [&](const auto& a){return a.tag() == mc->tag();}) != 0) {
                continue;
            }
            //assert(bimap->right.count(container->size() + i) == 0);
            //bimap->insert(boost::bimap<int, int>::value_type(cells[i].tag()(), container->size()));
            container->push_back(*mc);
        }
        //assert(container->size() == bimap->left.size());
    }

    /*void Mesh::add_walls(std::vector<MeshCell>& cells)
    {
        assert(wall_boundaries_.size() == wall_tag_index_map_.left.size());
        for (int i=0; i<cells.size(); ++i)
        {
            assert(std::count_if(wall_boundaries_.begin(), wall_boundaries_.end(), [&](const MeshCell& mc){return mc.tag() == cells[i].tag();}) == 0);
            assert(wall_tag_index_map_.left.count(cells[i].tag()()) == 0);
        }
        //int initial = wall_boundaries_.size();
        for (int i=0; i<cells.size(); ++i)
        {
            if (wall_tag_index_map_.left.count(cells[i].tag()()) != 0) {
                continue;
            }
            assert(wall_tag_index_map_.right.count(wall_boundaries_.size() + i) == 0);
            wall_tag_index_map_.insert(boost::bimap<int, int>::value_type(cells[i].tag()(), wall_boundaries_.size()));
            wall_boundaries_.push_back(cells[i]);
        }
        assert(wall_boundaries_.size() == wall_tag_index_map_.left.size());
    }*/

    /*void Mesh::add_dirichlets(std::vector<MeshCell>& cells)
    {
        assert(dirichlet_boundaries_.size() == dirichlet_tag_index_map_.left.size());
        for (const MeshCell& c: cells)
        {
            assert(std::count_if(dirichlet_boundaries_.begin(), dirichlet_boundaries_.end(), [&](const MeshCell& mc){return mc.tag() == c.tag();}) == 0);
            assert(dirichlet_tag_index_map_.left.count(c.tag()()) == 0);
        }
        //int initial = dirichlet_boundaries_.size();
        for (int i=0; i<cells.size(); ++i)
        {
            if (dirichlet_tag_index_map_.left.count(cells[i].tag()()) != 0) {
                continue;
            }
            assert(dirichlet_tag_index_map_.right.count(dirichlet_boundaries_.size() + i) == 0);
            dirichlet_tag_index_map_.insert(boost::bimap<int, int>::value_type(cells[i].tag()(), dirichlet_boundaries_.size()));
            dirichlet_boundaries_.push_back(cells[i]);
        }
        assert(dirichlet_boundaries_.size() == dirichlet_tag_index_map_.left.size());
    }*/

    void Mesh::add_wall_boundary(MeshCell&& mc)
    {
        if (std::find_if(wall_boundaries_.begin(), wall_boundaries_.end(), [&](const MeshCell& c){return c.tag() == mc.tag();}) != wall_boundaries_.end()) {
            return;
        }

        wall_boundaries_.push_back(mc);
        //wall_tag_index_map_.insert(boost::bimap<int, int>::value_type(mc.tag()(), wall_boundaries_.size() - 1));
    }

    void Mesh::add_wall_boundary(const MeshCell& mc)
    {
        // no points are added because interiors, therefore, points are added already.
        if (std::find_if(wall_boundaries_.begin(), wall_boundaries_.end(), [&](const MeshCell& c){return c.tag() == mc.tag();}) != wall_boundaries_.end()) {
            return;
        }

        wall_boundaries_.push_back(mc);
        //wall_tag_index_map_.insert(boost::bimap<int, int>::value_type(mc.tag()(), wall_boundaries_.size() - 1));
    }

    void Mesh::add_symmetry_boundary(MeshCell&& mc)
    {
        if (std::find_if(symmetry_boundaries_.begin(), symmetry_boundaries_.end(), [&](const MeshCell& c){return c.tag() == mc.tag();}) != symmetry_boundaries_.end()) {
            return;
        }

        symmetry_boundaries_.push_back(mc);
        //symmetry_tag_index_map_.insert(boost::bimap<int, int>::value_type(mc.tag()(), symmetry_boundaries_.size() - 1));
    }

    void Mesh::add_symmetry_boundary(const MeshCell& mc)
    {
        // no points are added because interiors, therefore, points are added already.
        if (std::find_if(symmetry_boundaries_.begin(), symmetry_boundaries_.end(), [&](const MeshCell& c){return c.tag() == mc.tag();}) != symmetry_boundaries_.end()) {
            return;
        }

        symmetry_boundaries_.push_back(mc);
        //symmetry_tag_index_map_.insert(boost::bimap<int, int>::value_type(mc.tag()(), symmetry_boundaries_.size() - 1));
    }

    void Mesh::add_dirichlet_boundary(MeshCell&& mc)
    {
        // no points are added because interiors, therefore, points are added already.
        if (std::find_if(dirichlet_boundaries_.begin(), dirichlet_boundaries_.end(), [&](const MeshCell& c){return c.tag() == mc.tag();}) != dirichlet_boundaries_.end()) {
            return;
        }
        dirichlet_boundaries_.push_back(mc);
        //dirichlet_tag_index_map_.insert(boost::bimap<int, int>::value_type(mc.tag()(), dirichlet_boundaries_.size() - 1));
    }

    void Mesh::add_dirichlet_boundary(const MeshCell& mc)
    {
        // no points are added because interiors, therefore, points are added already.
        if (std::find_if(dirichlet_boundaries_.begin(), dirichlet_boundaries_.end(), [&](const MeshCell& c){return c.tag() == mc.tag();}) != dirichlet_boundaries_.end()) {
            return;
        }

        dirichlet_boundaries_.push_back(mc);
        //dirichlet_tag_index_map_.insert(boost::bimap<int, int>::value_type(mc.tag()(), dirichlet_boundaries_.size() - 1));
    }

    void Mesh::add_farfield_boundary(const MeshCell& mc)
    {
        if (std::find_if(farfield_boundaries_.begin(), farfield_boundaries_.end(), [&](const MeshCell& c){return c.tag() == mc.tag();}) != farfield_boundaries_.end()) {
            return;
        }

        farfield_boundaries_.push_back(mc);
        //farfield_tag_index_map_.insert(boost::bimap<int, int>::value_type(mc.tag()(), farfield_boundaries_.size() - 1));
    }

    void Mesh::add_farfield_boundary(MeshCell&& mc)
    {
        if (std::find_if(farfield_boundaries_.begin(), farfield_boundaries_.end(), [&](const MeshCell& c){return c.tag() == mc.tag();}) != farfield_boundaries_.end()) {
            return;
        }

        farfield_boundaries_.push_back(mc);
        //farfield_tag_index_map_.insert(boost::bimap<int, int>::value_type(mc.tag()(), farfield_boundaries_.size() - 1));
    }

    void Mesh::add_empty_boundary(const MeshCell& mc)
    {
        // no points are added because interiors, therefore, points are added already.
        if (std::find_if(empty_boundaries_.begin(), empty_boundaries_.end(), [&](const MeshCell& c){return c.tag() == mc.tag();}) != empty_boundaries_.end()) {
            return;
        }

        empty_boundaries_.push_back(mc);
        //empty_tag_index_map_.insert(boost::bimap<int, int>::value_type(mc.tag()(), empty_boundaries_.size() - 1));
    }

    void Mesh::add_interog_boundary(const MeshCell& mc)
    {
        // no points are added because interiors, therefore, points are added already.
        if (std::find_if(interog_boundaries_.begin(), interog_boundaries_.end(), [&](const MeshCell& c){return c.tag() == mc.tag();}) != interog_boundaries_.end()) {
            return;
        }

        interog_boundaries_.push_back(mc);
        //interog_tag_index_map_.insert(boost::bimap<int, int>::value_type(mc.tag()(), interog_boundaries_.size() - 1));
    }

    void Mesh::add_interog_boundary(MeshCell&& mc)
    {
        if (std::find_if(interog_boundaries_.begin(), interog_boundaries_.end(), [&](const MeshCell& c){return c.tag() == mc.tag();}) != interog_boundaries_.end()) {
            return;
        }

        interog_boundaries_.push_back(mc);
        //interog_tag_index_map_.insert(boost::bimap<int, int>::value_type(mc.tag()(), interog_boundaries_.size() - 1));
    }

    void Mesh::add_empty_boundary(MeshCell&& mc)
    {
        // no points are added because interiors, therefore, points are added already.
        if (std::find_if(empty_boundaries_.begin(), empty_boundaries_.end(), [&](const MeshCell& c){return c.tag() == mc.tag();}) != empty_boundaries_.end()) {
            return;
        }

        empty_boundaries_.push_back(mc);
        //empty_tag_index_map_.insert(boost::bimap<int, int>::value_type(mc.tag()(), empty_boundaries_.size() - 1));
    }

    //void Mesh::add_wall_boundaries(const Mesh& wall_mesh)
    //{
    //for (const MeshCell mc: wall_mesh.cell())
    //{
    //add_wall_boundary(mc);
    //}
    //}

    //void Mesh::add_dirichlet_boundaries(const Mesh& dirichlet_mesh)
    //{
    //for (const MeshCell mc: dirichlet_mesh.cell())
    //{
    //add_dirichlet_boundary(mc);
    //}
    //}

    void Mesh::merge_dirichlet_to_interior(const Mesh& dirichlet_mesh)
    {
        for (const MeshCell mc: dirichlet_mesh.cell_)
        {
            auto pit = std::find_if(point_.begin(), point_.end(), [&](const MeshPoint& mp){return mp.tag() == mc.point(0).tag();});
            if (pit == point_.end()) continue;

            for (auto& pc: pit->parent_cell())
            {
                auto& cellpc = cell_p(pc);
                //auto& cellpc = *pc.addr();
                int count = std::count_if(cellpc.point().begin(), cellpc.point().end(), [&](const MeshPoint& mp){return mp.tag() == mc.point(1).tag();});
                if (count != 0)
                {
                    cellpc.add_dirichlet_boundary(mc.tag());
                    dirichlet_boundaries_.push_back(mc);
                    dirichlet_boundaries_.back().set_interior_boundary(pc);
                    //dirichlet_tag_index_map_.insert(boost::bimap<int, int>::value_type(mc.tag()(), dirichlet_boundaries_.size() - 1));
                    break;
                }
            }
        }
    }

    //void Mesh::set_cell_oga_status_nonresident(const Tag& ic)
    //{
    //assert(ic.isvalid());
    //assert(i >= 0);
    ////assert(i < cell_.size());
    ////cell_[i].set_oga_cell_type(OGA_cell_type_t::non_resident);
    //cell_p(ic).set_oga_cell_type(OGA_cell_type_t::non_resident);
    //}

    /*bool do_overlap_both_boundaries()
      {
      bool refine = true;

      while(refine)
      {
      refine = false;
      for (const HoleMapBin& b: hole_map.bin())
      {
      for (const MeshCell& mc: cell_)
      {
      if (mc.boundary_type() == boundary_t::dirichlet)
      {
      if (b.aabb().do_intersect(mc.polygon()))
      {
      overlap_dirichlet = true;
      break;
      }
      }
      }
      for (const MeshCell& mc: cell_)
      {
      if (mc.boundary_type() == boundary_t::wall)
      {
      if (b.aabb().do_intersect(mc.polygon()))
      {
      overlap_wall = true;
      break;
      }
      }
      }

      if (overlap_dirichlet && overlap_wall)
      {
      refine = true;
      break;
      }
      }
      }
      }*/

    //void Mesh::set_hole_aabb(double minx, double miny, double minz, double maxx, double maxy, double maxz)
    //{
        //hole_aabb_.set_bbox(minx, miny, minz, maxx, maxy, maxz);
    //}

    //void Mesh::set_hole_aabb(const AABB& aabb)
    //{
        //hole_aabb_.set_bbox(aabb.min(0), aabb.min(1), aabb.min(2), aabb.max(0), aabb.max(1), aabb.max(2));
    //}

    void Mesh::convert_candholes_to_holes()
    {
        for (MeshCell& mc: cell_)
        {
            if (mc.oga_cell_type() == OGA_cell_type_t::hole_candidate)
                mc.set_oga_cell_type(OGA_cell_type_t::hole);
        }
    }

    void Mesh::flood_fill(const Tag& start)
    {
        auto& mc = cell_p(start);
        mc.set_oga_cell_type(OGA_cell_type_t::hole);

        for (const auto& t: mc.pnei())
        {
            const MeshCell& nei = cell(t);
            //const MeshCell& nei = *t.const_addr();
            if (nei.oga_cell_type() == OGA_cell_type_t::hole) continue;

            flood_fill(t);
        }
    }

    bool Mesh::flood_fill(const Tag& start, const AABB& aabb, int dummyrank, int sptag, std::ofstream& out)
    {
        MeshCell& mc = cell_p(start);
        mc.set_oga_cell_type(OGA_cell_type_t::hole_candidate);

        {
            if (out.is_open()) {
                out << "mc.tag() " << mc.tag()() << " type: " << static_cast<int>(mc.oga_cell_type()) << std::endl;
            }
        }

        for (const auto& t: mc.pnei())
        {
            const MeshCell& nei = cell(t);
            //const MeshCell& nei = *t.const_addr();
            OGA_cell_type_t neitype = nei.oga_cell_type();

            if (neitype == OGA_cell_type_t::hole || neitype == OGA_cell_type_t::hole_candidate)
            {
                {
                    if (out.is_open()) {
                        out << "nei is hole or holecand. continuning." << std::endl;
                    }
                }
                continue;
            }

            /*if (neitype == OGA_cell_type_t::non_resident)
              {
              mc.set_oga_cell_type(OGA_cell_type_t::undefined);
              if (mc.tag()() == 31600 && mc.parent_mesh()() == 4)
              {
              std::cout << "22222222222222222222" << std::endl;
              }
              return false;
              }*/

            if (!aabb.do_intersect(AABB(nei.poly())))
            {
                mc.set_oga_cell_type(OGA_cell_type_t::undefined);
                //if (mc.parent_mesh()() == 3 && dummyrank == 2 && sptag == 5)
                {
                    if (out.is_open()) {
                        out << "hole aabb does not intersect aabb of nei. returning." << std::endl;
                    }
                }
                return false;
            }

            //if (mc.parent_mesh()() == 3 && dummyrank == 2 && sptag == 5)
            {
                if (out.is_open()) {
                    out << "hole aabb intersects aabb of nei. flooding for the neighbor." << std::endl;
                }
            }

            bool success = flood_fill(t, aabb, dummyrank, sptag, out);
            //flood_fill(t, aabb, dummyrank, verbose);
            if (!success)
            {
                mc.set_oga_cell_type(OGA_cell_type_t::undefined);
                if (out.is_open()) {
                    out << "no success with flooding. returning false." << std::endl;
                }
                return false;
            }

            if (out.is_open()) {
                out << "success with flooding." << std::endl;
            }
        }

        //if (mc.parent_mesh()() == 3 && dummyrank == 2 && sptag == 5)
        {
            if (out.is_open()) {
                out << "getting out of flood function. returning true." << std::endl;
            }
        }

        return true;
    }

    void Mesh::set_as_field(const Tag& t)
    {
        cell_p(t).set_oga_cell_type(OGA_cell_type_t::field);
    }

    void Mesh::set_as_hole(const Tag& t)
    {
        cell_p(t).set_oga_cell_type(OGA_cell_type_t::hole);
    }

    //void Mesh::set_as_mreceptor(const Tag& t)
    //{
        //cell_p(t).set_oga_cell_type(OGA_cell_type_t::mandat_receptor);
    //}

    void Mesh::set_as_hole(const std::set<Tag>& overlaps)
    {
        for (const Tag& t: overlaps)
            cell_p(t).set_oga_cell_type(OGA_cell_type_t::hole);
    }

    void Mesh::check_nei(const MeshCell& mc, const Polyhedron& pol, std::set<Tag>& overlaps, int cuttercelltag, bool verbose) const
    {
        for (const auto& t: mc.pnei())
        {
            const MeshCell nei = cell(t);
            //const MeshCell nei = *t.const_addr();
            //if (mc.tag()() == 30002 && mc.parent_mesh()() == 3 && cuttercelltag == 418)
            //{
                //verbose = true;
            //}
            //if (mc.tag()() == 30002 && mc.parent_mesh()() == 3 && cuttercelltag == 418)
                //std::cout << "kpkpkpkpkpkpk: " << mc.tag()() << " " << t() << " " << AABB(nei.poly()).do_intersect(AABB(pol)) << " " << cuttercelltag << " " << overlaps.count(t) << std::endl;
            //}
            if (overlaps.count(t) != 0) continue;
            //if (nei.is_ghost()) continue;

            //Segment s(pol.edge(0));
            //if (nei.polygon().do_intersect(pol))
            //if (nei.poly().do_intersect(pol))
            if (AABB(nei.poly()).do_intersect(AABB(pol))) // not accurate but should work for now.
            {
                overlaps.insert(t);
                check_nei(nei, pol, overlaps, cuttercelltag, verbose);
            }
        }
    }

    //size_t Mesh::n_ghost() const
    //{
    //return n_ghost_cell_;
    //}

    //std::vector<MeshCell>& Mesh::cell_p()
    //{
    //return cell_;
    //}

    //void Mesh::handle_donor_conflict(const std::deque<Mesh>& mesh)
    //{
    //// handle orphans.
    //// remove invalid candidate donors. this is why this is called before determine_best_donor() to avoid comparing invalid candidate donors.

    //for (MeshCell& c: cell_)
    //{
    //if (c.oga_cell_type() == OGA_cell_type_t::receptor || c.oga_cell_type() == OGA_cell_type_t::mandat_receptor)
    //{
    //if (c.oga_cell_type() == OGA_cell_type_t::receptor)
    //assert(!c.cand_donor().empty());
    //for (auto d = c.cand_donor().begin(); d != c.cand_donor().end();)
    //{
    //const Tag& cdm = d->first;
    //const Tag& cdc = d->second;

    //const Mesh* donor_mesh = NULL;
    //for (auto it=mesh.begin(); it!= mesh.end(); ++it)
    //{
    //if (it->tag() == cdm)
    //{
    //donor_mesh = &(*(it));
    //break;
    //}
    //}
    //assert(donor_mesh != NULL);

    //if (!donor_mesh->query(cdc))
    //{
    //std::cout << "c tag = " << c.tag()() << std::endl;
    //std::cout << "pm tag = " << tag_() << std::endl;
    //std::cout << "cdc = " << cdc() << std::endl;
    //std::cout << "donor mesh = " << donor_mesh->tag()() << std::endl;
    //}
    //assert(donor_mesh->query(cdc));


    //if (donor_mesh->cell(cdc).oga_cell_type() == OGA_cell_type_t::receptor || donor_mesh->cell(cdc).oga_cell_type() == OGA_cell_type_t::mandat_receptor || donor_mesh->cell(cdc).oga_cell_type() == OGA_cell_type_t::orphan)
    //{
    //c.remove_cand_donor(cdc, cdm);
    //}
    //else
    //{
    //++d;
    //}
    //}
    //if (c.cand_donor().empty())
    //{
    //c.set_oga_cell_type(OGA_cell_type_t::orphan);
    //}
    //}
    //}
    //}

    //void Mesh::insert_cells(const std::vector<MeshCell>& cells)
    //{
    //    cell_.reserve(cell_.size() + cells.size());
    //    cell_.insert(cell_.end(), cells.begin(), cells.end());
    //}

    /*void Mesh::set_status_of_cells(const SPCDonorInfo& spcdi)
    {
        for (MeshCell& mc: cell_)
        {
            auto cdi = spcdi.query(tag(), mc.tag());
            //if (cdi == nullptr) continue;
            if (cdi == nullptr)
            {
                //mc.set_oga_cell_type(OGA_cell_type_t::non_resident); 
                continue;
            }
            //assert(cdi != nullptr);

            if (cdi->final_oga_cell_type() == OGA_cell_type_t::receptor || cdi->final_oga_cell_type() == OGA_cell_type_t::mandat_receptor)
            {
                assert(cdi->donor_mesh().isvalid());
                assert(cdi->donor_cell().isvalid());
                mc.set_oga_cell_type(cdi->final_oga_cell_type(), cdi->donor_mesh(), cdi->donor_cell()); 
            }
            else
            {
                assert(cdi->final_oga_cell_type() != OGA_cell_type_t::undefined);
                mc.set_oga_cell_type(cdi->final_oga_cell_type()); 
            }
        }
    }*/


    void Mesh::sort_cells()
    {
        std::sort(cell_.begin(), cell_.end(), [](const MeshCell& left, const MeshCell& right){return left.tag() < right.tag();});
        //cell_.sort([](const MeshCell& left, const MeshCell& right){return left.tag() < right.tag();});
    }

    void Mesh::sort_points()
    {
        std::sort(point_.begin(), point_.end(), [](const MeshPoint& left, const MeshPoint& right){return left.tag() < right.tag();});
    }

    void Mesh::simple_merge(const Mesh& other_mesh, int rank)
    {
        assert(tag_ == other_mesh.tag());
        reserve_interior(cell_.size() + other_mesh.cell().size());
        cell_.insert(cell_.end(), other_mesh.cell().begin(), other_mesh.cell().end());
        sort_cells();

        for (const MeshCell& mc: other_mesh.wall_boundaries_)
        {
            int count = std::count_if(wall_boundaries_.begin(), wall_boundaries_.end(), [&](const MeshCell& wb){return wb.tag() == mc.tag();});
            //assert(count == 0);
            if (count != 0) continue;
            add_wall_boundary(mc);
        }
        for (const MeshCell& mc: other_mesh.symmetry_boundaries_)
        {
            int count = std::count_if(symmetry_boundaries_.begin(), symmetry_boundaries_.end(), [&](const MeshCell& wb){return wb.tag() == mc.tag();});
            //assert(count == 0);
            if (count != 0) continue;
            add_symmetry_boundary(mc);
        }
        for (const MeshCell& mc: other_mesh.dirichlet_boundaries_)
        {
            int count = std::count_if(dirichlet_boundaries_.begin(), dirichlet_boundaries_.end(), [&](const MeshCell& wb){return wb.tag() == mc.tag();});
            //assert(count == 0);
            if (count != 0) continue;
            add_dirichlet_boundary(mc);
        }
        for (const MeshCell& mc: other_mesh.farfield_boundaries_)
        {
            int count = std::count_if(farfield_boundaries_.begin(), farfield_boundaries_.end(), [&](const MeshCell& wb){return wb.tag() == mc.tag();});
            //assert(count == 0);
            if (count != 0) continue;
            add_farfield_boundary(mc);
        }
        for (const MeshCell& mc: other_mesh.empty_boundaries_)
        {
            int count = std::count_if(empty_boundaries_.begin(), empty_boundaries_.end(), [&](const MeshCell& wb){return wb.tag() == mc.tag();});
            //assert(count == 0);
            if (count != 0) continue;
            add_empty_boundary(mc);
        }
        for (const MeshCell& mc: other_mesh.interog_boundaries_)
        {
            int count = std::count_if(interog_boundaries_.begin(), interog_boundaries_.end(), [&](const MeshCell& wb){return wb.tag() == mc.tag();});
            //assert(count == 0);
            if (count != 0) continue;
            add_interog_boundary(mc);
        }
    }

    void Mesh::merge_batch(const Mesh& other_mesh, int rank)
    {
        assert(tag_ == other_mesh.tag());
        //reserve_interior(cell_.size() + other_mesh.cell().size());
        //cell_.reserve(cell_.size() + other_mesh.cell().size());
        cell_.insert(cell_.end(), other_mesh.cell_.begin(), other_mesh.cell_.end());
        //point_.reserve(point_.size() + other_mesh.point().size());
        point_.insert(point_.end(), other_mesh.point_.begin(), other_mesh.point_.end());
        //sort_cells();
        //shrink_cells();
        //shrink_points();
        for (const MeshCell& mc: other_mesh.wall_boundaries_)
        {
            int count = std::count_if(wall_boundaries_.begin(), wall_boundaries_.end(), [&](const MeshCell& wb){return wb.tag() == mc.tag();});
            //assert(count == 0);
            if (count != 0) continue;
            add_wall_boundary(mc);
        }
        for (const MeshCell& mc: other_mesh.symmetry_boundaries_)
        {
            int count = std::count_if(symmetry_boundaries_.begin(), symmetry_boundaries_.end(), [&](const MeshCell& wb){return wb.tag() == mc.tag();});
            //assert(count == 0);
            if (count != 0) continue;
            add_symmetry_boundary(mc);
        }
        for (const MeshCell& mc: other_mesh.dirichlet_boundaries_)
        {
            int count = std::count_if(dirichlet_boundaries_.begin(), dirichlet_boundaries_.end(), [&](const MeshCell& wb){return wb.tag() == mc.tag();});
            //assert(count == 0);
            if (count != 0) continue;
            add_dirichlet_boundary(mc);
        }
        for (const MeshCell& mc: other_mesh.farfield_boundaries_)
        {
            int count = std::count_if(farfield_boundaries_.begin(), farfield_boundaries_.end(), [&](const MeshCell& wb){return wb.tag() == mc.tag();});
            //assert(count == 0);
            if (count != 0) continue;
            add_farfield_boundary(mc);
        }
        for (const MeshCell& mc: other_mesh.empty_boundaries_)
        {
            int count = std::count_if(empty_boundaries_.begin(), empty_boundaries_.end(), [&](const MeshCell& wb){return wb.tag() == mc.tag();});
            //assert(count == 0);
            if (count != 0) continue;
            add_empty_boundary(mc);
        }
        for (const MeshCell& mc: other_mesh.interog_boundaries_)
        {
            int count = std::count_if(interog_boundaries_.begin(), interog_boundaries_.end(), [&](const MeshCell& wb){return wb.tag() == mc.tag();});
            //assert(count == 0);
            if (count != 0) continue;
            add_interog_boundary(mc);
        }
    }

    void Mesh::merge(const Mesh& other_mesh)
    {
        assert(tag_ == other_mesh.tag());
        //reserve_interior(cell_.size() + other_mesh.cell().size());
        for (const MeshCell& mc: other_mesh.cell())
        {
            //assert(!query(mc.tag()));
            //if (query(mc.tag())) continue;
            //auto it = std::find_if(cell_.begin(), cell_.end(), [&](const MeshCell& tmc){return tmc.tag() == mc.tag();});
            //if (it == cell_.end()) {
            bool exist;
            add_interior_cell_sorted(mc, exist);
            //add_interior_cell_nonsorted(mc);
            //}
        }
        //sort_cells();
        //shrink_cells();
        //shrink_points();
        for (const MeshCell& mc: other_mesh.wall_boundaries_)
        {
            int count = std::count_if(wall_boundaries_.begin(), wall_boundaries_.end(), [&](const MeshCell& wb){return wb.tag() == mc.tag();});
            //assert(count == 0);
            if (count != 0) continue;
            add_wall_boundary(mc);
        }
        for (const MeshCell& mc: other_mesh.symmetry_boundaries_)
        {
            int count = std::count_if(symmetry_boundaries_.begin(), symmetry_boundaries_.end(), [&](const MeshCell& wb){return wb.tag() == mc.tag();});
            //assert(count == 0);
            if (count != 0) continue;
            add_symmetry_boundary(mc);
        }
        for (const MeshCell& mc: other_mesh.dirichlet_boundaries_)
        {
            int count = std::count_if(dirichlet_boundaries_.begin(), dirichlet_boundaries_.end(), [&](const MeshCell& wb){return wb.tag() == mc.tag();});
            //assert(count == 0);
            if (count != 0) continue;
            add_dirichlet_boundary(mc);
        }
        for (const MeshCell& mc: other_mesh.farfield_boundaries_)
        {
            int count = std::count_if(farfield_boundaries_.begin(), farfield_boundaries_.end(), [&](const MeshCell& wb){return wb.tag() == mc.tag();});
            //assert(count == 0);
            if (count != 0) continue;
            add_farfield_boundary(mc);
        }
        for (const MeshCell& mc: other_mesh.empty_boundaries_)
        {
            int count = std::count_if(empty_boundaries_.begin(), empty_boundaries_.end(), [&](const MeshCell& wb){return wb.tag() == mc.tag();});
            //assert(count == 0);
            if (count != 0) continue;
            add_empty_boundary(mc);
        }
        for (const MeshCell& mc: other_mesh.interog_boundaries_)
        {
            int count = std::count_if(interog_boundaries_.begin(), interog_boundaries_.end(), [&](const MeshCell& wb){return wb.tag() == mc.tag();});
            //assert(count == 0);
            if (count != 0) continue;
            add_interog_boundary(mc);
        }
    }

    /*void Mesh::make_cell_tags()
      {
      cell_tag_.clear();
      cell_tag_.reserve(cell_.size());
      for (const MeshCell& mc: cell_)
      {
      cell_tag_.push_back(mc.tag()());
      }

      std::sort(cell_tag_.begin(), cell_tag_.end());
      }*/

    /*void Mesh::merge_no_check(const Mesh& other_mesh)
      {
      assert(tag_ == other_mesh.tag());
      reserve_interior(cell_.size() + other_mesh.cell().size());
      for (const MeshCell& mc: other_mesh.cell())
      {
      add_interior_cell_no_check(mc);
      }
      for (const MeshCell& mc: other_mesh.wall_boundaries_)
      {
      auto it = std::find_if(wall_boundaries_.begin(), wall_boundaries_.end(), [&](const MeshCell& wb){return wb.tag() == mc.tag();});
      if (it != wall_boundaries_.end()) continue;
      add_wall_boundary(mc);
      }
      for (const MeshCell& mc: other_mesh.dirichlet_boundaries_)
      {
      auto it = std::find_if(dirichlet_boundaries_.begin(), dirichlet_boundaries_.end(), [&](const MeshCell& wb){return wb.tag() == mc.tag();});
      if (it != dirichlet_boundaries_.end()) continue;
      add_dirichlet_boundary(mc);
      }
      remove_merge_dup_pts();
      remove_dup_cells();
      set_point_tag_index_map();
      set_cell_tag_index_map();
      shrink_points();
      shrink_cells();
      }*/

    void Mesh::remove_merge_duplicate_cells()
    {
        remove_merge_dup(cell_);
    }

    void Mesh::remove_duplicate_cells(Profiler* profiler, std::string s)
    {
        remove_dup(cell_);
        //if (profiler != nullptr) {profiler->start(s + "-rem-dup-sort");}
        //std::sort(cell_.begin(), cell_.end());
        //if (profiler != nullptr) {profiler->stop(s + "-rem-dup-sort");}

        //if (profiler != nullptr) {profiler->start(s + "-rem-dup-uniq");}
        //cell_.erase(std::unique(cell_.begin(), cell_.end()), cell_.end());
        //if (profiler != nullptr) {profiler->stop(s + "-rem-dup-uniq");}
    }

    void Mesh::remove_duplicate_points()
    {
        remove_dup(point_);
    }

    //void Mesh::receptor_to_field(std::deque<Mesh>& mesh)
    //{
    //for (MeshCell& mc: cell_)
    //{
    //if (mc.oga_cell_type() != OGA_cell_type_t::mandat_receptor) continue;
    ////assert(mc.donor_mesh().isvalid());
    ////assert(mc.donor_cell().isvalid());
    //assert(mc.donor().first.isvalid());
    //assert(mc.donor().second.isvalid());

    //for (Mesh& m: mesh)
    //{
    ////if (m.tag() != mc.donor_mesh()) continue;
    //if (m.tag() != mc.donor().first) continue;

    ////m.cell_p(mc.donor_cell()).set_oga_cell_type(OGA_cell_type_t::field);
    //m.cell_p(mc.donor().second).set_oga_cell_type(OGA_cell_type_t::field);
    //break;
    //}
    //}
    //}

    //void Mesh::fringe_to_field()
    //{
    //for (MeshCell& mc: cell_)
    //{
    ////if (mc.is_ghost()) continue;

    ////if (mc.polygon().edge().size() < mc.nei().size())
    ////{
    ////std::cout << mc.polygon().edge().size() << std::endl;
    ////std::cout << mc.nei().size() << std::endl;
    ////}
    ////assert(mc.polygon().edge().size() >= mc.nei().size());
    ////assert(!mc.pnei().empty());

    //if (mc.oga_cell_type() != OGA_cell_type_t::receptor) continue;
    //if (mc.oga_cell_type() == OGA_cell_type_t::mandat_receptor) continue;

    //int count = 0;
    //for (const Tag& n: mc.pnei())
    //{
    //auto it = cell_tag_index_map.left.find(n());
    //if (it == cell_tag_index_map.left.end())
    //{
    //std::cout << mc.pnei().size() << std::endl;
    //std::cout << n() << std::endl;
    //}
    //assert(it != cell_tag_index_map.left.end());

    ////if (cell(n).is_ghost()) continue;

    //if (cell(n).oga_cell_type() == OGA_cell_type_t::field) ++count;
    //}

    //if (count != 0 && count == mc.pnei().size())
    //mc.set_oga_cell_type(OGA_cell_type_t::field);
    //}
    //}

    const MeshPoint* Mesh::query_point(const Tag& ic) const
    {
        assert(ic.isvalid());

        auto it = std::lower_bound(point_.begin(), point_.end(), MeshPoint(ic, Tag()));
        if (it != point_.end() && it->tag() == ic)
        {
            return &*it;
        }

        return nullptr;
    }

    const MeshCell* Mesh::query_sorted(const Tag& ic) const
    {
        assert(ic.isvalid());

        MeshCell tmc;
        tmc.set_tag(ic);

        auto it = std::lower_bound(cell_.begin(), cell_.end(), tmc, [](const MeshCell& left, const MeshCell& right){return left.tag()() < right.tag()();});

        if (it != cell_.end())
        {
            if (it->tag() == ic) {
                return &(*it);
            }
        }

        return nullptr;
    }

    const MeshCell* Mesh::query(const Tag& ic) const
    {
        assert(ic.isvalid());

        auto it = std::lower_bound(cell_.begin(), cell_.end(), MeshCell(ic, Tag()));
        if (it != cell_.end() && it->tag() == ic)
        {
            return &*it;
        }

        return nullptr;
    }

    mcc::iterator Mesh::query_itp(const Tag& ic)
    {
        assert(ic.isvalid());

        auto it = std::lower_bound(cell_.begin(), cell_.end(), MeshCell(ic, Tag()));
        if (it != cell_.end() && it->tag() == ic)
        {
            return it;
        }

        return cell_.end();
    }

const MeshCell* Mesh::query_bou(const Tag& ic, BouType type) const
{
    if (type == BouType::wall)
    {
        //int count = wall_tag_index_map_.left.count(ic());
        //if (count != 0) {
            return wall_boundary(ic);
        //}
    }
    else if (type == BouType::symmetry)
    {
        //int count = wall_tag_index_map_.left.count(ic());
        //if (count != 0) {
            return symmetry_boundary(ic);
        //}
    }
    else if (type == BouType::dirichlet)
    {
        //int count = dirichlet_tag_index_map_.left.count(ic());
        //if (count != 0) {
            return dirichlet_boundary(ic);
        //}
    }
    else if (type == BouType::farfield)
    {
        //int count = farfield_tag_index_map_.left.count(ic());
        //if (count != 0) {
            return farfield_boundary(ic);
        //}
    }
    else if (type == BouType::empty)
    {
        //int count = empty_tag_index_map_.left.count(ic());
        //if (count != 0) {
            return empty_boundary(ic);
        //}
    }
    else if (type == BouType::interog)
    {
        //int count = interog_tag_index_map_.left.count(ic());
        //if (count != 0) {
            return interog_boundary(ic);
        //}
    }
    else
    {
        assert(false);
    }

    return nullptr;
}

            //std::vector<ADTPoint> Mesh::make_adtpoint() const
            //{
            //std::vector<ADTPoint> adtpts;
            //adtpts.reserve(cell_.size());

            //for (const MeshCell& mc: cell_)
            //{
            ////if (mc.is_ghost()) continue;
            //adtpts.push_back(ADTPoint(mc.geom_point(), mc.tag()()));
            //}

            //return adtpts;
            //}

            //size_t Mesh::n_wall_ghost() const
            //{
            //size_t  count = 0;
            //for (const MeshCell& mc: cell_)
            //{
            //if (mc.boundary_type() == boundary_t::wall)
            //++count;
            //}

            //return count;
            //}
            //size_t Mesh::n_interior_cell() const
            //{
            //return n_interior_cell_;
            //}

            /*void Mesh::determine_hole_aabb()
            {
                // currently only one hole could exist per mesh.

                double minx, miny, minz, maxx, maxy, maxz;
                double gminx = TAILOR_BIG_POS_NUM;
                double gminy = TAILOR_BIG_POS_NUM;
                double gminz = TAILOR_BIG_POS_NUM;
                double gmaxx = TAILOR_BIG_NEG_NUM;
                double gmaxy = TAILOR_BIG_NEG_NUM;
                double gmaxz = TAILOR_BIG_NEG_NUM;

                Vector3 _min(TAILOR_BIG_POS_NUM, TAILOR_BIG_POS_NUM, TAILOR_BIG_POS_NUM);
                Vector3 _max(TAILOR_BIG_NEG_NUM, TAILOR_BIG_NEG_NUM, TAILOR_BIG_NEG_NUM);

                bool exist = false;
                for (const MeshCell& mc: wall_boundaries_)
                    //for (const MeshCell& mc: cell_)
                {
                    //if (mc.boundary_type() == boundary_t::wall)
                    {
                        minx = std::min(mc.point(0).p().r(0), mc.point(1).p().r(0));
                        miny = std::min(mc.point(0).p().r(1), mc.point(1).p().r(1));
                        minz = std::min(mc.point(0).p().r(2), mc.point(1).p().r(2));
                        maxx = std::max(mc.point(0).p().r(0), mc.point(1).p().r(0));
                        maxy = std::max(mc.point(0).p().r(1), mc.point(1).p().r(1));
                        maxz = std::max(mc.point(0).p().r(2), mc.point(1).p().r(2));

                        gminx = std::min(gminx, minx);
                        gminy = std::min(gminy, miny);
                        gminz = std::min(gminz, minz);
                        gmaxx = std::max(gmaxx, maxx);
                        gmaxy = std::max(gmaxy, maxy);
                        gmaxz = std::max(gmaxz, maxz);

                        exist = true;
                    }
                }

                if (!exist)
                {
                    //hole_aabb_.set_min(_min);
                    //hole_aabb_.set_max(_max);
                    hole_aabb_.set_bbox(_min, _max);
                    return;
                }

                _min.set(gminx, gminy, gminz);
                _max.set(gmaxx, gmaxy, gmaxz);
                //hole_aabb_.set_min(_min);
                //hole_aabb_.set_max(_max);
                hole_aabb_.set_bbox(_min, _max);
                return;
            }*/

            //const AABB& Mesh::hole_aabb() const
            //{
                //return hole_aabb_;
            //}

            void Mesh::reset_all_cell_oga_status()
            {
                for (MeshCell& mc: cell_)
                {
                    mc.donor_ = Donor();
                    mc.receptor_ = Receptor();

                    if (mc.oga_cell_type() == OGA_cell_type_t::ghost) {
                        continue;
                    }
                    //if (mc.oga_cell_type() == OGA_cell_type_t::non_resident) {
                        //continue;
                    //}
                    mc.reset_oga_status();
                    //mc.reset_btype();
                }

                //for (MeshFace& mf: face_)
                //{
                    //mf.reset_btype(); // except boundaries and interiors.
                //}
            }

            void Mesh::remove_cell_from_cellhood(const Tag& ic)
            {
                assert(ic.isvalid());
                assert(query(ic));
                MeshCell& mc = cell_p(ic);

                for (auto& inei: mc.pnei())
                {
                    assert(inei.isvalid());
                    //if (inei.addr() == nullptr)
                    if (!query(inei))
                    {
                        std::cout << tag_() << std::endl;
                        std::cout << ic() << std::endl;
                        //std::cout << inei() << std::endl;
                        print_as_vtk("themachine.vtk");
                    }
                    assert(query(inei));
                    cell_p(inei).remove_neighbor(ic);
                    //inei.addr()->remove_neighbor(ic);
                    //std::find(cell(inei).pnei().begin(), cell(inei).pnei().end(), )cell(inei).pnei().find();
                }

                mc.remove_all_neighbors();
            }

            //void Mesh::deparent_cell_from_vertices(const Tag& ic)
            //{
            //    for (const MeshPoint& mp: cell_p(ic).point())
            //    {
            //        auto pp =  point_p(mp.tag());
            //        assert(pp != nullptr);
            //        //point_p(mp.tag()).remove_parent_cell(ic);
            //        pp->remove_parent_cell(ic);
            //    }
            //}

            unsigned int Mesh::point_index(const Tag& t) const
            {
                auto it = std::lower_bound(point_.begin(), point_.end(), MeshPoint(t, Tag()));
                assert(it != point_.end());
                assert(it->tag() == t);
                return std::distance(point_.begin(), it);
            }

            void Mesh::rotate(double angle, int axis, const Vector3& rot_point)
            {
                assert(!std::isnan(rot_point(0)));
                assert(!std::isnan(rot_point(1)));
                assert(!std::isnan(rot_point(2)));
                // rotate points.
                for (MeshPoint& point: point_)
                {
                    point.rotate_point(angle, axis, rot_point);
                }
                // rotate faces.
                //for (MeshFace& face: face_)
                //{
                    //face.rotate(angle, axis, rot_point);
                //}
                // rotate cell points.
                for (MeshCell& cell: cell_)
                {
                    cell.rotate_points(angle, axis, rot_point);
                }
                // rotate cell points.
                for (MeshCell& cell: wall_boundaries_)
                {
                    cell.rotate_points(angle, axis, rot_point);
                }
                for (MeshCell& cell: symmetry_boundaries_)
                {
                    cell.rotate_points(angle, axis, rot_point);
                }
                for (MeshCell& cell: dirichlet_boundaries_)
                {
                    cell.rotate_points(angle, axis, rot_point);
                }
                for (MeshCell& cell: farfield_boundaries_)
                {
                    cell.rotate_points(angle, axis, rot_point);
                }
                for (MeshCell& cell: empty_boundaries_)
                {
                    cell.rotate_points(angle, axis, rot_point);
                }
                for (MeshCell& cell: interog_boundaries_)
                {
                    cell.rotate_points(angle, axis, rot_point);
                }
            }

            void Mesh::move(const Vector3& v)
            {
                // this is a temporary function to mimic solver.
                // aim is to displace certain meshblocks.

                {
                    // move points.
                    for (MeshPoint& point: point_)
                    {
                        point.move_point(v);
                    }
                    // move faces.
                    //for (MeshFace& face: face_)
                    //{
                        //face.move(v);
                    //}
                    // move cell points.
                    for (MeshCell& cell: cell_)
                    {
                        cell.move_points(v);
                    }

                    // move wall points.
                    for (MeshCell& cell: wall_boundaries_)
                    {
                        cell.move_points(v);
                    }

                    for (MeshCell& cell: symmetry_boundaries_)
                    {
                        cell.move_points(v);
                    }

                    // move dirichlet points.
                    for (MeshCell& cell: dirichlet_boundaries_)
                    {
                        cell.move_points(v);
                    }
                    
                    for (MeshCell& cell: empty_boundaries_)
                    {
                        cell.move_points(v);
                    }

                    for (MeshCell& cell: farfield_boundaries_)
                    {
                        cell.move_points(v);
                    }

                    for (MeshCell& cell: interog_boundaries_)
                    {
                        cell.move_points(v);
                    }

                    //hole_aabb_.set_min(hole_aabb_.min()+v);
                    //hole_aabb_.set_max(hole_aabb_.max()+v);
                }
            }

            void Mesh::convert_undefined_to_field()
            {
                for (MeshCell& c: cell_)
                {
                    if (c.oga_cell_type() == OGA_cell_type_t::undefined)
                    {
                        c.set_oga_cell_type(OGA_cell_type_t::field);
                    }
                }
            }

            void Mesh::convert_receptor_to_hole()
            {
                for (MeshCell& c: cell_)
                {
                    for (const auto& f: c.pnei())
                    {
                        const auto& nei = cell(f);
                        if (nei.oga_cell_type() == OGA_cell_type_t::field)
                        {

                            assert(nei.oga_cell_type() != OGA_cell_type_t::hole);
                        }
                    }
                }

                for (MeshCell& c: cell_)
                {
                    if (c.oga_cell_type() == OGA_cell_type_t::receptor)
                    {
                        bool has_field_pnei = false;
                        for (const auto& t: c.pnei())
                        {
                            if (cell(t).oga_cell_type() == OGA_cell_type_t::field)
                            {
                                has_field_pnei = true;
                                break;
                            }
                        }
                        if (!has_field_pnei)
                        {
                            c.set_oga_cell_type(OGA_cell_type_t::hole);
                        }
                    }
                }
                
                for (MeshCell& c: cell_)
                {
                    for (const auto& f: c.pnei())
                    {
                            const auto& nei = cell(f);
                        if (nei.oga_cell_type() == OGA_cell_type_t::field)
                        {

                            assert(nei.oga_cell_type() != OGA_cell_type_t::hole);
                        }
                    }
                }
            }

            void Mesh::set_cell(const Tag& t, const MeshCell& c)
            {
                cell_p(t) = c;
            }

            const mpc& Mesh::point() const
            {
                return point_;
            }

            /*std::set<Point> Mesh::bou_raw_point(BouType btype) const
            {
                std::set<Point> s;
                if (btype == BouType::wall) {
                    s = bou_raw_point_(&wall_boundaries_);
                }
                else if (btype == BouType::dirichlet) {
                    s = bou_raw_point_(&dirichlet_boundaries_);
                }
                else if (btype == BouType::farfield) {
                    s = bou_raw_point_(&farfield_boundaries_);
                }
                else if (btype == BouType::empty) {
                    s = bou_raw_point_(&empty_boundaries_);
                }
                else if (btype == BouType::interog) {
                    s = bou_raw_point_(&interog_boundaries_);
                }
                else {
                    assert(false);
                }

                return s;
            }*/

            /*std::set<Point> Mesh::bou_raw_point_(const std::vector<MeshCell>* container) const
            {
                std::set<Point> rawpoints;
                if (container->empty()) {
                    return rawpoints;
                }

                for (const auto& mc: *container)
                {
                    for (const auto& p: mc.poly().vertices())
                    {
                        rawpoints.insert(p);
                    }
                }

                return rawpoints;
            }*/

            std::vector<Point> Mesh::rawpoint() const
            {
                std::vector<Point> rawpoints;
                if (point_.empty()) {
                    return rawpoints;
                }
                rawpoints.reserve(point_.size());
                for (const auto& p: point_)
                {
                    rawpoints.push_back(p.p());
                }

                return rawpoints;
            }

            const Tag& Mesh::tag() const
            {
                return tag_;
            }

            void Mesh::set_tag(const Tag& t)
            {
                tag_ = t;
            }

            const mcc& Mesh::cell() const
            {
                return cell_;
            }

            mcc& Mesh::cell_p()
            {
                return cell_;
            }

            /*std::vector<MeshCell>& Mesh::cell()
              {
              return cell_;
              }*/

            MeshPoint* Mesh::point_p(const Tag& ptag)
            {
                assert(ptag.isvalid());

                auto it = std::lower_bound(point_.begin(), point_.end(), MeshPoint(ptag, Tag()));
                if (it == point_.end())
                {
                    return nullptr;
                }
                //assert(it != point_.end());
                assert(it->tag() == ptag);

                return &(*it);
            }

            const MeshPoint* Mesh::point(const Tag& ptag) const
            {
                assert(ptag.isvalid());

                auto it = std::lower_bound(point_.begin(), point_.end(), MeshPoint(ptag, Tag()));
                if (it == point_.end())
                {
                    return nullptr;
                }
                //assert(it != point_.end());
                assert(it->tag() == ptag);

                return &(*it);
            }

            //MeshPoint& Mesh::point_p(const Tag& ptag)
            //{
            //    assert(ptag.isvalid());

            //    auto it = std::lower_bound(point_.begin(), point_.end(), MeshPoint(ptag, Tag()));
            //    assert(it != point_.end());
            //    assert(it->tag() == ptag);

            //    return *it;
            //}

            //const MeshPoint& Mesh::point(const Tag& ptag) const
            //{
            //    assert(ptag.isvalid());

            //    auto it = std::lower_bound(point_.begin(), point_.end(), MeshPoint(ptag, Tag()));
            //    assert(it != point_.end());
            //    assert(it->tag() == ptag);

            //    return *it;
            //}

            const MeshCell* Mesh::cell_ptr(const Tag& ctag) const
            {
                assert(ctag.isvalid());

                auto it = std::lower_bound(cell_.begin(), cell_.end(), MeshCell(ctag, Tag()));
                if (it == cell_.end())
                {
                    return nullptr;
                }
                else if (it->tag() != ctag)
                {
                    return nullptr;
                }

                return &(*it);
            }

            //MeshFace& Mesh::face_p(const FaceTag& ctag)
            //{
            //    assert(ctag.isvalid());

            //    auto it = std::find_if(face_.begin(), face_.end(), [&](const MeshFace& mf){return mf.tag() == ctag;});
            //    assert(it != face_.end());
            //    //if (it->tag() != ctag)
            //    //{
            //        //std::cout << "it->tag(): " << it->tag()() << std::endl;
            //        //std::cout << "ctag: " << ctag() << std::endl;
            //    //}
            //    assert(it->tag() == ctag);

            //    return *it;
            //}

            //const MeshFace& Mesh::face(const FaceTag& ctag) const
            //{
            //    assert(ctag.isvalid());

            //    auto it = std::find_if(face_.begin(), face_.end(), [&](const MeshFace& mf){return mf.tag() == ctag;});
            //    assert(it != face_.end());
            //    //if (it->tag() != ctag)
            //    //{
            //        //std::cout << "it->tag(): " << it->tag()() << std::endl;
            //        //std::cout << "ctag: " << ctag() << std::endl;
            //    //}
            //    assert(it->tag() == ctag);

            //    return *it;
            //}

            const MeshCell& Mesh::cell(const Tag& ctag) const
            {
                assert(ctag.isvalid());

                auto it = std::lower_bound(cell_.begin(), cell_.end(), MeshCell(ctag, Tag()));
                assert(it != cell_.end());
                if (it->tag() != ctag)
                {
                    std::cout << "it->tag(): " << it->tag()() << std::endl;
                    std::cout << "ctag: " << ctag() << std::endl;
                }
                assert(it->tag() == ctag);

                return *it;
            }

            MeshCell& Mesh::cell_p(const Tag& ctag)
            {
                assert(ctag.isvalid());

                auto it = std::lower_bound(cell_.begin(), cell_.end(), MeshCell(ctag, Tag()));
                assert(it != cell_.end());
                if (it->tag() != ctag)
                {
                    std::cout << "it->tag(): " << it->tag()() << std::endl;
                    std::cout << "ctag: " << ctag() << std::endl;
                }
                assert(it->tag() == ctag);

                return *it;
            }

            void Mesh::add_cell_only(MeshCell c, size_t size)
            {
                auto add = [&] ()
                {
                    /*c.set_parent_cell_of_vertices();
                      assert(tag().isvalid());
                      for (const MeshPoint& p: c.point())
                      {
                      assert(c.tag().isvalid());
                      point_p(p.tag()).add_parent_cell(c.tag());
                      }
                    // add to container.
                    cell_.push_back(std::move(c));
                    for (const MeshPoint& p: cell_.back().point())
                    {
                    assert(p.parent_cell().size() > 0);
                    }
                    cell_tag_index_map.insert(boost::bimap<int, int>::value_type(c.tag()(), cell_.size() - 1));*/

                    bool dummy;
                    //add_cell_sorted(c, dummy);
                    add_cell_nonsorted(c);
                };

                // reserve memory for cell_ if size is provided.
                if (size != 0)
                {
                    //if (cell_.capacity() != size)
                    //{
                        //cell_.reserve(size);
                    //}
                }

                assert(c.tag().isvalid());
                if (c.tag().isvalid())
                {
                    //auto it = std::find_if(cell_.begin(), cell_.end(), [&c](const MeshCell& _c){return _c.tag() == c.tag();});
                    //assert(it == cell_.end());
                    //if (it == cell_.end())
                    {
                        add();
                    }
                }
                else
                {
                    add();
                }
            }

            void Mesh::shrink_cells()
            {
                mcc(cell_).swap(cell_);
            }

            void Mesh::shrink_points()
            {
                mpc(point_).swap(point_);
            }

            void Mesh::reserve_interior_only(size_t size)
            {
                //cell_.reserve(size); // TODO
            }

            void Mesh::reserve_interior(size_t size)
            {
                //cell_.reserve(size); // TODO
                point_.reserve(size*4); // assuming quad.
            }
            void Mesh::reserve_wall(size_t size)
            {
                //wall_boundaries_.reserve(size); // TODO
            }
            void Mesh::reserve_dirichlet(size_t size)
            {
                //dirichlet_boundaries_.reserve(size); // TODO
            }
            void Mesh::reserve_farfield(size_t size)
            {
                //farfield_boundaries_.reserve(size); // TODO
            }
            void Mesh::reserve_empty(size_t size)
            {
                //empty_boundaries_.reserve(size); // TODO
            }

            /*vo]d Mesh::add_interior_cell(MeshCell&& c)
              {
              c.set_parent_mesh(tag());
              c.remove_parent_cells_of_vertices();
              cell_.push_back(c);
              cell_tag_index_map.insert(boost::bimap<int, int>::value_type(cell_.back().tag()(), cell_.size() - 1));
              cell_.back().set_parent_cell_of_vertices();

              for (const MeshPoint& p: cell_.back().point())
              {
              auto it = std::find_if(point_.begin(), point_.end(), [&p](const MeshPoint& mp){return mp.tag() == p.tag();});
              if (it == point_.end())
              { 
              point_.push_back(p);
            //auto itt = point_tag_index_map.left.find(p.tag()());
            //assert(itt == point_tag_index_map.left.end());
            point_tag_index_map.insert(boost::bimap<int, int>::value_type(p.tag()(), point_.size() - 1));
            }
            else
            {
            it->add_parent_cell(cell_.back().tag());
            }
            }
            }*/

            void Mesh::remove_dup_cells_and_points()
            {
                ////remove_merge_dup_pts();
                remove_merge_duplicate_cells(); // because there may be dup cells after add and push sp in remap. dups are due to being right on border of two or more bins.
                remove_merge_duplicate_points();

                ////set_point_tag_index_map();
                ////set_cell_tag_index_map();

                shrink_points();
                shrink_cells();
            }

            void Mesh::add_interior_nonsorted_addpoint(MeshCell mc)
            {
                mc.set_parent_mesh(tag());
                for (const MeshPoint& p: mc.point())
                {
                    point_.push_back(p);
                }
                cell_.push_back(std::move(mc));
            }

            void Mesh::add_interiors_nonsorted_nopoint(mcc& cells)
            {
                for (MeshCell& mc: cells)
                {
                    mc.set_parent_mesh(tag());
                    mc.remove_parent_cells_of_vertices();
                    mc.set_parent_cell_of_vertices();
                }

                cell_.insert(cell_.end(), cells.begin(), cells.end());
            }

            void Mesh::add_interior_cell_no_check(const MeshCell& c)
            {
                assert(tag().isvalid());
                assert(c.tag().isvalid());
                //assert(!query(c.tag()));

                cell_.push_back(c);
                cell_.back().set_parent_mesh(tag());
                cell_.back().remove_parent_cells_of_vertices();

                //assert(query(cell_.back().tag()));

                cell_.back().set_parent_cell_of_vertices();

                assert(!cell_.back().point().empty());

                point_.insert(point_.end(), cell_.back().point().begin(), cell_.back().point().end());

                for (const MeshPoint& p: cell_.back().point())
                {
                    assert(!p.parent_cell().empty());
                }

                for (const MeshPoint& _p: point_)
                {
                    for (const auto& _t: _p.parent_cell())
                    {
                        assert(_t.isvalid());
                    }
                }
            }

            /*void Mesh::set_cell_tag_index_map()
              {
              cell_tag_index_map.clear();
              for (const MeshCell& mc: cell_)
              {
              cell_tag_index_map.insert(boost::bimap<int, int>::value_type(mc.tag()(), cell_tag_index_map.size()));
              }
              }

              void Mesh::set_point_tag_index_map()
              {
              point_tag_index_map.clear();
              for (const MeshPoint& mp: point_)
              {
              point_tag_index_map.insert(boost::bimap<int, int>::value_type(mp.tag()(), point_tag_index_map.size()));
              }
              }*/

            /*void Mesh::add_interior_cell_find(const MeshCell& c)
              {
              assert(tag().isvalid());
              assert(c.tag().isvalid());
              assert(!query(c.tag()));

              cell_.push_back(c);
              cell_.back().set_parent_mesh(tag());
              cell_.back().remove_parent_cells_of_vertices();
              cell_tag_index_map.insert(boost::bimap<int, int>::value_type(c.tag()(), cell_.size() - 1));

              assert(query(cell_.back().tag()));

              cell_.back().set_parent_cell_of_vertices();

              assert(!cell_.back().point().empty());

              for (const MeshPoint& p: cell_.back().point())
              {
            //assert(p.parent_mesh() == c.parent_mesh());

            auto it = std::find_if(point_.begin(), point_.end(), [&p](const MeshPoint& mp){return mp.tag() == p.tag();});
            if (it == point_.end())
            { 
            assert(cell_.back().tag() == c.tag());
            assert(p.tag().isvalid());
            point_.push_back(p);
            //auto itt = point_tag_index_map.left.find(p.tag()()); 
            //assert(itt == point_tag_index_map.left.end());
            point_tag_index_map.insert(boost::bimap<int, int>::value_type(p.tag()(), point_.size() - 1));
            //for (const Tag& _t: point_.back().parent_cell())
            //{
            //assert(_t.isvalid());
            //}

            }
            else
            {
            it->add_parent_cell(cell_.back().tag());
            }
            }

            for (const MeshPoint& p: cell_.back().point())
            {
            assert(!p.parent_cell().empty());
            }

            for (const MeshPoint& _p: point_)
            {
            for (const Tag& _t: _p.parent_cell())
            {
            assert(_t.isvalid());
            }
            }
            }*/

            //std::vector<MeshCell>::iterator Mesh::add_cell_sorted(const MeshCell& c, bool& exist)
            mcc::iterator Mesh::add_cell_sorted(const MeshCell& c, bool& exist)
            {
                auto it = std::lower_bound(cell_.begin(), cell_.end(), c, [](const MeshCell& left, const MeshCell& right){return left.tag()() < right.tag()();});

                if (it == cell_.end())
                { 
                    cell_.push_back(c);
                    it = std::prev(cell_.end());
                }
                else if (it->tag() != c.tag())
                {
                    //assert(it->tag() != c.tag());
                    it = cell_.insert(it, c);
                }
                else
                {
                    exist = true;
                    return it;
                }

                it->set_parent_mesh(tag_);
                it->remove_parent_cells_of_vertices();
                it->set_parent_cell_of_vertices();

                assert(!it->point().empty());

                return it;
            }

            void Mesh::remove_parent_cells_of_all_points()
            {
                for (MeshPoint& mp: point_)
                {
                    mp.remove_parent_cells();
                }
            }

            void Mesh::remove_parent_cells_of_vertices_of_all_cells()
            {
                for (MeshCell& mc: cell_)
                {
                    mc.remove_parent_cells_of_vertices();
                }
            }

            void Mesh::set_parent_cell_of_vertices_of_all_cells()
            {
                for (MeshCell& mc: cell_)
                {
                    mc.set_parent_cell_of_vertices();
                }

                //for (const MeshPoint& mp: point())
                //{
                    //for (const auto& pc: mp.const_parent_cell())
                    //{
                        //if (mp.tag()() == 1569)
                        //{
                            //std::cout << "pc: " << pc.tag()() << std::endl;
                        //}

                        //assert(query(pc.tag()) != nullptr);
                        //if (pc.const_addr() == nullptr)
                        //{
                            //std::cout << "pc: " << pc.tag()() << std::endl;
                        //}
                        //assert(pc.const_addr() != nullptr);
                    //}
                //}
            }

            void Mesh::update_points_from_cell_vertices(int rank)
            {
                // do this after points are sorted and mesh cells vertices are ready.
                for (MeshCell& mc: cell_)
                {
                    for (const MeshPoint& mp: mc.point())
                    {
                        auto it = std::lower_bound(point_.begin(), point_.end(), mp, [](const MeshPoint& left, const MeshPoint& right){return left.tag()() < right.tag()();});
                        assert(it != point_.end());
                        if (it->tag() != mp.tag())
                        {
                            std::cout << "mesh: " << tag_() << std::endl;
                            std::cout << "mc: " << mc.tag()() << std::endl;
                            std::cout << "mp: " << mp.tag()() << std::endl;
                            std::cout << "rank: " << rank << std::endl;
                            std::cout << "mp parent size: " << mp.parent_cell().size() << std::endl;
                        }
                        assert(it->tag() == mp.tag());
                        it->add_parent_cell(mc.tag());
                        //std::cout << it->tag()() << std::endl;
                        //assert(mp.tag()() != 1569);
                        //assert(mp.tag()() != 1570);
                        //assert(tag()() != 6);
                        //assert(it->tag()() != 1570);

                        //for (const Tag& pc: it->parent_cell())
                        //{
                            //assert(query(pc) != nullptr);
                        //}
                    }
                }
            }

            mcc::iterator Mesh::add_cell_nonsorted(const MeshCell& c)
            {
                cell_.push_back(std::move(c));
                auto it = std::prev(cell_.end());

                assert(tag_.isvalid());
                it->set_parent_mesh(tag_);
                it->remove_parent_cells_of_vertices();
                it->set_parent_cell_of_vertices();

                assert(!it->point().empty());

                return it;
            }

            void Mesh::add_points(const mcc& cells, int rank)
            {
                for (const MeshCell& mc: cells)
                {
                    for (const MeshPoint& p: mc.point())
                    {
                        auto it = std::lower_bound(point_.begin(), point_.end(), p, [](const MeshPoint& left, const MeshPoint& right){return left.tag()() < right.tag()();});

                        if (it == point_.end())
                        { 
                            assert(p.tag().isvalid());
                            point_.push_back(p);
                        }
                        else if (it->tag() == p.tag())
                        {
                            it->add_parent_cell(mc.tag());
                        }
                        else
                        {
                            point_.insert(it, p);
                        }
                    }
                }

                std::sort(point_.begin(), point_.end());
            }

            mcc::iterator Mesh::add_interior_cell_sorted(const MeshCell& c, bool& exist)
            {
                assert(tag().isvalid());
                assert(c.tag().isvalid());
                //assert(!query(c.tag()));

                exist = false;
                auto mc = add_cell_sorted(c, exist);
                if (exist) {
                    return mc;
                }

                for (const MeshPoint& p: mc->point())
                {
                    //assert(p.parent_mesh() == c.parent_mesh());

                    //auto it = std::find_if(point_.begin(), point_.end(), [&p](const MeshPoint& mp){return mp.tag() == p.tag();});
                    // if point_ is not empty, must be already sorted.
                    auto it = std::lower_bound(point_.begin(), point_.end(), p, [](const MeshPoint& left, const MeshPoint& right){return left.tag()() < right.tag()();});
                    if (it == point_.end())
                    { 
                        assert(mc->tag() == c.tag());
                        assert(p.tag().isvalid());
                        point_.push_back(p);
                        //auto itt = point_tag_index_map.left.find(p.tag()()); 
                        //assert(itt == point_tag_index_map.left.end());
                        //point_tag_index_map.insert(boost::bimap<int, int>::value_type(p.tag()(), point_.size() - 1));
                        //for (const Tag& _t: point_.back().parent_cell())
                        //{
                        //assert(_t.isvalid());
                        //}

                    }
                    else if (it->tag() == p.tag())
                    {
                        it->add_parent_cell(mc->tag());
                    }
                    else
                    {
                        point_.insert(it, p);
                    }
                }

                for (const MeshPoint& p: cell_.back().point())
                {
                    assert(!p.parent_cell().empty());
                }

                //for (const MeshPoint& _p: point_)
                //{
                //    for (const auto& _t: _p.parent_cell())
                //    {
                //        assert(_t.isvalid());
                //    }
                //}

                return mc;
            }

            //void Mesh::add_meshpoint(const MeshCell& mc)
            void Mesh::add_meshpoint(MeshCell& mc, int rank)
            {
                for (const MeshPoint& p: mc.point())
                {
                    auto it = std::lower_bound(point_.begin(), point_.end(), p, [](const MeshPoint& left, const MeshPoint& right){return left.tag()() < right.tag()();});
                    if (it == point_.end())
                    { 
                        assert(p.tag().isvalid());
                        point_.push_back(p);

                    }
                    else if (it->tag() == p.tag())
                    {
                        it->add_parent_cell(mc.tag());
                    }
                    else
                    {
                        point_.insert(it, p);
                    }
                }
            }

            void Mesh::add_interior_cell_nonsorted(const MeshCell& c)
            {
                assert(tag().isvalid());
                assert(c.tag().isvalid());
                //assert(!query(c.tag()));

                auto mc = add_cell_nonsorted(c);

                for (const MeshPoint& p: mc->point())
                {
                    //assert(p.parent_mesh() == c.parent_mesh());

                    //auto it = std::find_if(point_.begin(), point_.end(), [&p](const MeshPoint& mp){return mp.tag() == p.tag();});
                    // if point_ is not empty, must be already sorted.
                    auto it = std::lower_bound(point_.begin(), point_.end(), p, [](const MeshPoint& left, const MeshPoint& right){return left.tag()() < right.tag()();});
                    if (it == point_.end())
                    { 
                        assert(mc->tag() == c.tag());
                        assert(p.tag().isvalid());
                        point_.push_back(p);
                        //auto itt = point_tag_index_map.left.find(p.tag()()); 
                        //assert(itt == point_tag_index_map.left.end());
                        //point_tag_index_map.insert(boost::bimap<int, int>::value_type(p.tag()(), point_.size() - 1));
                        //for (const Tag& _t: point_.back().parent_cell())
                        //{
                        //assert(_t.isvalid());
                        //}

                    }
                    else if (it->tag() == p.tag())
                    {
                        //if (it->oga_cell_type() == OGA_cell_type_t::non_resident && mc->oga_cell_type() != OGA_cell_type_t::non_resident) {
                            //it->set_oga_cell_type(OGA_cell_type_t::undefined);
                        //}
                        it->add_parent_cell(mc->tag());
                    }
                    else
                    {
                        point_.insert(it, p);
                    }
                }

                for (const MeshPoint& p: cell_.back().point())
                {
                    //assert(!p.parent_cell().empty());
                    assert(!p.parent_cell().empty());
                }

                for (const MeshPoint& _p: point_)
                {
                    for (const auto& _t: _p.parent_cell())
                    {
                        assert(_t.isvalid());
                    }
                }
            }

            void Mesh::add_interior_cell(const MeshCell& c)
            {
                assert(tag().isvalid());
                assert(c.tag().isvalid());
                assert(!query(c.tag()));

                //cell_.push_back(c);
                //cell_.back().set_parent_mesh(tag());
                //cell_.back().remove_parent_cells_of_vertices();
                //cell_tag_index_map.insert(boost::bimap<int, int>::value_type(c.tag()(), cell_.size() - 1));

                //assert(query(cell_.back().tag()));

                //cell_.back().set_parent_cell_of_vertices();

                //assert(!cell_.back().point().empty());

                bool dummy = false;
                auto mc = add_cell_sorted(c, dummy);

                for (const MeshPoint& p: mc->point())
                {
                    //assert(p.parent_mesh() == c.parent_mesh());

                    //auto it = std::find_if(point_.begin(), point_.end(), [&p](const MeshPoint& mp){return mp.tag() == p.tag();});
                    // if point_ is not empty, must be already sorted.
                    auto it = std::lower_bound(point_.begin(), point_.end(), p, [](const MeshPoint& left, const MeshPoint& right){return left.tag()() < right.tag()();});
                    if (it == point_.end())
                    { 
                        assert(mc->tag() == c.tag());
                        assert(p.tag().isvalid());
                        point_.push_back(p);
                        //auto itt = point_tag_index_map.left.find(p.tag()()); 
                        //assert(itt == point_tag_index_map.left.end());
                        //point_tag_index_map.insert(boost::bimap<int, int>::value_type(p.tag()(), point_.size() - 1));
                        //for (const Tag& _t: point_.back().parent_cell())
                        //{
                        //assert(_t.isvalid());
                        //}

                    }
                    else if (it->tag() == p.tag())
                    {
                        it->add_parent_cell(mc->tag());
                    }
                    else
                    {
                        point_.insert(it, p);
                    }
                }

                for (const MeshPoint& p: cell_.back().point())
                {
                    assert(!p.parent_cell().empty());
                }

                for (const MeshPoint& _p: point_)
                {
                    for (const auto& _t: _p.parent_cell())
                    {
                        assert(_t.isvalid());
                    }
                }
            }

            void Mesh::add_point(MeshPoint p, const Tag& t, size_t size)
            {
                // reserve memory for cell_ if size is provided.
                if (size != 0)
                {
                    //if (point_.capacity() != size)
                    //{
                        //point_.reserve(size);
                    //}
                }

                //auto it = std::find_if(point_.begin(), point_.end(), [&p](const MeshPoint& _p){return _p.tag() == p.tag();});
                //auto it = std::lower_bound(point_.begin(), point_.end(), p, [](const MeshPoint& left, const MeshPoint& right){return left.tag()() < right.tag()();});
                //assert(it == point_.end());
                //if (it == point_.end())
                {
                    // generate a new cell tag.
                    //p.set_tag(pointtag);
                    p.set_tag(t);
                    // set parent mesh.
                    //p.set_parent_mesh(tag_);
                    // add to container.
                    point_.push_back(std::move(p));
                    // insert to bimap.
                    //point_tag_index_map.insert(boost::bimap<int, int>::value_type(p.tag()(), point_.size() - 1));
                }
                /*else
                  {
                  p.set_tag(t);
                  p.set_parent_mesh(tag_);

                  assert(it->tag() != p.tag());
                  point_.insert(it, p);
                //point_tag_index_map.insert(boost::bimap<int, int>::value_type(p.tag()(), point_.size() - 1));
                }*/
            }

            bool Mesh::do_point_exist(const Tag& t) const
            {
                auto it = std::lower_bound(point_.begin(), point_.end(), MeshPoint(t, Tag()));
                if (it != point_.end() && it->tag() == t)
                {
                    return true;
                }

                return false;
            }

            void Mesh::remove_cell_boundary(const Tag& ic)
            {
                const MeshCell& mc = cell(ic);

                auto mark = [&](mcc& container, const Tag& t)
                {
                    auto it = std::find_if(container.begin(), container.end(), [&](const auto& a){return a.tag() == t;});
                    assert(it != container.end());
                    it->mark_to_be_erased();
                };

                auto erase = [&](mcc& container)
                {
                    auto itt = std::remove_if(container.begin(), container.end(), [&](const MeshCell& c){return c.erase() == true;});
                    container.erase(itt, container.end());
                };

                //for (const Tag& t: mc.wall_boundary()) {
                    //mark(wall_boundaries_, t);
                //}
                erase(wall_boundaries_);
                erase(symmetry_boundaries_);

                //for (const Tag& t: mc.dirichlet_boundary()) {
                    //mark(dirichlet_boundaries_, t);
                //}
                erase(dirichlet_boundaries_);

                //for (const Tag& t: mc.farfield_boundary()) {
                    //mark(farfield_boundaries_, t);
                //}
                erase(farfield_boundaries_);

                //for (const Tag& t: mc.empty_boundary()) {
                    //mark(empty_boundaries_, t);
                //}
                erase(empty_boundaries_);

                //for (const Tag& t: mc.interog_boundary()) {
                    //mark(interog_boundaries_, t);
                //}
                erase(interog_boundaries_);

                /*for (const Tag& t: mc.wall_boundary())
                {
                    //auto it = wall_tag_index_map_.left.find(t());
                    //if (it == wall_tag_index_map_.left.end()) {
                        //continue;
                    //}
                    //wall_boundaries_.erase(wall_boundaries_.begin() + it->second);
                    //wall_boundaries_.remove_if([&](const auto& a){return a.tag() == t;});
                    auto it = std::remove_if(wall_boundaries_.begin(), wall_boundaries_.end(), [&](const auto& a){return a.tag() == t;});
                    wall_boundaries_.erase(it, wall_boundaries_.end());

                    //int thres = it->second;
                    //wall_tag_index_map_.left.erase(t());
                    //for (int j=thres+1; j<=wall_boundaries_.size(); ++j)
                    //{
                    //    auto itt = wall_tag_index_map_.right.find(j);
                    //    if (itt == wall_tag_index_map_.right.end())
                    //    {
                    //        std::cout << "size map: " << wall_tag_index_map_.size() << std::endl;
                    //        std::cout << "size vec: " << wall_boundaries_.size() << std::endl;
                    //        std::cout << "j: " << j << std::endl;
                    //        std::cout << "thres: " << thres << std::endl;
                    //        for (auto yy=wall_tag_index_map_.right.begin(); yy != wall_tag_index_map_.right.end(); ++yy)
                    //        {
                    //            std::cout << "yy: " << yy->first << std::endl;
                    //        }
                    //    }
                    //    assert(itt != wall_tag_index_map_.right.end());
                    //    int rep = itt->first - 1;
                    //    bool successful_replace = wall_tag_index_map_.right.replace_key(itt, rep);
                    //    assert(successful_replace);
                    //    //if (wall_boundaries_[rep].tag()() != itt->second)
                    //    //{
                    //    //    std::cout << "wb size: " << wall_boundaries_.size() << std::endl;
                    //    //    std::cout << "map size: " << wall_tag_index_map_.size() << std::endl;
                    //    //    std::cout << "rep: " << wall_boundaries_[rep].tag()() << std::endl;
                    //    //    std::cout << "second: " << itt->second << std::endl;
                    //    //    std::cout << "j: " << j << std::endl;
                    //    //    std::cout << "thres: " << thres << std::endl;
                    //    //}
                    //    //assert(wall_boundaries_[rep].tag()() == itt->second);
                    //}
                    //assert(wall_tag_index_map_.left.count(t()) == 0);
                }

                for (const Tag& t: mc.dirichlet_boundary())
                {
                    //auto it = dirichlet_tag_index_map_.left.find(t());
                    //if (it == dirichlet_tag_index_map_.left.end()) {
                        //continue;
                    //}
                    //dirichlet_boundaries_.erase(dirichlet_boundaries_.begin() + it->second);
                    //dirichlet_boundaries_.remove_if([&](const auto& a){return a.tag() == t;});
                    auto it = std::remove_if(dirichlet_boundaries_.begin(), dirichlet_boundaries_.end(), [&](const auto& a){return a.tag() == t;});
                    dirichlet_boundaries_.erase(it, dirichlet_boundaries_.end());

                    //int thres = it->second;
                    //dirichlet_tag_index_map_.left.erase(t());
                    //for (int j=thres+1; j<=dirichlet_boundaries_.size(); ++j)
                    //{
                    //    auto itt = dirichlet_tag_index_map_.right.find(j);
                    //    if (itt == dirichlet_tag_index_map_.right.end())
                    //    {
                    //        std::cout << "size map: " << dirichlet_tag_index_map_.size() << std::endl;
                    //        std::cout << "size vec: " << dirichlet_boundaries_.size() << std::endl;
                    //        std::cout << "j: " << j << std::endl;
                    //        std::cout << "thres: " << thres << std::endl;
                    //        for (auto yy=dirichlet_tag_index_map_.right.begin(); yy != dirichlet_tag_index_map_.right.end(); ++yy)
                    //        {
                    //            std::cout << "yy: " << yy->first << std::endl;
                    //        }
                    //    }
                    //    assert(itt != dirichlet_tag_index_map_.right.end());
                    //    int rep = itt->first - 1;
                    //    bool successful_replace = dirichlet_tag_index_map_.right.replace_key(itt, rep);
                    //    assert(successful_replace);
                    //    //assert(dirichlet_boundaries_[rep].tag()() == itt->second);
                    //}
                    //assert(dirichlet_tag_index_map_.left.count(t()) == 0);
                }

                for (const Tag& t: mc.farfield_boundary())
                {
                    //auto it = farfield_tag_index_map_.left.find(t());
                    //if (it == farfield_tag_index_map_.left.end()) {
                        //continue;
                    //}
                    //farfield_boundaries_.erase(farfield_boundaries_.begin() + it->second);
                    //farfield_boundaries_.remove_if([&](const auto& a){return a.tag() == t;});
                    auto it = std::remove_if(farfield_boundaries_.begin(), farfield_boundaries_.end(), [&](const auto& a){return a.tag() == t;});
                    farfield_boundaries_.erase(it, farfield_boundaries_.end());

                    //int thres = it->second;
                    //farfield_tag_index_map_.left.erase(t());
                    //for (int j=thres+1; j<=farfield_boundaries_.size(); ++j)
                    //{
                    //    auto itt = farfield_tag_index_map_.right.find(j);
                    //    if (itt == farfield_tag_index_map_.right.end())
                    //    {
                    //        std::cout << "size map: " << farfield_tag_index_map_.size() << std::endl;
                    //        std::cout << "size vec: " << farfield_boundaries_.size() << std::endl;
                    //        std::cout << "j: " << j << std::endl;
                    //        std::cout << "thres: " << thres << std::endl;
                    //        for (auto yy=farfield_tag_index_map_.right.begin(); yy != farfield_tag_index_map_.right.end(); ++yy)
                    //        {
                    //            std::cout << "yy: " << yy->first << std::endl;
                    //        }
                    //    }
                    //    assert(itt != farfield_tag_index_map_.right.end());
                    //    int rep = itt->first - 1;
                    //    bool successful_replace = farfield_tag_index_map_.right.replace_key(itt, rep);
                    //    assert(successful_replace);
                    //    //assert(farfield_boundaries_[rep].tag()() == itt->second);
                    //}
                    //assert(farfield_tag_index_map_.left.count(t()) == 0);
                }
                
                for (const Tag& t: mc.empty_boundary())
                {
                    //auto it = empty_tag_index_map_.left.find(t());
                    //if (it == empty_tag_index_map_.left.end()) {
                        //continue;
                    //}
                    //empty_boundaries_.erase(empty_boundaries_.begin() + it->second);
                    //empty_boundaries_.remove_if([&](const auto& a){return a.tag() == t;});
                    auto it = std::remove_if(empty_boundaries_.begin(), empty_boundaries_.end(), [&](const auto& a){return a.tag() == t;});
                    empty_boundaries_.erase(it, empty_boundaries_.end());

                    //int thres = it->second;
                    //empty_tag_index_map_.left.erase(t());
                    //for (int j=thres+1; j<=empty_boundaries_.size(); ++j)
                    //{
                    //    auto itt = empty_tag_index_map_.right.find(j);
                    //    if (itt == empty_tag_index_map_.right.end())
                    //    {
                    //        std::cout << "size map: " << empty_tag_index_map_.size() << std::endl;
                    //        std::cout << "size vec: " << empty_boundaries_.size() << std::endl;
                    //        std::cout << "j: " << j << std::endl;
                    //        std::cout << "thres: " << thres << std::endl;
                    //        for (auto yy=empty_tag_index_map_.right.begin(); yy != empty_tag_index_map_.right.end(); ++yy)
                    //        {
                    //            std::cout << "yy: " << yy->first << std::endl;
                    //        }
                    //    }
                    //    assert(itt != empty_tag_index_map_.right.end());
                    //    int rep = itt->first - 1;
                    //    bool successful_replace = empty_tag_index_map_.right.replace_key(itt, rep);
                    //    assert(successful_replace);
                    //    //assert(empty_boundaries_[rep].tag()() == itt->second);
                    //}
                    //assert(empty_tag_index_map_.left.count(t()) == 0);
                }

                for (const Tag& t: mc.interog_boundary())
                {
                    //auto it = interog_tag_index_map_.left.find(t());
                    //if (it == interog_tag_index_map_.left.end()) {
                        //continue;
                    //}
                    //interog_boundaries_.erase(interog_boundaries_.begin() + it->second);
                    //interog_boundaries_.remove_if([&](const auto& a){return a.tag() == t;});
                    auto it = std::remove_if(interog_boundaries_.begin(), interog_boundaries_.end(), [&](const auto& a){return a.tag() == t;});
                    interog_boundaries_.erase(it, interog_boundaries_.end());

                    //int thres = it->second;
                    //interog_tag_index_map_.left.erase(t());
                    //for (int j=thres+1; j<=interog_boundaries_.size(); ++j)
                    //{
                    //    auto itt = interog_tag_index_map_.right.find(j);
                    //    if (itt == interog_tag_index_map_.right.end())
                    //    {
                    //        std::cout << "size map: " << interog_tag_index_map_.size() << std::endl;
                    //        std::cout << "size vec: " << interog_boundaries_.size() << std::endl;
                    //        std::cout << "j: " << j << std::endl;
                    //        std::cout << "thres: " << thres << std::endl;
                    //        for (auto yy=interog_tag_index_map_.right.begin(); yy != interog_tag_index_map_.right.end(); ++yy)
                    //        {
                    //            std::cout << "yy: " << yy->first << std::endl;
                    //        }
                    //    }
                    //    assert(itt != interog_tag_index_map_.right.end());
                    //    int rep = itt->first - 1;
                    //    bool successful_replace = interog_tag_index_map_.right.replace_key(itt, rep);
                    //    assert(successful_replace);
                    //    //assert(interog_boundaries_[rep].tag()() == itt->second);
                    //}
                    //assert(interog_tag_index_map_.left.count(t()) == 0);
                }*/
            }

            void Mesh::reset_erase_marks()
            {
                for (MeshCell& mc: cell_)
                {
                    mc.unmark_to_be_erased();
                }

                for (MeshCell& mc: wall_boundaries_)
                {
                    mc.unmark_to_be_erased();
                }

                for (MeshCell& mc: symmetry_boundaries_)
                {
                    mc.unmark_to_be_erased();
                }

                for (MeshCell& mc: dirichlet_boundaries_)
                {
                    mc.unmark_to_be_erased();
                }

                for (MeshCell& mc: farfield_boundaries_)
                {
                    mc.unmark_to_be_erased();
                }

                for (MeshCell& mc: empty_boundaries_)
                {
                    mc.unmark_to_be_erased();
                }

                for (MeshCell& mc: interog_boundaries_)
                {
                    mc.unmark_to_be_erased();
                }
            }

            //void Mesh::mark_to_be_erased(const MeshCell& mc)
            //{
                //cell_p(ic).mark_to_be_erased();
            //}

            void Mesh::mark_to_be_erased(const Tag& ic)
            {
                cell_p(ic).mark_to_be_erased();
            }

            void Mesh::erase_marked_cells(int rank)
            {
                for (MeshCell& mc: cell_)
                {
                    if (mc.erase())
                    {
                        assert(mc.oga_cell_type() != OGA_cell_type_t::ghost);
                        prepare_to_remove_cell(mc, rank);
                    }
                }
                auto it = std::remove_if(wall_boundaries_.begin(), wall_boundaries_.end(), [&](const MeshCell& mc){return mc.erase() == true;});
                wall_boundaries_.erase(it, wall_boundaries_.end());

                it = std::remove_if(symmetry_boundaries_.begin(), symmetry_boundaries_.end(), [&](const MeshCell& mc){return mc.erase() == true;});
                symmetry_boundaries_.erase(it, symmetry_boundaries_.end());

                it = std::remove_if(dirichlet_boundaries_.begin(), dirichlet_boundaries_.end(), [&](const MeshCell& mc){return mc.erase() == true;});
                dirichlet_boundaries_.erase(it, dirichlet_boundaries_.end());

                it = std::remove_if(farfield_boundaries_.begin(), farfield_boundaries_.end(), [&](const MeshCell& mc){return mc.erase() == true;});
                farfield_boundaries_.erase(it, farfield_boundaries_.end());

                it = std::remove_if(empty_boundaries_.begin(), empty_boundaries_.end(), [&](const MeshCell& mc){return mc.erase() == true;});
                empty_boundaries_.erase(it, empty_boundaries_.end());

                it = std::remove_if(interog_boundaries_.begin(), interog_boundaries_.end(), [&](const MeshCell& mc){return mc.erase() == true;});
                interog_boundaries_.erase(it, interog_boundaries_.end());

                it = std::remove_if(cell_.begin(), cell_.end(), [&](const MeshCell& mc){return mc.erase() == true;});
                cell_.erase(it, cell_.end());
            }

            void Mesh::disconnect_orphan_faces()
            {
                //for (auto it=face_.begin(); it != face_.end();)
                //{
                //    if (it->parent_cell().empty())
                //    {
                //        face_.erase(it);
                //    }
                //    else
                //    {
                //        ++it;
                //    }
                //}
//mesh: 1
            //2 mf: 
              //3 gmf: 4859 7570 10891
                //4 addr: 0x28e45220

                //face_.erase(std::remove_if(face_.begin(), face_.end(), [&](const auto& f){return f.parent_cell().empty();}), face_.end());

                //face_.remove_if([&](const auto& f){return f.parent_cell().empty();});


                /*for (auto& mf: face_)
                {
                    if (mf.is_boundary())
                    {
                        continue;
                    }

                    for (const Tag& pct: mf.parent_cell())
                    {
                        if (query(pct) == nullptr)
                        {
                            assert(false);
                            mf.remove_parent_cell(pct);
                            continue;
                        }

                        auto& pc = cell_p(pct);

                        auto itt = std::find_if(pc.face().begin(), pc.face().end(), [&](const auto& ff){return ff.tag() == mf.tag();});
                        if (itt == pc.face().end())
                        {
                            assert(false);
                            mf.remove_parent_cell(pct);
                        }
                    }
                }

                for (auto& mc: cell_)
                {
                    for (auto& mf: mc.face_)
                    {
                        if (mf.is_boundary())
                        {
                            continue;
                        }

                        for (const Tag& pc: mf.parent_cell())
                        {
                            if (pc == mc.tag()) {
                                continue;
                            }

                            if (query(pc) == nullptr)
                            {
                                assert(false);
                                mf.remove_parent_cell(pc);
                            }
                        }
                    }
                }*/
            }

            void Mesh::prepare_to_remove_ghost(const Tag& ic)
            {
                assert(ic.isvalid());

                auto mc = cell_p(ic);

                // remove the cell from parency of global faces.
                for (MeshFace& mf: mc.face_p())
                {
                    mf.remove_parent_cell(ic);
                }

                // remove the cell from parency of global points.
                for (const MeshPoint& mp: mc.point())
                {
                    auto p = point_p(mp.tag());
                    assert(p != nullptr);
                    p->remove_parent_cell(ic);
                }

                //deparent_cell_from_vertices(ic);
                mc.deparent_self_from_faces(); // not needed. for testing.

                // remove self from neighbors' faces.
                for (auto& inei: mc.pnei())
                {
                    auto& nei = cell_p(inei);
                    nei.deparent_cell_from_faces(ic);
                    nei.remove_neighbor(ic);
                }

                mc.remove_all_neighbors();
                //remove_cell_from_cellhood(ic); // needed to update neighbors.
                //remove_cell_boundary(ic);
            }

            void Mesh::prepare_to_remove_cell(MeshCell& mc, int rank)
            //void Mesh::prepare_to_remove_cell(const Tag& ic, int rank)
            {
                auto ic = mc.tag();
                //assert(ic.isvalid());

                // remove the cell from parency of global faces.
                //MeshCell& mc = cell_p(ic);
                //if (mc.oga_cell_type() != OGA_cell_type_t::ghost && mc.oga_cell_type() != OGA_cell_type_t::non_resident)
                {
                    //for (const MeshFace& mf: cell(ic).face())
                    for (MeshFace& mf: mc.face_p())
                    {
                        //assert(mf.parent_cell().size() == face(mf.tag()).parent_cell().size());

                        //auto it = mf.faceaddr();
                        //assert(it != nullptr);
                        //auto it = std::find_if(face_.begin(), face_.end(), [&](const MeshFace& f){return f.tag() == mf.tag();});
                        //assert(it != face_.end());
                        if (mf.is_boundary())
                        {
                            //it->remove_parent_cells();
                            mf.remove_parent_cells();
                        }
                        else
                        {
                            //it->remove_parent_cell(ic);
                            mf.remove_parent_cell(ic);
                        }
                    }
                }

                // remove the cell from parency of global points.
                for (const MeshPoint& mp: mc.point())
                {
                    //auto it = std::find_if(point_.begin(), point_.end(), [&](const MeshPoint& p){return p.tag() == mp.tag();});
                    //assert(it != point_.end());
                    //it->remove_parent_cell(ic);
                    auto p = point_p(mp.tag());
                    assert(p != nullptr);
                    p->remove_parent_cell(ic);
                }
                
                //deparent_cell_from_vertices(ic);
                //cell_p(ic).deparent_self_from_faces();
                //cell_p(ic).deparent_neis_from_faces(); // except boundary faces. // is this needed?

                // remove self from neighbors' faces.
                for (auto& inei: mc.pnei())
                {
                    auto& nei = cell_p(inei);
                    //auto& nei = *inei.addr();
                    nei.deparent_cell_from_faces(ic);
                    nei.remove_neighbor(ic);
                }

                //for (auto& mf: cell_p(ic).face())
                //{
                    //auto iter = std::find_if(face_.begin(), face_.end(), [&](const auto& f){return f.tag() == mf.tag();});
                    //assert(iter != face_.end());
                    //iter->remove_parent_cell(ic);
                //}

                mc.remove_all_neighbors();
                //remove_cell_from_cellhood(ic); // needed to update neighbors.
                //remove_cell_boundary(ic); //TODO

                //for (const auto& mf: mc.face())
                //{
                //    if (mf.faceaddr() != nullptr)
                //    {
                //        assert(mf.tag() == mf.faceaddr()->tag());
                //    }
                //}
            }

            void Mesh::remove_ghosts()
            {
                std::vector<bool> erase(cell_.size(), false);

                //for (auto mc = cell_.begin(); mc != cell_.end();)
                for (auto mc = cell_.begin(); mc != cell_.end(); ++mc)
                {
                    if (mc->oga_cell_type() == OGA_cell_type_t::ghost)
                    {
                        prepare_to_remove_ghost(mc->tag());
                        erase[std::distance(cell_.begin(), mc)] = true;
                        //mc = cell_.erase(mc);
                    }
                    //else
                    //{
                        //++mc;
                    //}
                }

                int i = 0;
                auto from = std::remove_if(cell_.begin(), cell_.end(), [&erase, &i](auto& mc){++i; return erase[i-1] == true;});
                cell_.erase(from, cell_.end());
            }

            void Mesh::remove_nonresidents()
            {
                std::vector<bool> erase(cell_.size(), false);

                //for (auto mc = cell_.begin(); mc != cell_.end();)
                for (auto mc = cell_.begin(); mc != cell_.end(); ++mc)
                {
                    if (mc->oga_cell_type() == OGA_cell_type_t::non_resident)
                    {
                        prepare_to_remove_cell(*mc, -1);
                        erase[std::distance(cell_.begin(), mc)] = true;
                        //mc = cell_.erase(mc);
                    }
                    //else
                    //{
                        //++mc;
                    //}
                }

                int i = 0;
                auto from = std::remove_if(cell_.begin(), cell_.end(), [&erase, &i](auto& mc){++i; return erase[i-1] == true;});
                cell_.erase(from, cell_.end());
            }

void Mesh::remove_cell(const Tag& ic)
{
    auto mc = query_itp(ic);
    assert(mc != cell_.end());
    prepare_to_remove_cell(*mc, -1);

    // erase cell.
    //auto it = cell_tag_index_map.left.find(ic());
    //assert(it != cell_tag_index_map.left.end());
    //cell_.erase(cell_.begin() + it->second);

    //auto it = std::lower_bound(cell_.begin(), cell_.end(), MeshCell(ic, Tag()));
                //assert(it != cell_.end() && it->tag() == ic);
                cell_.erase(mc);

                // remove from cell map.
                // example:
                // a b c d e f g h
                // 0 1 2 3 4 5 6 7
                // 0 1 2 x 4 5 6 7
                // thres+1 = 4, cell_.size()=7
                /*int thres = it->second;
                  cell_tag_index_map.left.erase(ic());
                  for (int j=thres+1; j<=cell_.size(); ++j)
                  {
                  auto itt = cell_tag_index_map.right.find(j);
                  assert(itt != cell_tag_index_map.right.end());
                  int rep = itt->first - 1;
                  bool successful_replace = cell_tag_index_map.right.replace_key(itt, rep);
                  assert(successful_replace);
                  assert(cell_[rep].tag()() == itt->second);
                //assert(query(itt->second));
                }
                assert(cell_tag_index_map.left.count(ic()) == 0);*/

            }

            void Mesh::remove_isolated_points()
            {
                for (MeshPoint& mp: point_)
                {
                    if (mp.parent_cell().empty())
                    {
                        //for (const auto& mc: cell_)
                        //{
                            //for (const auto& p: mc.point())
                            //{
                                //assert(p.tag() != mp.tag());
                            //}
                        //}

                        mp.mark_to_be_erased();
                    }
                }

                auto it = std::remove_if(point_.begin(), point_.end(), [&](const MeshPoint& mp){return mp.erase() == true;});
                point_.erase(it, point_.end());
            }

            void Mesh::remove_point(Tag ip)
            {
                assert(ip.isvalid());

                auto it = std::lower_bound(point_.begin(), point_.end(), MeshPoint(ip, Tag()));
                assert(it != point_.end() && it->tag() == ip);

                point_.erase(it);

                /*// remove from map.
                  int thres = it->second;
                  point_tag_index_map.left.erase(ip());
                  for (int j=thres+1; j<=point_.size(); ++j)
                  {
                  auto itt = point_tag_index_map.right.find(j);
                  assert(itt != point_tag_index_map.right.end()); 
                  int rep = itt->first - 1;
                  bool successful_replace = point_tag_index_map.right.replace_key(itt, rep);
                  assert(successful_replace);
                  }*/
            }

            size_t Mesh::calc_nhole_cell() const
            {
                size_t nhole_cell = 0;
                for (const MeshCell& mc: cell_)
                {
                    if (mc.oga_cell_type() == OGA_cell_type_t::hole)
                        ++nhole_cell;
                }

                return nhole_cell;
            }

int Mesh::priority() const
{
    return priority_;
}

            //Mesh::Mesh(const Tag& tag, const Tag& parent_mesh): tag_(tag), parent_mesh_(parent_mesh)
            Mesh::Mesh(const Tag& tag): tag_(tag), priority_(-1)
            {
                //if (file_name != "")
                //{
                //if (rank != 0)
                //{
                //file_name.append("_");
                //file_name.append(std::to_string(rank));
                //}
                //file_name.append(".msh");
                ////std::cout << "file name: " << file_name << std::endl;
                //read_mesh_GMSH(file_name);
                //}

                //hole_map_.set_nstripe(1, 1);

                //determine_hole_aabb();
            }

            //void Mesh::print_as_vtk(std::string file_name, const std::deque<std::pair<Tag, Tag>>& include) const
            //{    
            //int cell_list_size = 0;
            //std::ofstream out;    

            //out.open (file_name);

            //out << "# vtk DataFile Version 3.0" << std::endl;
            //out << "All in VTK format" << std::endl;
            //out << "ASCII" << std::endl;
            //out << "DATASET UNSTRUCTURED_GRID" << std::endl;
            //out << "POINTS " << point().size() << " float" << std::endl;

            //for (int i=0; i<point_.size(); ++i)
            //{
            //out << point_[i].p().r(0);
            //out << " ";
            //out << point_[i].p().r(1);
            //out << " ";
            //out << 0.;
            //out << std::endl;
            //}

            //// get cell list size.
            //int cell_size = 0;
            //for (int i=0; i<cell_.size(); ++i)    
            //{
            //int count = std::count_if(include.begin(), include.end(), [&](const std::pair<Tag, Tag>& pr){return true; (tag() == pr.first && cell_[i].tag() == pr.second);});
            //if (count == 0) continue;
            ////if (!cell_[i].is_ghost())
            //cell_list_size += (cell_[i].point().size() + 1);
            //++cell_size;
            //}

            //out << std::endl;    
            //out << "CELLS " << cell_size << " " << cell_list_size << std::endl;
            ////out << "CELLS " << n_interior_cell_ << " " << cell_list_size << std::endl;

            ////int counter = 0;
            ////for (int i=0; i<cell_.size(); ++i)    
            ////{        
            ////if (!cell_[i].is_ghost())
            ////++ counter;
            ////}

            //for (int i=0; i<cell_.size(); ++i)    
            //{        
            //int count = std::count_if(include.begin(), include.end(), [&](const std::pair<Tag, Tag>& pr){return (tag() == pr.first && cell_[i].tag() == pr.second);});
            //if (count == 0) continue;
            ////if (!cell_[i].is_ghost())
            //{
            //out << cell_[i].point().size();
            //out << " ";

            //for (int j=0; j<cell_[i].point().size(); ++j)
            //{
            //out << point_index(cell_[i].point(j).tag());
            //out << " ";
            //}

            //out << std::endl;
            //}
            //}    

            //out << "CELL_TYPES " << cell_size << std::endl;
            ////out << "CELL_TYPES " << n_interior_cell_ << std::endl;
            //for (int i=0; i<cell_.size(); ++i)
            //{
            //int count = std::count_if(include.begin(), include.end(), [&](const std::pair<Tag, Tag>& pr){return (tag() == pr.first && cell_[i].tag() == pr.second);});
            //if (count == 0) continue;

            ////if (!cell_[i].is_ghost())
            //{
            //if (cell_[i].polygon().vertex().size() == 3)
            ////if (cell_[i].polytope()->vertex().size() == 3)
            //{
            //out << 5;
            //}
            //else if (cell_[i].polygon().vertex().size() == 4)
            ////else if (cell_[i].polytope()->vertex().size() == 4)
            //{
            //out << 9;
            //}        
            //else if (cell_[i].polygon().vertex().size() == 2)
            ////else if (cell_[i].polytope()->vertex().size() == 2)
            //{
            //out << 3;
            //}        

            //out << std::endl;
            //}
            //}

            //out << "CELL_DATA " << cell_size << std::endl;
            ////out << "CELL_DATA " << n_interior_cell_ << std::endl;
            //out << "SCALARS " << "oga_cell_status " << "int " << "1" << std::endl;
            //out << "LOOKUP_TABLE default" << std::endl;    
            //for (int i=0; i<cell_.size(); ++i)
            //{
            //int count = std::count_if(include.begin(), include.end(), [&](const std::pair<Tag, Tag>& pr){return (tag() == pr.first && cell_[i].tag() == pr.second);});
            //if (count == 0) continue;
            ////if (!cell_[i].is_ghost())
            //out << cell_[i].oga_cell_type() << std::endl;
            //} 
            //out << "SCALARS " << "cell_tag " << "int " << "1" << std::endl;
            //out << "LOOKUP_TABLE default" << std::endl;    
            //for (int i=0; i<cell_.size(); ++i)
            //{
            //int count = std::count_if(include.begin(), include.end(), [&](const std::pair<Tag, Tag>& pr){return (tag() == pr.first && cell_[i].tag() == pr.second);});
            //if (count == 0) continue;
            ////if (!cell_[i].is_ghost())
            //out << cell_[i].tag()() << std::endl;
            //} 

            //out.close();
            //}

            void Mesh::print_interog_as_vtk(std::string file_name) const
            {
                int cell_list_size = 0;
                std::ofstream out;    

                out.open (file_name);

                out << "# vtk DataFile Version 3.0" << std::endl;
                out << "All in VTK format" << std::endl;
                out << "ASCII" << std::endl;
                out << "DATASET UNSTRUCTURED_GRID" << std::endl;
                out << "POINTS " << point().size() << " float" << std::endl;

                for (int i=0; i<point_.size(); ++i)
                {
                    out << point_[i].p().r(0);
                    out << " ";
                    out << point_[i].p().r(1);
                    out << " ";
                    out << point_[i].p().r(2);
                    out << std::endl;
                }

                // get cell list size.
                for (auto mc = interog_boundaries_.begin(); mc != interog_boundaries_.end(); ++mc)
                {
                    cell_list_size += (mc->point().size() + 1);
                }

                out << std::endl;    
                out << "CELLS " << interog_boundaries_.size() << " " << cell_list_size << std::endl;

                for (auto mc = interog_boundaries_.begin(); mc != interog_boundaries_.end(); ++mc)
                {        
                    out << mc->point().size();
                    out << " ";

                    for (int j=0; j<mc->point().size(); ++j)
                    {
                        out << point_index(mc->point(j).tag());
                        out << " ";
                    }

                    out << std::endl;
                }    

                out << "CELL_TYPES " << interog_boundaries_.size() << std::endl;
                for (auto mc = interog_boundaries_.begin(); mc != interog_boundaries_.end(); ++mc)
                {
                    //https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf

                    if (mc->poly().shape() == Shape::tri)
                    {
                        out << 5;
                    }
                    else if (mc->poly().shape() == Shape::quad)
                    {
                        out << 9;
                    }        
                    else
                    {
                        assert(false);
                    }

                    out << std::endl;
                }

                out << "CELL_DATA " << interog_boundaries_.size() << std::endl;
                out << "SCALARS " << "oga_cell_status " << "int " << "1" << std::endl;
                out << "LOOKUP_TABLE default" << std::endl;    
                for (auto mc = interog_boundaries_.begin(); mc != interog_boundaries_.end(); ++mc)
                {
                    out << 10 << std::endl;
                } 
                out << "SCALARS " << "cell_tag " << "int " << "1" << std::endl;
                out << "LOOKUP_TABLE default" << std::endl;    
                for (auto mc = interog_boundaries_.begin(); mc != interog_boundaries_.end(); ++mc)
                {
                    out << mc->tag()() << std::endl;
                } 
                out << "SCALARS " << "rho " << "float " << "1" << std::endl;
                out << "LOOKUP_TABLE default" << std::endl;
                for (auto mc = interog_boundaries_.begin(); mc != interog_boundaries_.end(); ++mc)
                {
                    out << mc->prim(0) << std::endl;
                }
                out << "SCALARS " << "u " << "float " << "1" << std::endl;
                out << "LOOKUP_TABLE default" << std::endl;
                for (auto mc = interog_boundaries_.begin(); mc != interog_boundaries_.end(); ++mc)
                {
                    out << mc->prim(1) << std::endl;
                }
                out << "SCALARS " << "v " << "float " << "1" << std::endl;
                out << "LOOKUP_TABLE default" << std::endl;
                for (auto mc = interog_boundaries_.begin(); mc != interog_boundaries_.end(); ++mc)
                {
                    out << mc->prim(2) << std::endl;
                }
                out << "SCALARS " << "w " << "float " << "1" << std::endl;
                out << "LOOKUP_TABLE default" << std::endl;
                for (auto mc = interog_boundaries_.begin(); mc != interog_boundaries_.end(); ++mc)
                {
                    out << mc->prim(3) << std::endl;
                }
                out << "SCALARS " << "p " << "float " << "1" << std::endl;
                out << "LOOKUP_TABLE default" << std::endl;
                for (auto mc = interog_boundaries_.begin(); mc != interog_boundaries_.end(); ++mc)
                {
                    out << mc->prim(4) << std::endl;
                }

                out.close();
            }

            void Mesh::print_wall_as_vtk(std::string file_name) const
            {
                int cell_list_size = 0;
                std::ofstream out;    

                out.open (file_name);

                out << "# vtk DataFile Version 3.0" << std::endl;
                out << "All in VTK format" << std::endl;
                out << "ASCII" << std::endl;
                out << "DATASET UNSTRUCTURED_GRID" << std::endl;
                out << "POINTS " << point().size() << " float" << std::endl;

                for (int i=0; i<point_.size(); ++i)
                {
                    out << point_[i].p().r(0);
                    out << " ";
                    out << point_[i].p().r(1);
                    out << " ";
                    out << point_[i].p().r(2);
                    out << std::endl;
                }

                // get cell list size.
                for (auto mc = wall_boundaries_.begin(); mc != wall_boundaries_.end(); ++mc)
                {
                    cell_list_size += (mc->point().size() + 1);
                }

                out << std::endl;    
                out << "CELLS " << wall_boundaries_.size() << " " << cell_list_size << std::endl;

                for (auto mc = wall_boundaries_.begin(); mc != wall_boundaries_.end(); ++mc)
                {        
                    out << mc->point().size();
                    out << " ";

                    for (int j=0; j<mc->point().size(); ++j)
                    {
                        out << point_index(mc->point(j).tag());
                        out << " ";
                    }

                    out << std::endl;
                }    

                out << "CELL_TYPES " << wall_boundaries_.size() << std::endl;
                for (auto mc = wall_boundaries_.begin(); mc != wall_boundaries_.end(); ++mc)
                {
                    //https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf

                    if (mc->poly().shape() == Shape::tri)
                    {
                        out << 5;
                    }
                    else if (mc->poly().shape() == Shape::quad)
                    {
                        out << 9;
                    }        
                    else
                    {
                        assert(false);
                    }

                    out << std::endl;
                }

                out << "CELL_DATA " << wall_boundaries_.size() << std::endl;
                out << "SCALARS " << "oga_cell_status " << "int " << "1" << std::endl;
                out << "LOOKUP_TABLE default" << std::endl;    
                for (auto mc = wall_boundaries_.begin(); mc != wall_boundaries_.end(); ++mc)
                {
                    out << 10 << std::endl;
                } 
                out << "SCALARS " << "cell_tag " << "int " << "1" << std::endl;
                out << "LOOKUP_TABLE default" << std::endl;    
                for (auto mc = wall_boundaries_.begin(); mc != wall_boundaries_.end(); ++mc)
                {
                    out << mc->tag()() << std::endl;
                } 
                out << "SCALARS " << "rho " << "float " << "1" << std::endl;
                out << "LOOKUP_TABLE default" << std::endl;
                for (auto mc = wall_boundaries_.begin(); mc != wall_boundaries_.end(); ++mc)
                {
                    out << mc->prim(0) << std::endl;
                }
                out << "SCALARS " << "u " << "float " << "1" << std::endl;
                out << "LOOKUP_TABLE default" << std::endl;
                for (auto mc = wall_boundaries_.begin(); mc != wall_boundaries_.end(); ++mc)
                {
                    out << mc->prim(1) << std::endl;
                }
                out << "SCALARS " << "v " << "float " << "1" << std::endl;
                out << "LOOKUP_TABLE default" << std::endl;
                for (auto mc = wall_boundaries_.begin(); mc != wall_boundaries_.end(); ++mc)
                {
                    out << mc->prim(2) << std::endl;
                }
                out << "SCALARS " << "w " << "float " << "1" << std::endl;
                out << "LOOKUP_TABLE default" << std::endl;
                for (auto mc = wall_boundaries_.begin(); mc != wall_boundaries_.end(); ++mc)
                {
                    out << mc->prim(3) << std::endl;
                }
                out << "SCALARS " << "p " << "float " << "1" << std::endl;
                out << "LOOKUP_TABLE default" << std::endl;
                for (auto mc = wall_boundaries_.begin(); mc != wall_boundaries_.end(); ++mc)
                {
                    out << mc->prim(4) << std::endl;
                }
                out << "NORMALS " << "n " << "float " << std::endl;
                for (auto mc = wall_boundaries_.begin(); mc != wall_boundaries_.end(); ++mc)
                {
                    const MeshCell& nei = cell(mc->interior_boundary());
                    auto mf = std::find_if(nei.face().begin(), nei.face().end(), [&](const auto& f){return f.tag() == mc->face()[0].tag();});
                    assert(mf != nei.face().end());
                    const auto& v0 = mf->face().vertex(0).r();
                    const auto& v1 = mf->face().vertex(1).r();
                    const auto& v2 = mf->face().vertex(2).r();
                    Vector3 cr = cross((v0-v1), (v2-v1));
                    if (len(cr) == 0.)
                    {
                        std::cout << "neiiiiiiiii: " << nei.tag()() << std::endl;
                    }
                    auto n = mf->face().normal();
                    out << n(0) << " " << n(1) << " " << n(2) << std::endl;
                }

                out.close();
            }

            void Mesh::print_as_vtk_geometry(std::string file_name) const
            {
                int cell_list_size = 0;
                std::ofstream out;    

                out.open (file_name);

                out << "# vtk DataFile Version 3.0" << std::endl;
                out << "All in VTK format" << std::endl;
                out << "ASCII" << std::endl;
                out << "DATASET UNSTRUCTURED_GRID" << std::endl;
                out << "POINTS " << point().size() << " float" << std::endl;

                for (int i=0; i<point_.size(); ++i)
                {
                    out << point_[i].p().r(0);
                    out << " ";
                    out << point_[i].p().r(1);
                    out << " ";
                    out << point_[i].p().r(2);
                    out << std::endl;
                }

                // get cell list size.
                for (auto mc = cell_.begin(); mc != cell_.end(); ++mc)
                {
                    cell_list_size += (mc->point().size() + 1);
                }

                out << std::endl;    
                out << "CELLS " << cell_.size() << " " << cell_list_size << std::endl;

                for (auto mc = cell_.begin(); mc != cell_.end(); ++mc)
                {        
                    out << mc->point().size();
                    out << " ";

                    for (int j=0; j<mc->point().size(); ++j)
                    {
                        out << point_index(mc->point(j).tag());
                        out << " ";
                    }

                    out << std::endl;
                }    

                out << "CELL_TYPES " << cell_.size() << std::endl;
                for (auto mc = cell_.begin(); mc != cell_.end(); ++mc)
                {
                    //https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf

                    if (mc->poly().shape() == Shape::tet)
                    {
                        out << 10;
                    }
                    else if (mc->poly().shape() == Shape::hex)
                    {
                        out << 12;
                    }        
                    else if (mc->poly().shape() == Shape::pri)
                    {
                        out << 13;
                    }        
                    else if (mc->poly().shape() == Shape::tri)
                    {
                        out << 5;
                    }        
                    else if (mc->poly().shape() == Shape::quad)
                    {
                        out << 9;
                    }        
                    else
                    {
                        assert(false);
                    }

                    out << std::endl;
                }

                out << "CELL_DATA " << cell().size() << std::endl;
                out << "SCALARS " << "oga_cell_status " << "int " << "1" << std::endl;
                out << "LOOKUP_TABLE default" << std::endl;    
                for (auto mc = cell_.begin(); mc != cell_.end(); ++mc)
                {
                    out << mc->oga_cell_type() << std::endl;
                } 
                out << "SCALARS " << "cell_tag " << "int " << "1" << std::endl;
                out << "LOOKUP_TABLE default" << std::endl;    
                for (auto mc = cell_.begin(); mc != cell_.end(); ++mc)
                {
                    out << mc->tag()() << std::endl;
                } 

                out.close();
            }

            void Mesh::print_as_vtk(std::string file_name) const
            {
                int cell_list_size = 0;
                std::ofstream out;    

                out.open (file_name);

                out << "# vtk DataFile Version 3.0" << std::endl;
                out << "All in VTK format" << std::endl;
                out << "ASCII" << std::endl;
                out << "DATASET UNSTRUCTURED_GRID" << std::endl;
                out << "POINTS " << point().size() << " float" << std::endl;

                for (int i=0; i<point_.size(); ++i)
                {
                    out << point_[i].p().r(0);
                    out << " ";
                    out << point_[i].p().r(1);
                    out << " ";
                    out << point_[i].p().r(2);
                    out << std::endl;
                }

                // get cell list size.
                for (auto mc = cell_.begin(); mc != cell_.end(); ++mc)
                {
                    cell_list_size += (mc->point().size() + 1);
                }

                out << std::endl;    
                out << "CELLS " << cell_.size() << " " << cell_list_size << std::endl;

                for (auto mc = cell_.begin(); mc != cell_.end(); ++mc)
                {        
                    out << mc->point().size();
                    out << " ";

                    for (int j=0; j<mc->point().size(); ++j)
                    {
                        out << point_index(mc->point(j).tag());
                        out << " ";
                    }

                    out << std::endl;
                }    

                out << "CELL_TYPES " << cell_.size() << std::endl;
                for (auto mc = cell_.begin(); mc != cell_.end(); ++mc)
                {
                    //https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf

                    if (mc->poly().shape() == Shape::tet)
                    {
                        out << 10;
                    }
                    else if (mc->poly().shape() == Shape::hex)
                    {
                        out << 12;
                    }        
                    else if (mc->poly().shape() == Shape::pri)
                    {
                        out << 13;
                    }        
                    else if (mc->poly().shape() == Shape::tri)
                    {
                        out << 5;
                    }        
                    else if (mc->poly().shape() == Shape::quad)
                    {
                        out << 9;
                    }        
                    else
                    {
                        assert(false);
                    }

                    out << std::endl;
                }

                out << "CELL_DATA " << cell().size() << std::endl;
                out << "SCALARS " << "oga_cell_status " << "int " << "1" << std::endl;
                out << "LOOKUP_TABLE default" << std::endl;    
                for (auto mc = cell_.begin(); mc != cell_.end(); ++mc)
                {
                    out << mc->oga_cell_type() << std::endl;
                } 
                out << "SCALARS " << "cell_tag " << "int " << "1" << std::endl;
                out << "LOOKUP_TABLE default" << std::endl;    
                for (auto mc = cell_.begin(); mc != cell_.end(); ++mc)
                {
                    out << mc->tag()() << std::endl;
                } 
                out << "SCALARS " << "btype " << "int " << "1" << std::endl;
                out << "LOOKUP_TABLE default" << std::endl;    
                for (auto mc = cell_.begin(); mc != cell_.end(); ++mc)
                {
                    out << static_cast<int>(mc->btype()) << std::endl;
                } 
                out << "SCALARS " << "rho " << "float " << "1" << std::endl;
                out << "LOOKUP_TABLE default" << std::endl;
                assert(!cell().empty());
                for (const MeshCell& mc: cell())
                {
                    out << mc.prim(0) << std::endl;
                }

                out << "VECTORS " << "u " << "float " << std::endl;
                assert(!cell().empty());
                for (const MeshCell& mc: cell())
                {
                    out << mc.prim(1) << " " << mc.prim(2) << " " << mc.prim(3) << std::endl;
                }

                out << "SCALARS " << "donor_mesh " << "int " << "1" << std::endl;
                out << "LOOKUP_TABLE default" << std::endl;    
                for (auto mc = cell_.begin(); mc != cell_.end(); ++mc)
                {
                    //if (mc->oga_cell_type() == OGA_cell_type_t::receptor || mc->oga_cell_type() == OGA_cell_type_t::mandat_receptor)
                    //{
                        out << mc->donor().mesh_tag_() << std::endl;
                    //}
                    //else
                    //{
                        //out << -1 << std::endl;
                    //}
                } 

                out << "SCALARS " << "donor_cell " << "int " << "1" << std::endl;
                out << "LOOKUP_TABLE default" << std::endl;    
                for (auto mc = cell_.begin(); mc != cell_.end(); ++mc)
                {
                    //if (mc->oga_cell_type() == OGA_cell_type_t::receptor || mc->oga_cell_type() == OGA_cell_type_t::mandat_receptor)
                    //{
                        out << mc->donor().cell_tag_() << std::endl;
                    //}
                    //else
                    //{
                        //out << -1 << std::endl;
                    //}
                } 

                out << "SCALARS " << "receptor_mesh " << "int " << "1" << std::endl;
                out << "LOOKUP_TABLE default" << std::endl;    
                for (auto mc = cell_.begin(); mc != cell_.end(); ++mc)
                {
                    //if (mc->oga_cell_type() == OGA_cell_type_t::receptor || mc->oga_cell_type() == OGA_cell_type_t::mandat_receptor)
                    //{
                        out << mc->receptor().mesh_tag_() << std::endl;
                    //}
                    //else
                    //{
                        //out << -1 << std::endl;
                    //}
                } 

                out << "SCALARS " << "receptor_cell " << "int " << "1" << std::endl;
                out << "LOOKUP_TABLE default" << std::endl;    
                for (auto mc = cell_.begin(); mc != cell_.end(); ++mc)
                {
                    //if (mc->oga_cell_type() == OGA_cell_type_t::receptor || mc->oga_cell_type() == OGA_cell_type_t::mandat_receptor)
                    //{
                        out << mc->receptor().cell_tag_() << std::endl;
                    //}
                    //else
                    //{
                        //out << -1 << std::endl;
                    //}
                } 

                //out << "SCALARS " << "u " << "float " << "1" << std::endl;
                //out << "LOOKUP_TABLE default" << std::endl;
                //assert(!cell().empty());
                //for (const MeshCell& mc: cell())
                //{
                //    out << mc.prim(1) << std::endl;
                //}

                //out << "SCALARS " << "v " << "float " << "1" << std::endl;
                //out << "LOOKUP_TABLE default" << std::endl;
                //assert(!cell().empty());
                //for (const MeshCell& mc: cell())
                //{
                //    out << mc.prim(2) << std::endl;
                //}

                //out << "SCALARS " << "w " << "float " << "1" << std::endl;
                //out << "LOOKUP_TABLE default" << std::endl;
                //assert(!cell().empty());
                //for (const MeshCell& mc: cell())
                //{
                //    out << mc.prim(3) << std::endl;
                //}

                out << "SCALARS " << "p " << "float " << "1" << std::endl;
                out << "LOOKUP_TABLE default" << std::endl;
                assert(!cell().empty());
                for (const MeshCell& mc: cell())
                {
                    out << mc.prim(4) << std::endl;
                }

                out << "SCALARS " << "q0 " << "float " << "1" << std::endl;
                out << "LOOKUP_TABLE default" << std::endl;
                assert(!cell().empty());
                for (const MeshCell& mc: cell())
                {
                    out << mc.cons_sp1_(0) << std::endl;
                }

                out << "SCALARS " << "q1 " << "float " << "1" << std::endl;
                out << "LOOKUP_TABLE default" << std::endl;
                assert(!cell().empty());
                for (const MeshCell& mc: cell())
                {
                    out << mc.cons_sp1_(1) << std::endl;
                }

                out << "SCALARS " << "q2 " << "float " << "1" << std::endl;
                out << "LOOKUP_TABLE default" << std::endl;
                assert(!cell().empty());
                for (const MeshCell& mc: cell())
                {
                    out << mc.cons_sp1_(2) << std::endl;
                }

                out << "SCALARS " << "q3 " << "float " << "1" << std::endl;
                out << "LOOKUP_TABLE default" << std::endl;
                assert(!cell().empty());
                for (const MeshCell& mc: cell())
                {
                    out << mc.cons_sp1_(3) << std::endl;
                }

                out << "SCALARS " << "q4 " << "float " << "1" << std::endl;
                out << "LOOKUP_TABLE default" << std::endl;
                assert(!cell().empty());
                for (const MeshCell& mc: cell())
                {
                    out << mc.cons_sp1_(4) << std::endl;
                }

                out << "SCALARS " << "R0 " << "float " << "1" << std::endl;
                out << "LOOKUP_TABLE default" << std::endl;
                for (const MeshCell& mc: cell())
                {
                    out << mc.R(0) << std::endl;
                    /*if (std::abs(mc.R(0)) >= 1e3 || std::isnan(mc.R(0)))
                    {
                        std::cout << "fn: " << file_name << std::endl;
                        std::cout << "tagggggggg: " << mc.tag()() << std::endl;
                        std::cout << "value: " << mc.R(0) << std::endl;
                        std::cout << "type: " << static_cast<int>(mc.oga_cell_type()) << std::endl;
                    }*/
                }

                out.close();
            }

            //void Mesh::set_cell_type(StencilWalkResult result, const Tag& my_cell_tag, const MeshCell* other_cell)
            void Mesh::set_cell_type(const Tag& my_cell_tag, const MeshCell* other_cell, const Mesh& other_mesh)
            {
                MeshCell& my_cell = cell_p(my_cell_tag);
                auto my_type = my_cell.oga_cell_type();

                //if (result == StencilWalkResult::inside_cell)
                {
                    assert(other_cell != nullptr);
                    //if (my_cell.near_boundary())
                    if (my_cell.near_interog())
                    {
                        //if (tag_() == 1 && my_cell.tag()() == 45386) {
                            //assert(false);
                        //}
                        //my_cell.set_oga_cell_type(OGA_cell_type_t::mandat_receptor, other_cell->parent_mesh(), other_cell->tag(), other_cell);
                        my_cell.set_oga_cell_type(OGA_cell_type_t::mandat_receptor);
                        my_cell.add_cand_donor(other_cell->parent_mesh(), other_cell->tag(), other_cell);
                        assert(!my_cell.cand_donor().empty());
                        return;
                    }

                    if (my_type == OGA_cell_type_t::receptor)
                    {
                        //my_cell.set_oga_cell_type(OGA_cell_type_t::receptor, other_cell->parent_mesh(), other_cell->tag(), other_cell);
                        my_cell.set_oga_cell_type(OGA_cell_type_t::receptor);
                        my_cell.add_cand_donor(other_cell->parent_mesh(), other_cell->tag(), other_cell);
                        return;
                    }

                    if (my_type == OGA_cell_type_t::mandat_receptor)
                    {
                        //my_cell.set_oga_cell_type(OGA_cell_type_t::mandat_receptor, other_cell->parent_mesh(), other_cell->tag(), other_cell);
                        my_cell.set_oga_cell_type(OGA_cell_type_t::mandat_receptor);
                        my_cell.add_cand_donor(other_cell->parent_mesh(), other_cell->tag(), other_cell);
                        return;
                    }

                    if (priority_ > other_mesh.priority())
                    {
                        my_cell.set_oga_cell_type(OGA_cell_type_t::field);
                        if (other_mesh.tag()() == 0)
                        {
                            my_cell.set_receptor(other_mesh.tag(), other_cell->tag(), other_cell);
                        }
                    }
                    else if (priority_ < other_mesh.priority())
                    {
                        my_cell.set_oga_cell_type(OGA_cell_type_t::receptor);
                        my_cell.add_cand_donor(other_cell->parent_mesh(), other_cell->tag(), other_cell);
                    }
                    else
                    {
                        if (std::abs(my_cell.poly().volume()) > std::abs(other_cell->poly().volume()))
                        {
                            my_cell.set_oga_cell_type(OGA_cell_type_t::receptor);
                            my_cell.add_cand_donor(other_cell->parent_mesh(), other_cell->tag(), other_cell);
                        }
                        else
                        {
                            my_cell.set_oga_cell_type(OGA_cell_type_t::field);
                            if (other_mesh.tag()() == 0)
                            {
                                assert(false);
                            }
                        }
                    }


                    if (my_cell.oga_cell_type() == OGA_cell_type_t::receptor)
                    {
                        assert(!my_cell.cand_donor().empty());
                    }
                }

                if (my_type == OGA_cell_type_t::receptor)
                {
                    assert(!my_cell.cand_donor().empty());
                }

                assert(my_cell.oga_cell_type() != OGA_cell_type_t::undefined);
            }

void Mesh::set_priority(int pri)
{
    assert(pri != -1);
    priority_ = pri;
}
        }
