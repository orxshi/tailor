#include "mesh.h"
#include "meshpoint.h"
#include "meshcell.h"

namespace Tailor
{
    FaceTag gen_face_tag(const MeshFace& mf)
    {
        std::vector<int> pts;
        for (const auto& p: mf.mesh_point())
        {
            pts.push_back(p());
        }
        auto ft = FaceTag(pts);
        return ft;
    }

    //void Mesh::update_connectivity(int rank)
    //{
    //    for (MeshCell& mc: cell_)
    //    {
    //        for (MeshFace& mf: mc.face_)
    //        {
    //            if (mf.is_boundary())
    //            {
    //                assert(mf.parent_cell().size() == 2);
    //                //assert(query(mf.left_cell()) != nullptr);
    //                //assert(query_bou(mf.right_cell(), mf.btype()) != nullptr);
    //            }
    //            else
    //            {
    //                for (auto& pc: mf.parent_cell())
    //                {
    //                    if (query(pc) == nullptr)
    //                    {
    //                        mf.remove_parent_cell(pc);
    //                        mf.set_btype(BouType::partition);
    //                    }
    //                }

    //                for (const auto& nei: mc.pnei())
    //                {
    //                    if (query(nei) == nullptr)
    //                    {
    //                        mc.remove_neighbor(nei);
    //                    }
    //                }
    //            }
    //        }
    //    }

    //    for (MeshCell& mc: cell_)
    //    {
    //        int nnei = 0;
    //        for (const MeshFace& mf: mc.face())
    //        {
    //            if (!mf.is_boundary())
    //            {
    //                ++nnei;
    //            }
    //        }
    //        if (mc.pnei().size() != nnei)
    //        {
    //            mc.set_btype(BouType::partition);
    //        }
    //    }

    //    for (const MeshCell& mc: cell_)
    //    {
    //        for (const MeshFace& mf: mc.face())
    //        {
    //            assert(mf.tag().isvalid());
    //        }
    //    }

    //    for (const MeshCell& mc: cell_)
    //    {
    //        for (const MeshFace& mf: mc.face())
    //        {
    //            assert(mf.parent_cell().size() <= 2);
    //            auto fiter = std::find_if(face_.begin(), face_.end(), [&](const auto& f){return f.tag() == mf.tag();});
    //            if (fiter == face_.end())
    //            {
    //                face_.push_back(mf);
    //            }
    //        }
    //    }
    //}

    /*void Mesh::updateface(MeshFace& mf, MeshFace& cf)
    {
        // this should be previously partition face.
        if (mf.btype() != BouType::undefined)
        {
            std::cout << "typeeeeeeeeeeeeeeee: " << static_cast<int>(mf.btype()) << std::endl;
        }
        assert(mf.btype() == BouType::undefined);
        mf.set_btype(BouType::interior);
        cf.set_btype(BouType::interior);
        //auto it = std::find_if(face().begin(), face().end(), [&](const MeshFace& f){return f.tag() == mf.tag();});
        cf.set_tag(mf.tag());
    }*/

    void Mesh::update_interior_face(MeshFace& mf, MeshFace& cf, MeshCell& mc, MeshCell& cc, MeshFace& gmf)
    {
        if (mf.parent_cell().size() == 2)
        {
            assert(false);
            assert(mf.btype() == BouType::interior);
            assert(cf.btype() == BouType::interior);
            assert(cf.parent_cell().size() == 2);
            return;
        }

        Tag mctag = mc.tag();
        Tag cctag = cc.tag();

        assert(mctag.isvalid());
        assert(cctag.isvalid());

        mf.set_btype(BouType::interior);
        cf.set_btype(BouType::interior);
        gmf.set_btype(BouType::interior);

        //mf.set_faceaddr(&gmf);
        //cf.set_faceaddr(&gmf);

        auto ft = gen_face_tag(mf);

        mf.set_tag(ft);
        cf.set_tag(ft);
        gmf.set_tag(ft);

        //auto it = std::find_if(face_.begin(), face_.end(), [&](const MeshFace& f){return f.tag() == ft;});
        //assert(it != face_.end());

        auto ittt = std::find_if(mf.parent_cell().begin(), mf.parent_cell().end(), [&](const Tag& pcc){return pcc == cctag;});
        if (ittt == mf.parent_cell().end())
        {
            mf.add_parent_cell(cctag);
        }

        ittt = std::find_if(cf.parent_cell().begin(), cf.parent_cell().end(), [&](const Tag& pcc){return pcc == mctag;});
        if (ittt == cf.parent_cell().end())
        {
            cf.add_parent_cell(mctag);
        }

        ittt = std::find_if(gmf.parent_cell().begin(), gmf.parent_cell().end(), [&](const Tag& pcc){return pcc == mctag;});
        if (ittt == gmf.parent_cell().end())
        {
            gmf.add_parent_cell(mctag);
        }

        ittt = std::find_if(gmf.parent_cell().begin(), gmf.parent_cell().end(), [&](const Tag& pcc){return pcc == cctag;});
        if (ittt == gmf.parent_cell().end())
        {
            gmf.add_parent_cell(cctag);
        }

        if (std::count(mc.pnei().begin(), mc.pnei().end(), cctag) == 0)
        {
            mc.add_pnei(cctag);
        }
        if (std::count(cc.pnei().begin(), cc.pnei().end(), mctag) == 0)
        {
            cc.add_pnei(mctag);
        }

        // for connect_cells this step is redundant.
        // for connect_after_exchange this step is necessary.
        for (const Tag& p: mf.mesh_point())
        {
            auto pp = point_p(p);
            assert(pp != nullptr);
            pp->add_parent_cell(cctag);
        }

        for (const Tag& p: cf.mesh_point())
        {
            auto pp = point_p(p);
            assert(pp != nullptr);
            pp->add_parent_cell(mctag);
        }

        assert(std::find_if(mc.face().begin(), mc.face().end(), [&](const MeshFace& ff){return ff.tag() == mf.tag();}) != mc.face().end());
        assert(std::find_if(cc.face().begin(), cc.face().end(), [&](const MeshFace& ff){return ff.tag() == mf.tag();}) != cc.face().end());
    }

    void Mesh::add_interior_face(MeshFace& mf, MeshFace& cf, MeshCell& mc, MeshCell& cc)
    {
        assert(mf.parent_cell().size() == 1);
        if (cf.parent_cell().size() != 1)
        {
            std::cout << "mc: " << mc.tag()() << std::endl;
            std::cout << "cc: " << cc.tag()() << std::endl;
            std::cout << "mf: " << mf.tag() << std::endl;
            std::cout << "cf: " << cf.tag() << std::endl;
            std::cout << "mf is bou: " << mf.is_boundary() << std::endl;
            std::cout << "cf is bou: " << cf.is_boundary() << std::endl;
            for (const auto& t: cf.parent_cell())
            {
                std::cout << "par: " << t() << std::endl;
            }
        }
        assert(cf.parent_cell().size() == 1);

        if (mf.parent_cell().size() == 2)
        {
            assert(false);
            assert(mf.btype() == BouType::interior);
            assert(cf.btype() == BouType::interior);
            assert(cf.parent_cell().size() == 2);
            return;
        }

        Tag mctag = mc.tag();
        Tag cctag = cc.tag();

        assert(mctag.isvalid());
        assert(cctag.isvalid());

        //{
        //auto ittt = std::find_if(cf.parent_cell().begin(), cf.parent_cell().end(), [&](const auto& pcc){return pcc == mctag;});
        //if (ittt != cf.parent_cell().end())
        //{
        //    std::cout << "mesh: " << tag_() << std::endl;
        //    std::cout << "mc: " << mc.tag()() << std::endl;
        //    std::cout << "cc: " << cc.tag()() << std::endl;
        //    std::cout << "mf: " << mf.tag() << std::endl;
        //    std::cout << "cf: " << cf.tag() << std::endl;
        //    std::cout << "cf btype: " << static_cast<int>(cf.btype()) << std::endl;
        //    std::cout << "mf btype: " << static_cast<int>(mf.btype()) << std::endl;
        //    std::cout << "cf pc size: " << cf.parent_cell().size() << std::endl;
        //    for (const auto& pc: cf.parent_cell())
        //    {
        //        std::cout << "cf pc: " << pc() << std::endl;
        //    }
        //}
        //assert(ittt == cf.parent_cell().end());
        //}

        mf.set_btype(BouType::interior);
        cf.set_btype(BouType::interior);

        auto ft = gen_face_tag(mf);

        mf.set_tag(ft);
        cf.set_tag(ft);

        //auto it = std::find_if(face_.begin(), face_.end(), [&](const MeshFace& f){return f.tag() == ft;});
        //assert(it == face_.end());

        auto ittt = std::find_if(mf.parent_cell().begin(), mf.parent_cell().end(), [&](const Tag& pcc){return pcc == cctag;});
        assert(ittt == mf.parent_cell().end());
        if (ittt == mf.parent_cell().end())
        {
            mf.add_parent_cell(cctag);
        }

        ittt = std::find_if(cf.parent_cell().begin(), cf.parent_cell().end(), [&](const Tag& pcc){return pcc == mctag;});
        if (ittt == cf.parent_cell().end())
        {
            cf.add_parent_cell(mctag);
            if (cf.parent_cell().size() == 2)
            {
                assert(cf.parent_cell(0) != cf.parent_cell(1));
            }
        }

        //if (mc.parent_mesh()() == 0)
        //{
            //if (mctag() == 18193 || cctag() == 18193)
            //{
                //if (mf.parent_cell().size() == 2)
                //{
                    //std::cout << mc.tag()() << std::endl;
                    //std::cout << cc.tag()() << std::endl;
                //}
                //assert(mf.parent_cell().size() != 2);
                //std::cout << mf.parent_cell().size() << std::endl;
            //}
        //}

        //if (mf.parent_cell().size() == 2)
        //{
            //face_.push_back(mf);
            //mf.set_faceaddr(&face_.back());
            //cf.set_faceaddr(&face_.back());
            //assert(mf.faceaddr() == &face_.back());
            //assert(mf.faceaddr() == &face(mf.tag()));
        //}

        if (std::count(mc.pnei().begin(), mc.pnei().end(), cctag) == 0)
        {
            mc.add_pnei(cctag);
        }
        if (std::count(cc.pnei().begin(), cc.pnei().end(), mctag) == 0)
        {
            cc.add_pnei(mctag);
        }

        // for connect_cells this step is redundant.
        // for connect_after_exchange this step is necessary.
        for (const Tag& p: mf.mesh_point())
        {
            auto pp = point_p(p);
            assert(pp != nullptr);
            pp->add_parent_cell(cctag);
        }

        for (const Tag& p: cf.mesh_point())
        {
            auto pp = point_p(p);
            assert(pp != nullptr);
            pp->add_parent_cell(mctag);
        }

        assert(std::find_if(mc.face().begin(), mc.face().end(), [&](const MeshFace& ff){return ff.tag() == mf.tag();}) != mc.face().end());
        assert(std::find_if(cc.face().begin(), cc.face().end(), [&](const MeshFace& ff){return ff.tag() == mf.tag();}) != cc.face().end());

        assert(mf.parent_cell().size() <= 2);
        assert(cf.parent_cell().size() <= 2);
    }

    void Mesh::add_partition_face(MeshFace& mf)
    {
        mf.set_btype(BouType::partition);
        // the rest is useless for connect_after_exchange().

        auto ft = gen_face_tag(mf);

        mf.set_tag(ft);
        //face_.push_back(mf);
        //mf.set_faceaddr(&face_.back());
        assert(mf.parent_cell().size() == 1);
    }

    void Mesh::change_bou(BouType btype, MeshCell& mc, const FaceTag& ftag)
    {
        const Array<Tag, 6>* container = nullptr;
        if (btype == BouType::wall) {
            container = &mc.wall_boundary();
        }
        else if (btype == BouType::dirichlet) {
            container = &mc.dirichlet_boundary();
        }
        else if (btype == BouType::farfield) {
            container = &mc.farfield_boundary();
        }
        else if (btype == BouType::empty) {
            container = &mc.empty_boundary();
        }

        for (const Tag& tt: *container)
        {
            auto bb = boundary(btype, tt);
            assert(bb != nullptr);

            bool facefound = false;
            for (const MeshFace& mf: mc.face())
            {
                assert(bb->face().size() == 1);
                if (mf.tag() == bb->face()[0].tag())
                {
                    assert(bb->face_p()[0].is_boundary());
                    if (bb->face_p()[0].tag() != ftag)
                    {
                        std::cout << "bou tag: " << bb->face_p()[0].tag() << std::endl;
                        std::cout << "ftag: " << ftag << std::endl;
                        std::cout << "mc: " << mc.tag()() << std::endl;
                        std::cout << "mesh: " << mc.parent_mesh()() << std::endl;
                        std::cout << "btype: " << static_cast<int>(btype) << std::endl;
                    }
                    assert(bb->face_p()[0].tag() == ftag);
                    //bb->face_p()[0].set_tag(ftag);
                    facefound = true;
                    break;
                }
            }
            if (!facefound)
            {
                std::cout << "mc: " << mc.tag()() << std::endl;
                std::cout << "btype: " << static_cast<int>(btype) << std::endl;
                std::cout << "container size: " << container->size() << std::endl;
                //std::cout << "bb: " << bb.face()[0].tag()() << std::endl;

                //for (const MeshFace& mf: mc.face())
                //{
                    //std::cout << "mf: " << mf.tag()() << " " << static_cast<int>(mf.btype()) << std::endl;
                //}

                for (const Tag& tt: *container)
                {
                    std::cout << "con: " << (tt)() << std::endl;
                }
            }
            assert(facefound);
        }
    }

    void Mesh::addface_bou(MeshFace& mf)
    {
        assert(mf.btype() != BouType::undefined);
        assert(mf.is_boundary());

        auto ft = gen_face_tag(mf);

        //mf.set_tag(FaceTag(mc.tag()(), rc(), mf.btype()));
        mf.set_tag(ft);

        //change_bou(BouType::wall, mc, mf.tag());
        //change_bou(BouType::dirichlet, mc, mf.tag());
        //change_bou(BouType::farfield, mc, mf.tag());
        //change_bou(BouType::empty, mc, mf.tag());

        //face_.push_back(mf);
        //mf.set_faceaddr(&face_.back());
    }

    void Mesh::opposing_cell_pc2(const MeshCell& mc, MeshFace& mf, MeshCell*& op_cell, MeshFace*& shared_face)
    {
        // returns opposing cell (if exists) which shares a face with mf.
        // returns shared face of opposing cell.
        // used when mf.parent_cell.size() = 2 so that potential opposing cell is already known.

        op_cell = nullptr;
        shared_face = nullptr;

        for (auto& t: mf.parent_cell())
        {
            if (t == mc.tag()) {
                continue;
            }

            assert(t.isvalid());
            //assert(query(t) != nullptr);

            op_cell = &cell_p(t);
            assert(op_cell != nullptr);

            shared_face = common_face(*op_cell, mc.tag());
            /*if (mc.tag()() == 998 && t() == 1011)
            {
                if (mf.tag()().size() == 3)
                {
                    if (mf.tag()()[0] == 4 && mf.tag()()[1] == 998 && mf.tag()()[2] == 1011)
                    {
                        std::cout << "null: " << (shared_face == nullptr) << std::endl;
                        assert(false);
                    }
                }
            }*/
            if (shared_face != nullptr)
            {
                assert(op_cell != nullptr);
                return;
            }
        }

        if (shared_face == nullptr) {
            op_cell = nullptr;
        }
    }

    void Mesh::opposing_cell_pc1(const MeshCell& mc, const MeshFace& mf, MeshCell*& op_cell, MeshFace*& shared_face)
    {
        // returns opposing cell (if exists) which shares a face with mf.
        // returns shared face of opposing cell.
        // used when mf.parent_cell.size() = 1 so that potential opposing cell cannot be found quickly.

        op_cell = nullptr;
        shared_face = nullptr;

        assert(!mf.mesh_point().empty());
        auto pp = point_p(mf.mesh_point(0));
        assert(pp != nullptr);
        //MeshPoint& p0 = point_p(mf.mesh_point(0));
        MeshPoint& p0 = *pp;
        assert(!p0.parent_cell().empty());

        for (auto& nei: p0.parent_cell())
        {
            //assert(query(nei) != nullptr);
            MeshCell& nc = cell_p(nei);
            //assert(nei.isvalid());
            //if (nei == mc.tag())
            if (nc.tag() == mc.tag())
            {
                continue;
            }

            if (is_potential_parent_cell(mf, nc.tag()))
            {
                shared_face = common_face(nc, mc.tag());
                op_cell = &nc;
                if (shared_face != nullptr) {
                    return;
                }
            }
        }

        if (shared_face == nullptr) {
            op_cell = nullptr;
        }
    }

    void Mesh::reset_mesh_connectivity()
    {
        for (MeshCell& mc: cell_)
        {
            mc.remove_all_neighbors();
            for (MeshFace& mf: mc.face_p())
            {
                if (!mf.is_boundary()) {
                    mf.remove_parent_cells();
                }
            }
        }

        remove_parent_cells_of_all_points();
        update_points_from_cell_vertices(-1);
    }

    void Mesh::connect_cells(std::function<bool(const vec3<double>&)> is_resi, int rank, Profiler* profiler, std::string name)
    {
        /*for (MeshCell& mc: cell_)
        {
            mc.remove_all_neighbors();
            for (MeshFace& mf: mc.face_p())
            {
                if (!mf.is_boundary()) {
                    mf.remove_parent_cells();
                }
            }
        }*/

        //for (MeshCell& mc: cell_)
        //{
        //    for (MeshFace& mf: mc.face_p())
        //    {
        //        assert(mf.parent_cell().size() != 1);
        //        if (mf.is_boundary())
        //        {
        //            if (mf.parent_cell().size() != 2)
        //            {
        //                std::cout << "rank: " << tag_() << std::endl;
        //                std::cout << "mesh tag: " << tag_() << std::endl;
        //                std::cout << "mc tag: " << mc.tag()() << std::endl;
        //                std::cout << "mf tag: " << mf.tag() << std::endl;
        //                std::cout << "mf btype: " << static_cast<int>(mf.btype()) << std::endl;
        //            }
        //            assert(mf.parent_cell().size() == 2);
        //        }
        //    }
        //}

        //if (profiler != nullptr) {profiler->start(name + "-connect-cells-pre");}
        //assert(face_.empty());

        for (MeshCell& mc: cell_)
        {
            for (MeshFace& mf: mc.face_p())
            {
                if (mf.is_boundary())
                {
                    continue;
                }

                assert(mf.parent_cell().empty());
                mf.add_parent_cell(mc.tag());
                assert(mf.parent_cell().size() == 1);
                //if (mf.tag() == FaceTag(std::vector{3279,3286,3340}))
                //{
                //    std::cout << "pc size: " << mf.parent_cell().size() << std::endl;
                //    for (const auto& t: mf.parent_cell())
                //    {
                //        std::cout << "pc: " << t() << std::endl;
                //    }
                //}
            }
        }

        for (MeshCell& mc: cell_)
        {
            assert(mc.tag().isvalid());

            for (MeshFace& mf: mc.face_p())
            {
                if (mf.is_boundary())
                {
                    assert(mf.parent_cell().size() == 2);
                    addface_bou(mf);
                    continue;
                }

                //assert(std::count(mf.parent_cell().begin(), mf.parent_cell().end(), mc.tag()) == 0);
                //assert(std::count_if(mf.parent_cell().begin(), mf.parent_cell().end(), [&](const auto& a){return a == mc.tag();}) == 0);
                //mf.add_parent_cell(mc.tag());

                if (mf.parent_cell().size() == 2)
                {
                    assert(mf.btype() != BouType::undefined);
                    continue;
                }

                assert(!mf.mesh_point().empty());
                //MeshPoint& p0 = point_p(mf.mesh_point(0));
                auto pp = point_p(mf.mesh_point(0));
                assert(pp != nullptr);
                MeshPoint& p0 = *pp;
                assert(!p0.parent_cell().empty());

                for (auto& nei: p0.parent_cell())
                {
                    //assert(query(nei) != nullptr);
                    MeshCell& nc = cell_p(nei);
                    assert(nei.isvalid());
                    if (nc.tag() == mc.tag())
                    {
                        continue;
                    }
                    
                    if (is_potential_parent_cell(mf, nc.tag()))
                    {
                        //if (std::count(mc.pnei().begin(), mc.pnei().end(), nei) == 0)
                        //{
                            //mc.add_pnei(nei);
                        //}
                        //if (std::count(nc.pnei().begin(), nc.pnei().end(), mc.tag()) == 0)
                        //{
                            //nc.add_pnei(mc.tag());
                        //}
                        //mf.add_parent_cell(nei);
                        //assert(mc.tag().isvalid());
                        //assert(query(nei));
                        //assert(mf.parent_cell().size() == 2);

                        MeshFace* cf = common_face(nc, mc.tag());
                        assert(cf != nullptr);
                        //cf->add_parent_cell(mc.tag());
                        //assert(query(mc.tag()));

                        add_interior_face(mf, *cf, mc, nc);

                        //if (mf.parent_cell().size() == 2)
                        //{
                            //if (mf.btype() != BouType::undefined)
                            //{
                                //std::cout << "type: " << static_cast<int>(mf.btype()) << std::endl;
                            //}
                            //assert(mf.btype() == BouType::undefined);
                            //add_interior_face(mf, *cf, mc, nc);
                        //}

                        break;
                    }
                }
            }
        }

        for (MeshCell& mc: cell_)
        {
            int nnei = 0;
            for (const MeshFace& mf: mc.face())
            {
                if (!mf.is_boundary())
                {
                    ++nnei;
                }
            }
            if (mc.pnei().size() != nnei)
            {
                //if (tag_() == 0 && mc.tag()() == 194)
                //{
                    //std::cout << "pnei size: " << mc.pnei().size() << std::endl;
                    //std::cout << "nnei: " << nnei << std::endl;
                    //assert(false);
                //}
                mc.set_btype(BouType::partition);
            }
            for (MeshFace& mf: mc.face_p())
            {
                if (mf.parent_cell().size() == 1)
                {
                    assert(!mf.is_boundary());
                    //assert(!mf.tag().isvalid());
                    add_partition_face(mf);
                    assert(mf.parent_cell().size() <= 1);
                }
            }
        }

        for (const MeshCell& mc: cell_)
        {
            for (const MeshFace& mf: mc.face())
            {
                assert(!mf.parent_cell().empty());

                if (mf.is_boundary())
                {
                    if (mf.parent_cell().size() != 2)
                    {
                        std::cout << mf.parent_cell().size() << std::endl;
                    }
                    assert(mf.parent_cell().size() == 2);
                }
                else
                {
                    if (mf.parent_cell().size() == 2)
                    {
                        assert(mf.btype() == BouType::interior);
                    }
                }
            }
        }

        //for (const MeshCell& mc: cell_)
        //{
        //    for (const MeshFace& mf: mc.face())
        //    {
        //        //if (mf.btype() == BouType::partition) {
        //            //continue;
        //        //}
        //        auto it = std::find_if(face_.begin(), face_.end(), [&](const MeshFace& f){return f.tag() == mf.tag();});
        //        if (it == face_.end())
        //        {
        //            std::cout << "mc: " << mc.tag()() << std::endl;
        //            //std::cout << "mf: " << mf.tag()() << std::endl;
        //            std::cout << "parent size: " << mf.const_parent_cell().size() << std::endl;
        //            std::cout << "is boundary: " << mf.is_boundary() << std::endl;
        //            std::cout << "btype: " << static_cast<int>(mc.btype()) << std::endl;
        //            std::cout << "mc.pnei().size(): " << mc.const_pnei().size() << std::endl;
        //            std::cout << "mc.poly().faces().size(): " << mc.poly().faces().size() << std::endl;
        //        }
        //        assert(it != face_.end());
        //    }
        //}

        //for (const auto& mf: face_)
        //{
        //    if (mf.btype() == BouType::undefined)
        //    {
        //        std::cout << "pc size: " << mf.parent_cell().size() << std::endl;
        //    }
        //    assert(mf.btype() != BouType::undefined);
        //}

        //for (const MeshFace& mf: face_)
        //{
        //    assert(mf.tag().isvalid());
        //}

        for (MeshCell& mc: cell_)
        {
            if (!is_resi(mc.poly().centroid()))
            {
                mc.set_oga_cell_type(OGA_cell_type_t::non_resident);
            }
        }
        //if (profiler != nullptr) {profiler->stop(name + "-connect-cells-pre");}

        for (const auto& mc: cell_)
        {
            for (const auto& mf: mc.face_)
            {
                if (mf.is_boundary()) {
                    continue;
                }
                if (mf.parent_cell(0) != mc.tag())
                {
                    std::cout << "mc: " << mc.tag()() << std::endl;
                    std::cout << "pc 0: " << mf.parent_cell(0)() << std::endl;
                    std::cout << "pc 1: " << mf.parent_cell(1)() << std::endl;
                    std::cout << "mf btype: " << static_cast<int>(mf.btype()) << std::endl;
                }
                assert(mf.parent_cell(0) == mc.tag());
                assert(mf.right_cell() == mc.tag());
            }
        }
    }

    void Mesh::connect_after_exchange(std::function<bool(const vec3<double>&)> is_resi, int rank, Profiler* profiler, std::string procorig)
    {
        std::string proc;
        proc = procorig;
        if (profiler != nullptr) {profiler->start(proc.append("-connect-after-exc-1"));}
        for (MeshCell& mc: cell_)
        {
            if (!is_resi(mc.poly().centroid()))
            {
                mc.set_oga_cell_type(OGA_cell_type_t::non_resident);
            }
        }
        proc = procorig;
        if (profiler != nullptr) {profiler->stop(proc.append("-connect-after-exc-1"));}

        proc = procorig;
        if (profiler != nullptr) {profiler->start(proc.append("-connect-after-exc-2"));}
        for (MeshCell& mc: cell_)
        {
            if (mc.pnei().size() == mc.poly().faces().size()) {
                continue;
            }

            for (MeshFace& mf: mc.face_p())
            {
                assert(mf.tag().isvalid());

                if (mf.is_boundary())
                {
                    if (mf.parent_cell().size() != 2)
                    {
                        std::cout << "pc size: " << mf.parent_cell().size() << std::endl;
                    }
                    assert(mf.parent_cell().size() == 2);
                    //auto it = std::find_if(face_.begin(), face_.end(), [&](const MeshFace& f){return f.tag() == mf.tag();});
                    //if (it == face_.end())
                    //{
                    //    addface_bou(mf);
                    //    assert(mf.parent_cell().size() == face_.back().parent_cell().size());
                    //}
                    continue;
                }

                if (mf.btype() == BouType::partition)
                {
                    MeshCell* op_cell = nullptr;
                    MeshFace* shared_face = nullptr;
                    opposing_cell_pc1(mc, mf, op_cell, shared_face);

                    if (op_cell != nullptr)
                    {
                        assert(shared_face != nullptr);
                        //auto it = std::find_if(face_.begin(), face_.end(), [&](const MeshFace& f){return f.tag() == mf.tag();});
                        //if (it == face_.end())
                        {
                            add_interior_face(mf, *shared_face, mc, *op_cell);
                            //assert(mf.parent_cell().size() == face_.back().parent_cell().size());
                            //assert(mf.faceaddr() == &face(mf.tag()));
                        }
                        //else
                        //{
                            //update_interior_face(mf, *shared_face, mc, *op_cell, *it);
                            //assert(mf.parent_cell().size() == it->parent_cell().size());
                            //assert(mf.faceaddr() == &face(mf.tag()));
                        //}
                    }
                    else
                    {
                        // will remain as partition face.
                        //auto it = std::find_if(face_.begin(), face_.end(), [&](const MeshFace& f){return f.tag() == mf.tag();});
                        //if (it == face_.end())
                        //{
                            add_partition_face(mf);
                            //assert(mf.parent_cell().size() == face_.back().parent_cell().size());
                            //assert(mf.faceaddr() == &face(mf.tag()));
                        //}
                    }

                    //auto fiter = std::find_if(face().begin(), face().end(), [&](const auto& ff){return ff.tag() == mf.tag();});
                    //assert(fiter != face().end());
                }
                else if (mf.btype() == BouType::interior)
                {
                    if (mf.parent_cell().size() == 1)
                    {
                        MeshCell* op_cell = nullptr;
                        MeshFace* shared_face = nullptr;
                        opposing_cell_pc1(mc, mf, op_cell, shared_face);

                        if (op_cell == nullptr)
                        {
                            //auto it = std::find_if(face_.begin(), face_.end(), [&](const MeshFace& f){return f.tag() == mf.tag();});
                            //if (it == face_.end())
                            //{
                                add_partition_face(mf);
                                //assert(mf.parent_cell().size() == face_.back().parent_cell().size());
                                //assert(mf.faceaddr() == &face(mf.tag()));
                            //}
                            //else
                            //{
                                //mf.set_btype(BouType::partition);
                                //it->set_btype(BouType::partition);
                                //if (mf.parent_cell().size() != it->parent_cell().size())
                                //{
                                    //std::cout << "mf parent size: " << mf.parent_cell().size() << std::endl;
                                    //std::cout << "gmf parent size: " << it->parent_cell().size() << std::endl;
                                //}
                                //assert(mf.tag() == it->tag());
                                //assert(mf.parent_cell().size() == it->parent_cell().size());
                                //assert(mf.faceaddr() == &face(mf.tag()));
                            //}
                        }
                        else
                        {
                            //auto it = std::find_if(face_.begin(), face_.end(), [&](const MeshFace& f){return f.tag() == mf.tag();});
                            //if (it != face_.end())
                            //{
                                //update_interior_face(mf, *shared_face, mc, *op_cell, *it);
                                //assert(mf.parent_cell().size() == it->parent_cell().size());
                                //assert(mf.faceaddr() == &face(mf.tag()));
                            //}
                            //else
                            //{
                                add_interior_face(mf, *shared_face, mc, *op_cell);
                                //assert(mf.parent_cell().size() == face_.back().parent_cell().size());
                                //assert(mf.faceaddr() == &face(mf.tag()));
                            //}
                        }

                        //auto fiter = std::find_if(face().begin(), face().end(), [&](const auto& ff){return ff.tag() == mf.tag();});
                        //assert(fiter != face().end());
                    }
                    else if (mf.parent_cell().size() == 2)
                    {
                        //auto it = std::find_if(face_.begin(), face_.end(), [&](const MeshFace& f){return f.tag() == mf.tag();});
                        //assert(it != face_.end());
                        continue;

                        MeshCell* op_cell = nullptr;
                        MeshFace* shared_face = nullptr;
                        opposing_cell_pc2(mc, mf, op_cell, shared_face);

                        assert(op_cell != nullptr);
                        assert(shared_face != nullptr);
                    }
                }
            }
        }
        proc = procorig;
        if (profiler != nullptr) {profiler->stop(proc.append("-connect-after-exc-2"));}

        //for (const MeshFace& mf: face())
        //{
            //if (mf.btype() == BouType::undefined)
            //{
            //    auto iter = std::find_if(cell_.begin(), cell_.end(), [&](const auto& c){return c.tag() == mf.parent_cell()[0];});
            //    assert(iter != cell_.end());
            //    //std::cout << "mf tag: " << mf.tag()() << std::endl;
            //    std::cout << "pc: " << mf.parent_cell()[0]() << std::endl;
            //    std::cout << "pc size: " << mf.parent_cell().size() << std::endl;
            //    std::cout << "is bou: " << mf.is_boundary() << std::endl;
            //    std::cout << "pc2 size: " << mf.parent_cell().size() << std::endl;
            //    std::cout << "btype " << static_cast<int>(mf.btype())<< std::endl;
            //    //for (const auto& f: cell(mf.parent_cell()[0]).face())
            //    //{
            //        //std::cout << "f btype " << f.tag()() << " " << static_cast<int>(f.btype())<< std::endl;
            //    //}
            //    //std::cout << "nei size: " << mc.pnei().size() << std::endl;
            //    //std::cout << "mc face size: " << mc.poly().faces().size() << std::endl;
            //}
            //assert(mf.btype() != BouType::undefined);
        //}

        //for (const auto& mf: face_)
        //{
        //    assert(mf.tag().isvalid());
        //}
        //for (const auto& mc: cell_)
        //{
        //    for (const auto& mf: mc.face())
        //    {
        //        assert(mf.tag().isvalid());
        //    }
        //}

        // delete hanging invalid faces.
        //for (auto mf = face_.begin(); mf != face_.end();)
        //{
            //if (!mf->tag().isvalid())
            //{
                //mf = face_.erase(mf);
            //}
            //else
            //{
                //++mf;
            //}
        //}

        proc = procorig;
        if (profiler != nullptr) {profiler->start(proc.append("-connect-after-exc-3"));}
        for (MeshCell& mc: cell_)
        {
            if (mc.btype() == BouType::partition)
            {
                int nnei = 0;
                for (const MeshFace& mf: mc.face())
                {
                    if (!mf.is_boundary())
                    {
                        ++nnei;
                    }
                }

                if (mc.pnei().size() == nnei)
                {
                    mc.set_btype(BouType::interior);
                }
            }
        }

        for (MeshCell& mc: cell_)
        {
            int nnei = 0;
            for (const MeshFace& mf: mc.face())
            {
                if (!mf.is_boundary())
                {
                    ++nnei;
                }
            }
            if (mc.pnei().size() != nnei)
            {
                mc.set_btype(BouType::partition);
            }
        }
        proc = procorig;
        if (profiler != nullptr) {profiler->stop(proc.append("-connect-after-exc-3"));}

        //for (const auto& mc: cell_)
        //{
        //    for (const auto& mf: mc.face())
        //    {
        //        auto fiter = std::find_if(face().begin(), face().end(), [&](const auto& ff){return ff.tag() == mf.tag();});
        //        if (fiter == face().end())
        //        {
        //            std::cout << "mf : " << mf.tag() << std::endl;
        //            std::cout << "mf btype : " << static_cast<int>(mf.btype()) << std::endl;
        //            std::cout << "mc : " << mc.tag()() << std::endl;
        //            std::cout << "mc pnei size : " << mc.const_pnei().size() << std::endl;
        //        }
        //        assert(fiter != face().end());
        //    }
        //}

        //for (const auto& mc: cell_)
        //{
        //    for (const Tag& bc: mc.farfield_boundary())
        //    {
        //        assert(query_bou(bc, BouType::farfield) != nullptr);
        //    }
        //}

        //for (const auto& mc: cell())
        //{
        //    for (const auto& mf: mc.face())
        //    {
        //        auto fi = std::find_if(face().begin(), face().end(), [&](const auto& ff){return ff.tag() == mf.tag();});
        //        assert(fi != face().end());
        //        if (mf.parent_cell().size() != face(mf.tag()).const_parent_cell().size())
        //        {
        //            std::cout << "mf: " << mf.tag() << std::endl;
        //            std::cout << "mf type: " << static_cast<int>(mf.btype()) << std::endl;
        //            std::cout << "mc type: " << static_cast<int>(mc.btype()) << std::endl;
        //            std::cout << "mc: " << mc.tag()() << std::endl;
        //            std::cout << "mf pc size: " << mf.const_parent_cell().size() << std::endl;
        //            std::cout << "gmf pc size: " << face(mf.tag()).const_parent_cell().size() << std::endl;
        //        }
        //        assert(mf.const_parent_cell().size() == face(mf.tag()).const_parent_cell().size());
        //    }
        //}

        for (const auto& mc: cell_)
        {
            for (const auto& mf: mc.face_)
            {
                if (mf.is_boundary()) {
                    continue;
                }
                if (mf.parent_cell(0) != mc.tag())
                {
                    std::cout << "mc: " << mc.tag()() << std::endl;
                    std::cout << "pc 0: " << mf.parent_cell(0)() << std::endl;
                    std::cout << "pc 1: " << mf.parent_cell(1)() << std::endl;
                    std::cout << "mf btype: " << static_cast<int>(mf.btype()) << std::endl;
                }
                assert(mf.parent_cell(0) == mc.tag());
            }
        }
    }

    /*void Mesh::connect_partition_cells(std::vector<Nei>& arrival_cell, int rank, std::function<bool(const vec3<double>&, int celltag)> is_resi)
    {
        for (const MeshFace& mf: face_)
        {
            if (mf.is_boundary()) {
                continue;
            }
            for (const Tag& t: mf.parent_cell())
            {
                bool found = false;
                assert(query(t) != nullptr);
                const auto& pc = cell(t);
                for (const MeshFace& cmf: pc.face())
                {
                    if (cmf.tag() == mf.tag())
                    {
                        found = true;
                        break;
                    }
                }
                if (!found)
                {
                    if (rank == 2)
                    {
                    std::cout << "pc: " << pc.tag()() << std::endl;
                    std::cout << "pc type: " << static_cast<int>(pc.oga_cell_type()) << std::endl;
                    std::cout << "mf btype: " << static_cast<int>(mf.btype()) << std::endl;
                    std::cout << "mc btype: " << static_cast<int>(pc.btype()) << std::endl;
                    for (const auto& nei: mf.parent_cell())
                    {
                        std::cout << "nei: " << nei() << std::endl;
                    }
                    print_as_vtk("two.vtk");
                    }
                }
                assert(found);
            }
        }

                for (const auto& mc: cell())
                {
                    for (const MeshFace& mf: mc.face())
                    {
                        auto it = std::find_if(face_.begin(), face_.end(), [&](const MeshFace& f){return f.tag() == mf.tag();});
                        if (it != face_.end()) {
                            int co = std::count(it->parent_cell().begin(), it->parent_cell().end(), mc.tag());
                            if (co == 0)
                            {
                                assert(mc.oga_cell_type() == OGA_cell_type_t::ghost);
                            }
                        }
                    }
                }

        for (MeshCell& mc: cell_)
        {
            if (!is_resi(mc.poly().centroid(), mc.tag()()))
            {
                mc.set_oga_cell_type(OGA_cell_type_t::non_resident);
            }
        }

        //if (mc.tag()() == 1047 && tag_() == 1)
        //{
            //assert(mc.oga_cell_type() == OGA_cell_type_t::non_resident);
        //}


        //for (Nei& nn: arrival_cell)
        //{
            //MeshCell& ac = nn.cell_;
            //ac.reset_face_tags();
        //}

        for (Nei& nn: arrival_cell)
        {
            MeshCell& ac = nn.cell_;
            ac.remove_all_neighbors();
            ac.remove_all_nonself_parents_of_faces();
        }
        
        std::vector<int> to_be_added;

        if (!arrival_cell.empty())
        {
            for (MeshCell& mc: cell_)
            {
                assert(mc.parent_mesh() == tag_);
                if (mc.parent_mesh() != tag_) {
                    continue;
                }

                if (mc.btype() == BouType::partition && mc.oga_cell_type() != OGA_cell_type_t::non_resident)
                {
                    for (MeshFace& mf: mc.face_)
                    {
                        if (mf.parent_cell().size() == 1)
                        {
                            bool acfound = false;
                            int actag;
                            int actype;
                            //for (const MeshCell& ac: arrival_cell)
                            for (Nei& nn: arrival_cell)
                            {
                                MeshCell& ac = nn.cell_;
                                if (ac.parent_mesh() != tag_) {
                                    continue;
                                }
                                actag = ac.tag()();
                                actype = static_cast<int>(ac.oga_cell_type());


                                assert(ac.tag().isvalid());
                                if (mc.tag() == ac.tag()) {
                                    continue;
                                }

                                for (MeshFace& acmf: ac.face_p())
                                {
                                    //if (!are_common_faces(mf, acmf)) {
                                    //continue;
                                    //}

                                    if (mf.tag() != acmf.tag()) {
                                        continue;
                                    }

                                    assert(acmf.tag().isvalid());
                                    assert(mf.parent_cell().size() <= 2);
                                    //assert(acmf.parent_cell().size() <= 2);
                                    assert(acmf.parent_cell().size() == 1);

                                    //if (mf.tag() != acmf.tag())
                                    //  {
                                    //  std::cout << "mc: " << mc.tag()() << std::endl;
                                    //  std::cout << "ac: " << ac.tag()() << std::endl;
                                    //  std::cout << "mftag: " << mf.tag()() << std::endl;
                                    //  std::cout << "acmftag: " << acmf.tag()() << std::endl;
                                    //  std::cout << "mf btype: " << static_cast<int>(mf.btype()) << std::endl;
                                    //  std::cout << "acmf btype: " << static_cast<int>(acmf.btype()) << std::endl;
                                    //  }

                                    assert(mf.tag().isvalid());
                                    assert(acmf.tag().isvalid());
                                    //assert(mf.tag() == acmf.tag());
                                    //if (rank == 7)
                                    //{
                                    //    if (acmf.btype() != BouType::interior)
                                    //    {
                                    //        std::cout << "gmf btype: " << static_cast<int>(face(mf.tag()).btype()) << std::endl;
                                    //        std::cout << "mf btype: " << static_cast<int>(mf.btype()) << std::endl;
                                    //        std::cout << "acmf btype: " << static_cast<int>(acmf.btype()) << std::endl;
                                    //        std::cout << "mesh: " << tag_() << std::endl;
                                    //        std::cout << "mc: " << mc.tag()() << std::endl;
                                    //        std::cout << "ac: " << ac.tag()() << std::endl;
                                    //        std::cout << "acmf: " << acmf.tag()() << std::endl;
                                    //        std::cout << "mf: " << mf.tag()() << std::endl;
                                    //        for (const auto& t: mf.face().vertex())
                                    //        {
                                    //            std::cout << "mf r: " << t.r(0) << " " << t.r(1) << " " << t.r(2) << std::endl;
                                    //        }
                                    //        for (const auto& t: acmf.face().vertex())
                                    //        {
                                    //            std::cout << "acmf r: " << t.r(0) << " " << t.r(1) << " " << t.r(2) << std::endl;
                                    //        }
                                    //        print_as_vtk_geometry("ranktwo.vtk");
                                    //    }
                                    //}
                                    //assert(acmf.btype() == BouType::interior);

                                    to_be_added.push_back(ac.tag()());

                                    assert(mf.tag() == acmf.tag());
                                    //mf.set_tag(acmf.tag()); // why?

                                    auto fiter = std::find_if(face_.begin(), face_.end(), [&](const auto& f){return f.tag() == mf.tag();});
                                    assert(fiter != face_.end());
                                    //if (fiter == face_.end())
                                    //{
                                    //face_.push_back(mf);
                                    //}

                                    MeshFace& gmf = face_p(mf.tag());

                                    //auto ex = std::find_if(mf.parent_cell().begin(), mf.parent_cell().end(), [&](const Tag& tt){return tt == ac.tag();});
                                    //assert(ex == mf.parent_cell().end());

                                    //ex = std::find_if(gmf.parent_cell().begin(), gmf.parent_cell().end(), [&](const Tag& tt){return tt == ac.tag();});
                                    //assert(ex == gmf.parent_cell().end());

                                    //auto exx = std::find_if(mc.pnei().begin(), mc.pnei().end(), [&](const Tag& tt){return tt == ac.tag();});
                                    //assert(exx == mc.pnei().end());

                                    assert(acmf.parent_cell().size() == 1);

                                    mf.add_parent_cell(ac.tag());
                                    gmf.add_parent_cell(ac.tag());
                                    acmf.add_parent_cell(mc.tag());
                                    mc.add_pnei(ac.tag());
                                    ac.add_pnei(mc.tag());

                                    //if (mc.tag()() == 43 && tag_() == 0) {
                                    //assert(mf.parent_cell().size() != 1);
                                    //assert(gmf.parent_cell().size() != 1);
                                    //std::cout << "themf: " << mf.tag() << std::endl;
                                    //}

                                    //for (const MeshFace& mf: face())
                                    //  {
                                    //  if (mf.parent_cell().size() == 1)
                                    //  {
                                    //  if (mf.parent_cell()[0]() == 43 && tag_() == 0)
                                    //  {
                                    //  if (cell(Tag(mf.parent_cell()[0])).oga_cell_type() != OGA_cell_type_t::non_resident && cell(Tag(mf.parent_cell()[0])).oga_cell_type() != OGA_cell_type_t::ghost)
                                    //  {
                                    //  std::cout << "mesh: " << tag()() << std::endl;
                                    //  std::cout << "mc: " << mf.parent_cell()[0]() << std::endl;
                                    //  std::cout << "oga type: " << static_cast<int>(cell(Tag(mf.parent_cell()[0])).oga_cell_type()) << std::endl;
                                    //  std::cout << "mf btype: " << static_cast<int>(mf.btype()) << std::endl;
                                    //  std::cout << "mf: " << mf << std::endl;
                                    //  }
                                    //  assert(cell(Tag(mf.parent_cell()[0])).oga_cell_type() == OGA_cell_type_t::non_resident || cell(Tag(mf.parent_cell()[0])).oga_cell_type() == OGA_cell_type_t::ghost);
                                    //  assert(!mf.is_boundary());
                                    //  continue;
                                    //  }
                                    //  }
                                    //  }

                                    assert(acmf.parent_cell().size() == 2);

                                    assert(ac.tag() != mc.tag());
                                    if (acmf.parent_cell().size() == 2)
                                    {
                                        assert(acmf.parent_cell()[0] != acmf.parent_cell()[1]);
                                    }

                                    //if (ac.tag()() == 1088 && mc.tag()() == 1077)
                                    //{
                                    //std::cout << "gmt tag: " << gmf.tag() << std::endl;
                                    //for (const Tag& tt: gmf.parent_cell())
                                    //{
                                    //std::cout << "pc: " << tt() << std::endl;
                                    //}
                                    //assert(false);
                                    //}
                                    //acmf.set_tag(mf.tag());
                                    if (acmf.parent_cell().size() > 2)
                                    {
                                        std::cout << "acmf pc size: " << acmf.parent_cell().size() << std::endl;
                                        for (const Tag& t: acmf.parent_cell())
                                        {
                                            std::cout << "pc: " << t() << std::endl;
                                        }
                                        std::cout << "mf pc size: " << mf.parent_cell().size() << std::endl;
                                        std::cout << "mc pnei size: " << mc.pnei().size() << std::endl;
                                        std::cout << "ac pnei size: " << ac.pnei().size() << std::endl;
                                    }
                                    assert(acmf.parent_cell().size() <= 2);
                                    //assert(ac.pnei().size() >= 4);
                                    acfound = true;
                                    break;
                                }
                                if (acfound) {
                                    break;
                                }
                            }
                            if (rank == 13)
                            {
                                if (!acfound)
                                {
                                    // make a mesh from single cell.
                                    {
                                        std::cout << rank << " mesh: " << tag_() << std::endl;
                                        std::cout << rank << " mc: " << mc.tag()() << std::endl;
                                        std::cout << rank << " ac: " << actag << std::endl;
                                        std::cout << rank << " mc oga: " << static_cast<int>(mc.oga_cell_type()) << std::endl;
                                        std::cout << rank << " mf is bou: " << mf.is_boundary() << std::endl;
                                        std::cout << rank << " mf btype: " << static_cast<int>(mf.btype()) << std::endl;
                                        std::cout << rank << " ac type: " << actype << std::endl;
                                        std::cout << rank << " mf: " << mf.tag() << std::endl;
                                        for (const auto& pp: mf.face().vertex())
                                        {
                                            std::cout << rank << " pp: " << pp.r(0) << " " << pp.r(1) << " " << pp.r(2) << std::endl;

                                        }
                                        std::string fn = "acnotfound-";
                                        std::string fn2 = "tempmesh-";
                                        fn.append(std::to_string(rank));
                                        fn2.append(std::to_string(rank));
                                        fn.append(".vtk");
                                        fn2.append(".vtk");
                                        assert(!cell_.empty());
                                        print_as_vtk_geometry(fn);
                                        Mesh tempmesh;
                                        tempmesh.set_tag(Tag(9));
                                        assert(!arrival_cell.empty());
                                        for (const auto& cc: arrival_cell)
                                        {
                                            const MeshCell& ac = cc.cell();
                                            //std::cout << "acc: " << ac.parent_mesh()() << "-" << ac.tag()() << ",";
                                            std::cout << ac.tag()() << ",";
                                            if (ac.parent_mesh() != tag_) {
                                                continue;
                                            }
                                            assert(ac.tag().isvalid());
                                            if (tempmesh.query(ac.tag()) != nullptr)
                                            {
                                                for (const Nei& nn: arrival_cell)
                                                {
                                                    const MeshCell& acac = nn.cell();
                                                    std::cout << "acac: " << acac.parent_mesh()() << " " << acac.tag()() << std::endl;
                                                }
                                            }
                                            assert(tempmesh.query(ac.tag()) == nullptr);
                                            tempmesh.add_interior_cell(ac);
                                        }
                                        std::cout << std::endl;
                                        assert(!tempmesh.cell().empty());
                                        tempmesh.print_as_vtk_geometry(fn2);
                                    }
                                }
                            }
                            assert(acfound);
                        }
                    }
                }
            }
        }


        std::sort(to_be_added.begin(), to_be_added.end(), [&](int a, int b){return a < b;});
        to_be_added.erase(std::unique(to_be_added.begin(), to_be_added.end()), to_be_added.end());
        
        int count = 0;
        for (const Nei& nn: arrival_cell)
        {
            const MeshCell& ac = nn.cell();
            //if (ac.tag() != tag_) {
            if (ac.parent_mesh() != tag_) {
                continue;
            }

            if (query(ac.tag()) == nullptr)
            {
                ++count;
            }
        }

        //for (const MeshCell& ac: arrival_cell)
        for (const Nei& nn: arrival_cell)
        {
            const MeshCell& ac = nn.cell();
            //if (ac.tag() != tag_) {
            if (ac.parent_mesh() != tag_) {
                continue;
            }

            //if (rank == 2) {
                //assert(ac.tag()() != 317);
            //}
            auto iac = std::find_if(to_be_added.begin(), to_be_added.end(), [&](int i){return i == ac.tag()();});

            if (iac == to_be_added.end()) {
                continue;
            }

            //if (rank == 2) {
                //assert(ac.tag()() != 317);
            //}

            bool exist;
            auto iter = add_interior_cell_sorted(ac, exist);
            //assert(ac.tag()() != 242);
            //if (exist)
            //{
                //std::cout << "ac.tag()(): " << ac.tag()() << std::endl;
                //std::string fn = "exist-";
                //fn.append(std::to_string(rank));
                //fn.append(".vtk");
                //print_as_vtk(fn);
            //}
            assert(!exist);
            assert(query(ac.tag()) != nullptr);

            //iter->set_oga_cell_type(OGA_cell_type_t::non_resident);
            iter->set_oga_cell_type(OGA_cell_type_t::ghost);
            iter->set_btype(BouType::undefined);

            //if (rank == 2) {
                //assert(query(Tag(317)) == nullptr);
            //}

        }

        //for (MeshCell& mc: cell_)
        //{
        //    //if (mc.tag()() == 1047 && tag_() == 1)
        //    //{
        //        //std::cout << "rankconnect: " << rank << std::endl;
        //        //std::cout << "isresi: " << is_resi(mc.poly().centroid()) << std::endl;
        //        //std::cout << "oga: " << static_cast<int>(mc.oga_cell_type()) << std::endl;
        //    //}
        //    if (mc.oga_cell_type() == OGA_cell_type_t::ghost)
        //    {
        //        continue;
        //    }
        //    if (!is_resi(mc.poly().centroid())) // duplicate and redundant???
        //    //if (!sp_aabb.do_intersect(mc.poly().centroid()))
        //    {
        //        mc.set_oga_cell_type(OGA_cell_type_t::non_resident);
        //    }
        //}

        //for (const MeshFace& mf: face_)
        //{
        //        for (const Tag& pc: mf.parent_cell())
        //        {
        //            int count = 0;
        //            if (query(pc) == nullptr)
        //            {
        //            }
        //            if (cell(pc).oga_cell_type() == OGA_cell_type_t::non_resident)
        //            {
        //                ++count;
        //            }
        //        }
        //        assert(count != 0);
        //}

        //if (rank == 7)
        if (!arrival_cell.empty())
        {
            //if (to_be_added.size() != arrival_cell.size())
            if (to_be_added.size() != count)
            {
                std::cout << "to_be_added.size(): " << to_be_added.size() << std::endl;
                std::cout << "count: " << count << std::endl;
                //for (const MeshCell& ac: arrival_cell)
                //for (const Nei& cc: arrival_cell)
                //{
                    //const MeshCell& ac = cc.cell();
                    //if (ac.tag() != tag_) {
                        //continue;
                    //}
                    //std::cout << ac.tag()() << ",";
                //}
                std::cout << std::endl;
                //for (int i: to_be_added)
                //{
                    //int co = std::count_if(arrival_cell.begin(), arrival_cell.end(), [&](const MeshCell& c){return c.tag()() == i;});
                    //if (co == 0) {
                        //std::cout << i << ",";
                    //}
                //}
                //std::cout << std::endl;

                std::cout << "to be added: ";
                for (int i: to_be_added)
                {
                    std::cout << i << ",";
                }
                std::cout << std::endl;

                //for (int i: non_resi)
                //{
                    //std::cout << i << ",";
                //}
                //std::cout << std::endl;

                print_as_vtk("noteqqq.vtk");

            }
            assert(to_be_added.size() == count);
        }

        //if (rank == 1)
        {
        for (const MeshFace& mf: face_)
        {
            if (mf.is_boundary()) {
                continue;
            }
            if (mf.btype() == BouType::partition) {
                continue;
            }
            for (const Tag& t: mf.parent_cell())
            {
                bool found = false;
                if (query(t) == nullptr)
                {
                    std::cout << rank << " btype: " << static_cast<int>(mf.btype()) << std::endl;
                    std::cout << rank << " pc size: " << mf.parent_cell().size() << std::endl;
                    std::cout << rank << " pc 0: " << mf.parent_cell()[0]() << std::endl;
                    std::cout << rank << " pc 1: " << mf.parent_cell()[1]() << std::endl;
                    std::cout << rank << " t: " << t() << std::endl;
                    print_as_vtk("bb.vtk");
                }
                assert(query(t) != nullptr);
                const auto& pc = cell(t);
                if (pc.oga_cell_type() == OGA_cell_type_t::non_resident) {
                    continue;
                }
                for (const MeshFace& cmf: pc.face())
                {
                    if (cmf.tag() == mf.tag())
                    {
                        found = true;
                        break;
                    }
                }
                if (!found)
                {
                    std::cout << "mf tag size : " << mf.tag()().size() << std::endl;
                    //std::cout << "mf: " << mf.tag()()[0] << " " << mf.tag()()[1] << " " << mf.tag()()[2] << std::endl;
                    std::cout << "typeeeeeeee: " << static_cast<int>(pc.oga_cell_type()) << std::endl;
                }
                assert(found);
            }
        }
        }
        //for (const MeshCell& mc: cell_)
        //{
        //    for (const MeshFace& mf: mc.face())
        //    {
        //        auto it = std::find_if(face_.begin(), face_.end(), [&](const MeshFace& f){return f.tag() == mf.tag();});
        //        if (it == face_.end())
        //        {
        //            std::cout << "rank: " << rank << std::endl;
        //            std::cout << "mc: " << mc.tag()() << std::endl;
        //            std::cout << "mf: " << mf.tag()() << std::endl;
        //            std::cout << "face_ size: " << face_.size() << std::endl;
        //            std::cout << "mc ogt: " << static_cast<int>(mc.oga_cell_type()) << std::endl;
        //        }
        //        assert(it != face_.end());
        //    }
        //}

        for (MeshCell& mc: cell_)
        {
            assert(mc.oga_cell_type() != OGA_cell_type_t::field);
        }

        for (const MeshFace& mf: face())
        {
            if (mf.parent_cell().size() == 1)
            {
                const auto& mc = cell(Tag(mf.parent_cell()[0]));
                if (mc.oga_cell_type() != OGA_cell_type_t::non_resident && mc.oga_cell_type() != OGA_cell_type_t::ghost)
                {
                    //if (mf.parent_cell()[0]() == 43)
                    {
                    std::cout << "mesh: " << tag()() << std::endl;
                    std::cout << "mc: " << mc.tag()() << std::endl;
                    std::cout << "oga type: " << static_cast<int>(mc.oga_cell_type()) << std::endl;
                    std::cout << "mf btype: " << static_cast<int>(mf.btype()) << std::endl;
                    std::cout << "mc btype: " << static_cast<int>(mc.btype()) << std::endl;
                    std::cout << "mf: " << mf.tag() << std::endl;
                    std::cout << "rank: " << rank << std::endl;
                    //std::cout << "mc is resi: " << is_resi(mc.poly().centroid()) << std::endl;
                    print_as_vtk("dodo.vtk");
                    }
                }
                assert(mc.oga_cell_type() == OGA_cell_type_t::non_resident || mc.oga_cell_type() == OGA_cell_type_t::ghost);
                assert(!mf.is_boundary());
            }
        }
    }*/

    void Mesh::connect_partition_cells(ArrCon<Nei>& arrival_cell, int rank, std::function<bool(const vec3<double>&, int celltag)> is_resi, Profiler* profiler)
    {
        //for (const auto& mf: face_)
        //{
        //    if (mf.tag() == FaceTag(std::vector<int>{374, 408, 481}))
        //    {
        //        std::cout << "pc size: " << mf.const_parent_cell().size() << std::endl;
        //        for (const auto& pcc: mf.const_parent_cell())
        //        {
        //            std::cout << "pc: " << pcc.tag()() << std::endl;
        //        }
        //        assert(false);
        //    }
        //}
        //profiler->start("sol-cpc-1");
        //for (MeshCell& mc: cell_)
        //{
            //if (!is_resi(mc.poly().centroid(), mc.tag()()))
            //{
                //mc.set_oga_cell_type(OGA_cell_type_t::non_resident);
            //}
        //}

        for (Nei& nn: arrival_cell)
        {
            MeshCell& ac = nn.cell_;
            ac.remove_all_neighbors();
            ac.remove_all_nonself_parents_of_faces();
            //for (auto& mf: ac.face_p())
            //{
                //mf.set_faceaddr(nullptr);
            //}
        }
        //profiler->stop("sol-cpc-1");
        
        std::vector<int> to_be_added;

        //profiler->start("sol-cpc-2");
        if (!arrival_cell.empty())
        {
            for (MeshCell& mc: cell_)
            {
                assert(mc.parent_mesh() == tag_);
                if (mc.parent_mesh() != tag_) {
                    continue;
                }

                if (mc.btype() == BouType::partition && mc.oga_cell_type() != OGA_cell_type_t::non_resident)
                {
                    for (MeshFace& mf: mc.face_)
                    {
                        if (mf.parent_cell().size() == 1)
                        {
                            bool acfound = false;
                            int actag = -1;
                            int actype = -1;
                            for (Nei& nn: arrival_cell)
                            {
                                MeshCell& ac = nn.cell_;
                                if (ac.parent_mesh() != tag_) {
                                    continue;
                                }
                                actag = ac.tag()();
                                actype = static_cast<int>(ac.oga_cell_type());


                                assert(ac.tag().isvalid());
                                if (mc.tag() == ac.tag()) {
                                    continue;
                                }


                                    //if (tag_() == 0 && mc.tag()() == 16189)
                                    //{
                                    //    std::cout << "fffffffffffffff: " << mf.tag() << std::endl;
                                    //    std::cout << "ac: " << ac.tag()() << std::endl;
                                    //    //std::cout << "acmf: " << acmf.tag() << std::endl;
                                    //    //if (ac.tag()() == 17707)
                                    //    {
                                    //        for (const auto& mf: mc.face())
                                    //        {
                                    //            std::cout << "mfff face: " << mf.tag() << std::endl;
                                    //        }
                                    //        for (const auto& mf: ac.face())
                                    //        {
                                    //            std::cout << "nccc face: " << mf.tag() << std::endl;
                                    //        }
                                    //        //assert(false);
                                    //    }
                                    //}

                                for (MeshFace& acmf: ac.face_p())
                                {
                                    //if (tag_() == 0 && mc.tag()() == 16189)
                                    //{
                                    //    //std::cout << "mf: " << mf.tag() << std::endl;
                                    //    //std::cout << "acmf: " << acmf.tag() << std::endl;
                                    //    if (ac.tag()() == 17707)
                                    //    {
                                    //        //for (const auto& mf: mc.face())
                                    //        //{
                                    //            //std::cout << "mf face: " << mf.tag() << std::endl;
                                    //        //}
                                    //        //for (const auto& mf: ac.face())
                                    //        //{
                                    //            //std::cout << "nc face: " << mf.tag() << std::endl;
                                    //        //}
                                    //        //assert(false);
                                    //    }
                                    //}
                                    if (mf.tag() != acmf.tag()) {
                                        continue;
                                    }



                                    to_be_added.push_back(ac.tag()());

                                    //auto fiter = std::find_if(face_.begin(), face_.end(), [&](const auto& f){return f.tag() == mf.tag();});
                                    //auto fiter = std::find_if(face_.begin(), face_.end(), [&](const auto& f){return f.tag() == mf.tag();});

                                    //MeshFace& gmf = face_p(mf.tag());
                                    //MeshFace& gmf = *mf.faceaddr();

                                    //bool exist;
                                    //auto iter = add_interior_cell_sorted(ac, exist);

                                    //iter->set_oga_cell_type(OGA_cell_type_t::ghost);
                                    //iter->set_btype(BouType::undefined);

                                    const auto& tam = mc.tag();
                                    const auto& tac = ac.tag();

                                    //mf.add_parent_cell(ac.tag());
                                    mf.add_parent_cell(tac);
                                    if (mf.parent_cell().size() > 2)
                                    {
                                        std::cout << "rank: " << rank << std::endl;
                                        std::cout << "mf tag: " << mf.tag() << std::endl;
                                        std::cout << "mf btype: " << static_cast<int>(mf.btype()) << std::endl;
                                        std::cout << "pc size: " << mf.parent_cell().size() << std::endl;
                                        for (const auto& pcc: mf.parent_cell())
                                        {
                                            std::cout << "pc: " << pcc() << std::endl;
                                        }
                                    }
                                    assert(mf.parent_cell().size() <= 2);
                                    //gmf.add_parent_cell(ac.tag());
                                    //gmf.add_parent_cell(tac);
                                    //assert(gmf.parent_cell().size() <= 2);
                                    //acmf.add_parent_cell(mc.tag());
                                    acmf.add_parent_cell(tam);
                                    assert(acmf.parent_cell().size() <= 2);
                                    //mc.add_pnei(ac.tag());
                                    mc.add_pnei(tac);
                                    //ac.add_pnei(mc.tag());
                                    ac.add_pnei(tam);

                                    //acmf.set_faceaddr(&gmf);

                                    acfound = true;
                                    assert(mf.parent_cell().size() == 2);
                                    break;
                                }
                                if (acfound) {
                                    break;
                                }
                            }
                            //if (rank == 7)
                            {
                                if (!acfound)
                                {
                                    // make a mesh from single cell.
                                    {
                                        std::cout << rank << " mesh: " << tag_() << std::endl;
                                        std::cout << rank << " mc: " << mc.tag()() << std::endl;
                                        std::cout << rank << " ac: " << actag << std::endl;
                                        std::cout << rank << " mc oga: " << static_cast<int>(mc.oga_cell_type()) << std::endl;
                                        std::cout << rank << " mf is bou: " << mf.is_boundary() << std::endl;
                                        std::cout << rank << " mf btype: " << static_cast<int>(mf.btype()) << std::endl;
                                        std::cout << rank << " ac type: " << actype << std::endl;
                                        std::cout << rank << " mf: " << mf.tag() << std::endl;
                                        for (const auto& pp: mf.face().vertex())
                                        {
                                            std::cout << rank << " pp: " << pp.r(0) << " " << pp.r(1) << " " << pp.r(2) << std::endl;

                                        }
                                        std::string fn = "acnotfound-";
                                        std::string fn2 = "tempmesh-";
                                        fn.append(std::to_string(rank));
                                        fn2.append(std::to_string(rank));
                                        fn.append(".vtk");
                                        fn2.append(".vtk");
                                        assert(!cell_.empty());
                                        print_as_vtk_geometry(fn);
                                        Mesh tempmesh;
                                        tempmesh.set_tag(Tag(9));
                                        assert(!arrival_cell.empty());
                                        for (const auto& cc: arrival_cell)
                                        {
                                            const MeshCell& ac = cc.cell();
                                            //std::cout << "acc: " << ac.parent_mesh()() << "-" << ac.tag()() << ",";
                                            std::cout << ac.tag()() << ",";
                                            if (ac.parent_mesh() != tag_) {
                                                continue;
                                            }
                                            assert(ac.tag().isvalid());
                                            //if (tempmesh.query(ac.tag()) != nullptr)
                                            //{
                                            //    for (const Nei& nn: arrival_cell)
                                            //    {
                                            //        const MeshCell& acac = nn.cell();
                                            //        std::cout << "acac: " << acac.parent_mesh()() << " " << acac.tag()() << std::endl;
                                            //    }
                                            //}
                                            //assert(tempmesh.query(ac.tag()) == nullptr);
                                            tempmesh.add_interior_cell(ac);
                                        }
                                        std::cout << std::endl;
                                        assert(!tempmesh.cell().empty());
                                        tempmesh.print_as_vtk_geometry(fn2);
                                    }
                                }
                            assert(acfound);
                            }
                        }
                    }
                }
            }
        }
        //profiler->stop("sol-cpc-2");

        //profiler->start("sol-cpc-3");
        std::sort(to_be_added.begin(), to_be_added.end(), [&](int a, int b){return a < b;});
        to_be_added.erase(std::unique(to_be_added.begin(), to_be_added.end()), to_be_added.end());
        //profiler->stop("sol-cpc-3");
        
        //profiler->start("sol-cpc-4");
        for (const Nei& nn: arrival_cell)
        {
            const MeshCell& ac = nn.cell();
            if (ac.parent_mesh() != tag_) {
                continue;
            }

            auto iac = std::find_if(to_be_added.begin(), to_be_added.end(), [&](int i){return i == ac.tag()();});
            if (iac == to_be_added.end()) {
                continue;
            }

            bool exist;
            auto iter = add_interior_cell_sorted(ac, exist);

            iter->set_oga_cell_type(OGA_cell_type_t::ghost);
            iter->set_btype(BouType::undefined);
        }
        //profiler->stop("sol-cpc-4");

        //for (const auto& mc: cell_)
        //{
            //for (const auto& mf: mc.face())
            //{
                //if (mf.btype() == BouType::interior)
                //{
                    //assert(mf.parent_cell().size() == 2);
                //}
            //}
        //}
        for (const auto& mc: cell_)
        {
            for (const auto& mf: mc.face_)
            {
                if (mf.is_boundary()) {
                    continue;
                }
                if (mf.parent_cell(0) != mc.tag())
                {
                    std::cout << "mc: " << mc.tag()() << std::endl;
                    std::cout << "pc 0: " << mf.parent_cell(0)() << std::endl;
                    std::cout << "pc 1: " << mf.parent_cell(1)() << std::endl;
                    std::cout << "mf btype: " << static_cast<int>(mf.btype()) << std::endl;
                }
                assert(mf.parent_cell(0) == mc.tag());
            }
        }
    }
}
