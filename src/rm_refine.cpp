#include "regular_mesh.h"

namespace Tailor
{
    const RegularMesh* RegularMesh::refine_adaptive(const std::deque<Mesh>& meshes, const BinRMTag& heaviest_bt, int rank, int& new_rmtag, int& new_bt, bool pseudo3D, RegType regtype)
    {
        //auto _rmtag = rmtag_address_map_.find(heaviest_bt.rmtag());
        //assert(_rmtag != rmtag_address_map_.end());

        //auto _rm = _rmtag->second;
        //assert(_rm != nullptr);

        //Bin& heavy_bin = _rm->bin_p(heaviest_bt.bintag());
        Bin& heavy_bin = bin_p(heaviest_bt.bintag());
        new_rmtag = get_new_rmtag();
        new_bt = get_new_bt();
        //if (new_bt >= 20)
        //{
            //std::cout << "new bt: " << new_bt << std::endl;
            //std::cout << "bin size: " << bin_.size() << std::endl;
            //std::cout << "rm tag: " << tag()() << std::endl;
        //}
        //assert(new_bt < 20);
        heavy_bin.init_rm(new_bt, Tag(new_rmtag), nstripe_, pseudo3D);
        assert(aabb_.face_cross_len() == false);
        assert(heavy_bin.aabb().face_cross_len() == false);
        assert(heavy_bin.rm()->aabb().face_cross_len() == false);
        {
            if (!heavy_bin.mesh_tag_index_map_res().left.empty())
            {
                assert(!heavy_bin.cell().empty());
            }
        }
        insert_to_rmtag_address_map(Tag(new_rmtag), heavy_bin.rm());
        insert_bin_addresses(heavy_bin.rm());
        assert(!bintag_address_map_.empty());
        //for (const auto& b: heavy_bin.rm()->bin()) {
            //insert_to_bintag_address_map(b.tag(), &b);
        //}

        /*if (rank == 0)
        {
            for (const Bin& b: heavy_bin.rm()->bin())
            {
                std::cout << "rm size: " << size() << std::endl;
                std::cout << "heavy bin: " << tag_() << " " << heavy_bin.tag()() << std::endl;
                std::cout << "aabb_.min(0): " << aabb_.min(0) << std::endl;
                std::cout << "aabb_.min(1): " << aabb_.min(1) << std::endl;
                std::cout << "aabb_.min(2): " << aabb_.min(2) << std::endl;
                std::cout << "aabb_.max(0): " << aabb_.max(0) << std::endl;
                std::cout << "aabb_.max(1): " << aabb_.max(1) << std::endl;
                std::cout << "aabb_.max(2): " << aabb_.max(2) << std::endl;
                std::cout << "heavy_bin->aabb().min(0): " << heavy_bin.aabb().min(0) << std::endl;
                std::cout << "heavy_bin->aabb().min(1): " << heavy_bin.aabb().min(1) << std::endl;
                std::cout << "heavy_bin->aabb().min(2): " << heavy_bin.aabb().min(2) << std::endl;
                std::cout << "heavy_bin->aabb().max(0): " << heavy_bin.aabb().max(0) << std::endl;
                std::cout << "heavy_bin->aabb().max(1): " << heavy_bin.aabb().max(1) << std::endl;
                std::cout << "heavy_bin->aabb().max(2): " << heavy_bin.aabb().max(2) << std::endl;
                std::cout << "heavy_bin->rm()->aabb().min(0): " << heavy_bin.rm()->aabb().min(0) << std::endl;
                std::cout << "heavy_bin->rm()->aabb().min(1): " << heavy_bin.rm()->aabb().min(1) << std::endl;
                std::cout << "heavy_bin->rm()->aabb().min(2): " << heavy_bin.rm()->aabb().min(2) << std::endl;
                std::cout << "heavy_bin->rm()->aabb().max(0): " << heavy_bin.rm()->aabb().max(0) << std::endl;
                std::cout << "heavy_bin->rm()->aabb().max(1): " << heavy_bin.rm()->aabb().max(1) << std::endl;
                std::cout << "heavy_bin->rm()->aabb().max(2): " << heavy_bin.rm()->aabb().max(2) << std::endl;
                std::cout << "b->aabb().min(0): " << b.aabb().min(0) << std::endl;
                std::cout << "b->aabb().min(1): " << b.aabb().min(1) << std::endl;
                std::cout << "b->aabb().min(2): " << b.aabb().min(2) << std::endl;
                std::cout << "b->aabb().max(0): " << b.aabb().max(0) << std::endl;
                std::cout << "b->aabb().max(1): " << b.aabb().max(1) << std::endl;
                std::cout << "b->aabb().max(2): " << b.aabb().max(2) << std::endl;
                print();
                for (const auto& ff: b.aabb().faces())
                {
                    std::cout << "face vertices" << std::endl;
                    for (const auto& v: ff.vertex())
                    {
                        std::cout << v.r(0) << " " << v.r(1) << " " << v.r(2) << std::endl;
                    }
                }
                assert(b.aabb().face_cross_len() == false);
            }
        }*/

        for (const Bin& b: heavy_bin.rm()->bin())
        {
            //if (b.cell().empty())
            //{
                //std::cout << b.parent_rm()() << " " << b.tag()() << std::endl;
            //}
            //assert(!b.cell().empty());
            if (b.aabb().face_cross_len() != false)
            {
                std::cout << "rm suze: " << size() << std::endl;
            }
            assert(b.aabb().face_cross_len() == false);
        }

        heavy_bin.register_cells_to_rm(meshes, pseudo3D, rank, regtype);

            if (!heavy_bin.mesh_tag_index_map_res().left.empty())
            {
                assert(!heavy_bin.cell().empty());
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

            /*for (const auto& b: bin_)
            {
                if (b.tag()() == 17)
                {
                    for (const auto& bc: b.4))).cell())
                    {
                        if (bc.mesh()() == 1 && bc.cell()() == 274)
                        {
                            std::cout << "zzzzzzzzzzzzzzz: " << rank << std::endl;
                            assert(false);
                        }
                    }
                }
            }*/


        return heavy_bin.rm();
    }
    
    void RegularMesh::refine(const std::vector<Mesh>& meshes, int rank, bool pseudo3D)
    {
        // Define new step length.
        Vector3Int new_nstripe;
        if (pseudo3D)
        {
            new_nstripe = Vector3Int(2 * nstripe(0), 2 * nstripe(1), 1);
        }
        else
        {
            new_nstripe = Vector3Int(2 * nstripe(0), 2 * nstripe(1), 2 * nstripe(2));
        }
        Vector3 newh = h()/2.;
        // Create new bins.
        std::vector<Bin> newbin;
        newbin.resize(new_nstripe(0)*new_nstripe(1)*new_nstripe(2));

        // set indices of new bins.
        for (int depth=0; depth<nstripe(2); ++depth)
        {
            for (int row=0; row<nstripe(1); ++row)
            {
                for (int col=0; col<nstripe(0); ++col)
                {
                    // indices of lower left quadrant.
                    int ll_col = 2. * col;
                    int ll_row = 2. * row;
                    int ll_depth = 2. * depth;
                    newbin[ll_col + ll_row*new_nstripe(0) + ll_depth*new_nstripe(0)*new_nstripe(1)].set_index(RegularMeshIndex(ll_row, ll_col, ll_depth));

                    int col_ = ll_col + 1;
                    int row_ = ll_row;
                    int depth_ = ll_depth;
                    newbin[col_ + row_*new_nstripe(0) + depth_*new_nstripe(0)*new_nstripe(1)].set_index(RegularMeshIndex(row_, col_, depth_));

                    col_ = ll_col;
                    row_ = ll_row + 1;
                    depth_ = ll_depth;
                    newbin[col_ + row_*new_nstripe(0) + depth_*new_nstripe(0)*new_nstripe(1)].set_index(RegularMeshIndex(row_, col_, depth_));

                    col_ = ll_col + 1;
                    row_ = ll_row + 1;
                    depth_ = ll_depth;
                    newbin[col_ + row_*new_nstripe(0) + depth_*new_nstripe(0)*new_nstripe(1)].set_index(RegularMeshIndex(row_, col_, depth_));

                    if (!pseudo3D)
                    {
                        col_ = ll_col;
                        row_ = ll_row;
                        depth_ = ll_depth + 1;
                        newbin[col_ + row_*new_nstripe(0) + depth_*new_nstripe(0)*new_nstripe(1)].set_index(RegularMeshIndex(row_, col_, depth_));

                        col_ = ll_col + 1;
                        row_ = ll_row;
                        depth_ = ll_depth + 1;
                        newbin[col_ + row_*new_nstripe(0) + depth_*new_nstripe(0)*new_nstripe(1)].set_index(RegularMeshIndex(row_, col_, depth_));

                        col_ = ll_col;
                        row_ = ll_row + 1;
                        depth_ = ll_depth + 1;
                        newbin[col_ + row_*new_nstripe(0) + depth_*new_nstripe(0)*new_nstripe(1)].set_index(RegularMeshIndex(row_, col_, depth_));

                        col_ = ll_col + 1;
                        row_ = ll_row + 1;
                        depth_ = ll_depth + 1;
                        newbin[col_ + row_*new_nstripe(0) + depth_*new_nstripe(0)*new_nstripe(1)].set_index(RegularMeshIndex(row_, col_, depth_));
                    }
                }
            }
        }

        // iterate through cells of unrefined bin and register them to new bins.
        // register cell to all overlapping bins in AABB.
        for (int depth=0; depth<nstripe(2); ++depth)
        {
            for (int row=0; row<nstripe(1); ++row)
            {
                for (int col=0; col<nstripe(0); ++col)
                {
                    RegularMeshIndex indorg(row, col, depth);
                    Bin& binorg = bin_p(indorg); // reference to original bin.
                    int ll_col = 2. * col;
                    int ll_row = 2. * row;
                    int ll_depth;
                    if (pseudo3D) {
                        ll_depth= depth;
                    }
                    else {
                        ll_depth= 2. * depth;
                    }
                    int col_ = ll_col + 1;
                    int row_ = ll_row + 1;
                    int depth_;
                    if (pseudo3D) {
                        depth_ = ll_depth;
                    }
                    else {
                        depth_ = ll_depth + 1;
                    }
                    RegularMeshIndex ind(ll_row, ll_col, ll_depth); // lower-left index of new bin.
                    // loop through new bins which are subset of unrefined bin.
                    for (int k=ll_depth; k<=depth_; ++k)
                    {
                        for (int i=ll_row; i<=row_; ++i)
                        {
                            for (int j=ll_col; j<=col_; ++j)
                            {
                                ind.set(i, j, k);
                                Bin& _bin = newbin[ind.j()+ind.i()*new_nstripe(0)+ind.k()*new_nstripe(0)*nstripe(1)]; // reference to new bin.
                                //_bin.resize_mesh_load(binorg.mesh_load().size());
                                // make a quad to represent new bin.
                                double d0 = aabb().min(0) + ind.j() * newh(0);
                                double d1 = aabb().min(1) + ind.i() * newh(1);
                                double d2 = aabb().min(2) + ind.k() * newh(2);
                                if (pseudo3D) {
                                    assert(ind.k() == 0);
                                }
                                Vector3 llc(d0, d1, d2);
                                Vector3 ur = llc + newh;
                                if (pseudo3D) {
                                    ur(2) = llc(2) + h(2);
                                }
                                _bin.set_aabb(AABB(llc, ur));
                                //std::vector<Point> pts = {Point(llc), Point(llc(0)+newh(0), llc(1)), Point(llc+newh), Point(llc(0), llc(1)+newh(1))};
                                //Polygon _quad(pts);
                                //AABB quad(_quad);
                                //Polytope quad(pts);
                                // loop through bincells of unrefined bin.
                                for (int t=0; t<binorg.cell().size(); ++t)
                                {
                                    Tag mt = binorg.cell(t).mesh();
                                    Tag ct = binorg.cell(t).cell();
                                    for (auto m=meshes.begin(); m!=meshes.end(); ++m)
                                    {
                                        if (mt == m->tag())
                                        {
                                            if (_bin.aabb().do_intersect(AABB(m->cell(ct).poly())))
                                                //if (quad.do_intersect(*m->cell(ct).polytope()))
                                            {
                                                //BinCell bc;
                                                //bc.set_cell(ct);
                                                //bc.set_mesh(mt);
                                                assert(ct.isvalid());
                                                assert(mt.isvalid());
                                                //_bin.add_bincell(bc);
                                                _bin.add_bincell(ct, mt, rank);

                                                //if (!m->cell(ct).is_ghost())
                                                {
                                                    bool point_in_bin = _bin.aabb().do_intersect(m->cell(ct).poly().centroid());
                                                    if (point_in_bin)
                                                    {
                                                        _bin.increment_mesh_load(mt);
                                                    }
                                                    else
                                                    {
                                                        _bin.insert_to_mesh_tag_index_map_all(mt);
                                                    }
                                                }
                                            }
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                }
                }
            }
    }

        set_nstripe(new_nstripe);
        set_h(newh);
        bin_ = newbin;
        assert(bin_.size() == newbin.size());
    }

    /*void RegularMesh::refine(const std::vector<Mesh>& meshes)
      {
    // Define new step length.
    unsigned int new_nstripe = 2. * nstripe;
    vec3 newh = h/2.;
    // Create new bins.
    std::vector<Bin> newbin;
    newbin.resize(new_nstripe*new_nstripe, Bin(bin[0].mesh_load.size()));

    for (int row=0; row<nstripe; ++row)
    {
    for (int col=0; col<nstripe; ++col)
    {
    // indices of lower left quadrant.
    int ll_col = 2. * col;
    int ll_row = 2. * row;
    newbin[ll_col + ll_row*new_nstripe].row = ll_row;
    newbin[ll_col + ll_row*new_nstripe].col = ll_col;

    //
    int col_ = ll_col + 1;
    int row_ = ll_row;
    newbin[col_ + row_*new_nstripe].row = row_;
    newbin[col_ + row_*new_nstripe].col = col_;

    //
    col_ = ll_col;
    row_ = ll_row + 1;
    newbin[col_ + row_*new_nstripe].row = row_;
    newbin[col_ + row_*new_nstripe].col = col_;

    //
    col_ = ll_col + 1;
    row_ = ll_row + 1;
    newbin[col_ + row_*new_nstripe].row = row_;
    newbin[col_ + row_*new_nstripe].col = col_;
    }
    }

    // iterate through cells of unrefined bin and register them to new bins.
    // register cell to all overlapping bins in AABB.
    for (int row=0; row<nstripe; ++row)
    {
    for (int col=0; col<nstripe; ++col)
    {
    Index indorg(row, col, 0);
    Bin& binorg = (*this)(indorg);
    int ll_col = 2. * col;
    int ll_row = 2. * row;
    int col_ = ll_col + 1;
    int row_ = ll_row + 1;
    Index ind(ll_row, ll_col, 0);
    for (ind.i=ll_row; ind.i<=row_; ++ind.i)
    {
    for (ind.j=ll_col; ind.j<=col_; ++ind.j)
    {
    Bin& bin_ = newbin[ind.j+ind.i*new_nstripe];
    vec3 llc;
    llc[0] = aabb.min[0] + ind.j * newh[0];
    llc[1] = aabb.min[1] + ind.i * newh[1];
    std::vector<Point> pts = {Point(llc), Point(llc[0]+newh[0], llc[1]), Point(llc+newh), Point(llc[0], llc[1]+newh[1])};
    Polygon quad(pts);
    for (int t=0; t<binorg.cell.size(); ++t)
    {
    int mt = binorg.cell[t].mesh;
    int ct = binorg.cell[t].cell;
    if (quad.do_intersect(meshes[mt].cell[ct].polygon))
    {
    BinCell bc;
    bc.cell = ct;
    bc.mesh = mt;
    //bin_.cell_tag.push_back(ct);
    //bin_.mesh_tag.push_back(mt);
    // make point in bin intersection test to determine centroid owner.
    vec3 p = meshes[mt].cell[ct].polygon.centroid().r;
    bool point_in_bin = false;
    if (p[0] >= llc[0] && p[0] <= llc[0]+newh[0])
    {
        if (p[1] >= llc[1] && p[1] <= llc[1]+newh[1])
        {
            point_in_bin = true;
        }
    }
    //bin_.centroid_owner.push_back(point_in_bin);
    bc.resident = point_in_bin;
    bin_.cell.push_back(bc);
    if (point_in_bin)
    {
        ++bin_.mesh_load[mt];
    }
}
}
}
}
}
}

nstripe = new_nstripe;
h = newh;
bin = newbin;
}*/
}
