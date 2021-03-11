#include "regular_mesh.h"
#include <iomanip>

namespace Tailor
{
    void RegularMesh::get_adjacency_from_graph(const Graph& graph)
    {
        auto pair = graph.id_to_hex().find(tag_());
        assert(pair != graph.id_to_hex().end());

        auto hex = pair->second;
        assert(hex != nullptr);

        for (auto& bin: bin_)
        {
            if (bin.rm() == nullptr)
            {
                auto vtx = std::find_if(hex->vertex_.begin(), hex->vertex_.end(), [&](const auto& v){return v.id_ == bin.tag()();});
                assert(vtx != hex->vertex_.end());

                for (const auto nei: vtx->nei_)
                {
                    assert(nei != nullptr);

                    auto rmtag = nei->parent_hex_;
                    auto bintag = nei->id_;

                    bin.add_nei(rmtag, bintag);
                }
            }
            else
            {
                bin.rm()->get_adjacency_from_graph(graph);
            }
        }
    }

    void make_adt_recur(const RegularMesh& rm, ADT& adt)
    {
        for (const auto& bin: rm.bin())
        {
            if (bin.rm() == nullptr)
            {
                const auto& vertices = bin.aabb().vertices();
                adt.insert(ADTPoint(vertices.begin(), vertices.end(), bin.tag()()));
            }
            else
            {
                make_adt_recur(*bin.rm(), adt);
            }
        }
    }
    void make_adt(const RegularMesh& rm, ADT& adt, int rank)
    {
        std::vector<double> dim(TAILOR_ADT_VAR);
        const auto& aabb = rm.aabb();
        dim[0] = aabb.min(0);
        dim[1] = aabb.max(0);
        dim[2] = aabb.min(1);
        dim[3] = aabb.max(1);
        dim[4] = aabb.min(2);
        dim[5] = aabb.max(2);

        adt = ADT(dim);

        assert(rm.aabb().min(0) == adt.aabb().min(0));
        assert(rm.aabb().min(1) == adt.aabb().min(1));
        assert(rm.aabb().min(2) == adt.aabb().min(2));
        assert(rm.aabb().max(0) == adt.aabb().max(0));
        assert(rm.aabb().max(1) == adt.aabb().max(1));
        assert(rm.aabb().max(2) == adt.aabb().max(2));

        //std::cout << "dim(0): " << dim[0] << std::endl;
        //std::cout << "dim(1): " << dim[1] << std::endl;
        //std::cout << "dim(2): " << dim[2] << std::endl;
        //std::cout << "dim(3): " << dim[3] << std::endl;
        //std::cout << "dim(4): " << dim[4] << std::endl;
        //std::cout << "dim(5): " << dim[5] << std::endl;

        for (const auto& bin: rm.bin())
        {
            if (bin.rm() == nullptr)
            {
                const auto& vertices = bin.aabb().vertices();
                adt.insert(ADTPoint(vertices.begin(), vertices.end(), bin.tag()()));
                
                //for (auto i: vertices)
                    //std::cout << "r(2): " << i.r(2) << std::endl;
            }
            else
            {
                make_adt_recur(*bin.rm(), adt);
            }
        }

        assert(rm.aabb().min(0) == adt.aabb().min(0));
        assert(rm.aabb().min(1) == adt.aabb().min(1));
        assert(rm.aabb().min(2) == adt.aabb().min(2));
        assert(rm.aabb().max(0) == adt.aabb().max(0));
        assert(rm.aabb().max(1) == adt.aabb().max(1));
        assert(rm.aabb().max(2) == adt.aabb().max(2));
    }

    void RegularMesh::flood(const Tag& start, const AABB& aabb)
    {
        auto current = start;
        assert(start.isvalid());
        assert(start() < bin_.size());
        Bin& bin = bin_p(Tag(start));
        bin.set_load(-1);

        //for (int t: pneis_of_bin(current))
        for (const auto& t: pneis_of_bin(current()))
        {
            if (!t.isvalid())
            {
                continue;
            }

            assert(t() < bin_.size());
            Bin& nei = bin_p(t);

            if (nei.load() == -2 || nei.load() == -1)
            {
                continue;
            }

            if (!aabb.do_intersect(nei.aabb()))
            {
                continue;
            }

            flood(t, aabb);
        }
    }

    void RegularMesh::move(const Vector3& v)
    {
        aabb_.move_points(v);
        for (Bin& b: bin_)
        {
            b.move(v);
        }
    }

    void RegularMesh::rotate(double ang, int axis, const Vector3& rot_point)
    {
        aabb_.rotate_points(ang, axis, rot_point);
        for (Bin& b: bin_)
        {
            b.rotate(ang, axis, rot_point);
        }
    }

    bool RegularMesh::operator<(const RegularMesh& other) const
    {
        return tag_ < other.tag();
    }

    bool RegularMesh::operator==(const RegularMesh& r) const
    {
        return tag_ == r.tag();
    }

    void RegularMesh::clear_cells()
    {
        for (Bin& b: bin_)
        {
            b.clear_cells();
        }
    }

    void RegularMesh::calc_step_length()
    {
        double d0 = (aabb().max(0) - aabb().min(0)) / nstripe_(0);
        double d1 = (aabb().max(1) - aabb().min(1)) / nstripe_(1);
        double d2 = (aabb().max(2) - aabb().min(2)) / nstripe_(2);
        set_h(d0, d1, d2);
    }

    //void RegularMesh::show_bin_size() const
    //{
        //for (const Bin& b: bin_)
        //{
            //if (b.rm() == nullptr)
            //{
                ////std::cout << "bin tag = " << b.tag()() << "   " << b.cell().size() << std::endl;
                ////std::cout << "bin tag = " << b.tag()() << "   " << b.load_basic() << std::endl;
                //std::cout << "bin tag = " << b.tag()() << "   " << b.mesh_load()[0] << std::endl;
            //}
            //else
            //{
                //b.rm()->show_bin_size();
            //}
        //}
    //}
    Tag closest_non_empty_bin(const RegularMeshIndex& ind, const RegularMesh& rm, const Mesh& m)
    {
        auto non_empty = [&](int i, int j, int k)
        {
            if (i >= 0 && j >= 0 && k >= 0 && i < rm.nstripe(0) && j < rm.nstripe(1) && k < rm.nstripe(2))
            {
                RegularMeshIndex new_ind(i, j, k);
                if (!rm.bin(new_ind).cell().empty()) return true;
            }

            return false;
        };

        auto has_non_ghost = [](const Bin& b, const Mesh& m)
        {
            assert(!b.cell().empty());
            for (const BinCell& _bcc: b.cell())
            {
                //if (!m.cell(_bcc.cell()).is_ghost())
                {
                    assert(_bcc.cell().isvalid());
                    return _bcc.cell();
                }
            }

            return Tag(-1);
        };

        auto valid = [&](int i, int j, int k, Tag& t)
        {
            if (non_empty(i, j, k))
            {
                t = has_non_ghost(rm.bin(i, j, k), m);
                if (t.isvalid())
                {
                    return true;
                }
            }

            return false;
        };

        int orig_i = ind.i();
        int orig_j = ind.j();
        int orig_k = ind.k();

        int i_w = orig_i, i_e = orig_i, i_n = orig_i, i_s = orig_i, i_ne = orig_i, i_nw = orig_i, i_se = orig_i, i_sw = orig_i;
        int j_w = orig_j, j_e = orig_j, j_n = orig_j, j_s = orig_j, j_ne = orig_j, j_nw = orig_j, j_se = orig_j, j_sw = orig_j;
        int k_w = orig_k, k_e = orig_k, k_n = orig_k, k_s = orig_k, k_ne = orig_k, k_nw = orig_k, k_se = orig_k, k_sw = orig_k;

        int i_bw = orig_i, i_be = orig_i, i_bn = orig_i, i_bs = orig_i, i_bne = orig_i, i_bnw = orig_i, i_bse = orig_i, i_bsw = orig_i;
        int j_bw = orig_j, j_be = orig_j, j_bn = orig_j, j_bs = orig_j, j_bne = orig_j, j_bnw = orig_j, j_bse = orig_j, j_bsw = orig_j;
        int k_bw = orig_k, k_be = orig_k, k_bn = orig_k, k_bs = orig_k, k_bne = orig_k, k_bnw = orig_k, k_bse = orig_k, k_bsw = orig_k;

        int i_fw = orig_i, i_fe = orig_i, i_fn = orig_i, i_fs = orig_i, i_fne = orig_i, i_fnw = orig_i, i_fse = orig_i, i_fsw = orig_i;
        int j_fw = orig_j, j_fe = orig_j, j_fn = orig_j, j_fs = orig_j, j_fne = orig_j, j_fnw = orig_j, j_fse = orig_j, j_fsw = orig_j;
        int k_fw = orig_k, k_fe = orig_k, k_fn = orig_k, k_fs = orig_k, k_fne = orig_k, k_fnw = orig_k, k_fse = orig_k, k_fsw = orig_k;

        Tag t;
        
        while(true)
        {
            // west.
            j_w -= 1;
            if (valid(i_w, j_w, k_w, t)) return t;
            // south-west.
            i_sw -= 1;
            j_sw -= 1;
            if (valid(i_sw, j_sw, k_sw, t)) return t;
            // south.
            i_s -= 1;
            if (valid(i_s, j_s, k_s, t)) return t;
            // south-east.
            i_se -= 1;
            j_se += 1;
            if (valid(i_se, j_se, k_se, t)) return t;
            // east.
            j_e += 1;
            if (valid(i_e, j_e, k_e, t)) return t;
            // north-east.
            i_ne += 1;
            j_ne += 1;
            if (valid(i_ne, j_ne, k_ne, t)) return t;
            // north.
            i_n += 1;
            if (valid(i_n, j_n, k_n, t)) return t;
            // north-west.
            i_nw += 1;
            j_nw -= 1;
            if (valid(i_nw, j_nw, k_nw, t)) return t;


            // back-west.
            j_bw -= 1;
            k_bw -= 1;
            if (valid(i_bw, j_bw, k_bw, t)) return t;
            // back-south-west.
            i_bsw -= 1;
            j_bsw -= 1;
            k_bsw -= 1;
            if (valid(i_bsw, j_bsw, k_bsw, t)) return t;
            // back-south.
            i_bs -= 1;
            k_bs -= 1;
            if (valid(i_bs, j_bs, k_bs, t)) return t;
            // back-south-east.
            i_bse -= 1;
            j_bse += 1;
            k_bse += 1;
            if (valid(i_bse, j_bse, k_bse, t)) return t;
            // back-east.
            j_be += 1;
            k_be += 1;
            if (valid(i_be, j_be, k_be, t)) return t;
            // back-north-east.
            i_bne += 1;
            j_bne += 1;
            k_bne += 1;
            if (valid(i_bne, j_bne, k_bne, t)) return t;
            // back-north.
            i_bn += 1;
            k_bn += 1;
            if (valid(i_bn, j_bn, k_bn, t)) return t;
            // back-north-west.
            i_bnw += 1;
            j_bnw -= 1;
            k_bnw -= 1;
            if (valid(i_bnw, j_bnw, k_bnw, t)) return t;
            

            // front-west.
            j_fw -= 1;
            k_fw -= 1;
            if (valid(i_fw, j_fw, k_fw, t)) return t;
            // front-south-west.
            i_fsw -= 1;
            j_fsw -= 1;
            k_fsw -= 1;
            if (valid(i_fsw, j_fsw, k_fsw, t)) return t;
            // front-south.
            i_fs -= 1;
            k_fs -= 1;
            if (valid(i_fs, j_fs, k_fs, t)) return t;
            // front-south-east.
            i_fse -= 1;
            j_fse += 1;
            k_fse += 1;
            if (valid(i_fse, j_fse, k_fse, t)) return t;
            // front-east.
            j_fe += 1;
            k_fe += 1;
            if (valid(i_fe, j_fe, k_fe, t)) return t;
            // front-north-east.
            i_fne += 1;
            j_fne += 1;
            k_fne += 1;
            if (valid(i_fne, j_fne, k_fne, t)) return t;
            // front-north.
            i_fn += 1;
            k_fn += 1;
            if (valid(i_fn, j_fn, k_fn, t)) return t;
            // front-north-west.
            i_fnw += 1;
            j_fnw -= 1;
            k_fnw -= 1;
            if (valid(i_fnw, j_fnw, k_fnw, t)) return t;


            if (j_w <= 0 && i_sw <= 0 && j_sw <= 0 && i_s <= 0 && i_se <= 0 && j_se > rm.nstripe(1) && j_e > rm.nstripe(1) && i_ne > rm.nstripe(0) && j_ne > rm.nstripe(1) && i_n > rm.nstripe(0) && i_nw > rm.nstripe(0) && j_nw <= 0) return Tag(-1);
        }
    }

    void RegularMesh::max_bin_tag(int& s) const
    {
        for (const Bin& b: bin_)
        {
            if (b.rm() == nullptr)
            {
                s = std::max(s, b.tag()());
            }
            else
            {
                b.rm()->max_bin_tag(s);
            }
        }
    }

    void RegularMesh::size(size_t& s) const
    {
        for (const Bin& b: bin_)
        {
            if (b.rm() == nullptr)
            {
                ++s;
            }
            else
            {
                b.rm()->size(s);
            }
        }
    }

    int RegularMesh::max_bin_tag() const
    {
        int s = -1;
        max_bin_tag(s);
        return s;
    }

    size_t RegularMesh::size() const
    {
        size_t s = 0;
        size(s);
        return s;
    }

    BinRMTag RegularMesh::index_to_bin(size_t& s, size_t i) const
    {
        assert(!bin_.empty());

        for (const Bin& b: bin_)
        {
            if (b.rm() == nullptr)
            {
                if (s == i)
                {
                    return BinRMTag(b.tag(), tag());
                }
                ++s;
            }
            else
            {
                return b.rm()->index_to_bin(s, i);
            }
        }
    }

    BinRMTag RegularMesh::index_to_bin(size_t i) const
    {
        size_t s = 0;
        return index_to_bin(s, i);
    }

    bool RegularMesh::extend_aabb(const AABB& other_aabb)
    {
        return aabb_.extend(other_aabb);
    }

    RegularMesh& RegularMesh::operator=(const RegularMesh& other)
    {
        tag_ = other.tag_;
        nstripe_ = other.nstripe_;
        h_ = other.h_;
        aabb_ = other.aabb_;
        bin_ = other.bin_;
        global_min_ = other.global_min_;
        global_max_ = other.global_max_;
        rmtag_address_map_ = other.rmtag_address_map_;
        bintag_address_map_ = other.bintag_address_map_;

        if (tag_() == 0)
            update_address();

        return *this;
    }

    RegularMesh::RegularMesh(const RegularMesh& other):
        tag_(other.tag_),
        nstripe_(other.nstripe_),
        h_(other.h_),
        aabb_(other.aabb_),
        bin_(other.bin_),
        global_min_(other.global_min_),
        global_max_(other.global_max_),
        rmtag_address_map_(other.rmtag_address_map_),
        bintag_address_map_(other.bintag_address_map_)
    {
        if (tag_() == 0)
            update_address();
    }

    void RegularMesh::insert_bin_addresses(RegularMesh* rm)
    {
        assert(tag_() == 0);

        for (auto& b: rm->bin_) {
            assert(&b != nullptr);
            insert_to_bintag_address_map(b.tag(), &b);
        }
    }

    void RegularMesh::insert_to_bintag_address_map(const Tag& t, Bin* bin)
    {
        assert(bin != nullptr);
        bintag_address_map_.insert(std::make_pair(t, bin));
    }

    void RegularMesh::insert_to_rmtag_address_map(const Tag& t, RegularMesh* rmp)
    {
        rmtag_address_map_.insert(std::make_pair(t, rmp));
        //auto it = rmtag_address_map_.find(t);
        //assert(it != rmtag_address_map_.end());
    }

    RegularMesh::RegularMesh()
    {
        tag_ = Tag(0);
        //rmtag_address_map_.insert(std::make_pair(0, this));
        rmtag_address_map_.insert(std::make_pair(tag_, this));
    }

    RegularMesh::RegularMesh(const Tag& rmtag)
    {
        tag_ = rmtag;
        //rmtag_address_map_.insert(boost::bimap<Tag, RegularMesh*>::value_type(rmtag(), this));
    }

    void RegularMesh::set_tag(const Tag& t)
    {
        tag_ = t;
    }

    BinRMTag::BinRMTag(const Tag& bintag, const Tag& rmtag): bintag_(bintag), rmtag_(rmtag) {}

    bool BinRMTag::operator<(const BinRMTag& other) const
    {
        if (rmtag_ != other.rmtag())
        {
            return rmtag_ < other.rmtag();
        }

        return bintag_ < other.bintag();
    }

    bool BinRMTag::operator==(const BinRMTag& other) const
    {
        if (rmtag_ == other.rmtag())
        {
            if (bintag_ == other.bintag())
                return true;
        }

        return false;
    }

    const Tag& BinRMTag::bintag() const
    {
        return bintag_;
    }

    const Tag& BinRMTag::rmtag() const
    {
        return rmtag_;
    }

    //boost::shared_ptr<RegularMesh> BinTag::rm() const
    /*RegularMesh* BinTag::rm() const
    {
        return rm_;
    }*/

    void RegularMesh::get_bin_tags(std::vector<BinRMTag>& tag) const
    {
        for (const Bin& b: bin_)
        {
            if (b.rm() == nullptr)
                //tag.push_back(BinTag(b.tag(), std::shared_ptr<RegularMesh>(f())));
                tag.push_back(BinRMTag(b.tag(), tag_));
            else
                b.rm()->get_bin_tags(tag);
        }
    }

    //void RegularMesh::get_bin_loads(std::vector<int>& load) const
    //{
        //for (const Bin& b: bin_)
        //{
            //if (b.rm() == nullptr)
                //load.push_back(b.load());
            //else
                //b.rm()->get_bin_loads(load);
        //}
    //}

    int RegularMesh::get_new_rmtag() const
    {
        auto it = std::max_element(rmtag_address_map_.begin(), rmtag_address_map_.end(), [](const auto& b0, const auto& b1){return b0 < b1;});

        return (it->first() + 1);
    }

    int RegularMesh::get_new_bt() const
    {
        assert(tag_() == 0);
        int t = -1;
        for (const auto& _rm: rmtag_address_map_)
        {
            int _t = _rm.second->get_new_bt_p();
            t = std::max(t, _t);
        }

        return t;
    }

    int RegularMesh::get_new_bt_p() const
    {
        auto it = std::max_element(bin_.begin(), bin_.end(), [](const Bin& b0, const Bin& b1){return b0.tag()() < b1.tag()();});
        assert(it != bin_.end());

        return (it->tag()() + 1);
    }

    void RegularMesh::calc_h()
    {
        double h0 = (aabb_.max(0) - aabb_.min(0)) / nstripe_(0);
        double h1 = (aabb_.max(1) - aabb_.min(1)) / nstripe_(1);
        double h2 = (aabb_.max(2) - aabb_.min(2)) / nstripe_(2);

        set_h(h0, h1, h2);
    }

    void RegularMesh::clear_bincells()
    {
        for (Bin& b: bin_)
            b.clear_cells();
    }

    int get_binindex(unsigned int index, RelativePosition rp, const Vector3Int& nstripe, unsigned int binsize)
    {
        assert(rp != RelativePosition::undefined);
        assert(index >= 0);
        assert(index < binsize);

        auto index_have_west = [&]() {return (index % nstripe(0)) != 0;};
        auto index_have_east = [&]() {return ((index+1) % nstripe(0)) != 0;};
        auto index_have_north = [&]() {return index < (nstripe(0)*nstripe(1)-nstripe(0));};
        auto index_have_south = [&]() {return index >= nstripe(0);};
        auto index_have_front = [&]() {return (index > nstripe(0)*nstripe(1));};
        auto index_have_back = [&]() {return (index < nstripe(0)*nstripe(1)*nstripe(2));};

        switch (rp)
        {
            case RelativePosition::ccc:
                return index;
            case RelativePosition::cwc:
                if (index_have_west())
                    return index - 1;
            case RelativePosition::cec:
                if (index_have_east())
                    return index + 1;
            case RelativePosition::scc:
                if (index_have_south())
                    return index - nstripe(0);
            case RelativePosition::swc:
                if (index_have_south() && index_have_west())
                    return index - nstripe(0) - 1;
            case RelativePosition::sec:
                if (index_have_south() && index_have_east())
                    return index - nstripe(0) + 1;
            case RelativePosition::ncc:
                if (index_have_north())
                    return index + nstripe(0);
            case RelativePosition::nwc:
                if (index_have_north() && index_have_west())
                    return index + nstripe(0) - 1;
            case RelativePosition::nec:
                if (index_have_north() && index_have_east())
                    return index + nstripe(0) + 1;


            case RelativePosition::cca:
                if (index_have_back())
                    return index + nstripe(0)*nstripe(1);
            case RelativePosition::cwa:
                if (index_have_west() && index_have_back())
                    return index - 1 + nstripe(0)*nstripe(1);
            case RelativePosition::cea:
                if (index_have_east() && index_have_back())
                    return index + 1 + nstripe(0)*nstripe(1);
            case RelativePosition::sca:
                if (index_have_south() && index_have_back())
                    return index - nstripe(0) + nstripe(0)*nstripe(1);
            case RelativePosition::swa:
                if (index_have_south() && index_have_west() && index_have_back())
                    return index - nstripe(0) - 1 + nstripe(0)*nstripe(1);
            case RelativePosition::sea:
                if (index_have_south() && index_have_east() && index_have_back())
                    return index - nstripe(0) + 1 + nstripe(0)*nstripe(1);
            case RelativePosition::nca:
                if (index_have_north() && index_have_back())
                    return index + nstripe(0) + nstripe(0)*nstripe(1);
            case RelativePosition::nwa:
                if (index_have_north() && index_have_west() && index_have_back())
                    return index + nstripe(0) - 1 + nstripe(0)*nstripe(1);
            case RelativePosition::nea:
                if (index_have_north() && index_have_east() && index_have_back())
                    return index + nstripe(0) + 1 + nstripe(0)*nstripe(1);


            case RelativePosition::ccf:
                if (index_have_back())
                    return index - nstripe(0)*nstripe(1);
            case RelativePosition::cwf:
                if (index_have_west() && index_have_back())
                    return index - 1 - nstripe(0)*nstripe(1);
            case RelativePosition::cef:
                if (index_have_east() && index_have_back())
                    return index + 1 - nstripe(0)*nstripe(1);
            case RelativePosition::scf:
                if (index_have_south() && index_have_back())
                    return index - nstripe(0) - nstripe(0)*nstripe(1);
            case RelativePosition::swf:
                if (index_have_south() && index_have_west() && index_have_back())
                    return index - nstripe(0) - 1 - nstripe(0)*nstripe(1);
            case RelativePosition::sef:
                if (index_have_south() && index_have_east() && index_have_back())
                    return index - nstripe(0) + 1 - nstripe(0)*nstripe(1);
            case RelativePosition::ncf:
                if (index_have_north() && index_have_back())
                    return index + nstripe(0) - nstripe(0)*nstripe(1);
            case RelativePosition::nwf:
                if (index_have_north() && index_have_west() && index_have_back())
                    return index + nstripe(0) - 1 - nstripe(0)*nstripe(1);
            case RelativePosition::nef:
                if (index_have_north() && index_have_east() && index_have_back())
                    return index + nstripe(0) + 1 - nstripe(0)*nstripe(1);
        }

        return -1;
    }

    std::vector<Tag> RegularMesh::pneis_of_bin(unsigned int index)
    {
        assert(index >= 0);
        assert(index < bin_.size());

        auto index_have_west = [&]() {return (index % nstripe(0)) != 0;};
        auto index_have_east = [&]() {return ((index+1) % nstripe(0)) != 0;};
        auto index_have_north = [&]() {return index < (nstripe(0)*nstripe(1)-nstripe(0));};
        auto index_have_south = [&]() {return index >= nstripe(0);};
        auto index_have_front = [&]() {return (index >= nstripe(0)*nstripe(1));};
        auto index_have_back = [&]() {return (index < nstripe(0)*nstripe(1));};

        std::vector<Tag> nei(6);

        if (index_have_west())
        {
            nei[0] = Tag(index - 1);
            assert(nei[0]() < bin_.size());
        }
        if (index_have_east())
        {
            nei[1] = Tag(index + 1);
            assert(nei[1]() < bin_.size());
        }
        if (index_have_south())
        {
            nei[2] = Tag(index - nstripe_(0));
            assert(nei[2]() < bin_.size());
        }
        if (index_have_north())
        {
            nei[3] = Tag(index + nstripe_(0));
            assert(nei[3]() < bin_.size());
        }
        if (index_have_back())
        {
            nei[4] = Tag(index + nstripe_(0)*nstripe(1));
            if (nei[4]() >= bin_.size())
            {
                std::cout << "nstripe(0): " << nstripe(0) << std::endl;
                std::cout << "nstripe(1): " << nstripe(1) << std::endl;
                std::cout << "nstripe(2): " << nstripe(2) << std::endl;
                std::cout << "index: " << index << std::endl;
            }
            assert(nei[4]() < bin_.size());
        }
        if (index_have_front())
        {
            nei[5] = Tag(index - nstripe_(0)*nstripe(1));
            assert(nei[5]() < bin_.size());
        }

        return nei;
    }

    //int RegularMesh::get_binindex(unsigned int index, RelativePosition rp)
    //{
        //return get_binindex(index, rp, nstripe_, bin_.size())
    //}

    /*bool RegularMesh::query_cell(const Tag& im, const Tag& ic) const
    {
        for (const Bin& _bin: bin_)
        {
            bool exist = _bin.query_bincell(im, ic);
            if (exist)
                return true;
        }

        return false;
    }*/

    RegularMeshIndex RegularMesh::global_index(const RegularMeshIndex& i) const
    {
        return RegularMeshIndex(i + global_min_);
    }

    /*Bin& RegularMesh::bin_p(const Tag& t)
    {
        assert(t.isvalid());
        auto it = std::lower_bound(bin_.begin(), bin_.end(), Bin(t));
        assert(it != bin_.end() && it->tag() == t);
        return *it;
    }*/

    Bin& RegularMesh::bin_p(const Tag& t)
    {
        auto it = bintag_address_map_.find(t);
        if (it == bintag_address_map_.end())
        {
            std::cout << "t: " << t() << std::endl;
            std::cout << "bintag_add size: " << bintag_address_map_.size() << std::endl;
            std::cout << "bin size: " << bin_.size() << std::endl;
        }
        assert(it != bintag_address_map_.end());
        assert(it->second != nullptr);
        return *(it->second);
    }

    const Bin& RegularMesh::bin(const Tag& t) const
    {
        auto it = bintag_address_map_.find(t);
        if (it == bintag_address_map_.end())
        {
            std::cout << "t: " << t() << std::endl;
            std::cout << "size: " << bintag_address_map_.size() << std::endl;
        }
        assert(it != bintag_address_map_.end());
        assert(it->second != nullptr);
        return *(it->second);
    }

    /*const Bin& RegularMesh::bin_local(const Tag& t) const
    {
        assert(t.isvalid());
        auto it = std::lower_bound(bin_.begin(), bin_.end(), Bin(t));
        if (it == bin_.end())
        {
            std::cout << "t: " << t() << std::endl;
            std::cout << "bin size: " << bin_.size() << std::endl;
            for (const Bin& b: bin_)
            {
                std::cout << "bin tag: " << b.tag()() << std::endl;
            }
        }
        assert(it != bin_.end());
        if (it->tag() != t)
        {
            std::cout << "t: " << t() << std::endl;
            std::cout << "bin size: " << bin_.size() << std::endl;
        }
        assert(it->tag() == t);
        return *it;
    }*/

    /*const Bin& RegularMesh::bin(const BinRMTag& bt) const
    {
        //assert(tag_() == 0);
        auto it = rmtag_address_map_.find(bt.rmtag());
        if (it == rmtag_address_map_.end())
        {
            std::cout << "rm map size = " << rmtag_address_map_.size() << std::endl;
            std::cout << "bin size = " << bin_.size() << std::endl;
            std::cout << "bintag = " << bt.bintag()() << std::endl;
            std::cout << "rmtag = " << bt.rmtag()() << std::endl;
            std::cout << "this rm tag = " << tag_() << std::endl;
            for (const auto& a: rmtag_address_map_)
            {
                std::cout << "inside = " << a.first() << std::endl;
            }
            
        }
        assert(it != rmtag_address_map_.end());
        assert(it->second != nullptr);
        return it->second->bin(bt.bintag());
    }*/

    /*Bin& RegularMesh::bin_p(const BinRMTag& bt)
    {
        assert(tag_() == 0);
        auto it = rmtag_address_map_.find(bt.rmtag());
        assert(it != rmtag_address_map_.end());
        return it->second->bin_p(bt.bintag());
    }*/

    const Tag& RegularMesh::tag() const
    {
        return tag_;
    }

    const std::map<Tag, RegularMesh*>& RegularMesh::rmtag_address_map() const
    {
        return rmtag_address_map_;
    }

    const std::map<Tag, Bin*>& RegularMesh::bintag_address_map() const
    {
        return bintag_address_map_;
    }

    void RegularMesh::set_rmtag_address_map(const std::map<Tag, RegularMesh*>& map)
    {
        rmtag_address_map_ = map;
    }

    void RegularMesh::update_address()
    {
        //assert(tag_() == 0);
        update_address(rmtag_address_map_, bintag_address_map_);
    }

    void RegularMesh::update_address(std::map<Tag, RegularMesh*>& map, std::map<Tag, Bin*>& binmap)
    {
        if (map.empty()) {
            return;
        }

        auto it = map.find(tag_);
        if (it == map.end())
        {
            std::cout << tag_() << std::endl;
            std::cout << map.size() << std::endl;
            for (const auto& ii: map)
            {
                std::cout << ii.first() << " " << ii.second << std::endl;
            }
        }
        assert(it != map.end());
        it->second = this;
        for (auto& bin: bin_)
        {
            auto bit = binmap.find(bin.tag());
            bit->second = &bin;
        }
        //bool successful_replace = map.replace_data(it, this);
        //assert(successful_replace);

        for (const Bin& b: bin_)
        {
            if (b.rm() != nullptr)
                b.rm()->update_address(map, binmap);
        }
    }

    void RegularMesh::merge(const RegularMesh& other)
    {
        for (const Bin& _bin: other.bin())
        {
            Bin& mybin = bin_p(_bin.tag());

            mybin.merge(_bin);
        }
    }

    bool RegularMesh::validate(const RegularMeshIndex& global_ind, RegularMeshIndex& local_ind) const
    {
        if (global_ind.i() >= global_min_.i() && global_ind.j() >= global_min_.j() && global_ind.k() >= global_min_.k())
        {
            if (global_ind.i() <= global_max_.i() && global_ind.j() <= global_max_.j() && global_ind.k() <= global_max_.k())  
            {
                int i0 = global_ind.i() - global_min_.i();
                int i1 = global_ind.j() - global_min_.j();
                int i2 = global_ind.k() - global_min_.k();
                local_ind.set(i0, i1, i2);
                return true;
            }
        }

        return false;
    }

    void RegularMesh::set_global_min(const RegularMeshIndex& i)
    {
        global_min_ = i;
    }
    void RegularMesh::set_global_min(int i, int j, int k)
    {
        global_min_.set(i, j, k);
    }
    void RegularMesh::set_global_max(const RegularMeshIndex& i)
    {
        global_max_ = i;
    }
    void RegularMesh::set_global_max(int i, int j, int k)
    {
        global_max_.set(i, j, k);
    }

    const RegularMeshIndex& RegularMesh::global_min() const
    {
        return global_min_;
    }

    const RegularMeshIndex& RegularMesh::global_max() const
    {
        return global_max_;
    }

    bool RegularMesh::is_resident(const Vector3& cnt) const
    {
        /*std::cout << "sssssssssssssss" << std::endl;
          std::cout << "cnt(0) = " << cnt(0) << std::endl;
          std::cout << "cnt(1) = " << cnt(1) << std::endl;
          std::cout << "aabb.min(0) = " << aabb_.min(0) << std::endl;
          std::cout << "aabb.min(1) = " << aabb_.min(1) << std::endl;
          std::cout << "aabb.max(0) = " << aabb_.max(0) << std::endl;
          std::cout << "aabb.max(1) = " << aabb_.max(1) << std::endl;*/

        return aabb_.do_intersect(cnt);

        /*for (int i=0; i<TAILOR_N_DIM; ++i)
        {
            if (cnt(i) < aabb_.min(i))
                return false;

            if (cnt(i) > aabb_.max(i))
            {
                //std::cout << "invalidddd 2" << std::endl;
                return false;
            }
        }

        return true;*/
    }

    /*void RegularMesh::set_binload(const RegularMeshIndex& ind, double l)
    {
        bin_p(ind).set_load(l);
    }*/

    const std::vector<Bin>& RegularMesh::bin() const
    {
        return bin_;
    }

    const Bin& RegularMesh::bin(int row, int col, int depth) const
    {
        assert(!bin_.empty());
        if (row < 0 || col < 0)
        {
            std::cout << "row: " << row << std::endl;
            std::cout << "col: " << col << std::endl;
        }
        assert(row >= 0);
        assert(col >= 0);
        assert(row < nstripe(1));
        assert(col < nstripe(0));
        assert(depth < nstripe(2));

        assert(col + row*nstripe(0) + depth*nstripe(0)*nstripe(1) >= 0);

        if (col + row*nstripe(0) + depth*nstripe(0)*nstripe(1) >= bin_.size())
        {
            std::cout << "row: " << row << std::endl;
            std::cout << "col: " << col << std::endl;
            std::cout << "depth: " << depth  << std::endl;
            std::cout << "bin size: " << bin_.size()  << std::endl;
            std::cout << "ind: " << col + row*nstripe(0) + depth*nstripe(0)*nstripe(1) << std::endl;
            std::cout << "nstripe 0: " << nstripe(0) << std::endl;
            std::cout << "nstripe 1: " << nstripe(1) << std::endl;
            std::cout << "nstripe 2: " << nstripe(2) << std::endl;
        }

        assert(col + row*nstripe(0) + depth*nstripe(0)*nstripe(1) < bin_.size());

        return bin_[col + row*nstripe(0) + depth*nstripe(0)*nstripe(1)];
    }

    const Bin& RegularMesh::bin(const RegularMeshIndex& ind) const
    {
        return bin(ind.i(), ind.j(), ind.k());
    }

    Bin& RegularMesh::bin_p(int row, int col, int depth)
    {
        assert(row >= 0);
        assert(col >= 0);
        if (row >= nstripe(1))
        {
            std::cout << "row = " << row << std::endl;
            std::cout << "col = " << col << std::endl;
            std::cout << "nstripe(1) = " << nstripe(1) << std::endl;
        }
        assert(row < nstripe(1));
        assert(col < nstripe(0));
        assert(col + row*nstripe(0) < bin_.size());

        return bin_[col + row*nstripe(0) + depth*nstripe(0)*nstripe(1)];
    }

    Bin& RegularMesh::bin_p(const RegularMeshIndex& ind)
    {
        return bin_p(ind.i(), ind.j(), ind.k());
    }

    //void RegularMesh::set_bincell(const RegularMeshIndex& ind, const std::vector<BinCell>& bc)
    void RegularMesh::set_bincell(const RegularMeshIndex& ind, const std::deque<BinCell>& bc)
    {
        bin_p(ind).set_cell(bc);
    }

    const Bin& RegularMesh::bin(int i) const
    {
        assert(i >= 0);
        assert(i < bin_.size());
        return bin_[i];
    }

    void RegularMesh::set_aabb(const AABB& aabb)
    {
        aabb_ = aabb;
    }

    //void RegularMesh::set_aabb_min(const Vector3& m)
    //{
        //aabb_.set_min(m);
    //}
    //void RegularMesh::set_aabb_max(const Vector3& m)
    //{
        //aabb_.set_max(m);
    //}

    void RegularMesh::set_nstripe(int x, int y, int z)
    {
        //if(x <= 1 || y <= 1 || z <= 1)
        //{
            //std::cout << "x = " << x << std::endl;
            //std::cout << "y = " << y << std::endl;
            //std::cout << "z = " << z << std::endl;
        //}
        //assert(x > 1);
        //assert(y > 1);
        //assert(z > 1);
        nstripe_.set(x, y, z);
    }

    void RegularMesh::set_nstripe(const Vector3Int& ns)
    {
        set_nstripe(ns(0), ns(1), ns(2));
    }

    void RegularMesh::set_h(const Vector3& v)
    {
        set_h(v(0), v(1), v(2));
    }

    void RegularMesh::set_h(double v0, double v1, double v2)
    {
        h_.set(v0, v1, v2);
    }

    double RegularMesh::h(int i) const
    {
        return h_(i);
    }

    const AABB& RegularMesh::aabb() const
    {
        return aabb_;
    }

    const Vector3& RegularMesh::h() const
    {
        return h_;
    }

    const Vector3Int& RegularMesh::nstripe() const
    {
        return nstripe_;
    }

    int RegularMesh::nstripe(int i) const
    {
        return nstripe_(i); 
    }

    bool RegularMeshIndex::sorter(const RegularMeshIndex& lhs, const RegularMeshIndex& rhs)
    {
        if (lhs.i() != rhs.i())
        {
            return lhs.j() < rhs.j();
        }
        return lhs.i() < rhs.i();
    }
    bool RegularMeshIndex::operator==(const RegularMeshIndex& other) const
    {
        if (i() == other.i() && j() == other.j())
        {
            return true;
        }

        return false;
    }
    /*bool RegularMeshIndex::operator!=(const RegularMeshIndex& other) const
      {
      return !(*this == other);
      }*/

    RelativePosition get_relative_position(int index, int targetindex, const Vector3Int& nstripe, size_t maxbinsize, bool& exact)
    {
        // Determines relative position of target index (ti) wrt index (in).
        /* |    |    |    |
           ----------------
           |    | in |    |
           ----------------
           |    |    | ti |
           */

        // Relative positions.
        /* | nw | nc | ne |
           ----------------
           | cw | cc | ce |
           ----------------
           | sw | sc | se |
           */

        // Example: in=2, ti=1, ns=2, max=4 --> se.
        /* | 2 | 3 |
           ---------
           | 0 | 1 |
           */

        exact = true;

        auto index_have_west = [&]() {return (index % nstripe(0)) != 0;};
        auto index_have_east = [&]() {return ((index+1) % nstripe(0)) != 0;};
        auto index_have_north = [&]() {return index < (nstripe(0)*nstripe(1)-nstripe(0));};
        auto index_have_south = [&]() {return index >= nstripe(0);};
        auto index_have_front = [&]() {return (index > nstripe(0)*nstripe(1));};
        auto index_have_back = [&]() {return (index < nstripe(0)*nstripe(1)*nstripe(2));};

        assert(index >= 0);
        assert(index < maxbinsize);

        if (index == targetindex) {return RelativePosition::ccc;}
        if (index == targetindex+1 && index_have_west()) {return RelativePosition::cwc;}
        if (index == targetindex-1 && index_have_east()) {return RelativePosition::cec;}
        if (index == targetindex-nstripe(0) && index_have_north()) {return RelativePosition::ncc;}
        if (index == targetindex-nstripe(0)+1 && index_have_west() && index_have_north()) {return RelativePosition::nwc;}
        if (index == targetindex-nstripe(0)-1 && index_have_east() && index_have_north()) {return RelativePosition::nec;}
        if (index == targetindex+nstripe(0) && index_have_south()) {return RelativePosition::scc;}
        if (index == targetindex+nstripe(0)+1 && index_have_west() && index_have_south()) {return RelativePosition::swc;}
        if (index == targetindex+nstripe(0)-1 && index_have_east() && index_have_south()) {return RelativePosition::sec;}

        if (index == targetindex+nstripe(0)*nstripe(1) && index_have_front()) {return RelativePosition::ccf;}
        if (index == targetindex+1+nstripe(0)*nstripe(1) && index_have_west() && index_have_front()) {return RelativePosition::cwf;}
        if (index == targetindex-1+nstripe(0)*nstripe(1) && index_have_east() && index_have_front()) {return RelativePosition::cef;}
        if (index == targetindex-nstripe(0)+nstripe(0)*nstripe(1) && index_have_north() && index_have_front()) {return RelativePosition::ncf;}
        if (index == targetindex-nstripe(0)+1+nstripe(0)*nstripe(1) && index_have_west() && index_have_north() && index_have_front()) {return RelativePosition::nwf;}
        if (index == targetindex-nstripe(0)-1+nstripe(0)*nstripe(1) && index_have_east() && index_have_north() && index_have_front()) {return RelativePosition::nef;}
        if (index == targetindex+nstripe(0)+nstripe(0)*nstripe(1) && index_have_south() && index_have_front()) {return RelativePosition::scf;}
        if (index == targetindex+nstripe(0)+1+nstripe(0)*nstripe(1) && index_have_west() && index_have_south() && index_have_front()) {return RelativePosition::swf;}
        if (index == targetindex+nstripe(0)-1+nstripe(0)*nstripe(1) && index_have_east() && index_have_south() && index_have_front()) {return RelativePosition::sef;}

        if (index == targetindex-nstripe(0)*nstripe(1) && index_have_back()) {return RelativePosition::cca;}
        if (index == targetindex+1-nstripe(0)*nstripe(1) && index_have_west() && index_have_back()) {return RelativePosition::cwa;}
        if (index == targetindex-1-nstripe(0)*nstripe(1) && index_have_east() && index_have_back()) {return RelativePosition::cea;}
        if (index == targetindex-nstripe(0)-nstripe(0)*nstripe(1) && index_have_north() && index_have_back()) {return RelativePosition::nca;}
        if (index == targetindex-nstripe(0)+1-nstripe(0)*nstripe(1) && index_have_west() && index_have_north() && index_have_back()) {return RelativePosition::nwa;}
        if (index == targetindex-nstripe(0)-1-nstripe(0)*nstripe(1) && index_have_east() && index_have_north() && index_have_back()) {return RelativePosition::nea;}
        if (index == targetindex+nstripe(0)-nstripe(0)*nstripe(1) && index_have_south() && index_have_back()) {return RelativePosition::sca;}
        if (index == targetindex+nstripe(0)+1-nstripe(0)*nstripe(1) && index_have_west() && index_have_south() && index_have_back()) {return RelativePosition::swa;}
        if (index == targetindex+nstripe(0)-1-nstripe(0)*nstripe(1) && index_have_east() && index_have_south() && index_have_back()) {return RelativePosition::sea;}

        exact = false;

        /* | 12 | 13 | 14 | 15 |
           ---------------------
           | 8  | 9  | in | 11 |
           ---------------------
           | 4  | tr | 6  | 7  |
           ---------------------
           | ti | 1  | 2  | 3  |
           */

       /* unsigned int inrow = index % nstripe(0);
        unsigned int tirow = targetindex % nstripe(0);

        if (tirow < inrow)
        {
            if (index == targetindex + nstripe(0))
                return RelativePosition::sc;
            if (index > targetindex + nstripe(0))
                return RelativePosition::sw;
            if (index < targetindex + nstripe(0))
                return RelativePosition::se;
        }
        else if (tirow > inrow)
        {
            if (index == targetindex + nstripe(0))
                return RelativePosition::nc;
            if (index > targetindex + nstripe(0))
                return RelativePosition::nw;
            if (index < targetindex + nstripe(0))
                return RelativePosition::ne;
        }
        else
        {
            assert(index != targetindex + nstripe(0));

            if (index > targetindex + nstripe(0))
                return RelativePosition::cw;
            if (index < targetindex + nstripe(0))
                return RelativePosition::ce;
        }*/

        return RelativePosition::undefined;
    }

    bool in_range(const RegularMeshIndex& index, const RegularMeshIndex& global_min, const RegularMeshIndex& global_max)
    {
        if (index.i() < global_min.i())
        {
            return true;
        }

        if (index.j() < global_min.j())
        {
            return true;
        }

        if (index.i() > global_max.i())
        {
            return true;
        }

        if (index.j() > global_max.j())
        {
            return true;
        }

        return false;
    }

    unsigned int mat_to_vec_index(const RegularMeshIndex& index, const Vector3Int& nstripe)
    {
        return index.j() + index.i() * nstripe(0);
    }

    Vector3 RegularMesh::centroid(const RegularMeshIndex& index) const
    {
        assert(index.i() < 0);
        assert(index.j() < 0);
        assert(index.i() >= nstripe_(1));
        assert(index.j() >= nstripe_(0));

        Vector3 c;

        double d0 = aabb().min(0) + index.j() * h(0) + h(0)/2.;
        double d1 = aabb().min(1) + index.i() * h(1) + h(1)/2.;
        double d2 = aabb().min(2) + index.k() * h(2) + h(2)/2.;
        c.set(d0, d1, d2);

        return c;
    }

    Vector3 llcoor(const RegularMeshIndex& index, const Vector3& aabb_min, const Vector3& h)
    {
        double d0 = aabb_min(0) + index.j() * h(0);
        double d1 = aabb_min(1) + index.i() * h(1);
        double d2 = aabb_min(2) + index.k() * h(2);
        return Vector3(d0, d1, d2);
    }

    Vector3 RegularMesh::llcoor(const RegularMeshIndex& index) const
    {
        double d0 = aabb().min(0) + index.j() * h(0);
        double d1 = aabb().min(1) + index.i() * h(1);
        double d2 = aabb().min(2) + index.k() * h(2);
        return Vector3(d0, d1, d2);
    }

    void RegularMesh::set_props(const Mesh& mesh, const AABB& forced_min_aabb, int _nstripe)
    {
        /*// set AABB.
          AABB _aabb(mesh.point());
          set_aabb_min(_aabb.min());
          set_aabb_max(_aabb.max());

          if (aabb_.min(0) > forced_min_aabb.min(0))
          {
          set_aabb_min(vec3(forced_min_aabb.min(0), aabb_.min(1)));
          }
          if (aabb_.min(1) > forced_min_aabb.min(1))
          {
          set_aabb_min(vec3(aabb_.min(0), forced_min_aabb.min(1)));
          }
          if (aabb_.max(0) < forced_min_aabb.max(0))
          {
          set_aabb_max(vec3(forced_min_aabb.max(0), aabb_.max(1)));
          }
          if (aabb_.max(1) < forced_min_aabb.max(1))
          {
          set_aabb_max(vec3(aabb_.max(0), forced_min_aabb.max(1)));
          }*/

        set_aabb(AABB(forced_min_aabb.min(), forced_min_aabb.max()));

        // set initial nstripe.
        if (_nstripe == 0)
        {
            int _size = mesh.cell().size();
            assert(_size != 0);
            //int ns = std::floor(std::sqrt(_size));
            int ns = std::floor(std::pow(_size, 1./3.));
            //set_nstripe(std::floor(std::sqrt(_size)));
            if (ns < 2) {
                ns = 2;
            }
            //set_nstripe(ns, ns, ns);
            //set_nstripe(ns, ns, 1);
            set_nstripe(1, 1, 1);
        }
        else
        {
            set_nstripe(_nstripe, _nstripe, _nstripe);
        }

        // calculate step length.
        double d0 = (aabb().max(0) - aabb().min(0)) / nstripe(0);
        double d1 = (aabb().max(1) - aabb().min(1)) / nstripe(1);
        double d2 = (aabb().max(2) - aabb().min(2)) / nstripe(2);
        set_h(d0, d1, d2);
    }

    /*void RegularMesh::set_props(const MeshBlock& meshblock, unsigned int nstripe_)
      {
    // set AABB.
    AABB aabb_(meshblock);
    aabb.min = aabb_.min;
    aabb.max = aabb_.max;

    // set initial nstripe.
    if (nstripe_ == 0)
    {
    int size = meshblock.cell.size();
    nstripe = std::floor(std::sqrt(size));
    }
    else
    {
    nstripe = nstripe_;
    }

    // calculate step length.
    h[0] = (aabb.max[0] - aabb.min[0]) / nstripe;
    h[1] = (aabb.max[1] - aabb.min[1]) / nstripe;
    }

    void RegularMesh::set_props(const std::vector<MeshBlock>& meshblock, unsigned int nstripe_)
    {
    // set AABB.
    vec3 min(BIG_POS_NUM, BIG_POS_NUM);
    vec3 max(BIG_NEG_NUM, BIG_NEG_NUM);
    for (const MeshBlock& mb: meshblock)
    {
    AABB aabb_(mb);
    min[0] = std::min(min[0], aabb_.min[0]);
    min[1] = std::min(min[1], aabb_.min[1]);
    max[0] = std::max(max[0], aabb_.max[0]);
    max[1] = std::max(max[1], aabb_.max[1]);
    }
    aabb.min = min;
    aabb.max = max;

    // set initial nstripe.
    if (nstripe_ == 0)
    {
    std::size_t size = 0;
    for (const MeshBlock& mb: meshblock)
    {
    size = std::max(size, mb.cell.size());
    }
    nstripe = std::floor(std::sqrt(size));
    }
    else
    {
    nstripe = nstripe_;
    }

    // calculate step length.
    h[0] = (aabb.max[0] - aabb.min[0]) / nstripe;
    h[1] = (aabb.max[1] - aabb.min[1]) / nstripe;
    }*/


    void RegularMesh::insert_bins(int mintag, int mesh_load_size, int rank)
    {
        bin_.clear();
        // insert bins with given indices.
        bin_.reserve(nstripe(0)*nstripe(1)*nstripe(2));
        RegularMeshIndex _index(0, 0, 0);
        int current_tag = mintag;
        for (int k=0; k<nstripe(2); ++k)
        {
            for (int i=0; i<nstripe(1); ++i)
            {
                for (int j=0; j<nstripe(0); ++j)
                {
                    _index.set(i, j, k);
                    Vector3 ll = llcoor(_index);
                    Vector3 ur = ll + h();

                    assert(tag_.isvalid());

                    assert(current_tag != -1);
                    Bin temp_bin(Tag(current_tag), tag_, _index, AABB(ll, ur), mesh_load_size);
                    auto it = std::lower_bound(bin_.begin(), bin_.end(), temp_bin);
                    if (it == bin_.end())
                    {
                        bin_.push_back(temp_bin);
                    }
                    else
                    {
                        assert(it->tag()() != current_tag);
                        bin_.insert(it, temp_bin);
                    }
                    //bin_.push_back(temp_bin);

                    /*assert(bin_.back().rm() == nullptr);
                      if (bin_.back().aabb().min(0) == bin_.back().aabb().max(0))
                      {
                      std::cout << "h(0) = " << h(0) << std::endl;
                      std::cout << "h(1) = " << h(1) << std::endl;
                      }
                      assert(bin_.back().aabb().min(0) != bin_.back().aabb().max(0));
                      assert(h(0) != 0);
                      assert(h(1) != 0);
                      assert(bin_.back().aabb().min(0) != bin_.back().aabb().max(0));*/

                    //bintag_index_map_.insert(boost::bimap<Tag, int>::value_type(current_tag, bin_.size() - 1));
                    ++current_tag;
                }
            }
    }
        if (bin_.size() != nstripe(0)*nstripe(1)*nstripe(2))
        {
            std::cout << "bin_.size() = " << bin_.size() << std::endl;
            std::cout << "nstripe(0) = " << nstripe(0) << std::endl;
            std::cout << "nstripe(1) = " << nstripe(1) << std::endl;
        }
        assert(bin_.size() == nstripe(0)*nstripe(1)*nstripe(2));
        /*if(bintag_index_map_.left.empty())
        {
            std::cout << "bin_.size() = " << bin_.size() << std::endl;
            std::cout << "nstripe(0) = " << nstripe(0) << std::endl;
            std::cout << "nstripe(1) = " << nstripe(1) << std::endl;
        }
        assert(!bintag_index_map_.left.empty());*/
        for (const Bin& b: bin_)
        {
            assert(b.cell().empty());
            assert(b.mesh_tag_index_map_res().empty());
            if (b.aabb().min(0) == b.aabb().max(0))
            {
                std::cout << "bin_.size() = " << bin_.size() << std::endl;
                std::cout << "h(0) = " << h(0) << std::endl;
                std::cout << "h(1) = " << h(1) << std::endl;
                std::cout << "h(2) = " << h(2) << std::endl;
            }
            assert(b.aabb().min(0) != b.aabb().max(0));
        }
    }

    void RegularMesh::calc_aabb_max()
    {
        double max_0 = aabb_.min(0) + nstripe(0) * h(0);
        double max_1 = aabb_.min(1) + nstripe(1) * h(1);
        double max_2 = aabb_.min(2) + nstripe(2) * h(2);
        Vector3 max_aabb(max_0, max_1, max_2);
        //set_aabb_max(max_aabb);
        set_aabb(AABB(aabb_.min(), max_aabb));
    }

    /*bool get_index_regular(const Vector3& cv, std::vector<RegularMeshIndex>& index, const Vector3& aabb_min, const Vector3& h, const vec3<int>& nstripe, bool closest)
    {
        RegularMesh rm;
        rm.set_h(h);
        rm.set_nstripe(nstripe);
        //rm.set_aabb_min(aabb_min);
        double max_0 = aabb_min(0) + nstripe(0) * h(0);
        double max_1 = aabb_min(1) + nstripe(1) * h(1);
        double max_2 = aabb_min(2) + nstripe(2) * h(2);
        Vector3 max_aabb(max_0, max_1, max_2);
        //rm.set_aabb_max(max_aabb);
        rm.set_aabb(AABB(aabb_min, max_aabb));

        return rm.get_index_regular(cv, index, closest);
    }*/

    RegularMeshIndex RegularMesh::get_index_regular_unique(const Vector3& cv) const
    {
        RegularMeshIndex _index;

        for (int i=0; i<TAILOR_N_DIM; ++i)
        {
            if (TAILOR_ZERO < aabb_.min(i) - cv(i))
            {
                std::cout << "aabb.min(0): " << aabb_.min(0) << std::endl;
                std::cout << "aabb.min(1): " << aabb_.min(1) << std::endl;
                std::cout << "aabb.min(2): " << aabb_.min(2) << std::endl;

                std::cout << "cv(0): " << cv(0) << std::endl;
                std::cout << "cv(1): " << cv(1) << std::endl;
                std::cout << "cv(2): " << cv(2) << std::endl;

                assert(false);

                return _index;
            }



            if (cv(i) - aabb_.max(i) > TAILOR_ZERO)
            {
                std::cout << "aabb.max(0): " << aabb_.max(0) << std::endl;
                std::cout << "aabb.max(1): " << aabb_.max(1) << std::endl;
                std::cout << "aabb.max(2): " << aabb_.max(2) << std::endl;

                std::cout << "cv(0): " << cv(0) << std::endl;
                std::cout << "cv(1): " << cv(1) << std::endl;
                std::cout << "cv(2): " << cv(2) << std::endl;

                assert(false);

                return _index;
            }
        }
        
        double _d0 = (cv(0) - aabb_.min(0)) / h(0);
        double _d1 = (cv(1) - aabb_.min(1)) / h(1);
        double _d2 = (cv(2) - aabb_.min(2)) / h(2);

        //std::cout << "_d0: " << _d0 << std::endl;
        //std::cout << "_d1: " << _d1 << std::endl;
        //std::cout << "_d2: " << _d2 << std::endl;

        /*if (std::abs(_d0 - std::round(_d0)) < TAILOR_ZERO) {
            _d0 = std::round(_d0);
        }
        if (std::abs(_d1 - std::round(_d1)) < TAILOR_ZERO) {
            _d1 = std::round(_d1);
        }
        if (std::abs(_d2 - std::round(_d2)) < TAILOR_ZERO) {
            _d2 = std::round(_d2);
        }*/

        //assert(_d0 >= 0.);
        //assert(_d1 >= 0.);
        //assert(_d2 >= 0.);

        if (_d2 < 0) {
            _d2 = 0;
        }
        if (_d1 < 0) {
            _d1 = 0;
        }
        if (_d0 < 0) {
            _d0 = 0;
        }

        _index.set(std::floor(_d1), std::floor(_d0), std::floor(_d2));

        if (_index.j() == nstripe_(0)) {
            _index.decj();
        }
        if (_index.i() == nstripe_(1)) {
            _index.deci();
        }
        if (_index.k() == nstripe_(2)) {
            _index.deck();
        }

        assert(_index.isvalid());

        for (int i=0; i<TAILOR_N_DIM; ++i)
        {
            assert(TAILOR_ZERO >= bin(_index).aabb().min(i) - cv(i));
            assert (cv(i) - aabb_.max(i) <= TAILOR_ZERO);
        }

        return _index;
    }

    bool RegularMesh::get_index_regular(const Vector3& cv, std::vector<RegularMeshIndex>& index, bool closest, int rank, int celltag) const
    {
        if (closest)
        {
            Vector3 newcv(cv);
            // cv might be out of bounds of rm. if closest is true return closest index to cv.
            if (newcv(0) < aabb_.min(0))
                newcv.set(aabb_.min(0), newcv(1), newcv(2));
            if (newcv(1) < aabb_.min(1))
                newcv.set(newcv(0), aabb_.min(1), newcv(2));
            if (newcv(2) < aabb_.min(2))
                newcv.set(newcv(0), newcv(1), aabb_.min(2));

            if (newcv(0) > aabb_.max(0))
                newcv.set(aabb_.max(0), newcv(1), newcv(2));
            if (newcv(1) > aabb_.max(1))
                newcv.set(newcv(0), aabb_.max(1), newcv(2));
            if (newcv(2) > aabb_.max(2))
                newcv.set(newcv(0), newcv(1), aabb_.max(2));

            assert(newcv(0) >= aabb_.min(0));
            assert(newcv(1) >= aabb_.min(1));
            assert(newcv(0) <= aabb_.max(0));
            assert(newcv(1) <= aabb_.max(1));

            return get_index_regular(newcv, index, false, rank, celltag);
            //return;
        }
        else
        {
            //if (cv(0) < aabb_.min(0)) return;
            //if (cv(1) < aabb_.min(1)) return;
            //if (cv(0) > aabb_.max(0)) return;
            //if (cv(1) > aabb_.max(1)) return;
        }

        /*std::cout << "checking if inside aabb" << std::endl;
        std::cout << "aabb_.min(0): " << aabb_.min(0) << std::endl;
        std::cout << "aabb_.min(1): " << aabb_.min(1) << std::endl;
        std::cout << "aabb_.min(2): " << aabb_.min(2) << std::endl;
        std::cout << "aabb_.max(0): " << aabb_.max(0) << std::endl;
        std::cout << "aabb_.max(1): " << aabb_.max(1) << std::endl;
        std::cout << "aabb_.max(2): " << aabb_.max(2) << std::endl;*/

        if (TAILOR_ZERO < aabb_.min(0) - cv(0)) return false;
        if (TAILOR_ZERO < aabb_.min(1) - cv(1)) return false;
        if (TAILOR_ZERO < aabb_.min(2) - cv(2)) return false;
        if (cv(0) - aabb_.max(0) > TAILOR_ZERO) return false;
        if (cv(1) - aabb_.max(1) > TAILOR_ZERO) return false;
        if (cv(2) - aabb_.max(2) > TAILOR_ZERO) return false;

        //std::cout << "inside aabb" << std::endl;

        bool remx_is_zero = false;
        bool remy_is_zero = false;
        bool remz_is_zero = false;

        double _d0 = (cv(0) - aabb_.min(0)) / h(0);
        double _d1 = (cv(1) - aabb_.min(1)) / h(1);
        double _d2 = (cv(2) - aabb_.min(2)) / h(2);
        if (std::abs(_d0 - std::round(_d0)) < TAILOR_ZERO)
        {
            _d0 = std::round(_d0);
            remx_is_zero = true;
        }
        if (std::abs(_d1 - std::round(_d1)) < TAILOR_ZERO)
        {
            _d1 = std::round(_d1);
            remy_is_zero = true;
        }
        if (std::abs(_d2 - std::round(_d2)) < TAILOR_ZERO)
        {
            _d2 = std::round(_d2);
            remz_is_zero = true;
        }
        if (_d2 < 0) {
            _d2 = 0;
        }
        if (_d1 < 0) {
            _d1 = 0;
        }
        if (_d0 < 0) {
            _d0 = 0;
        }


        /*if (std::floor(_d2) >= nstripe(2))
        {
            std::cout << "floor d0: " << std::floor(_d0) << std::endl;
            std::cout << "floor d1: " << std::floor(_d1) << std::endl;
            std::cout << "floor d2: " << std::floor(_d2) << std::endl;

            std::cout << "nstripe 0: " << nstripe(0) << std::endl;
            std::cout << "nstripe 1: " << nstripe(1) << std::endl;
            std::cout << "nstripe 2: " << nstripe(2) << std::endl;

            std::cout << "_d0: " << _d0 << std::endl;
            std::cout << "_d1: " << _d1 << std::endl;
            std::cout << "_d2: " << _d2 << std::endl;

            std::cout << "cv(0): " << cv(0) << std::endl;
            std::cout << "cv(1): " << cv(1) << std::endl;
            std::cout << "cv(2): " << cv(2) << std::endl;

            std::cout << "h(0): " << h(0) << std::endl;
            std::cout << "h(1): " << h(1) << std::endl;
            std::cout << "h(2): " << h(2) << std::endl;

            std::cout << "aabb.min(0): " << aabb_.min(0) << std::endl;
            std::cout << "aabb.min(1): " << aabb_.min(1) << std::endl;
            std::cout << "aabb.min(2): " << aabb_.min(2) << std::endl;
        }*/

        //assert(std::floor(_d0) < nstripe(0));
        //assert(std::floor(_d1) < nstripe(1));
        //assert(std::floor(_d2) < nstripe(2));

        const RegularMeshIndex _index(std::floor(_d1), std::floor(_d0), std::floor(_d2));


        // the point is not on the last vertical and horizontal grid edge.
        // upper right.
        if (_index.j() < nstripe_(0) && _index.i() < nstripe_(1) && _index.k() < nstripe_(2))
        {
            assert(_index.i() >= 0 && _index.j() >= 0);
            index.push_back(_index);
            assert(index.back().j() < nstripe_(0));
            assert(index.back().i() < nstripe_(1));
            assert(index.back().k() < nstripe_(2));
        }
        // the point is not on the first but other vertical edge.
        // upper left.
        if (_index.j() > 0 && remx_is_zero)
        {
            // the point is not on the last horizontal edge.
            if (_index.i() < nstripe_(1) && _index.k() < nstripe_(2))
            {
                assert(_index.i() >= 0 && _index.j()-1 >= 0);
                index.push_back(RegularMeshIndex(_index.i(), _index.j()-1, _index.k()));
                assert(index.back().j() < nstripe_(0));
                assert(index.back().i() < nstripe_(1));
                assert(index.back().k() < nstripe_(2));
            }
        }
        // the point is not on the first but other horizontal edge.
        // lower right.
        if (_index.i() > 0 && remy_is_zero)
        {
            // the point is not on the last vertical edge.
            if (_index.j() < nstripe_(0) && _index.k() < nstripe_(2))
            {
                assert(_index.i()-1 >= 0 && _index.j() >= 0);
                index.push_back(RegularMeshIndex(_index.i()-1, _index.j(), _index.k()));
                assert(index.back().j() < nstripe_(0));
                assert(index.back().i() < nstripe_(1));
                assert(index.back().k() < nstripe_(2));
            }
        }
        // the point is not on the first but other lateral edge.
        if (_index.k() > 0 && remz_is_zero)
        {
            // the point is not on the last vertical edge.
            if (_index.i() < nstripe_(1) && _index.j() < nstripe_(0))
            {
                //assert(_index.i()-1 >= 0 && _index.j() >= 0);
                index.push_back(RegularMeshIndex(_index.i(), _index.j(), _index.k()-1));
                assert(index.back().j() < nstripe_(0));
                assert(index.back().i() < nstripe_(1));
                assert(index.back().k() < nstripe_(2));
            }
        }
        // the point is not on the first vertical and horizontal edges but on other h. and v. edges.
        // lower left.
        if (_index.i() > 0 && _index.j() > 0 && _index.k() > 0 && remx_is_zero && remy_is_zero && remz_is_zero)
        {
            assert(_index.i()-1 >= 0 && _index.j()-1 >= 0 && _index.k()-1 >= 0);
            index.push_back(RegularMeshIndex(_index.i()-1, _index.j()-1, _index.k()-1));
            assert(index.back().j() < nstripe_(0));
            assert(index.back().i() < nstripe_(1));
            assert(index.back().k() < nstripe_(2));
        }
        if (_index.i() > 0 && _index.j() > 0 && remx_is_zero && remy_is_zero)
        {
            assert(_index.i()-1 >= 0 && _index.j()-1 >= 0);
            if (_index.k() < nstripe_(2))
            {
                index.push_back(RegularMeshIndex(_index.i()-1, _index.j()-1, _index.k()));
                assert(index.back().j() < nstripe_(0));
                assert(index.back().i() < nstripe_(1));
                assert(index.back().k() < nstripe_(2));
            }
        }
        if (_index.i() > 0 && _index.k() > 0 && remx_is_zero && remz_is_zero)
        {
            assert(_index.i()-1 >= 0 && _index.k()-1 >= 0);
            if (_index.j() < nstripe_(0))
            {
                index.push_back(RegularMeshIndex(_index.i()-1, _index.j(), _index.k()-1));
                assert(index.back().j() < nstripe_(0));
                assert(index.back().i() < nstripe_(1));
                assert(index.back().k() < nstripe_(2));
            }
        }
        if (_index.j() > 0 && _index.k() > 0 && remy_is_zero && remz_is_zero)
        {
            assert(_index.j()-1 >= 0 && _index.k()-1 >= 0);
            if (_index.i() < nstripe_(1))
            {
                index.push_back(RegularMeshIndex(_index.i(), _index.j()-1, _index.k()-1));
                assert(index.back().j() < nstripe_(0));
                assert(index.back().i() < nstripe_(1));
                assert(index.back().k() < nstripe_(2));
            }
        }

        /*std::cout << "centroid(0): " << cv(0) << std::endl;
        std::cout << "centroid(1): " << cv(1) << std::endl;
        std::cout << "centroid(2): " << cv(2) << std::endl;
        std::cout << "remxiszero: " << remx_is_zero << std::endl;
        std::cout << "remyiszero: " << remy_is_zero << std::endl;
        std::cout << "remziszero: " << remz_is_zero << std::endl;
        std::cout << "index size: " << index.size() << std::endl;*/

        /*for (const auto& ii: index)
        {
            std::cout << "i: " << ii.i() << std::endl;
            std::cout << "j: " << ii.j() << std::endl;
            std::cout << "k: " << ii.k() << std::endl;
            //std::cout << "ns(0): " << nstripe_(0) << std::endl;
            //std::cout << "ns(1): " << nstripe_(1) << std::endl;
            //std::cout << "ns(2): " << nstripe_(2) << std::endl;
            //std::cout << "cond 1: " << (_index.i() > 0 && remy_is_zero) << std::endl;
            //std::cout << "cond 2: " << (_index.j() < nstripe_(0) && _index.k() < nstripe_(2)) << std::endl;
            //std::cout << "rank: " << rank << std::endl;
            //std::cout << "_min(0): " << aabb_.min(0) << std::endl;
            //std::cout << "_min(1): " << aabb_.min(1) << std::endl;
            //std::cout << "_min(2): " << aabb_.min(2) << std::endl;
            //std::cout << "_max(0): " << aabb_.max(0) << std::endl;
            //std::cout << "_max(1): " << aabb_.max(1) << std::endl;
            //std::cout << "_max(2): " << aabb_.max(2) << std::endl;
        }*/
        //assert(!index.empty());*/

        /*if (celltag == 49 && rank == 1)
        {
            std::cout << "floor d0: " << std::floor(_d0) << std::endl;
            std::cout << "floor d1: " << std::floor(_d1) << std::endl;
            std::cout << "floor d2: " << std::floor(_d2) << std::endl;

            std::cout << "nstripe 0: " << nstripe(0) << std::endl;
            std::cout << "nstripe 1: " << nstripe(1) << std::endl;
            std::cout << "nstripe 2: " << nstripe(2) << std::endl;

            std::cout << "_d0: " << _d0 << std::endl;
            std::cout << "_d1: " << _d1 << std::endl;
            std::cout << "_d2: " << _d2 << std::endl;

            std::cout << "cv(0): " << cv(0) << std::endl;
            std::cout << "cv(1): " << cv(1) << std::endl;
            std::cout << "cv(2): " << cv(2) << std::endl;

            std::cout << "h(0): " << h(0) << std::endl;
            std::cout << "h(1): " << h(1) << std::endl;
            std::cout << "h(2): " << h(2) << std::endl;

            std::cout << "aabb.min(0): " << aabb_.min(0) << std::endl;
            std::cout << "aabb.min(1): " << aabb_.min(1) << std::endl;
            std::cout << "aabb.min(2): " << aabb_.min(2) << std::endl;

            std::cout << "i: " << _index.i() << std::endl;
            std::cout << "j: " << _index.j() << std::endl;
            std::cout << "k: " << _index.k() << std::endl;

            std::cout << "remx_is_zero: " << remx_is_zero << std::endl;
            std::cout << "remy_is_zero: " << remy_is_zero << std::endl;
            std::cout << "remz_is_zero: " << remz_is_zero << std::endl;

            std::cout << "cond1: " << (_index.j() > 0 && remx_is_zero) << std::endl;
            std::cout << "cond2: " << (_index.i() < nstripe_(1) && _index.k() < nstripe_(2)) << std::endl;

            std::cout << "index size: " << index.size() << std::endl;
        }*/

        return true;
    }

    void RegularMesh::get_bintag_adaptive_unique_(const Vector3& cv, std::vector<BinRMTag>& tag, int rank, int celltag) const
    {
        //auto _index = get_index_regular_unique(cv);
        std::vector<RegularMeshIndex> _index;
        bool res = get_index_regular(cv, _index, false, rank, celltag);

        //if (!_index.isvalid()) {
        //    //std::cout << "i: " << _index.i() << std::endl;
        //    //std::cout << "j: " << _index.j() << std::endl;
        //    //std::cout << "k: " << _index.k() << std::endl;
        //    std::cout << "ccccv(0): " << cv(0) << std::endl;
        //    std::cout << "ccccv(1): " << cv(1) << std::endl;
        //    std::cout << "ccccv(2): " << cv(2) << std::endl;
        //    //std::cout << "ns(0): " << nstripe_(0) << std::endl;
        //    //std::cout << "ns(1): " << nstripe_(1) << std::endl;
        //    //std::cout << "ns(2): " << nstripe_(2) << std::endl;
        //    assert(false);
        //    return;
        //}

        //std::cout << "ccv(0): " << cv(0) << std::endl;
        //std::cout << "ccv(1): " << cv(1) << std::endl;
        //std::cout << "ccv(2): " << cv(2) << std::endl;

        //if (celltag == 365177 && rank == 8)
        //{
            //std::cout << "index size: " << _index.size() << std::endl;
        //}
        //if (_index.empty())
        //{
        //    //std::cout << "i: " << _index.i() << std::endl;
        //    //std::cout << "j: " << _index.j() << std::endl;
        //    //std::cout << "k: " << _index.k() << std::endl;
        //    std::cout << "rankk: " << rank << std::endl;
        //    std::cout << "celltag: " << celltag << std::endl;
        //    std::cout << "ccccv(0): " << cv(0) << std::endl;
        //    std::cout << "ccccv(1): " << cv(1) << std::endl;
        //    std::cout << "ccccv(2): " << cv(2) << std::endl;
        //    std::cout << "ns(0): " << nstripe_(0) << std::endl;
        //    std::cout << "ns(1): " << nstripe_(1) << std::endl;
        //    std::cout << "ns(2): " << nstripe_(2) << std::endl;
        //    std::cout << "aabb_min(0): " << aabb_.min(0) << std::endl;
        //    std::cout << "aabb_min(1): " << aabb_.min(1) << std::endl;
        //    std::cout << "aabb_min(2): " << aabb_.min(2) << std::endl;
        //    std::cout << "aabb_max(0): " << aabb_.max(0) << std::endl;
        //    std::cout << "aabb_max(1): " << aabb_.max(1) << std::endl;
        //    std::cout << "aabb_max(2): " << aabb_.max(2) << std::endl;
        //    print();
        //}
        //assert(!_index.empty());

        if (_index.empty())
        {
            std::cout << "rankk: " << rank << std::endl;
            std::cout << "celltag: " << celltag << std::endl;
            std::cout << "ccccv(0): " << cv(0) << std::endl;
            std::cout << "ccccv(1): " << cv(1) << std::endl;
            std::cout << "ccccv(2): " << cv(2) << std::endl;
            std::cout << "h(0): " << h_(0) << std::endl;
            std::cout << "h(1): " << h_(1) << std::endl;
            std::cout << "h(2): " << h_(2) << std::endl;
            std::cout << "nstripe(0): " << nstripe_(0) << std::endl;
            std::cout << "nstripe(1): " << nstripe_(1) << std::endl;
            std::cout << "nstripe(2): " << nstripe_(2) << std::endl;
            std::cout << "min(0)" << aabb_.min(0) << std::endl;
            std::cout << "min(1)" << aabb_.min(1) << std::endl;
            std::cout << "min(2)" << aabb_.min(2) << std::endl;
            std::cout << "max(0)" << aabb_.max(0) << std::endl;
            std::cout << "max(1)" << aabb_.max(1) << std::endl;
            std::cout << "max(2)" << aabb_.max(2) << std::endl;
            std::cout << "index size: " << _index.size() << std::endl;
            std::cout << "_index_.size(): " << _index.size() << std::endl;
            for (const auto& ii: _index)
            {
                std::cout << "i: " << ii.i() << std::endl;
                std::cout << "j: " << ii.j() << std::endl;
                std::cout << "k: " << ii.k() << std::endl;
            }
        }
        assert(!_index.empty());

        //_index.erase(std::remove_if(_index.begin(), _index.end(), [&](const auto& i){return !bin(i).aabb().do_intersect(cv);}), _index.end());

        //if (_index.size() != 1)
        //{
            //_index.erase(std::next(_index.begin()), _index.end());
        //}

        //assert(_index.size() == 1);

        for (const auto& i: _index)
        {
            assert(i.isvalid());
            //auto subrm = bin(_index).rm();
            auto subrm = bin(i).rm();
            //std::cout << "bin_min(0): " << bin(i).aabb().min(0) << std::endl;
            //std::cout << "bin_min(1): " << bin(i).aabb().min(1) << std::endl;
            //std::cout << "bin_min(2): " << bin(i).aabb().min(2) << std::endl;
            //std::cout << "bin_max(0): " << bin(i).aabb().max(0) << std::endl;
            //std::cout << "bin_max(1): " << bin(i).aabb().max(1) << std::endl;
            //std::cout << "bin_max(2): " << bin(i).aabb().max(2) << std::endl;
            //std::cout << "rankk: " << rank << std::endl;
            //std::cout << "subrm nullptr: " << (subrm == nullptr) << std::endl;
            if (subrm == nullptr)
            {
                //std::cout << "inter: " << bin(i).aabb().do_intersect(cv) << std::endl;
                //bool verbose = false;
                //if (celltag == 179 && rank == 1 && bin(i).tag()() == 50)
                //if (celltag == 365177 && rank == 8 && bin(i).tag()() == 18)
                //if (celltag == 365177 && rank == 8)
                //{
                //    std::cout << "bintag: " << bin(i).tag()() << std::endl;

                //    std::cout << "bin_min(0): " << bin(i).aabb().min(0) << std::endl;
                //    std::cout << "bin_min(1): " << bin(i).aabb().min(1) << std::endl;
                //    std::cout << "bin_min(2): " << bin(i).aabb().min(2) << std::endl;
                //    std::cout << "bin_max(0): " << bin(i).aabb().max(0) << std::endl;
                //    std::cout << "bin_max(1): " << bin(i).aabb().max(1) << std::endl;
                //    std::cout << "bin_max(2): " << bin(i).aabb().max(2) << std::endl;

                //    std::cout << "ccccv(0): " << cv(0) << std::endl;
                //    std::cout << "ccccv(1): " << cv(1) << std::endl;
                //    std::cout << "ccccv(2): " << cv(2) << std::endl;

                //    std::cout << "intersect: " << bin(i).aabb().do_intersect(cv, true) << std::endl;

                //    //verbose = true;
                //}
                //if (rank == 5)
                //{
                //    if (!bin(i).aabb().do_intersect(cv))
                //    {
                //        bin(i).aabb().do_intersect(cv, true);
                //        std::cout << "rank: " << rank << std::endl;
                //        std::cout << "bin_min(0): " << bin(i).aabb().min(0) << std::endl;
                //        std::cout << "bin_min(1): " << bin(i).aabb().min(1) << std::endl;
                //        std::cout << "bin_min(2): " << bin(i).aabb().min(2) << std::endl;
                //        std::cout << "bin_max(0): " << bin(i).aabb().max(0) << std::endl;
                //        std::cout << "bin_max(1): " << bin(i).aabb().max(1) << std::endl;
                //        std::cout << "bin_max(2): " << bin(i).aabb().max(2) << std::endl;

                //        std::cout << "ccccv(0): " << cv(0) << std::endl;
                //        std::cout << "ccccv(1): " << cv(1) << std::endl;
                //        std::cout << "ccccv(2): " << cv(2) << std::endl;

                //        std::cout << "bin tag: " << bin(i).tag()() << std::endl;
                //    }
                //}
                //assert(bin(i).aabb().do_intersect(cv));
                if (bin(i).aabb().do_intersect(cv))
                {
                    tag.push_back(BinRMTag(bin(i).tag(), tag_));
                }
                //tag = BinRMTag(bin(_index).tag(), tag_);
                //assert(tag.isvalid());
                //return;
            }
            else {
                //if (celltag == 365177 && rank == 8)
                //{
                    //std::cout << "checking subrm of bin: " << bin(i).tag()() << std::endl;
                //}
                subrm->get_bintag_adaptive_unique_(cv, tag, rank, celltag);
                //std::cout << "cheked subrm" << std::endl;
            }
        }
    }

    void RegularMesh::get_bintag_adaptive_unique(const Vector3& cv, BinRMTag& tag, int rank, int celltag) const
    {
        std::vector<BinRMTag> ttag;
        get_bintag_adaptive_unique_(cv, ttag, rank, celltag);
        // ttag size might not be 1. e.g. cv(0) = 1e-14. since aabb.do_intersect() involves tolerance multiple bins could be overlapped.
        //if (ttag.size() != 1)
        //{
        //    std::cout << cv(0) << " " << cv(1) << " " << cv(2) << std::endl;
        //    // print aabbs too!
        //    for (const auto& t: ttag)
        //    {
        //        std::cout << "min(0) :" << " " << bin(t).aabb().min(0) << std::endl;
        //        std::cout << "min(1) :" << " " << bin(t).aabb().min(1) << std::endl;
        //        std::cout << "min(2) :" << " " << bin(t).aabb().min(2) << std::endl;
        //        std::cout << "max(0) :" << " " << bin(t).aabb().max(0) << std::endl;
        //        std::cout << "max(1) :" << " " << bin(t).aabb().max(1) << std::endl;
        //        std::cout << "max(2) :" << " " << bin(t).aabb().max(2) << std::endl;
        //        std::cout << t.rmtag()() << " " << t.bintag()() << std::endl;
        //        print();
        //    }
        //}
        //assert(ttag.size() == 1);

        if (ttag.size() < 1)
        {
            if (rank == 8)
            {
                std::cout << "celltag: " << celltag << std::endl;
                std::cout << cv(0) << " " << cv(1) << " " << cv(2) << std::endl;
                std::cout << "ttag size: " << ttag.size() << std::endl;
            }
        }
        assert(ttag.size() >= 1);
        tag = ttag[0];
    }

    /*bool RegularMesh::get_bintag_adaptive(const Vector3& cv, std::vector<BinRMTag>& tag, bool closest) const
    {
        assert(false); // want to make sure that only get_bintag_adaptive_unique is used.
        std::vector<RegularMeshIndex> _index;
        bool res = get_index_regular(cv, _index, closest, -1);
        if (res == false)
            return res;

        for (const auto& i: _index)
        {
            auto subrm = bin(i).rm();
            if (subrm == nullptr)
            {
                tag.push_back(BinRMTag(bin(i).tag(), tag_));
                continue;
            }

            subrm->get_bintag_adaptive(cv, tag, closest);
        }

        //std::sort(tag.begin(), tag.end(), &Tag::sorter);
        //tag.erase(std::unique(tag.begin(), tag.end()), tag.end());
        //assert(tag.size() <= 1);
        
        return true;
    }*/

    /*bool RegularMesh::get_bintag_regular(const Vector3& cv, std::vector<BinRMTag>& tag, bool closest) const
    {
        std::vector<RegularMeshIndex> _index;
        bool res = get_index_regular(cv, _index, closest);
        if (res == false)
            return res;

        for (const auto& i: _index)
            tag.push_back(bin(i).tag());

        return true;
    }*/

    /*bool RegularMesh::get_index_adaptive(const Vector3& cv, std::vector<RegularMeshIndex>& index, bool closest) const
    {
        std::vector<RegularMeshIndex> _index;
        bool res = get_index_regular(cv, _index, closest);
        if (res == false)
            return res;

        for (const auto& i: _index)
        {
            auto subrm = bin(i).rm();
            if (subrm == nullptr)
            {
                index.push_back(i);
                continue;
            }

            subrm->get_index_adaptive(cv, index, closest);
        }

        return true;
    }*/

    void RegularMesh::get_bintag_adaptive(const AABB& aabb, std::vector<BinRMTag>& tag, int& nlevel, int celltag) const
    {
        nlevel = 0;
        get_bintag_adaptive_(aabb, tag, nlevel, celltag);
        std::sort(tag.begin(), tag.end(), &BinRMTag::sorter);
        tag.erase(std::unique(tag.begin(), tag.end()), tag.end());
        //for (const auto& t: tag)
        //{
            //if (t.bintag()() == 74)
            //{
                //if (celltag == 986)
                //{
                    //assert(false);
                //}
            //}
        //}
    }

    void RegularMesh::get_bintag_adaptive_(const AABB& aabb, std::vector<BinRMTag>& tag, int& nlevel, int celltag) const
    {
        ++nlevel;
        std::vector<RegularMeshIndex> _index;
        //_index.reserve(64);
        std::vector<Vector3> p;
        p.reserve(8);
        p.push_back(Vector3(aabb.min(0), aabb.min(1), aabb.min(2))); // xmin, ymin, zmin
        p.push_back(Vector3(aabb.min(0), aabb.max(1), aabb.min(2))); // xmin, ymax, zmin
        p.push_back(Vector3(aabb.max(0), aabb.min(1), aabb.min(2))); // xmax, ymin, zmin
        p.push_back(Vector3(aabb.max(0), aabb.max(1), aabb.min(2))); // xmax, ymax, zmin

        p.push_back(Vector3(aabb.min(0), aabb.min(1), aabb.max(2))); // xmin, ymin, zmax
        p.push_back(Vector3(aabb.min(0), aabb.max(1), aabb.max(2))); // xmin, ymax, zmax
        p.push_back(Vector3(aabb.max(0), aabb.min(1), aabb.max(2))); // xmax, ymin, zmax
        p.push_back(Vector3(aabb.max(0), aabb.max(1), aabb.max(2))); // xmax, ymax, zmax


        for (const auto& _p: p) 
        {
            std::vector<RegularMeshIndex> __index;
            //__index.reserve(8);
            //get_index_regular(_p, __index, true);
            get_index_regular(_p, __index, true, -1, celltag); // closest=true is needed when working with aabb.
            _index.insert(_index.end(), __index.begin(), __index.end());
        }

        int mini = TAILOR_BIG_POS_NUM;
        int minj = TAILOR_BIG_POS_NUM;
        int mink = TAILOR_BIG_POS_NUM;
        int maxi = TAILOR_BIG_NEG_NUM;
        int maxj = TAILOR_BIG_NEG_NUM;
        int maxk = TAILOR_BIG_NEG_NUM;

        for (const RegularMeshIndex& i: _index)
        {
            const Bin& _bin = bin(i);

            mini = std::min(mini, _bin.index().i());
            minj = std::min(minj, _bin.index().j());
            mink = std::min(mink, _bin.index().k());
            maxi = std::max(maxi, _bin.index().i());
            maxj = std::max(maxj, _bin.index().j());
            maxk = std::max(maxk, _bin.index().k());
        }

        if (_index.empty()) {
            return;
        }

        _index.clear();
        for (int k=mink; k<=maxk; ++k)
        {
            for (int j=minj; j<=maxj; ++j)
            {
                for (int i=mini; i<=maxi; ++i)
                {
                    _index.push_back(RegularMeshIndex(i, j, k));
                }
            }
        }

        for (const auto& i: _index)
        {
            auto subrm = bin(i).rm();

            if (subrm == nullptr)
            {
                if (bin(i).aabb().do_intersect(aabb))
                {
                    tag.push_back(BinRMTag(bin(i).tag(), tag_));
                }
                continue;
            }

            subrm->get_bintag_adaptive_(aabb, tag, nlevel, celltag);
        }

        //std::sort(tag.begin(), tag.end(), &BinRMTag::sorter);
        //tag.erase(std::unique(tag.begin(), tag.end()), tag.end());

        //return true;
    }

    bool BinRMTag::sorter(const BinRMTag& lhs, const BinRMTag& rhs)
    {
        if (lhs.rmtag() < rhs.rmtag())
        {
            return lhs.bintag() < rhs.bintag();
        }

        return lhs.bintag() < rhs.bintag();
    }

    /*bool RegularMesh::get_bintag_regular(const AABB& aabb, std::vector<Tag>& tag) const
    {
        std::vector<RegularMeshIndex> _index;
        std::vector<Vector3> p;
        p.push_back(Vector3(aabb.min(0), aabb.min(1))); // xmin, ymin
        p.push_back(Vector3(aabb.min(0), aabb.max(1))); // xmin, ymax
        p.push_back(Vector3(aabb.max(0), aabb.min(1))); // xmax, ymin
        p.push_back(Vector3(aabb.max(0), aabb.max(1))); // xmax, ymax

        for (const auto& _p: p) 
        {
            std::vector<RegularMeshIndex> __index;
            get_index_regular(_p, __index, false);
            _index.insert(_index.end(), __index.begin(), __index.end());
        }

        if (_index.empty())
            return false;

        for (const auto& i: _index)
            tag.push_back(bin(i).tag());

        std::sort(tag.begin(), tag.end(), &Tag::sorter);
        tag.erase(std::unique(tag.begin(), tag.end()), tag.end());

        return true;
    }*/

    /*bool RegularMesh::get_index_adaptive(const AABB& aabb, std::vector<RegularMeshIndex>& index, bool closest) const
    {
        std::vector<RegularMeshIndex> _index;
        std::vector<Vector3> p;
        p.push_back(aabb.min(0), aabb.max(0));
        p.push_back(aabb.min(0), aabb.max(1));
        p.push_back(aabb.min(1), aabb.max(0));
        p.push_back(aabb.min(1), aabb.max(1));

        for (const auto& _p: p) 
        {
            std::vector<RegularMeshIndex> __index;
            get_index_regular(_p, __index, closest);
            _index.insert(_index.end(), __index.begin(), __index.end());
        }

        if (_index.empty())
            return false;

        for (const auto& i: _index)
        {
            auto subrm = bin(i).rm();
            if (subrm == nullptr)
            {
                index.push_back(i);
                continue;
            }

            subrm->get_index_adaptive(cv, index, closest);
        }

        return true;
    }*/

    /*void RegularMesh::get_index(const vec3& cv, std::vector<Index>& index)
      {
      bool remx_is_zero = false;
      bool remy_is_zero = false;

      double j_ = (cv[0] - aabb.min[0]) / h[0];
      double i_ = (cv[1] - aabb.min[1]) / h[1];

      int j = std::floor(j_);
      int i = std::floor(i_);

    // check if the point is on vertical grid edge.
    if (std::abs(std::remainder(cv[0] - aabb.min[0], h[0])) < 1e-12)
    {
    remx_is_zero = true;
    }
    // check if the point is on horizontal grid edge.
    if (std::abs(std::remainder(cv[1] - aabb.min[1], h[1])) < 1e-12)
    {
    remy_is_zero = true;
    }

    // the point is not on the last vertical and horizontal grid edge.
    // upper right.
    if (j < nstripe && i < nstripe)
    {
    index.push_back(Index(i, j));
    }
    // the point is not on the first but other vertical edge.
    // upper left.
    if (j > 0 && remx_is_zero)
    {
    // the point is not on the last horizontal edge.
    if (i < nstripe)
    {
    index.push_back(Index(i, j-1));
    }
    }
    // the point is not on the first but other horizontal edge.
    // lower right.
    if (i > 0 && remy_is_zero)
    {
    // the point is not on the last vertical edge.
    if (j < nstripe)
    {
    index.push_back(Index(i-1, j));
    }
    }
    // the point is not on the first vertical and horizontal edges but on other h. and v. edges.
    // lower left.
    if (i > 0 && j > 0 && remx_is_zero && remy_is_zero)
    {
    index.push_back(Index(i-1, j-1));
    }
    }*/

    bool RegularMesh::point_inside_bin(const RegularMeshIndex& ind, const Vector3& p)
    {
        Vector3 ll = llcoor(ind);

        if (p(0) >= ll(0) && p(0) <= ll(0)+h(0))
        {
            if (p(1) >= ll(1) && p(1) <= ll(1)+h(1))
            {
                return true;
            }
        }

        return false;
    }

    void RegularMesh::print(std::string aux) const
    {
        std::string s = "rm_";
        s.append(aux);
        s.append("_");
        s.append(std::to_string(tag_()));
        s.append(".vtk");
        print_(s);

        for (const Bin& b: bin_)
        {
            if (b.rm() == nullptr) {
                continue;
            }
            b.rm()->print(aux);
        }
    }

    void RegularMesh::print_(std::string file_name) const
    {    
        std::ofstream out;    
        out.open (file_name);

        out << "# vtk DataFile Version 3.0" << std::endl;
        out << "Regular Mesh" << std::endl;
        out << "ASCII" << std::endl;
        out << "DATASET STRUCTURED_GRID" << std::endl;
        out << "DIMENSIONS " << nstripe(0)+1 << " " << nstripe(1)+1 << " " << nstripe(2)+1 << std::endl;
        out << "POINTS " << (nstripe(0)+1)*(nstripe(1)+1)*(nstripe(2)+1) << " float" <<  std::endl;

        for (int k=0; k<nstripe(2)+1; ++k)
        {
            for (int i=0; i<nstripe(1)+1; ++i)
            {
                for (int j=0; j<nstripe(0)+1; ++j)
                {
                    out << aabb().min(0) + j*h(0);
                    out << " ";
                    out << aabb().min(1) + i*h(1);
                    out << " ";
                    out << aabb().min(2) + k*h(2);
                    out << std::endl;
                }
            }
        }

        out << "CELL_DATA " << nstripe(0)*nstripe(1)*nstripe(2) << std::endl;
        //out << "SCALARS " << "load " << "float " << "1" << std::endl;
        //out << "LOOKUP_TABLE default" << std::endl;    
        //for (int k=0; k<nstripe(2); ++k)
        //{
            //for (int i=0; i<nstripe(1); ++i)
            //{
                //for (int j=0; j<nstripe(0); ++j)
                //{
                    ////out << bin(i,j).load(mesh);
                    //out << bin(i,j,k).load_without_calc();
                    //out << std::endl;
                //}
            //}
        //}

        out << "SCALARS " << "tag " << "int " << "1" << std::endl;
        out << "LOOKUP_TABLE default" << std::endl;    
        for (const Bin& b: bin_)
        {
            out << b.tag()();
            out << std::endl;
        }

        out << "SCALARS " << "parentrm " << "int " << "1" << std::endl;
        out << "LOOKUP_TABLE default" << std::endl;    
        for (const Bin& b: bin_)
        {
            out << b.parent_rm()();
            out << std::endl;
        }

        out << "SCALARS " << "load " << "int " << "1" << std::endl;
        out << "LOOKUP_TABLE default" << std::endl;    
        for (const Bin& b: bin_)
        {
            out << b.load();
            out << std::endl;
        }

        out.close();
    }


    void RegularMesh::info() const
    {
        for (const Bin& b: bin_)
        {
            for (const BinCell& bc: b.cell())
            {
                std::cout << "bin: " << &b - bin_.data() << " mesh: " << bc.mesh()() << " cell: " << bc.cell()() << std::endl;
            }
        }
    }
}
