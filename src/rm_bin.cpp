#include "rm_bin.h"
#include "regular_mesh.h"

namespace Tailor
{
    void Bin::add_nei(int rmtag, int bintag)
    {
        nei_.push_back(BinRMTag(Tag(bintag), Tag(rmtag)));
    }

    const std::deque<BinRMTag>& Bin::nei() const
    {
        return nei_;
    }

    double Bin::load() const
    {
        return load_;
    }

    void Bin::set_load(double l)
    {
        load_ = l;
    }

    void Bin::move(const vec3<double>& v)
    {
        aabb_.move_points(v);
        if (rm_ != nullptr)
        {
            rm_->move(v);
        }
    }

    void Bin::rotate(double ang, int axis, const vec3<double>& rot_point)
    {
        aabb_.rotate_points(ang, axis, rot_point);
        if (rm_ != nullptr)
        {
            rm_->rotate(ang, axis, rot_point);
        }
    }

    bool Bin::operator<(const Bin& other) const
    {
        return tag_ < other.tag_;
    }

    bool Bin::operator==(const Bin& other) const
    {
        return tag_ == other.tag_;
    }

    void Bin::clear_cells()
    {
        //std::vector<BinCell>().swap(cell_);
        cell_.clear(); // I suspect that swapping is slow due to reallocation. lets try clear now. this should affect second-remap-remove.
        mesh_tag_index_map_all_.clear();
        mesh_tag_index_map_res_.clear();

        if (rm_ != nullptr)
        {
            rm_->clear_cells();
        }
    }

    double Bin::load_without_calc() const
    {
        return load_;
    }

    /*std::deque<std::vector<Point>> Bin::mesh_pts(const std::deque<Mesh>& mesh) const
    {
        //std::deque<std::vector<MeshPoint>> mpts(mesh.size());
        std::deque<std::vector<Point>> pts(mesh.size());

        for (const BinCell& bc: cell_)
        {
            auto m = std::find_if(mesh.begin(), mesh.end(), [&](const Mesh& tm){return tm.tag() == bc.mesh();});
            assert(!mesh.empty());
            if (m == mesh.end())
            {
                std::cout << "m.tag(): " << bc.mesh()() << std::endl;
                std::cout << "mesh size: " << mesh.size() << std::endl;

            }
            assert(m != mesh.end());
            int i = std::distance(mesh.begin(), m);

            const MeshCell& mc = m->cell(bc.cell());
            if (!is_resident(mc.poly().centroid())) {
                continue;
            }

            std::vector<Point>& pt = pts[i];

            for (const MeshPoint& mp: mc.point())
            {
                pt.push_back(mp.p());
            }
        }

        //for (auto& a: mpts)
        //{
            //remove_dup(a);
        //}

        //for (int i=0; i<mpts.size(); ++i)
        //{
            //for (const auto& a: mpts[i])
            //{
                //pts[i].push_back(a.p());
            //}
        //}
        
        return pts;
    }*/

    std::deque<AABB> Bin::mesh_aabb(const std::deque<Mesh>& mesh) const
    {
        std::deque<AABB> aabb(mesh.size());

        if (cell_.empty()) {
            return aabb;
        }

        for (AABB& ab: aabb)
        {
            vec3<double>minn(TAILOR_BIG_POS_NUM, TAILOR_BIG_POS_NUM, TAILOR_BIG_POS_NUM);
            vec3<double>maxx(TAILOR_BIG_NEG_NUM, TAILOR_BIG_NEG_NUM, TAILOR_BIG_NEG_NUM);
            ab.set_bbox(minn, maxx);
        }

        for (const BinCell& bc: cell_)
        {
            auto m = std::find_if(mesh.begin(), mesh.end(), [&](const Mesh& tm){return tm.tag() == bc.mesh();});
            assert(!mesh.empty());
            if (m == mesh.end())
            {
                std::cout << "m.tag(): " << bc.mesh()() << std::endl;
                std::cout << "mesh size: " << mesh.size() << std::endl;

            }
            assert(m != mesh.end());
            int i = std::distance(mesh.begin(), m);

            const MeshCell& mc = m->cell(bc.cell());
            if (!is_resident(mc.poly().centroid())) {
                continue;
            }
            AABB& ab = aabb[i];

            vec3<double> min_, max_;
            min_.set(ab.min(0), ab.min(1), ab.min(2));
            max_.set(ab.max(0), ab.max(1), ab.max(2));

            for (const MeshPoint& mp: mc.point())
            {
                double vx = mp.p().r(0);
                double vy = mp.p().r(1);
                double vz = mp.p().r(2);

                min_.set_x(std::min(min_(0), vx));
                min_.set_y(std::min(min_(1), vy));
                min_.set_z(std::min(min_(2), vz));
                max_.set_x(std::max(max_(0), vx));
                max_.set_y(std::max(max_(1), vy));
                max_.set_z(std::max(max_(2), vz));
            }

            //ab.set_min(min_);
            //ab.set_max(max_);
            ab.set_bbox(min_, max_);
        }

        for (AABB& ab: aabb)
        {
            ab.set_vertices_from_bbox();
            ab.set_faces();
        }

        for (const AABB& ab: aabb)
        {
            assert(!ab.vertices().empty());
            assert(!ab.faces().empty());
        }

        return aabb;
    }

    const Tag& Bin::parent_rm() const
    {
        return parent_rm_;
    }

    void BinCell::set_proc(int proc)
    {
        proc_ = proc;
    }

    int BinCell::proc() const
    {
        return proc_;
    }

    Bin::Bin(const Tag& t): tag_(t), rm_(nullptr), load_(0.)
    {
    }

    Bin::Bin(const Tag& t, const Tag& parent_rm, const RegularMeshIndex& index, const AABB& aabb, int mesh_load_size): aabb_(aabb), tag_(t), rm_(nullptr), parent_rm_(parent_rm), load_(0.)
    {
        assert(parent_rm.isvalid());
        set_index(index);
        assert(mesh_load_size >= 0);
        //if (mesh_load_size > 0)
        //{
            //mesh_load_.resize(mesh_load_size);
        //}
    }

    Bin::~Bin()
    {
        //std::cout << "deleting" << std::endl;
        if (rm_ != nullptr)
        {
            //std::cout << "rrrrrrrrrrrrrrrrrrrrrrrrrrrr" << std::endl;
            delete rm_;
        }
        //std::cout << "deleted" << std::endl;
    }

    Bin& Bin::operator=(const Bin& other)
    {
        tag_ = other.tag_;
        parent_rm_ = other.parent_rm_;
        rmindex_ = other.rmindex_;
        mesh_load_ = other.mesh_load_;
        cell_ = other.cell_;
        mesh_tag_index_map_all_ = other.mesh_tag_index_map_all_;
        mesh_tag_index_map_res_ = other.mesh_tag_index_map_res_;
        aabb_ = other.aabb_;
        load_ = other.load_;

        rm_ = nullptr;
        if (other.rm_ != nullptr)
            rm_ = new RegularMesh(*other.rm_);
        else
            assert(rm_ == nullptr);

        return *this;
    }

    Bin::Bin(const Bin& other):
        tag_(other.tag_),
        parent_rm_(other.parent_rm_),
        rmindex_(other.rmindex_),
        mesh_load_(other.mesh_load_),
        cell_(other.cell_),
        mesh_tag_index_map_all_(other.mesh_tag_index_map_all_),
        mesh_tag_index_map_res_(other.mesh_tag_index_map_res_),
        aabb_(other.aabb_),
        load_(other.load_)
    {
        rm_ = nullptr;
        if (other.rm_ != nullptr)
            rm_ = new RegularMesh(*other.rm_);
        else
            assert(rm_ == nullptr);
    }

    const boost::bimap<int, int>& Bin::mesh_tag_index_map_res() const
    {
        return mesh_tag_index_map_res_;
    }
    const boost::bimap<int, int>& Bin::mesh_tag_index_map_all() const
    {
        return mesh_tag_index_map_all_;
    }
    const Tag& Bin::tag() const
    {
        return tag_;
    }

    void Bin::register_cells_to_rm(const std::deque<Mesh>& meshes, bool pseudo3D, int rank, RegType regtype)
    {
        for (const BinCell& bc: cell_)
        {
            Tag mt = bc.mesh();
            Tag ct = bc.cell();

            for (auto m=meshes.begin(); m!=meshes.end(); ++m)
            {
                if (mt == m->tag()) {
                    assert(m->query(ct));
                }
            }
        }

        assert(aabb_.face_cross_len() == false);
        for (const Bin& b: rm_->bin())
        {
            assert(b.aabb().face_cross_len() == false);
        }

        rm_->register_bincells(cell_, meshes, pseudo3D, rank, regtype);

        if (!mesh_tag_index_map_res_.left.empty())
        {
            assert(!cell_.empty());
        }

        //std::vector<BinCell>().swap(cell_);
        cell_.clear();
        mesh_tag_index_map_res_.clear();
        mesh_tag_index_map_all_.clear();

        if (!mesh_tag_index_map_res_.left.empty())
        {
            assert(!cell_.empty());
        }
    }

    void Bin::init_rm(const Bin& other)
    {
        //rm_ = boost::make_shared<RegularMesh>(*other.rm());
        //rm_ = new RegularMesh();
        rm_ = new RegularMesh(*other.rm());
        //*rm_ = *other.rm();
        //rm_->tag_ = other.rm().tag();
        //rm_->nstripe_ = other.rm().nstripe();
        //rm_->h_ = other.rm().h();
        //rm_->aabb_ = other.rm().aabb();
        //rm_->bin_.resize(other.rm()->bin().size());
        //for (int i=0; i<rm_->bin_.size(); ++i)
        //{
            //rm_->bin_[i].set_tag(other.rm().bin(i).tag());
        //}
    }

    void Bin::init_rm(int mintag, const Tag& rmtag, const vec3<int>& nstripe, bool pseudo3D)
    {
        //rm_ = boost::make_shared<RegularMesh>();
        rm_ = new RegularMesh();
        rm_->set_tag(rmtag);
        rm_->set_aabb(aabb_);
        assert(std::abs(rm_->aabb().min(0) - rm_->aabb().max(0)) > TAILOR_ZERO);
        if (pseudo3D) {
            rm_->set_nstripe(nstripe(0), nstripe(1), 1);
        }
        else {
            rm_->set_nstripe(nstripe);
        }
        rm_->calc_h();
        if (rm_->h(0) <= TAILOR_ZERO)
        {
            std::cout << rm_->h(0) << std::endl;
            std::cout << "min0: " << aabb_.min(0) << std::endl;
            std::cout << "min1: " << aabb_.min(1) << std::endl;
            std::cout << "max0: " << aabb_.max(0) << std::endl;
            std::cout << "max1: " << aabb_.max(1) << std::endl;
        }
        assert(rm_->h(0) > TAILOR_ZERO);
        assert(rm_->h(1) > TAILOR_ZERO);
        rm_->insert_bins(mintag, 0, 0);
        for (const Bin& b: rm_->bin())
        {
            assert(b.cell().empty());
            assert(b.aabb().min(0) != b.aabb().max(0));
        }
    }

    //boost::shared_ptr<RegularMesh> Bin::rm() const
    RegularMesh* Bin::rm() const
    {
        return rm_;
    }

    bool Bin::is_resident(const vec3<double>& cnt) const
    {
        return aabb_.do_intersect(cnt);
    }

    void Bin::set_aabb(const AABB& aabb)
    {
        aabb_ = aabb;
    }

    //void Bin::prepare_transfer_info_with_mesh(int dest, std::vector<MeshTransferInfo>& mti, std::vector<int>& mti_proc, std::vector<int>& dest_rank, const std::deque<Mesh>& mesh) const
    //{
        //assert(parent_rm_.isvalid());
        //assert(!cell_.empty());
        //for (int i=0; i<cell_.size(); ++i)
        //{
            //const Tag& mt = cell(i).mesh();
            //const Tag& ct = cell(i).cell();
            //int pc = cell(i).proc();
            //auto itt = std::find_if(mesh.begin(), mesh.end(), [&](const Mesh& m){return m.tag() == mt;});
            //assert(itt != mesh.end());
            //if (!is_resident(itt->cell(ct).polygon().centroid())) {
                //continue;
            //}

            //auto it = std::find_if(mti.begin(), mti.end(), [&] (const MeshTransferInfo& mti_) {return mti_.source() == pc;});
            //if (it == mti.end())
            //{
                //mti.push_back(MeshTransferInfo(pc, dest, mt, ct, BinRMTag(tag_, parent_rm_)));
                //mti_proc.push_back(pc);
                //if (pc != dest)
                //{
                    //++dest_rank[dest];
                //}
                //assert(pc != 0);
            //}
            //else
            //{
                ////if (pc != dest && !it->bin_exist(tag_))
                //if (pc != dest && !it->dest_exist(dest))
                //{
                    //++dest_rank[dest];
                //}
                //it->add(dest, mt, ct, BinRMTag(tag_, parent_rm_));
            //}
        //}

        //assert(!mti.empty());
    //}
    
    /*void Bin::prepare_transfer_info(int dest, std::vector<MeshTransferInfo>& mti, std::vector<int>& mti_proc, std::vector<int>& dest_rank, int nbin) const
    {
        assert(parent_rm_.isvalid());
        if (cell_.empty()) {
            return;
        }
        //assert(!cell_.empty());
        for (int i=0; i<cell_.size(); ++i)
        {
            const Tag& mt = cell(i).mesh();
            const Tag& ct = cell(i).cell();
            int pc = cell(i).proc();

            auto it = std::find_if(mti.begin(), mti.end(), [&] (const MeshTransferInfo& mti_) {return mti_.source() == pc;});
            if (it == mti.end())
            {
                mti.push_back(MeshTransferInfo(pc, dest, mt, ct, BinRMTag(tag_, parent_rm_)));
                mti_proc.push_back(pc);
                if (pc != dest)
                {
                    ++dest_rank[dest];
                }
                assert(pc != 0);
            }
            else
            {
                //if (pc != dest && !it->bin_exist(tag_))
                if (pc != dest && !it->dest_exist(dest))
                {
                    ++dest_rank[dest];
                }
                it->add(dest, mt, ct, BinRMTag(tag_, parent_rm_));
            }

            if (dest_rank[dest] >= nbin)
            {
                std::cout << "nbin: " << nbin << std::endl;
                std::cout << "dest_rank[dest]: " << dest_rank[dest] << std::endl;
                std::cout << "dest size: " << it->dest().size() << std::endl;
                std::cout << "dest: " << dest << std::endl;
                for (const auto& dd: it->dest())
                {
                    std::cout << "rank: " << dd.rank_ << std::endl;
                }
            }
            assert(dest_rank[dest] < nbin);
        }

        assert(!mti.empty());
    }*/

    void Bin::group_mesh_res(const std::deque<Mesh>& mesh, std::deque<Mesh>& mb, int rank) const
    {
        //std::deque<Mesh> mb(mesh_tag_index_map_.left.size());
        //assert(!mesh.empty());
        //assert(mesh_tag_index_map_.empty());
        /*if (rank == 3)
        {
            if (tag_() == 17)
            {
            std::cout << "map zize: " << mesh_tag_index_map_.size() << std::endl;
            assert(false);
            }
        }*/
        mb.resize(mesh_tag_index_map_res_.left.size());
                        //if(mb.size() != 5)
                            //std::cout << "mb size= " << mb.size() << std::endl;
                        //assert(mb.size() == 5);
        int i = 0;
        for (const auto itt: mesh_tag_index_map_res_.left)
        {
            mb[i].set_tag(Tag(itt.first));
            //mb[i].set_parent_mesh(itt.first);
            ++i;
        }

        //for (const Mesh& m: mb)
        //{
            //assert(m.cell().size() == 0);
            //assert(m.n_interior_cell() == 0);
            //assert(m.n_ghost() == 0);
        //}

        if (!mesh_tag_index_map_res_.left.empty()) {
            assert(!cell_.empty());
        }

        for (int i=0; i<cell_.size(); ++i)
        {
            const Tag& mt = cell(i).mesh();
            const Tag& ct = cell(i).cell();

            int j = 0;
            for (const auto itt: mesh_tag_index_map_res_.left) // why use this? // this causes only bin-resident cells to be grouped. better solution is to add non-bin-resident cells into mesh_tag_index_map as well. The reason they were't was due to calculate loads based on residents only. We need to distinguish properly.
            {
                if (itt.first == mt())
                {
                    auto it = std::find_if(mesh.begin(), mesh.end(), [&mt](const Mesh& _m){return _m.tag() == mt;});
                    assert(it != mesh.end());
                    //assert(it->query(ct) != nullptr);
                    //assert(it->cell(ct).pnei().empty());

                    const MeshCell& mc = MeshCell(it->cell(ct));
                    mb[j].add_interior_cell(mc);

                    for (const Tag& t: mc.wall_boundary()) {
                        auto bb = it->wall_boundary(t);
                        assert(bb != nullptr);
                        mb[j].add_wall_boundary(*bb);
                    }
                    for (const Tag& t: mc.dirichlet_boundary()) {
                        auto bb = it->dirichlet_boundary(t);
                        assert(bb != nullptr);
                        mb[j].add_dirichlet_boundary(*bb);
                    }
                    for (const Tag& t: mc.farfield_boundary()) {
                        auto bb = it->farfield_boundary(t);
                        assert(bb != nullptr);
                        mb[j].add_farfield_boundary(*bb);
                    }
                    for (const Tag& t: mc.empty_boundary()) {
                        auto bb = it->empty_boundary(t);
                        assert(bb != nullptr);
                        mb[j].add_empty_boundary(*bb);
                    }
                    for (const Tag& t: mc.interog_boundary()) {
                        auto bb = it->interog_boundary(t);
                        assert(bb != nullptr);
                        mb[j].add_interog_boundary(*bb);
                    }
                    break;
                }
                ++j;
            }
            /*if (!found)
            {
                std::cout << "sss: " << mesh_tag_index_map_.size() << std::endl;
                std::cout << "mt: " << mt() << std::endl;
                std::cout << "ct: " << ct() << std::endl;
                std::cout << "bin tag: " << tag_() << std::endl;
            }
            assert(found);*/
        }

        //for (const Mesh& m: mb)
        //{
            //assert(m.cell().size() == (m.n_interior_cell()) + m.n_ghost());
        //}

        //return mb;
    }

    void Bin::group_mesh_all(const std::deque<Mesh>& mesh, std::deque<Mesh>& mb, int rank) const
    {
        //std::deque<Mesh> mb(mesh_tag_index_map_.left.size());
        //assert(!mesh.empty());
        //assert(mesh_tag_index_map_.empty());
        /*if (rank == 3)
        {
            if (tag_() == 17)
            {
            std::cout << "map zize: " << mesh_tag_index_map_.size() << std::endl;
            assert(false);
            }
        }*/
        mb.resize(mesh_tag_index_map_all_.left.size());
                        //if(mb.size() != 5)
                            //std::cout << "mb size= " << mb.size() << std::endl;
                        //assert(mb.size() == 5);
        int i = 0;
        for (const auto itt: mesh_tag_index_map_all_.left)
        {
            mb[i].set_tag(Tag(itt.first));
            //mb[i].set_parent_mesh(itt.first);
            ++i;
        }

        //for (const Mesh& m: mb)
        //{
            //assert(m.cell().size() == 0);
            //assert(m.n_interior_cell() == 0);
            //assert(m.n_ghost() == 0);
        //}

        if (!mesh_tag_index_map_all_.left.empty()) {
            assert(!cell_.empty());
        }

        for (int i=0; i<cell_.size(); ++i)
        {
            const Tag& mt = cell(i).mesh();
            const Tag& ct = cell(i).cell();

            int j = 0;
            for (const auto itt: mesh_tag_index_map_all_.left) // why use this? // this causes only bin-resident cells to be grouped. better solution is to add non-bin-resident cells into mesh_tag_index_map as well. The reason they were't was due to calculate loads based on residents only. We need to distinguish properly.
            {
                if (itt.first == mt())
                {
                    auto it = std::find_if(mesh.begin(), mesh.end(), [&mt](const Mesh& _m){return _m.tag() == mt;});
                    assert(it != mesh.end());
                    assert(it->query(ct) != nullptr);
                    //assert(it->cell(ct).pnei().empty());

                    const MeshCell& mc = MeshCell(it->cell(ct));

                    //if (mc.oga_cell_type() == OGA_cell_type_t::non_resident || mc.oga_cell_type() == OGA_cell_type_t::ghost)
                    //if (mc.oga_cell_type() == OGA_cell_type_t::ghost)
                    //{
                        //continue;
                    //}

                    mb[j].add_interior_cell(mc);

                    for (const Tag& t: mc.wall_boundary()) {
                        auto bb = it->wall_boundary(t);
                        if (bb == nullptr)
                        {
                            std::cout << "oga type: " << static_cast<int>(mc.oga_cell_type()) << std::endl;
                        }
                        assert(bb != nullptr);
                        mb[j].add_wall_boundary(*bb);
                    }
                    for (const Tag& t: mc.dirichlet_boundary()) {
                        auto bb = it->dirichlet_boundary(t);
                        assert(bb != nullptr);
                        mb[j].add_dirichlet_boundary(*bb);
                    }
                    for (const Tag& t: mc.farfield_boundary()) {
                        auto bb = it->farfield_boundary(t);
                        assert(bb != nullptr);
                        mb[j].add_farfield_boundary(*bb);
                    }
                    for (const Tag& t: mc.empty_boundary()) {
                        auto bb = it->empty_boundary(t);
                        assert(bb != nullptr);
                        mb[j].add_empty_boundary(*bb);
                    }
                    for (const Tag& t: mc.interog_boundary()) {
                        auto bb = it->interog_boundary(t);
                        assert(bb != nullptr);
                        mb[j].add_interog_boundary(*bb);
                    }
                    break;
                }
                ++j;
            }
            /*if (!found)
            {
                std::cout << "sss: " << mesh_tag_index_map_.size() << std::endl;
                std::cout << "mt: " << mt() << std::endl;
                std::cout << "ct: " << ct() << std::endl;
                std::cout << "bin tag: " << tag_() << std::endl;
            }
            assert(found);*/
        }

        //for (const Mesh& m: mb)
        //{
            //assert(m.cell().size() == (m.n_interior_cell()) + m.n_ghost());
        //}

        //return mb;
    }

    const BinCell& Bin::cell(int i) const
    {
        return cell_[i];
    }

    bool BinCell::operator<(const BinCell& other) const
    {
        if (mesh_() == other.mesh()())
        {
            return cell_() < other.cell()();
        }

        return mesh_() < other.mesh()();
    }

    /*void Bin::append_to_cell(const std::vector<BinCell>& other)
    {
        cell_.insert(cell_.end(), other.begin(), other.end());
        for (int i=other.size()-1; i>=0; --i)
            bincell_index_map.insert(boost::bimap<BinCell, int>::value_type(cell_[cell_.size() - 1 - i], cell_.size() - i - 1));
    }*/

    void Bin::merge(const Bin& other)
    {
        if (rm_ == nullptr)
        {
            if (!mesh_tag_index_map_res_.left.empty())
            {
                assert(!cell_.empty());
            }
            //for (const auto& itt: other.mesh_tag_index_map_.left)
            //{
                //int other_tag = itt.first;
                //increment_mesh_load(Tag(other_tag), other.mesh_load(other_tag));
                //assert(!std::isnan(other.load_without_calc()));
                //load_ += other.load_without_calc();
                //auto it = mesh_tag_index_map_.left.find(other_tag);
                //if (it == mesh_tag_index_map_.left.end())
                    //mesh_tag_index_map_.insert(boost::bimap<int, int>::value_type(other_tag, mesh_tag_index_map_.size()));
            //}
            
            for (const auto& itt: other.mesh_tag_index_map_res_.left)
            {
                int other_tag = itt.first;
                increment_mesh_load(Tag(other_tag), other.mesh_load(Tag(other_tag)));
                //assert(!std::isnan(other.load_without_calc()));
                //load_ += other.load_without_calc();
                //auto it = mesh_tag_index_map_.left.find(other_tag);
                //if (it == mesh_tag_index_map_.left.end())
                    //mesh_tag_index_map_.insert(boost::bimap<int, int>::value_type(other_tag, mesh_tag_index_map_.size()));
            }

            cell_.insert(cell_.end(), other.cell().begin(), other.cell().end());

            //if (!other.mesh_tag_index_map_.left.empty())
            //{
                //assert(!mesh_tag_index_map_.left.empty());
                //assert(!cell_.empty());
            //}

            //if (cell_.empty())
            //{
                //assert(mesh_tag_index_map_.left.empty());
                //assert(other.mesh_tag_index_map_.left.empty());
            //}

            for (int i=other.cell().size()-1; i>=0; --i)
            {
                //bincell_index_map.insert(boost::bimap<BinCell, int>::value_type(cell_[cell_.size() - 1 - i], cell_.size() - i - 1));
                //bincell_index_map.insert(std::pair<BinCell, int>(cell_[cell_.size() - 1 - i], cell_.size() - i - 1)); // must present if.
            }
        }
        else
        {
            assert(other.rm() != nullptr);
            rm_->merge(*other.rm());
        }
    }

    const std::vector<double>& Bin::mesh_load() const
    {
        return mesh_load_;
    }

    void Bin::resize_mesh_load(size_t size)
    {
        assert(mesh_load_.empty());
        mesh_load_.resize(size);
    }

    const RegularMeshIndex& Bin::index() const
    {
        return rmindex_;
    }

    /*int BinCell::residency() const
      {
      return resident_;
      }*/

    const AABB& Bin::aabb() const
    {
        return aabb_;
    }

    /*Tag Bin::generate_celltag() const
    {
        int i = TAILOR_BIG_NEG_NUM;
        Tag t;

        if (cell_.empty())
        {
            t.set(0);
            return t;
        }
        for (const BinCell& c: cell_)
        {
            i = std::max(i, c.tag()()); 
        }

        t.set(i+1);

        return t;
    }*/

    void Bin::add_bincell(const Tag& ctag, const Tag& mtag, int proc)
    {
        BinCell c;
        c.set_cell(ctag);
        c.set_mesh(mtag);
        c.set_proc(proc);
        //assert(!query_bincell(c.mesh(), c.cell()));
        cell_.push_back(c);
        // insert to bimap.
        //cell_tag_index_map.insert(boost::bimap<int, int>::value_type(c.tag()(), cell_.size() - 1));
        //bincell_index_map.insert(boost::bimap<BinCell, int>::value_type(cell_.back(), cell_.size() - 1));
        //bincell_index_map.insert(std::pair<BinCell, int>(cell_.back(), cell_.size() - 1)); // must present if map exists.
    }

    void Bin::reserve(size_t size)
    {
        //cell_.reserve(size);
    }

    void Bin::shrink()
    {
        //std::vector<BinCell>(cell_).swap(cell_);
    }

    void Bin::copy_bincell(const BinCell& c)
    {
        //assert(!query_bincell(c.mesh(), c.cell()));
        cell_.push_back(c);
        //bincell_index_map.insert(boost::bimap<BinCell, int>::value_type(cell_.back(), cell_.size() - 1));
        //bincell_index_map.insert(std::pair<BinCell, int>(cell_.back(), cell_.size() - 1)); for testing uncommented.
    }

    //const boost::bimap<int, int>& Bin::mesh_list() const
    //{
        //return mesh_tag_index_map_;
    //}

    int Bin::mesh_load(const Tag& t) const
    {
        assert(t.isvalid());

        auto it = mesh_tag_index_map_res_.left.find(t());
        assert(it != mesh_tag_index_map_res_.left.end());
        return mesh_load_[it->second];
    }

    void Bin::insert_to_mesh_tag_index_map_res(const Tag& t)
    {
        assert(t.isvalid());

        auto it = mesh_tag_index_map_res_.left.find(t());

        if (it == mesh_tag_index_map_res_.left.end())
        {
            mesh_tag_index_map_res_.insert(boost::bimap<int, int>::value_type(t(), mesh_tag_index_map_res_.left.size()));
        }
    }

    void Bin::insert_to_mesh_tag_index_map_all(const Tag& t)
    {
        assert(t.isvalid());

        auto it = mesh_tag_index_map_all_.left.find(t());

        if (it == mesh_tag_index_map_all_.left.end())
        {
            mesh_tag_index_map_all_.insert(boost::bimap<int, int>::value_type(t(), mesh_tag_index_map_all_.left.size()));
        }
    }

    void Bin::increment_mesh_load(const Tag& t, int load)
    {
        insert_to_mesh_tag_index_map_res(t);
        insert_to_mesh_tag_index_map_all(t);

        if (mesh_tag_index_map_res_.size() > mesh_load_.size())
        {
            mesh_load_.resize(mesh_tag_index_map_res_.left.size());
        }

        auto it = mesh_tag_index_map_res_.left.find(t());
        assert(it != mesh_tag_index_map_res_.left.end());
        assert(it->second >= 0);
        assert(it->second < mesh_load_.size());
        if (load == -1) {
            ++mesh_load_[it->second];
        }
        else {
            mesh_load_[it->second] += load;
        }
    }

    void BinCell::set_mesh(const Tag& tag)
    {
        mesh_ = tag;
    }

    /*void BinCell::set_tag(const Tag& tag)
    {
        tag_ = tag;
    }*/

    void BinCell::set_cell(const Tag& tag)
    {
        cell_ = tag;
    }
    /*void BinCell::set_residency(int i)
      {
      assert(i >= 0);
      assert(i < 2);
      resident_ = i;
      }*/
    const Tag& BinCell::mesh() const
    {
        return mesh_;
    }
    /*const Tag& BinCell::tag() const
    {
        return tag_;
    }*/

    const Tag& BinCell::cell() const
    {
        return cell_;
    }

    //void Bin::set_cell(const std::vector<BinCell>& c)
    void Bin::set_cell(const std::deque<BinCell>& c)
    {
        cell_ = c;
    }

    //const std::vector<BinCell>& Bin::cell() const
    const std::deque<BinCell>& Bin::cell() const
    {
        return cell_;
    }

    /*bool Bin::query_bincell(const Tag& im, const Tag& ic) const
    {
        //auto itm = mesh_tag_index_map.left.find(im());
        //auto itc = cell_tag_index_map.left.find(ic());
        BinCell _bc;
        _bc.set_mesh(im);
        _bc.set_cell(ic);
        //auto it = bincell_index_map.left.find(_bc);
        auto it = bincell_index_map.find(_bc);

        //if (it == bincell_index_map.left.end())
        if (it == bincell_index_map.end())
            return false;

        //if (itm == mesh_tag_index_map.left.end())
        //return false;
        //if (itc == cell_tag_index_map.left.end())
        //return false;
        //if (itm->second != itc->second)
        //return false;

        return true;
    }*/

    void Bin::set_index(const RegularMeshIndex& ind)
    {
        rmindex_ = ind;
    }
}
