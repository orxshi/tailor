#include "spc.h"

namespace Tailor
{
    void SpatialPartitionContainer::set_comm(boost::mpi::communicator* comm)
    {
        assert(comm != nullptr);
        comm_ = comm;

        for (auto& sp: sp_)
        {
            sp.set_comm(comm_);
        }
    }

    void SpatialPartitionContainer::set_profiler(Profiler* prof)
    {
        assert(prof != nullptr);
        profiler_ = prof;

        for (auto& sp: sp_)
        {
            sp.set_profiler(prof);
        }
    }

    int SpatialPartitionContainer::nresi() const
    {
        int n = 0;
        for (const auto& sp: sp_)
        {
            for (const auto& m: sp.mesh())
            {
                for (const auto& mc: m.cell())
                {
                    Tag t;
                    bool isresi = is_resident(mc.poly().centroid(), -1, t);
                    if (isresi)
                    {
                        ++n;
                    }
                }
            }
        }

        int nn;
        boost::mpi::all_reduce(*comm_, n, nn, std::plus<double>());

        return nn;
    }

    void SpatialPartitionContainer::reset_mesh_connectivity()
    {
        for (auto& sp: sp_)
        {
            sp.reset_mesh_connectivity();
        }
    }

    double cell_imbalance(const SpatialPartitionContainer& spc, double& max, double& min)
    {
        std::vector<double> load_global;

        boost::mpi::all_gather(*spc.comm(), solver_load(spc.sp().front()), load_global);
        spc.comm()->barrier();

        auto result = std::minmax_element(load_global.begin(), load_global.end());
        double minload = *result.first;
        double maxload = *result.second;
        max = maxload;
        min = minload;
        double fl = ((maxload - minload) / maxload) * 100.;
        return fl;
    }

    double field_imbalance(const SpatialPartitionContainer& spc)
    {
        std::vector<double> load_global;

        boost::mpi::all_gather(*spc.comm(), field_load(spc.sp().front()), load_global);
        spc.comm()->barrier();

        auto result = std::minmax_element(load_global.begin(), load_global.end());
        double minload = *result.first;
        double maxload = *result.second;
        double fl = ((maxload - minload) / maxload) * 100.;
        return fl;
    }

    double imbalance(const SpatialPartitionContainer& spc)
    {
        {
            std::vector<double> load_global;

            boost::mpi::all_gather(*spc.comm(), solver_load(spc.sp().front()), load_global);
            spc.comm()->barrier();

            auto result = std::minmax_element(load_global.begin(), load_global.end());
            double minload = *result.first;
            double maxload = *result.second;
            //tl = ((maxload - minload) / maxload) * 100.;
            double cratio = maxload / minload;
            return cratio;
            //double tl = (maxload - minload);
        }
        //{
        //    std::vector<double> load_global;

        //    boost::mpi::all_gather(*spc.comm(), field_load(spc.sp().front()), load_global);
        //    spc.comm()->barrier();

        //    auto result = std::minmax_element(load_global.begin(), load_global.end());
        //    double minload = *result.first;
        //    double maxload = *result.second;
        //    fl = ((maxload - minload) / maxload) * 100.;
        //    //double fl = (maxload - minload);
        //}
    }

    void get_coef_(const Mesh& mesh, AeroCoef& local_coef, const AeroCoefPara& aero_para, const SpatialPartitionContainer* spc)
    {
        const auto& wall = mesh.wall_boundaries();

        double moment_length = aero_para.moment_length; 

        Vector3 F(0., 0., 0.);
        Vector3 M(0., 0., 0.);
        std::vector<std::tuple<Vector3, double>> P;
        P.reserve(wall.size());

        for (auto mc = wall.begin(); mc != wall.end(); ++mc)
        {
            const MeshCell& in = mesh.cell(mc->interior_boundary());

            assert(mc->face().size() == 1);
            const auto& mf = mc->face()[0];

            auto cnt = mf.face().centroid();

            Tag temp;
            if (!spc->is_resident(in.poly().centroid(), -1, temp)) {
                continue;
            }

            P.push_back(std::make_tuple(cnt, (mc->prim(4) - aero_para.p_ref) / aero_para.dpp()));
            F = F + mc->face()[0].face().normal() * mc->prim(4) * std::abs(mc->face()[0].face().signed_area());
            M = M + cross(cnt - moment_length, F);
        }

        local_coef.get_pressure_coef(P);
        local_coef.get_force_coef(F, aero_para);
        local_coef.get_moment_coef(M, aero_para);
    }

    void SpatialPartitionContainer::get_coef(const std::vector<AeroCoefPara>& aero_para, int iter) const
    {
        int mss = mesh_system_size();
        for (int i = 0; i < mss; ++i)
        {
            std::vector<double> local, global;
            local.resize(7, 0.);
            global.resize(7, 0.);

            AeroCoef local_coef;

            auto mesh = std::find_if(sp_.front().mesh().begin(), sp_.front().mesh().end(), [&](const auto& m){return m.tag()() == i;});
            if (mesh != sp_.front().mesh().end())
            {
                get_coef_(*mesh, local_coef, aero_para[i], this);
            }

            local[0] = local_coef.f(0);
            local[1] = local_coef.f(1);
            local[2] = local_coef.f(2);
            local[3] = local_coef.thrust;
            local[4] = local_coef.m(0);
            local[5] = local_coef.m(1);
            local[6] = local_coef.m(2);

            boost::mpi::all_reduce(*comm_, local.data(), 7, global.data(), std::plus<double>());

            if (comm_->rank() == 0)
            {
                std::string fn = "coef-";
                fn.append(std::to_string(i));
                fn.append(".dat");
                std::ofstream out;
                out.open(fn, std::fstream::app);

                out << iter; // mesh tag
                out << " ";
                out << global[0]; // cN
                out << " "; 
                out << global[1]; // cA
                out << " "; 
                out << global[2]; // cY
                out << " "; 
                out << global[3]; // cT
                out << " "; 
                out << global[4]; // cn
                out << " "; 
                out << global[5]; // cm
                out << " "; 
                out << global[6]; // cl
                out << std::endl; 
                out.close();

                if (mesh != sp_.front().mesh().end())
                {
                    std::string fn = "pres_coef-";
                    fn.append(std::to_string(comm_->rank()));
                    fn.append("-");
                    fn.append(std::to_string(i));
                    fn.append("-");
                    fn.append(std::to_string(iter));
                    fn.append(".dat");

                    std::ofstream out;
                    out.open(fn);

                    for (const auto& P: local_coef.p)
                    {
                        auto [cnt, p] = P;

                        out << cnt(0);
                        out << " "; 
                        out << cnt(1);
                        out << " "; 
                        out << cnt(2);
                        out << " "; 
                        out << p;
                        out << std::endl;
                    }

                    out.close();
                }
            }
        }
    }

    void SpatialPartitionContainer::make_global_adt()
    {
        make_adt(global_rm_, global_adt_, comm_->rank());

        assert(global_rm_.aabb().min(0) == global_adt_.aabb().min(0));
        assert(global_rm_.aabb().min(1) == global_adt_.aabb().min(1));
        assert(global_rm_.aabb().min(2) == global_adt_.aabb().min(2));
        assert(global_rm_.aabb().max(0) == global_adt_.aabb().max(0));
        assert(global_rm_.aabb().max(1) == global_adt_.aabb().max(1));
        assert(global_rm_.aabb().max(2) == global_adt_.aabb().max(2));
    }

    int SpatialPartitionContainer::mesh_system_size() const
    {
        int lmax_tag = TAILOR_BIG_NEG_NUM;
        for (const auto& sp: sp_)
        {
            for (const auto& m: sp.mesh())
            {
                lmax_tag = std::max(lmax_tag, m.tag()());
            }
        }

        int mesh_system_size = TAILOR_BIG_NEG_NUM;
        //std::cout << "mss: " << comm_->rank() << " " << lmax_tag << " " << mesh_system_size << std::endl;
        boost::mpi::all_reduce(*comm_, lmax_tag, mesh_system_size, boost::mpi::maximum<int>());
        return (mesh_system_size + 1); // because tags start from zero.
    }

    std::vector<AABB> get_hole_aabb(boost::mpi::communicator& comm, const SpatialPartitionContainer& spc)
    {
        int mesh_system_size = spc.mesh_system_size();

        std::vector<AABB> laabb(mesh_system_size);

        for (const auto& sp: spc.sp())
        {
            for (const auto& mesh: sp.mesh())
            {
                Vector3 min, max;
                minmax(mesh.wall_boundaries(), min, max);
                if (mesh.tag()() >= laabb.size())
                {
                    std::cout << "mesh tag: " << mesh.tag()()  << std::endl;
                    std::cout << "laabb size: " << laabb.size()  << std::endl;
                }
                assert(mesh.tag()() < laabb.size());
                laabb[mesh.tag()()] = AABB(min, max);
            }
        }

        std::vector<AABB> hole_aabb;
        boost::mpi::all_gather(comm, laabb.data(), laabb.size(), hole_aabb);

        //for (int i=0; i<mesh_system_size; ++i)
        //{
        //    for (int j=1; j<mesh_system_size; ++j)
        //    {
        //        hole_aabb[i].extend(hole_aabb[j * mesh_system_size]);
        //    }
        //}

        for (int i=0; i<mesh_system_size; ++i)
        {
            for (int j=1; j<comm.size(); ++j)
            {
                hole_aabb[i].extend(hole_aabb[j * mesh_system_size + i]);
            }
        }

        return hole_aabb;
    }

    const boost::mpi::communicator* SpatialPartitionContainer::comm() const
    {
        return comm_;
    }

    void SpatialPartitionContainer::sort_sp_meshes()
    {
        for (auto& sp: sp_)
        {
            sp.sort_mesh();
        }
    }

    void SpatialPartitionContainer::remove_ghosts()
    {
        for (auto& sp: sp_)
        {
            sp.remove_ghosts();
        }
    }
    void SpatialPartitionContainer::remove_nonresidents()
    {
        for (auto& sp: sp_)
        {
            sp.remove_nonresidents();
        }
    }
    void SpatialPartitionContainer::update_connectivity()
    {
        for (auto& sp: sp_)
        {
            sp.update_connectivity(comm_->rank());
        }
    }
    void SpatialPartitionContainer::reset_all_oga_status()
    {
        for (auto& sp: sp_)
        {
            sp.reset_all_oga_status();
        }
    }

    /*bool SpatialPartitionContainer::brmt_rank(int bintag, BinRMTag& brmt, int& rank) const
    {
        auto it = bintag_proc_map_.find(bintag);
        if (it != bintag_proc_map_.end())
        //for (const auto& pair: bintag_proc_map_)
        {
            //if (pair.first.bintag()() == bintag)
            //{
                //brmt = pair.first;
                //rank = pair.second;
                rank = it->second;
                return true;
            //}
        }

        return false;
    }*/

    std::map<Tag, int> SpatialPartitionContainer::search_rm(const AABB& aabb, const Vector3& cnt, std::pair<Tag, int>& resbin, bool verbose) const
    {
        auto adtpoint = ADTPoint(aabb.vertices().begin(), aabb.vertices().end());
        auto res = global_adt_.search(adtpoint, verbose);

        std::map<Tag, int> map;

        for (int r: res)
        {
            auto it = bintag_proc_map_.find(r);
            if (it != bintag_proc_map_.end())
            {
                auto pair = std::make_pair(Tag(it->first), it->second);
                map.insert(pair);
                assert(!global_rm_.bin().empty());
                if (global_rm_.bin(Tag(it->first)).aabb().do_intersect(cnt)) {
                    resbin = pair;
                }
            }

        }

        return map;
    }

    //std::map<Tag, int> SpatialPartitionContainer::search_rm(const AABB& aabb, bool verbose) const
    //{
    //    //std::vector<Point> pts;
    //    //pts.reserve(aabb.vertices().size());
    //    //for (const auto& v: aabb.vertices())
    //    //{
    //    //    pts.push_back(v);
    //    //}
    //    //auto adtpoint = ADTPoint(pts.begin(), pts.end());
    //    auto adtpoint = ADTPoint(aabb.vertices().begin(), aabb.vertices().end());
    //    auto res = global_adt_.search(adtpoint, verbose);

    //    //std::map<BinRMTag, int> map;
    //    std::map<Tag, int> map;

    //    for (int r: res)
    //    {
    //        //BinRMTag brmt;
    //        //int rank;
    //        //bool resu = brmt_rank(r, brmt, rank);
    //        //assert(resu);
    //        auto it = bintag_proc_map_.find(r);
    //        if (it != bintag_proc_map_.end()) {
    //            map.insert(std::make_pair(Tag(it->first), it->second));
    //            //map.insert(*it);
    //        }
    //    }

    //    return map;
    //}

    bool SpatialPartitionContainer::search_rm(const Vector3& cnt, Tag& bintag, int& rank, bool verbose) const
    {
        auto adtpoint = ADTPoint(cnt);
        std::vector<int> res = global_adt_.search(adtpoint, verbose);

        if (verbose) {
            std::cout << "res size: " << res.size() << std::endl;
        }

        if (res.empty()) {
            return false;
        }

        //assert(res.size() == 1);
        auto t = res[0];

        auto it = bintag_proc_map_.find(t);
        assert(it != bintag_proc_map_.end());
        //if (it != bintag_proc_map_.end()) {
            bintag = Tag(it->first);
            rank = it->second;
            assert(bintag.isvalid());
            return true;
        //}

        //return false;

        //return brmt_rank(res[0], brmt, rank);
    }

    //bool SpatialPartitionContainer::is_resident(const Vector3& cnt, int celltag, BinRMTag& brmt) const
    //bool SpatialPartitionContainer::is_resident(const Vector3& cnt, int celltag, Tag& brmt, int& dest_rank) const
    bool SpatialPartitionContainer::is_resident(const Vector3& cnt, int celltag, Tag& brmt) const
    {
        if (comm_->size() == 1) {
            return true;
        }
        // use if merge_bins is true.
        //std::vector<BinRMTag> ttag;
        //BinRMTag tag;
        //global_rm_.get_bintag_adaptive(cnt, tag, temp);
        //global_rm_.get_bintag_adaptive_unique(cnt, tag, comm_->rank(), celltag);
        //auto adtpoint = ADTPoint(cnt);
        //std::vector<int> res = global_adt_.search(adtpoint);
        //assert(res.size() == 1);
        int dest_rank;
        bool resu = search_rm(cnt, brmt, dest_rank); 
        if (!resu)
        {
            //if (comm_->rank() == 9)
            {
                std::cout << "cnt: " << cnt(0) << " " << cnt(1) << " " << cnt(2) << std::endl;
                std::cout << "min(0): " << global_adt_.aabb().min(0) << std::endl;
                std::cout << "min(1): " << global_adt_.aabb().min(1) << std::endl;
                std::cout << "min(2): " << global_adt_.aabb().min(2) << std::endl;
                std::cout << "max(0): " << global_adt_.aabb().max(0) << std::endl;
                std::cout << "max(1): " << global_adt_.aabb().max(1) << std::endl;
                std::cout << "max(2): " << global_adt_.aabb().max(2) << std::endl;
            }
        }
        assert(resu);
        //assert(tag.isvalid());

        /*for (const auto& brmt: tag)
          {
          auto it = bintag_proc_map_.find(brmt);
          assert(it != bintag_proc_map_.end());
          if (it->second == comm_->rank())
          {
          if (global_rm_.bin(it->first).aabb().do_intersect(cnt)) {
          return true;
          }
          }
          }*/

        //if (mt() == 1 && ct() == 91)
        //{
        //}

        //for (const auto& pair: bintag_proc_map_)
        //{
            //if (pair.first.bintag()() == res[0])
            //{
                //dest_rank = pair.second;
                //bin = pair.first.bintag()();
                //break;
            //}
        //}
        //assert(dest_rank != -1);

        //auto it = bintag_proc_map_.find(tag);
        //assert(it != bintag_proc_map_.end());
        if (dest_rank != comm_->rank()) {
        //if (it->second != comm_->rank()) {
            return false;
        }

        //bin = it->first.bintag();
        return true;
        //return global_rm_.bin(it->first).aabb().do_intersect(cnt); // do I need this? did not I already checked this with ...unique?
    }

    void SpatialPartitionContainer::connect_partition_cells(ArrCon<Nei>& arrival)
    {
        std::function<bool(const Vector3&, int)> is_resi = [this](const Vector3& cnt, int celltag)
        {
            Tag spt;
            return is_resident(cnt, celltag, spt);
        };

        for (SpatialPartition& sp: sp_)
        {
            bool found = false;
            //for (auto& stor: storage)
            {
                //if (sp.tag()() == stor.tag())
                {
                    //sp.connect_partition_cells(stor.arrival_p(), is_resi, profiler_);
                    sp.connect_partition_cells(arrival, is_resi, profiler_);
                    found = true;
                    break;
                }
            }
            assert(found);
        }
    }

    //void SpatialPartitionContainer::update_ghost_primitives(const std::vector<Storage<Var>>& stor)
    //void SpatialPartitionContainer::update_ghost_primitives(const RecvCon<Var>& stor)
    //{
        //assert(sp_.size() == 1);
        //assert(stor.size() == 1);
        //auto& sp = sp_.front();
        //auto& s = stor.front();
        //assert(sp.tag()() == s.tag());
        //sp.update_ghost_primitives(s.arrival());
    //}

    void SpatialPartitionContainer::oga_interpolate(const ArrCon<Var>& arrival)
    {
        assert(sp_.size() == 1);
        //assert(stor.size() == 1);
        auto& sp = sp_.front();
        //auto& s = stor.front();
        //assert(sp.tag()() == s.tag());
        sp.oga_interpolate(arrival);
    }

    void SpatialPartitionContainer::convert_to_fortran() const
    {
        for (const SpatialPartition& sp: sp_)
        {
            sp.convert_to_fortran();
        }
    }

    //void SpatialPartitionContainer::solve(Solver& solver, const Vector3& vel, std::function<void(CellExchanger& cell_exchanger)> exchange_ghosts)
    //void SpatialPartitionContainer::solve(Solver& solver, const Vector3& vel, std::function<void()> exchange_ghosts)
    //{
        //bool called_exc_ghosts;
        //for (int i=0; i<sp_.size(); ++i)
        //{
            //if (i == 0) {
                //called_exc_ghosts = false;
            //}
            //else {
                //called_exc_ghosts = true;
            //}
            //sp_[i].solve(solver, vel, exchange_ghosts, called_exc_ghosts);
        //}
    //}

    /*void SpatialPartitionContainer::init_interior()
    {
        for (SpatialPartition& sp: sp_)
        {
            sp.init_interior();
        }
    }

    void SpatialPartitionContainer::init_sod()
    {
        for (SpatialPartition& sp: sp_)
        {
            sp.init_sod();
        }
    }

    void SpatialPartitionContainer::init_dirichlet()
    {
        for (SpatialPartition& sp: sp_)
        {
            sp.init_dirichlet();
        }
    }*/

    void SpatialPartitionContainer::init()
    {
        for (SpatialPartition& sp: sp_)
        {
            sp.init();
        }
    }

    void SpatialPartitionContainer::init_sod()
    {
        for (SpatialPartition& sp: sp_)
        {
            sp.init_sod();
        }
    }

    void SpatialPartitionContainer::connect_after_exchange(Profiler* profiler, std::string proc)
    {
        std::function<bool(const Vector3&)> is_resi = [this](const Vector3& cnt)
        {
            Tag spt;
            return is_resident(cnt, -1, spt);
        };

        for (SpatialPartition& sp: sp_)
        {
            sp.connect_after_exchange(is_resi, profiler, proc);
        }
    }

    void SpatialPartitionContainer::update_prev_cand_donor()
    {
        for (SpatialPartition& sp: sp_)
        {
            sp.update_prev_cand_donor();
        }
    }

    bool duplicate_exist(const SpatialPartitionContainer& spc)
    {
        for (const SpatialPartition& sp: spc.sp())
        {
            for (const Mesh& mesh: sp.mesh())
            {
                for (const MeshCell& mc: mesh.cell())
                {
                    for (const MeshPoint& mp: mc.point())
                    {
                        assert(!std::isnan(mp.p().r(0)));
                        assert(!std::isnan(mp.p().r(1)));
                        assert(!std::isnan(mp.p().r(2)));
                    }
                }
                {
                    auto copy = mesh.cell();
                    std::sort(copy.begin(), copy.end(), [&](const MeshCell& left, const MeshCell& right){return left.tag() < right.tag();});
                    //copy.sort([&](const MeshCell& left, const MeshCell& right){return left.tag() < right.tag();});
                    auto it = std::adjacent_find(copy.begin(), copy.end());
                    if (it != copy.end())
                    {
                        return true;
                    }
                }

                {
                    auto copy = mesh.wall_boundaries();
                    std::sort(copy.begin(), copy.end(), [&](const MeshCell& left, const MeshCell& right){return left.tag() < right.tag();});
                    auto it = std::adjacent_find(copy.begin(), copy.end());
                    if (it != copy.end())
                    {
                        return true;
                    }
                }
            }
        }

        return false;
    }

    //const std::map<BinRMTag, int>& SpatialPartitionContainer::bintag_proc_map() const
    const std::map<int, int>& SpatialPartitionContainer::bintag_proc_map() const
    {
        return bintag_proc_map_;
    }

    const RegularMesh& SpatialPartitionContainer::global_rm() const
    {
        return global_rm_;
    }

    RegularMesh& SpatialPartitionContainer::global_rm()
    {
        return global_rm_;
    }

    const ADT& SpatialPartitionContainer::global_adt() const
    {
        return global_adt_;
    }

    void SpatialPartitionContainer::set_global_rm(const RegularMesh& rm)
    {
        global_rm_ = rm;
    }

    void SpatialPartitionContainer::make_mesh_aabbs()
    {
        for (SpatialPartition& s: sp_)
        {
            s.make_mesh_aabb();
        }
    }

    void SpatialPartitionContainer::merge_sp_meshes(std::deque<Mesh>& mesh) const
    {
        //assert(!sp_.empty());

        for (const SpatialPartition& s: sp_)
        {
            //s.merge_meshes_no_check(mesh);
            s.merge_meshes(mesh);
        }

        for (Mesh& m: mesh)
        {
            m.shrink_points();
            m.sort_cells();
            //m.set_cell_tag_index_map(); 
            //m.set_point_tag_index_map();
        }

        //for (Mesh& m: mesh)
        //{
            //m.remove_dup_cells_and_points();
        //}
    }

    void SpatialPartitionContainer::remap_uniproc(const std::deque<Mesh>& mesh)
    {
        SpatialPartition sp;
        sp.set_tag(Tag(0));

        AABB aabb;

        for (const Mesh& m: mesh)
        {
            sp.add_mesh(m);
            aabb.extend(AABB(m.rawpoint()));
        }

        aabb.set_faces();
        aabb.set_vertices_from_bbox();
        aabb.set_faces();

        sp.set_aabb(aabb);

        add_sp(sp);
    }

    void SpatialPartitionContainer::profstart(std::string fun)
    {
        if (profiler_ == nullptr) {
            return;
        }

        profiler_->start(fun);
    }

    void SpatialPartitionContainer::profstop(std::string fun)
    {
        if (profiler_ == nullptr) {
            return;
        }

        profiler_->stop(fun);
    }

    void SpatialPartitionContainer::remove_dup_cells_and_points()
    {
        for (SpatialPartition& s: sp_) {
            s.remove_dup_cells_and_points();
        }
    }

    void SpatialPartitionContainer::add_sp(const std::vector<std::vector<BinRMTag>>& dist, const ArrCon<SpatialPartition>& arrival)
    {
        //assert(!dist.empty());
        //for (int proc=0; proc<dist.size(); ++proc)
        {
            //if (proc != comm_->rank()) {
                //continue;
            //}

            //for (const auto& brmt: dist[proc])
            {
                //auto arrival = std::find_if(storage.begin(), storage.end(), [&](const auto& st){return st.tag() == brmt.bintag()();});
                //assert(arrival != storage.end());
                //}
                //for (int i=0; i<sp.size(); ++i)
                //if (stor->arrival().empty())
                //{
                    //std::cout << "bintag: " << brmt.bintag()() << std::endl;
                    //std::cout << "rank: " << comm_->rank() << std::endl;
                //}
                // sp might be inside hole. so arrival might be empty.
                //if (stor->arrival().empty()) {
                //if (arrival == storage.end()) {
                    //continue;
                //}

                //assert(!stor->arrival().empty());
                //for (int i=0; i<stor->arrival().size(); ++i)
                for (int i=0; i<arrival.size(); ++i)
                {
                    //const auto& other_sp = sp[i];
                    //const auto& other_sp = stor->arrival()[i];
                    const auto& other_sp = arrival[i];

                    auto exist_sp = std::find_if(sp_.begin(), sp_.end(), [&](const auto& s){return s.tag() == other_sp.tag();});
                    if (exist_sp == sp_.end())
                    {
                        sp_.push_back(other_sp);
                        sp_.back().set_comm(comm_); 
                        if (profiler_ != nullptr) {
                            sp_.back().set_profiler(profiler_);
                        }
                        //cond[i] = true;
                        assert(!sp_.back().aabb().faces().empty());
                    } 
                    else
                    {
                        exist_sp->merge(other_sp);
                        assert(!exist_sp->aabb().faces().empty());
                    }
                }
            }
        }
    }

    /*void SpatialPartitionContainer::push_sp(const std::vector<std::vector<BinRMTag>>& dist, const std::vector<Storage<SpatialPartition>>& storage, std::vector<bool>& cond)
    {
        // This function it time consuming for low number of meshes with high number of cells. Each time a cell is added to mesh, we iterate through points of mesh. So these monster meshes which usually reside in low-rank processors, cause time consumption.
        // I assumed that point insertion is the problem because we insert cells in batch with reserved space.
        // Besides point insertion, sorting of cells could also be a consumer. But sorting complexity is O(nlogn). I am not sure if partitioned sort takes same time as all in one sort. After all, n is nearly the same for all processors. 

        for (int proc=0; proc<dist.size(); ++proc)
        {
            if (proc != comm_->rank()) {
                continue;
            }

            for (const auto& brmt: dist[proc])
            {
                auto stor = std::find_if(storage.begin(), storage.end(), [&](const auto& st){return st.tag() == brmt.bintag()();});
                assert(stor != storage.end());

                //for (int i=0; i<sp.size(); ++i)
                for (int i=0; i<stor->arrival().size(); ++i)
                {
                    if (cond[i] == true) {
                        continue;
                    }
                    //auto& other_sp = sp[i];
                    const auto& other_sp = stor->arrival()[i];
                    auto exist_sp = std::find_if(sp_.begin(), sp_.end(), [&](const auto& s){return s.tag() == other_sp.tag();});

                    exist_sp->merge(other_sp);
                    assert(!exist_sp->aabb().faces().empty());

                    for (const Mesh& mesh: exist_sp->mesh())
                    {
                        for (const MeshCell& mc: mesh.cell())
                        {
                            for (const Tag& t: mc.empty_boundary())
                            {
                                auto it = std::find_if(mesh.empty_boundaries().begin(), mesh.empty_boundaries().end(), [&](const MeshCell& tmc){return tmc.tag() == t;});
                                assert(it != mesh.empty_boundaries().end());
                            }
                        }
                    }
                }
            }
        }
    }*/

    /*void SpatialPartitionContainer::make_outline(std::vector<std::vector<Vector3>> allpts)
    {
        Outline outline;
        if (!master())
        {
            outline.set_tag(world_.rank());
            if (!allpts.empty())
            {
                outline.build(allpts, world_.rank());
            }
        }

        boost::mpi::all_gather(world_, outline, outline_);

        if (!master())
        {
            outline_.erase(outline_.begin());
            for (Outline& ol: outline_) {
                ol.build_polygon(world_.rank());
            }
        }
    }*/

    void SpatialPartitionContainer::reduce_sp(const std::vector<std::vector<BinRMTag>>& dist, const ArrCon<SpatialPartition>& incoming, std::string name)
    {
        //assert(dist.size() == comm_->size());
        //for (int proc=0; proc<dist.size(); ++proc)
        {
            //if (proc != comm_->rank()) {
                //continue;
            //}

            //for (const auto& brmt: dist[proc])
            {
                //auto stor = std::find_if(storage.begin(), storage.end(), [&](const auto& st){return st.tag() == brmt.bintag()();});
                //auto incoming = std::find_if(storage.begin(), storage.end(), [&](const auto& st){return st.tag() == brmt.bintag()();});
                //assert(stor != storage.end());
                //assert(incoming != storage.end());
                //auto incoming = stor->arrival();

                assert(!incoming.empty());

                //if (!incoming.empty())
                {
                    auto inc = incoming.begin();

                    //if (profiler_ != nullptr) {profiler_->start(name+"-reduce-empty");}
                    if (sp_.empty())
                    {
                        sp_.push_back(incoming[0]);
                        sp_.back().set_comm(comm_);
                        sp_.back().set_tag(Tag(comm_->rank())); // I want sp tags to be unique to avoid problems with mpi.
                        inc = std::next(inc);
                    }
                    //if (profiler_ != nullptr) {profiler_->stop(name+"-reduce-empty");}

                    //if (profiler_ != nullptr) {profiler_->start(name+"-reduce-merge");}
                    for (; inc != incoming.end(); ++inc)
                    {
                        assert(sp_.size() > 0);
                        sp_[0].add_merge_mesh_leave_dups(*inc);
                        sp_[0].extend_aabb(inc->aabb());
                    }
                    //if (profiler_ != nullptr) {profiler_->stop(name+"-reduce-merge");}

                    //if (profiler_ != nullptr) {profiler_->start(name+"-reduce-dup");}
                    sp_[0].remove_dups(name, nullptr);
                    //if (profiler_ != nullptr) {profiler_->stop(name+"-reduce-dup");}
                }
            }
        }
    }

    /*void SpatialPartitionContainer::reduce_sp()
    {
        if (master()) {
            return;
        }

        assert(!sp_.empty());

        SpatialPartition sp = sp_[0];

        for (int i=1; i<sp_.size(); ++i)
        {
            sp.set_tag(sp_[i].tag()); // otherwise merge_mesh() won't let the merge.
            sp.merge_mesh(sp_[i]);
            sp.extend_aabb(sp_[i].aabb());
        }

        sp_.clear();
        sp.set_tag(Tag(0));
        sp_.push_back(std::move(sp));
    }*/

    void SpatialPartitionContainer::merge_sp(const std::deque<SpatialPartition>& sp)
    {
        for (const auto& other_sp: sp)
        {
            auto exist_sp = std::find_if(sp_.begin(), sp_.end(), [&](const auto& s){return s.tag() == other_sp.tag();});
            if (exist_sp != sp_.end())
            {
                exist_sp->merge_mesh(other_sp);
            } 
        }
    }

    /*void SpatialPartitionContainer::merge_sp(const std::deque<SpatialPartition>& sp)
    {
        for (const auto& other_sp: sp)
        {
            auto exist_sp = std::find_if(sp_.begin(), sp_.end(), [&](const auto& s){return s.tag() == other_sp.tag();});
            if (exist_sp != sp_.end())
            {
                exist_sp->merge(other_sp);
            } 
        }
    }*/

    void SpatialPartitionContainer::add_sp(SpatialPartition& other_sp)
    {
        auto exist_sp = std::find_if(sp_.begin(), sp_.end(), [&](const auto& s){return s.tag() == other_sp.tag();});
        if (exist_sp == sp_.end())
        {
            sp_.push_back(other_sp);
            sp_.back().set_comm(comm_); 
            if (profiler_ != nullptr) {
                sp_.back().set_profiler(profiler_);
            }
        } 
        else
        {
            exist_sp->merge(other_sp);
        }

    }

    double SpatialPartitionContainer::dev() const
    {
        //assert(master());
        return dev_;
    }

    /*void SpatialPartitionContainer::output_load(int iteration)
    {
        std::vector<double> _loads;
        gather(world_, load_, _loads, MASTER);

        if (!master()) {
            return;
        }

        double dev;
        double dev_ave;
        auto result = std::minmax_element(_loads.begin()+1, _loads.end());
        double minload = *result.first;
        double maxload = *result.second;

        double average = std::accumulate(_loads.begin()+1, _loads.end(), 0.) / (_loads.size() - 1);
        dev_ave = std::max(maxload-average, average-minload);
        dev_ave *= 100. / average;

        //dev = static_cast<double>(maxload - minload) * 100. / static_cast<double>(maxload);
        dev = (maxload - minload) * 100. / maxload;
        //for (double ll: _loads)
        //{
            //std::cout << "spc load: " << ll << std::endl;
        //}
        dev_ = dev;

        std::fstream out;
        out.open("dev-vs-iter.csv", std::fstream::out | std::fstream::app);
        if (iteration == 0)
            out << "iter,dev\n";
        out << iteration << "," << dev << "\n";
        out.close();

        out.open("devave-vs-iter.csv", std::fstream::out | std::fstream::app);
        if (iteration == 0)
            out << "iter,dev\n";
        out << iteration << "," << dev_ave << "\n";
        out.close();

        out.open("load_iter.csv", std::fstream::out | std::fstream::app);
        if (iteration == 0)
        {
            out << "iter,";
            for (int i=1; i<world_.size(); ++i)
            {
                out << "proc" << i;
                if (i != world_.size() - 1)
                    out << ",";
            }
            out << "\n";
        }
        out << iteration << ",";
        for (int i=1; i<world_.size(); ++i)
        {
            out << _loads[i];
            if (i != world_.size() - 1)
                out << ",";
        }
        out << "\n";
        out.close();
    }*/

    SpatialPartitionContainer::SpatialPartitionContainer(): verbose_(false), nmesh_(0), profiler_(nullptr), comm_(nullptr), tl_(0.), fl_(0.), dev_(0.)
    {
    }

    //SpatialPartitionContainer::SpatialPartitionContainer(boost::mpi::communicator* comm, Profiler* profiler, bool verbose, const std::map<BinRMTag, int>& bintag_proc_map, int nmesh): comm_(comm), profiler_(profiler), verbose_(verbose)
    SpatialPartitionContainer::SpatialPartitionContainer(boost::mpi::communicator* comm, Profiler* profiler, bool verbose, const std::map<int, int>& bintag_proc_map, int nmesh): comm_(comm), profiler_(profiler), verbose_(verbose), tl_(0.), fl_(0.), dev_(0.)
    {
        /*if (master())
        {
            for (const Bin& b: lm.rm().bin())
            {
                if (!b.mesh_tag_index_map().left.empty())
                {
                    assert(!b.cell().empty());
                }
            }
        }*/

        //global_nstripe_ = lm.rm().nstripe();
        //global_steplength_ = lm.rm().h();
        //global_aabb_min_ = lm.rm().aabb().min();
        bintag_proc_map_ = bintag_proc_map;
        nmesh_ = nmesh;
        //send_proc_win = MPI_WIN_NULL;
        /*MPI_Win_allocate(world_.size()*sizeof(int), sizeof(int), MPI_INFO_NULL, world_, &send_size, &send_size_win);
        MPI_Win_allocate(world_.size()*sizeof(int), sizeof(int), MPI_INFO_NULL, world_, &send_proc_2, &send_proc_win);
        MPI_Win_allocate(world_.size()*sizeof(int), sizeof(int), MPI_INFO_NULL, world_, &send_proc_wall_, &send_proc_win_wall_);
        MPI_Win_allocate(world_.size()*sizeof(int), sizeof(int), MPI_INFO_NULL, world_, &send_proc_outer_, &send_proc_win_outer_);*/
        //boost::mpi::group world_group = world_.group();
        //std::vector<int> excluded_proc = {0};
        //boost::mpi::group worker_group = world_group.exclude(excluded_proc.begin(), excluded_proc.end());
        //worker_comm = boost::mpi::communicator(world_, worker_group);
    }

    SpatialPartitionContainer::SpatialPartitionContainer(boost::mpi::communicator* comm, Profiler* profiler, bool verbose, int nmesh): comm_(comm), profiler_(profiler), verbose_(verbose), nmesh_(nmesh), tl_(0.), fl_(0.), dev_(0.)
    {
        //global_nstripe_ = lm.rm().nstripe();
        //global_steplength_ = lm.rm().h();
        //global_aabb_min_ = lm.rm().aabb().min();
        //bintag_proc_map_ = lm.bintag_proc_map();
        //nmesh_ = lm.nmesh();
        //send_proc_win = MPI_WIN_NULL;
        //send_size_win = MPI_WIN_NULL;
        //MPI_Win_allocate(world_.size()*sizeof(int), sizeof(int), MPI_INFO_NULL, world_, &send_proc_2, &send_proc_win);
        //MPI_Win_allocate(world_.size()*sizeof(int), sizeof(int), MPI_INFO_NULL, world_, &send_proc_wall_, &send_proc_win_wall_);
        //MPI_Win_allocate(world_.size()*sizeof(int), sizeof(int), MPI_INFO_NULL, world_, &send_proc_outer_, &send_proc_win_outer_);
        //boost::mpi::group world_group = world_.group();
        //std::vector<int> excluded_proc = {0};
        //boost::mpi::group worker_group = world_group.exclude(excluded_proc.begin(), excluded_proc.end());
        //worker_comm = boost::mpi::communicator(world_, worker_group);
    }

    /*void SpatialPartitionContainer::info() const
    {
        if (!master())
        {
            if (verbose_) {
                std::cout << "rank: " << comm_->rank() << " | spc has " << sp_.size() << " sp." << std::endl;
            }

            for (const SpatialPartition& _sp: sp_)
            {
                _sp.info();
            }
        }
    }*/

    const std::vector<SpatialPartition>& SpatialPartitionContainer::sp() const
    {
        return sp_;
    }


    //SpatialPartitionContainer::~SpatialPartitionContainer()
    //{
        /*if (send_proc_win != MPI_WIN_NULL)
        {
            MPI_Win_free(&send_proc_win);
            MPI_Win_free(&send_proc_win_wall_);
            MPI_Win_free(&send_proc_win_outer_);
            MPI_Win_free(&send_size_win);
        }*/
    //}


    void SpatialPartitionContainer::rotate_meshblocks(const Tag& _parent_mesh, double ang, int axis, const Vector3& rot_axis)
    {
        // this is a temporary function to mimic solver.
        // aim is to displace certain meshblocks.

            //std::cout << rot_axis(0) << " " << rot_axis(1) << " " << rot_axis(2) << std::endl;
            assert(!std::isnan(rot_axis(0)));
            assert(!std::isnan(rot_axis(1)));
            assert(!std::isnan(rot_axis(2)));

        for (SpatialPartition& _sp: sp_)
        {
            _sp.rotate_mesh(_parent_mesh, ang, axis, rot_axis);
        }
    }

    void SpatialPartitionContainer::move_meshblocks(const Tag& _parent_mesh, const Vector3& v)
    {
        // this is a temporary function to mimic solver.
        // aim is to displace certain meshblocks.

        for (SpatialPartition& _sp: sp_)
        {
            _sp.move_mesh(_parent_mesh, v);
        }
    }

    void SpatialPartitionContainer::connect_cells(std::string name)
    {
        std::function<bool(const Vector3&)> is_resi = [this](const Vector3& cnt)
        {
            Tag spt;
            return is_resident(cnt, -1, spt);
        };

        for (SpatialPartition& _sp: sp_) {
            _sp.connect_cells(is_resi, profiler_, name);
        }
    }

    /*bool SpatialPartitionContainer::master() const 
    {
        if (comm_->rank() == MASTER)
        {
            return true;
        }

        return false;
    }*/

    void SpatialPartitionContainer::make_regular_maps_uniproc()
    {
        for (SpatialPartition& _sp: sp_)
        {
            _sp.make_regular_maps_for_mesh_uniproc(profiler_);
        }
    }

    void SpatialPartitionContainer::make_regular_maps()
    {
        for (SpatialPartition& _sp: sp_)
        {
            _sp.make_regular_maps_for_mesh();
            //for (const auto& _rm: _sp.rm())
            //{
                //std::cout << "rank: " << world_.rank() << " rm size: " << _rm.bin().size() << std::endl;
            //}
        }
    }

    /*void SpatialPartitionContainer::print_all_meshes_in_partitions_uniproc()
    {
        int nmesh = 0;
        for (SpatialPartition& _sp: sp_)
        {
            nmesh += _sp.mesh().size();
            for (int i=0; i<_sp.mesh().size(); ++i)
            {
                std::string s = "spc_rank_";
                s.append(std::to_string(comm_->rank()));
                s.append("_sp_");
                s.append(std::to_string(&_sp - sp_.data()));
                s.append("_mb_");
                //s.append(std::to_string(i));
                s.append(std::to_string(_sp.mesh()[i].tag()()));
                s.append(".vtk");
                _sp.mesh()[i].print_as_vtk(s);
            }
        }
    }*/

    /*void SpatialPartitionContainer::print_all_meshes_in_partitions()
    {
        if (master()) return;

        int nmesh = 0;
        for (SpatialPartition& _sp: sp_)
        {
            nmesh += _sp.mesh().size();
            for (int i=0; i<_sp.mesh().size(); ++i)
            {
                std::string s = "spc_rank_";
                s.append(std::to_string(comm_->rank()));
                s.append("_sp_");
                s.append(std::to_string(_sp.tag()()));
                s.append("_mb_");
                s.append(std::to_string(_sp.mesh()[i].tag()()));
                s.append(".vtk");
                _sp.mesh()[i].print_as_vtk(s);

                s = "spc_wall_rank_";
                s.append(std::to_string(comm_->rank()));
                s.append("_sp_");
                s.append(std::to_string(_sp.tag()()));
                s.append("_mb_");
                s.append(std::to_string(_sp.mesh()[i].tag()()));
                s.append(".vtk");
                _sp.mesh()[i].print_wall_as_vtk(s);
            }
        }
    }*/
}
