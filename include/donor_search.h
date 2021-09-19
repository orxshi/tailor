#ifndef DONOR_SEARCH_H
#define	DONOR_SEARCH_H

#include "spc.h"
//#include "hole_profiler.h"
//#include "direct_cut.h"
#include "holemap.h"
#include "donor_info2.h"
#include "tcr_exchanger.h"

namespace Tailor
{
    enum class DonorSearchAlgo
    {
        mixed = 0,
        alladt = 1,
    };

    class DonorSearcher
    {
        public:

            DonorSearcher(boost::mpi::communicator* comm_, SpatialPartitionContainer* spc, Profiler* profiler, DonorSearchAlgo dsa, bool mergebins, bool verbose, bool print_donor_info);
            void donor_search(int nassemble, bool pseudo3D);
            void handle_donor_conflict(const ArrCon<DonorInfo2>& nonresi);
            void determine_best_donor(const ArrCon<DonorInfo2>& nonresi);
            void convert_orphan_to_field(const ArrCon<DonorInfo2>& nonresi);
            void convert_undefined_to_field();
            void receptor_to_field(const ArrCon<DonorInfo2>& nonresi);
            void determine_orphan(const ArrCon<DonorInfo2>& nonresi);
            void convert_receptor_to_hole();
            void check_donor_validity();
            void increase_overlap_thickness(int nlayer);

        private:

            bool print_info_;
            bool mergebins_;
            DonorSearchAlgo dsa_;
            bool verbose_;
            bool uniproc_;
            Profiler* profiler_;
            boost::mpi::communicator* comm_;
            SpatialPartitionContainer* spc_;

            //void apply_stencil_walk_to_bincells(const Bin& active_bin, const RegularMesh& active_rm, const RegularMesh& passive_rm, Mesh& active_mesh, const Mesh& passive_mesh, /*const ADT& passive_wall_adt, const ADT& passive_dirichlet_adt,*/ HoleCutType hole_cut_type, size_t& nsw, int dummyrank, const std::deque<Mesh>& wall, const std::deque<Mesh>& dirichlet, const HoleMap& holemap, bool use_past_seeds, int& totaliter, int& nseed);
            //bool master() const;
            void read_priority();
            //void donor_search_rm(SpatialPartition& sp, const std::deque<Mesh>& wall, const std::deque<Mesh>& dirichlet, const std::deque<HoleMap>& holemap, size_t& nsw, std::vector<double>& ds, int& totaliter, int& nseed);
            void donor_search_mixed(SpatialPartition& sp, const std::deque<HoleMap>& holemap, int& adthit, int nassemble, int& adtsearch, int& swsearch, int& min_iter, int& max_iter, const std::vector<AABB>& hole_aabb);
            //void donor_search_adt(SpatialPartition& sp, const std::deque<HoleMap>& holemap, size_t& adtmiss, size_t& adthit, std::vector<double>& ds);
            void convert_undefined_to_nonresident(SpatialPartition& sp);
            void reset_oga_status(SpatialPartition& sp);
            //std::deque<ADT> make_adt(SpatialPartition& sp);
            //std::deque<AABB> make_aabb(SpatialPartition& sp);
            //void handle_wall_contained_meshes(SpatialPartition& sp, const std::deque<HoleMap>& holemap);
            void handle_wall_contained_meshes(SpatialPartition& sp, const std::deque<HoleMap>& holemap, const std::vector<AABB>& hole_aabb);
            void donor_search_adt_mreceptor(SpatialPartition& sp, const std::deque<HoleMap>& holemap, const std::vector<AABB>& hole_aabb);
            void determine_nonresidents(SpatialPartition& sp);
    };

    ADT make_cell_adt(const Mesh& mesh);
}

#endif
