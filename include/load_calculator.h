#ifndef LOAD_CALCULATOR_H
#define LOAD_CALCULATOR_H

//#include "spc.h"
#include "regular_mesh.h"
#include <boost/program_options.hpp>
//#include <CGAL/Polyline_simplification_2/simplify.h>
//#include <CGAL/convex_hull_2.h>

//namespace PS = CGAL::Polyline_simplification_2;
//typedef PS::Stop_below_count_ratio_threshold StopRatio;
//typedef PS::Stop_below_count_threshold StopBelow;
//typedef PS::Squared_distance_cost            Cost;

namespace Tailor
{
    enum class AreaRep
    {
        aabb = 0,
        concave = 1,
        convex = 2,
    };

    enum class LoadEstim
    {
        area = 0,
        hybrid = 1,
        solver = 2,
        type1 = 3,
    };

    class LoadCalculator
    {
        public:

            LoadCalculator(LoadEstim type, boost::mpi::communicator* comm_);
            //bool load(const SpatialPartitionContainer& spc, double refine_tol, int iteration) const;
            //std::vector<BinRMTag> sorted_bin_tags(const RegularMesh& rm, std::vector<double>& load_global, std::vector<double>& nonsorted_load_global, std::vector<double>& load_global_r, const std::deque<Mesh>& mesh);
            std::vector<BinRMTag> sorted_bin_tags(const RegularMesh& rm, std::vector<double>& load_global, std::vector<double>& nonsorted_load_global, const std::deque<Mesh>& mesh);
            //void gather_load(const RegularMesh& rm, std::vector<double>& load_global, std::vector<double>& load_global_r, const std::deque<Mesh>& mesh);
            void gather_load(const RegularMesh& rm, std::vector<double>& load_global, const std::deque<Mesh>& mesh);
            double load(const Bin& bin, const std::deque<Mesh>& mesh) const;
            double load_r(const Bin& bin, const std::deque<Mesh>& mesh) const;

        private:
            AreaRep arearep_;
            LoadEstim type_;
            boost::mpi::communicator* comm_;

            void read_settings();
            //double area_load_alpha(const Bin& bin, const std::deque<Mesh>& mesh) const;
            //double area_load_convex(const Bin& bin, const std::deque<Mesh>& mesh) const;
            double area_load(const Bin& bin, const std::deque<Mesh>& mesh) const;
            double type1_load(const Bin& bin, const std::deque<Mesh>& mesh) const;
            //double area_load(const SpatialPartition& sp, const Outline& outline) const;
            //double area_load(const SpatialPartition& sp) const;
            //double type1_load(const SpatialPartition& sp) const;
            //double area_load_alpha(const SpatialPartition& sp, const Outline& outline) const;
            //double area_load_alpha(const SpatialPartition& sp) const;
            //double area_load_convex(const SpatialPartition& sp, const Outline& outline) const;
            //double area_load_convex(const SpatialPartition& sp) const;
            double hybrid_load(const Bin& bin) const;
            //double hybrid_load(const SpatialPartition& sp, const Outline& outline) const;
            //double hybrid_load(const SpatialPartition& sp) const;
            double solver_load(const Bin& bin) const;
            //double solver_load(const SpatialPartition& sp, const Outline& outline) const;
            //double solver_load(const SpatialPartition& sp) const;
            double minmesh_load(const Bin& bin) const;
            //double minmesh_load(const SpatialPartition& sp, const Outline& outline) const;
            //double load(const SpatialPartition& sp, const Outline& outline) const;
            //double load(const SpatialPartition& sp) const;
            //void gather_load(const RegularMesh& rm, std::vector<double>& load_local, std::vector<double>& load_local_r, int& j, const std::deque<Mesh>& mesh);
            void gather_load(const RegularMesh& rm, std::vector<double>& load_local, int& j, const std::deque<Mesh>& mesh);
    };

    //void print_regmesh(std::string file_name, const std::map<BinRMTag, int>& bintag_proc_map, const LoadCalculator& lc, const RegularMesh& rm, const std::deque<Mesh>& mesh);
    void print_regmesh(std::string file_name, const std::map<int, int>& bintag_proc_map, const LoadCalculator& lc, const RegularMesh& rm, const std::deque<Mesh>& mesh);
}

#endif
