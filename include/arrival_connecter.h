#ifndef ARRIVAL_CONNECTER_H
#define ARRIVAL_CONNECTER_H

// isn't it better to put these into connectivity.cpp?

#include "cell_exchanger.h"

namespace Tailor
{
    class ArrivalConnecter
    {
        public:
            //ArrivalConnecter(SpatialPartitionContainer* spc, std::vector<Storage<Cell>>* storage, Profiler* profiler);
            ArrivalConnecter(SpatialPartitionContainer* spc, ArrCon<Cell>* arrival, Profiler* profiler);
            void add_arrival_cells();

        private:
            //std::vector<Storage<Cell>>* storage_;
            ArrCon<Cell>* arrival_;
            Profiler* profiler_;
            SpatialPartitionContainer* spc_;
            std::vector<bool> erase_interior_;
            //std::vector<bool> erase_wall_;
            //std::vector<bool> erase_dirichlet_;
            //std::vector<std::vector<MeshCell>> arr_interior_;
            //std::vector<std::vector<MeshCell>> arr_wall_;
            //std::vector<std::vector<MeshCell>> arr_dirichlet_;
            //std::vector<std::vector<MeshCell>> arr_farfield_;
            //std::vector<std::vector<MeshCell>> arr_empty_;
            //std::vector<std::vector<MeshCell>> arr_interog_;
            std::vector<mcc> arr_interior_;
            std::vector<mcc> arr_wall_;
            std::vector<mcc> arr_symmetry_;
            std::vector<mcc> arr_dirichlet_;
            std::vector<mcc> arr_farfield_;
            std::vector<mcc> arr_empty_;
            std::vector<mcc> arr_interog_;

            void add_arrival_cells(SpatialPartition& sp);
            void add_arrival_bou(SpatialPartition& sp, BouType bou);
            //void determine_added_cells(SpatialPartition& sp, Storage<Cell>&);
            void determine_added_cells(SpatialPartition& sp, ArrCon<Cell>&);
            //void determine_added_bou(SpatialPartition& sp, const Storage<Cell>&, std::string bou);
            void determine_added_bou(SpatialPartition& sp, ArrCon<Cell>& storage, BouType bou);
            void check_cell_exist2(const Mesh& m, const MeshCell& ac, int i);
    };
}

#endif
