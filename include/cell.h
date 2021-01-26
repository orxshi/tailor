#ifndef CELL_H
#define CELL_H

#include "meshcell.h"

namespace Tailor
{
    struct Cell
    {
        int source_rank_;
        int source_tag_;
        MeshCell cell_;
        std::list<MeshCell> wall_;
        std::list<MeshCell> dirichlet_;
        std::list<MeshCell> farfield_;
        std::list<MeshCell> empty_;
        std::list<MeshCell> interog_;
        bool overlap_;

        Cell(int source_rank, int source_tag, const MeshCell& cell, std::list<MeshCell>&& wall, std::list<MeshCell>&& dirichlet, std::list<MeshCell>&& farfield, std::list<MeshCell>&& empty, std::list<MeshCell>&& interog);
        Cell(): source_rank_(-1), source_tag_(-1), overlap_(false) {}

        const MeshCell& cell() const;

        template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            ar & source_rank_;
            ar & source_tag_;
            ar & cell_;
            ar & wall_;
            ar & farfield_;
            ar & dirichlet_;
            ar & interog_;
            ar & empty_;
        }

        friend class boost::serialization::access;
    };
}

#endif
