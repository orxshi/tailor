#ifndef VAR_H
#define VAR_H

#include "array.h"

namespace Tailor
{
    struct Var
    {
        int source_rank_;
        int source_tag_;
        std::pair<int, int> mesh_cell_;
        vararray var_;

        Var(int source_rank, int source_tag, const vararray& var, int mesh, int cell);
        Var() = default; // for serialization (easy way)

        template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            ar & source_rank_;
            ar & source_tag_;
            ar & var_;
            ar & mesh_cell_;
        }

        friend class boost::serialization::access;
    };
}

#endif
