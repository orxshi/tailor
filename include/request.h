#ifndef REQUEST_H
#define REQUEST_H

#include <boost/mpi.hpp>

namespace Tailor
{
    struct Request
    {
        int source_rank_;
        int source_tag_;
        //std::vector<std::pair<int, int>> mesh_cell_;
        std::pair<int, int> mesh_cell_;

        Request(int mesh, int cell, int source_rank, int source_tag);
        Request() = default; // for serialization (easy way)

        template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            ar & source_rank_;
            ar & source_tag_;
            ar & mesh_cell_;
        }

        friend class boost::serialization::access;
    };
}

#endif
