#ifndef NEI_H
#define NEI_H

#include "meshcell.h"

namespace Tailor
{
    struct Nei
    {
        int source_rank_;
        int source_tag_;
        //std::vector<MeshCell> cell_;
        MeshCell cell_;

        Nei(int source_rank, int source_tag, const MeshCell& cell);
        Nei() = default; // for serialization (easy way)

        //const std::vector<MeshCell>& cell() const;
        const MeshCell& cell() const;
        //void remove_dup();
        bool operator==(const Nei& other) const;
        bool operator<(const Nei& other) const;

        template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            ar & source_rank_;
            ar & source_tag_;
            ar & cell_;
        }

        friend class boost::serialization::access;
    };
}

#endif
