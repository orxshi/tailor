#ifndef TAILOR_DONOR_H
#define TAILOR_DONOR_H

#include "tag.h"

namespace Tailor
{
    class MeshCell;

    struct Donor
    {
        Tag mesh_tag_;
        Tag cell_tag_;
        const MeshCell* addr_; // to retreive centroid. needed when mapping donor info from assembler partition to solver's. meaningful only in assembler where overlapping meshes are grouped together.

        Donor(const Tag& mesh_tag, const Tag& cell_tag, const MeshCell* addr);
        Donor();

        /*friend class boost::serialization::access;
        template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            ar & mesh_tag_;
            ar & cell_tag_;
            ar & addr_;
        }*/
    };
}

#endif
