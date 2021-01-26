#ifndef DONOR_INFO2_H
#define DONOR_INFO2_H

#include <boost/mpi.hpp>

namespace Tailor
{
    // DonorInfo2 is used to transfer donor info from assembler to solver.
    struct DonorInfo2
    {
        int source_rank_;
        int source_tag_;
        int meshtag_;
        int celltag_;
        OGA_cell_type_t oga_cell_type_;
        int donor_mesh_;
        int donor_cell_;
        int donor_rank_; // to be used by solver to retreive info from donor_rank_ in solver partition.

        DonorInfo2(int meshtag, int celltag, OGA_cell_type_t type, int donor_mesh, int donor_cell, int donor_rank, int source_rank, int source_tag);
        DonorInfo2(int meshtag, int celltag, OGA_cell_type_t type, int source_rank, int source_tag);
        DonorInfo2() = default;

        template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            ar & meshtag_;
            ar & celltag_;
            ar & oga_cell_type_;
            ar & donor_mesh_;
            ar & donor_cell_;
            ar & donor_rank_;
            ar & source_rank_;
            ar & source_tag_;
        }

        friend class boost::serialization::access;
    };
}

#endif
