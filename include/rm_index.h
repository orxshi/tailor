#ifndef REGULAR_MESH_INDEXH
#define REGULAR_MESH_INDEXH

#include <boost/mpi.hpp>

namespace Tailor
{
    class RegularMeshIndex
    {
        int i_, j_, k_;
        friend class boost::serialization::access;

        public:

        RegularMeshIndex(int i, int j, int k);
        RegularMeshIndex(): RegularMeshIndex(-1, -1, -1) {}

        bool isvalid() const;
        const int i() const;
        const int j() const;
        const int k() const;
        void set(int i, int j, int k);
        void inci();
        void incj();
        void inck();
        void deci();
        void decj();
        void deck();
        static bool sorter(const RegularMeshIndex& lhs, const RegularMeshIndex& rhs);
        bool operator==(const RegularMeshIndex& other) const;
        RegularMeshIndex operator+(const RegularMeshIndex& other) const;
        //bool operator!=(const RegularMeshIndex& other) const;
        template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            ar & i_;
            ar & j_;
            ar & k_;
        }
    };
}

#endif
