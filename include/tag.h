#ifndef TAG_H
#define	TAG_H

#include <boost/serialization/access.hpp>
#include <cassert>
#include <iostream>
//#include <pybind11/pybind11.h>

//namespace py = pybind11;

namespace Tailor
{
    class Tag
    {
        int tag_;
        friend class boost::serialization::access;

        public:

        explicit Tag(int tag): tag_(tag) {};
        Tag(): Tag(-1) {};
        operator int() const = delete;
        int operator()() const;
        //Tag& operator=(const Tag& t);
        bool operator==(const Tag& other) const;
        bool operator<(const Tag& other) const;
        bool operator!=(const Tag& other) const;
        bool isvalid() const;
        void set(int t);
        static bool sorter(const Tag& lhs, const Tag& rhs);
        template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            ar & tag_;
        }
    };

    class BinRMTag
    {
        Tag bintag_;
        //boost::shared_ptr<RegularMesh> rm_;
        //RegularMesh* rm_;
        Tag rmtag_;
        friend class boost::serialization::access;

        public:

        BinRMTag() = default;
        //BinTag(const Tag& t, RegularMesh* r);
        BinRMTag(const Tag& bintag, const Tag& rmtag);

        bool isvalid() const;
        static bool sorter(const BinRMTag& lhs, const BinRMTag& rhs);
        bool operator<(const BinRMTag& other) const;
        bool operator==(const BinRMTag& other) const;
        const Tag& bintag() const;
        const Tag& rmtag() const;
        //boost::shared_ptr<RegularMesh> rm() const;
        //RegularMesh* rm() const;
        void set_bin_tag(const Tag& t);
        void set_rm_tag(const Tag& t);
        template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            ar & bintag_;
            ar & rmtag_;
        }
    };

}

#endif
