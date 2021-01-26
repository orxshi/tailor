#ifndef FACETAG_H
#define	FACETAG_H

#include <boost/serialization/access.hpp>
#include <cassert>
#include <iostream>
#include <vector>
#include <algorithm>

namespace Tailor
{
    class FaceTag
    {
        std::vector<int> tag_;
        friend class boost::serialization::access;

        public:

        FaceTag(std::vector<int> pts);
        FaceTag() = default;
        //int operator()() const;
        bool operator==(const FaceTag& other) const;
        bool operator<(const FaceTag& other) const;
        bool operator!=(const FaceTag& other) const;
        const std::vector<int>& operator()() const;
        bool isvalid() const;
        void set(std::vector<int> pts);
        template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            ar & tag_;
        }
    };

    std::ostream& operator<<(std::ostream& os, const FaceTag& res);
}

#endif
