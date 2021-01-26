#include "facetag.h"

namespace Tailor
{
    FaceTag::FaceTag(std::vector<int> pts)
    {
        set(pts);
    }

    bool FaceTag::isvalid() const
    {
        if (tag_.empty())
        {
            return false;
        }

        return true;
    }

    void FaceTag::set(std::vector<int> pts)
    {
        std::sort(pts.begin(), pts.end());
        tag_ = pts;
    }

    bool FaceTag::operator==(const FaceTag& other) const
    {
        if (tag_.size() != other.tag_.size()) {
            return false;
        }

        for (int i=0; i<tag_.size(); ++i)
        {
            if (tag_[i] != other.tag_[i]) {
                return false;
            }
        }

        return true;
    }

    /*bool FaceTag::operator<(const FaceTag& other) const
    {
        if (tag_[0] != other.tag_[0])
        {
            return tag_[0] < other.tag_[0];
        }
        else
        {
            if (tag_[1] != other.tag_[1])
            {
                return tag_[1] < other.tag_[1];
            }
            else
            {
                return tag_[2] < other.tag_[2];
            }
        }
    }*/

    bool FaceTag::operator!=(const FaceTag& other) const
    {
        return !(*this == other);
    }

    const std::vector<int>& FaceTag::operator()() const
    {
        return tag_;
    }

    std::ostream& operator<<(std::ostream& os, const FaceTag& ft)
    {
        for (int i: ft())
        {
            os << i << " ";
        }

        return os;
    }
}
