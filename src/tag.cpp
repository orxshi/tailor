#include "tag.h"

/*PYBIND11_MODULE(Tag, m)
{
    py::class_<Tailor::Tag>(m, "Tag")
        .def(py::init<int>())
        .def("set", &Tailor::Tag::set)
        .def(py::self);
}*/

namespace Tailor
{
    void BinRMTag::set_bin_tag(const Tag& t)
    {
       bintag_ = t; 
    }

    void BinRMTag::set_rm_tag(const Tag& t)
    {
       rmtag_ = t; 
    }

    bool Tag::sorter(const Tag& lhs, const Tag& rhs)
    {
        return lhs() < rhs();
    }

    bool BinRMTag::isvalid() const
    {
        if (!bintag_.isvalid() || !rmtag_.isvalid())
        {
            return false;
        }

        return true;
    }

    bool Tag::isvalid() const
    {
        if (tag_ == -1)
        {
            return false;
        }

        return true;
    }

    void Tag::set(int t)
    {
        assert(t != -1);
        tag_ = t;
    }

    /*Tag& Tag::operator=(const Tag& t)
      {
      assert(t.isvalid());
      tag_ = t();
      return *this;
      }*/

    bool Tag::operator==(const Tag& other) const
    {
        return (tag_ == other.tag_);
    }

    bool Tag::operator<(const Tag& other) const
    {
        return tag_ < other.tag_;
    }

    bool Tag::operator!=(const Tag& other) const
    {
        return !(*this == other);
    }

    int Tag::operator()() const
    {
        return tag_;
    }
}
