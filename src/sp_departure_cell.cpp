#include "sp_departure_cell.hpp"

namespace Tailor
{
    const Tag& DepartureCell::cell_tag() const
    {
        return icell_;
    }
    const Tag& DepartureCell::mesh_tag() const
    {
        return imesh_;
    }
    void DepartureCell::set_cell_tag(const Tag& t)
    {
        icell_ = t;
    }
    void DepartureCell::set_mesh_tag(const Tag& t)
    {
        imesh_ = t;
    }
}
