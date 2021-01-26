#include "sp.h"

namespace Tailor
{
    void SpatialPartition::connect_cells(std::function<bool(const vec3<double>&)> is_resi, Profiler* profiler, std::string name)
    {
        for (Mesh& m: mesh_)
        {
            m.connect_cells(is_resi, comm_->rank(), profiler, name);
        }
    }
}

