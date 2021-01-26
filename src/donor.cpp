#include "donor.h"

namespace Tailor
{
    Donor::Donor(): addr_(nullptr)
    {
    }

    Donor::Donor(const Tag& mesh_tag, const Tag& cell_tag, const MeshCell* addr): mesh_tag_(mesh_tag), cell_tag_(cell_tag), addr_(addr)
    {
    }
}
