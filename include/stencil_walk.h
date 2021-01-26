#ifndef STENCIL_WALK_H
#define	STENCIL_WALK_H

#include "mesh.h"
#include "adt.h"

namespace Tailor
{
    BouType find_intersected_face(const Mesh& mesh, const MeshCell*& current_cell, Tag& prev_nei, Tag& prev_prev_nei, const MeshCell*& closest_cell, const MeshFace*& inter_mf, const MeshCell*& nei, const Segment& start_to_target, bool verbose);
    StencilWalkResult stencil_walk(const Mesh& donor_mesh, const Point& target, const MeshCell& starting_cell, const MeshCell*& closest_cell, int dummyrank, int& iter, bool verbose);

    std::ostream& operator<<(std::ostream& os, const StencilWalkResult& res);
}

#endif
