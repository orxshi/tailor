#ifndef READ_MESH_H
#define	READ_MESH_H

#include "mesh.h"

namespace Tailor
{
    enum class gmsh
    {
        tri = 2,
        quad = 3,
        tet = 4,
        hex = 5,
        pri = 6,
    };

    void read_mesh_GMSH(Mesh& mesh, std::string file_name);
    void read_cells_single(Mesh& mesh, std::string file_name);
    void read_interior_cells(Mesh& mesh, std::string file_name, int rank, bool uniproc);
    void read_symmetry_bou(Mesh& mesh, std::string file_name);
    void read_dirichlet_bou(Mesh& mesh, std::string file_name);
    void read_wall_bou(Mesh& mesh, std::string file_name);
    void read_farfield_bou(Mesh& mesh, std::string file_name);
    void read_empty_bou(Mesh& mesh, std::string file_name);
}

#endif
