//#include "meshpoint.h"
//#include "meshcell.h"
#include "read_mesh.h"
//#include <fstream>
//#include <limits>
//#include <iostream>

namespace Tailor
{
    std::ifstream& go_to_beg_of_line(std::ifstream& file, int num)
    {
        file.seekg(std::ios::beg);
        for(int i=0; i<num-1; ++i)
        {
            file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
        }

        return file;
    }

    void read_mesh_GMSH(Mesh& mesh, std::string file_name)
    {
        std::string temps;
        int tempi;
        int line_number = 1;
        int tag, igs, n_tags, geo, n_part, n_total, n_point, part_tag;
        gmsh gs;
        BouType phys = BouType::undefined;

        std::ifstream in;
        in.open (file_name);
        if (!in.is_open())
        {
            std::cout << "file " << file_name << " could not be opened." << std::endl;
            return;
        }
        assert(in.is_open());

        in >> temps; // mesh format.
        in >> temps; // version.
        assert(temps == "2.2"); // msh2 format.
        in >> temps; // version.
        in >> temps; // version.
        in >> temps; // end mesh format.
        in >> temps; // nodes.
        in >> n_point; // number under "$Nodes" is the number of points.
        line_number += 5;    

        //point_.reserve(n_point);
        for (int i=0; i<n_point; ++i)
        {
            int ptag;
            double x, y, z;
            in >> ptag;
            in >> x;
            in >> y;
            in >> z;
            ++line_number;
            //assert(ptag != 1449);

            //MeshPoint mp;
            //mp.set_tag(ptag);
            //mp.set_parent_mesh(tag_);
            //point_.push_back(std::move(mp));
            mesh.add_point(MeshPoint(x, y, z), Tag(ptag), n_point);
        }

        //mesh.set_point_tag_index_map();
        //mesh.shrink_points(); // uncomment
        mesh.sort_points(); // uncomment

        in >> temps; // end of nodes.
        in >> temps; // elements.
        in >> n_total; // the number under "$Elements" is total number of elements which includes boundary faces and cells.
        std::cout << "Total number of Elements: " << n_total << std::endl;
        line_number += 3;

        // read elements.
        int n_bface = 0;
        for (int e=0; e<n_total; ++e)
        {
            in >> tag; // tag of boundary face or cell.
            in >> igs; // geometric shape of boundary face or cell.        
            gs = static_cast<gmsh>(igs);
            in >> n_tags; // number of GMSH tags.

            int tag_count = 0;

            // read GMSH tags.
            while (n_tags > 0)
            {
                // read physical number.
                in >> tempi;
                ++ tag_count;
                if (tag_count == n_tags) break;
                in >> geo; // geometrical number.
                ++ tag_count;
                if (tag_count == n_tags) break;
                in >> n_part; // number of partitions to which element belongs.
                for (int i=0; i<n_part; ++i)
                {
                    in >> temps;
                }            

                break;
            }

            if (gs == gmsh::tri)
            {
                in >> tempi;
                in >> tempi;
                in >> tempi;
            }
            else if (gs == gmsh::quad)
            {
                in >> tempi;
                in >> tempi;
                in >> tempi;
                in >> tempi;
            }
            else
            {
                break;
            }        

            ++n_bface;
        }

        std::cout << "Total number of boundary elements: " << n_bface << std::endl;

        go_to_beg_of_line(in, line_number);

        /*if (n_bface == 0)
        {
            cell_.reserve(n_bface);
        }
        else
        {
            assert(n_bface == n_total);
            cell_.reserve(n_total);
        }*/

        for (int e=0; e<n_total; ++e)
        {
            in >> tag; // read tag of boundary face or cell.
            in >> igs; // read geometric shape of boundary face or cell.        
            gs = static_cast<gmsh>(igs);
            in >> n_tags; // read number of GMSH tags.

            int tag_count = 0;

            while (n_tags > 0)
            {
                // read physical number.
                in >> tempi;
                phys = static_cast<BouType>(tempi);
                ++ tag_count;
                if (tag_count == n_tags) break;
                in >> geo; // read geometrical number.
                ++ tag_count;
                if (tag_count == n_tags) break;
                in >> n_part; // read number of partitions to which element belongs.
                //assert(n_part == 1);
                for (int i=0; i<n_part; ++i)
                {
                    in >> part_tag;
                }            

                break;
            }

            if (e < n_bface)
            {
                Shape shape = Shape::undef;
                int n_vertex;
                switch (gs)
                {
                    case gmsh::tri:					    
                        n_vertex = 3;
                        shape = Shape::tri;
                        break;
                    case gmsh::quad:			
                        n_vertex = 4;
                        shape = Shape::quad;
                        break;
                    case gmsh::tet:
                        n_vertex = 4;						
                        shape = Shape::tet;
                        break;
                    case gmsh::hex:
                        n_vertex = 8;						
                        shape = Shape::hex;
                        break;
                    case gmsh::pri:
                        n_vertex = 6;						
                        shape = Shape::pri;
                        break;
                    default:
                        std::cout << "file: " << file_name << std::endl;
                        std::cout << "igs: " << igs << std::endl;
                        std::cout << "invalid n_vertex" << std::endl;
                        std::cout << "igs = " << igs << std::endl;
                        std::cout << "tag = " << tag << std::endl;
                        std::cout << "n_tags = " << n_tags << std::endl;
                        exit(0);
                }
                std::vector<int> vtx;
                vtx.reserve(n_vertex);
                for (int i=0; i<n_vertex; ++i)
                {
                    in >> tempi;
                    vtx.push_back(tempi);					
                    if (i != 0)
                    {
                        assert(vtx[i] != vtx[i-1]);
                    }
                }

                std::vector<MeshPoint> mpts;

                for (int z=0; z<vtx.size(); ++z)
                {
                    //mpts.push_back(mesh.point(Tag(vtx[z])));
                    auto pp = mesh.point(Tag(vtx[z]));
                    assert(pp != nullptr);
                    mpts.push_back(*pp);
                }

                auto mc = MeshCell(Tag(tag), mesh.tag(), mpts, phys, shape);
                //mesh.add_cell_only(mc);
                mesh.add_element(mc);

                //for (int i=0; i<mesh.cell().back().point().size(); ++i)
                //{
                    //if (i != 0)
                    //{
                        //assert(mesh.cell().back().point(i).tag() != mesh.cell().back().point(i-1).tag());
                    //}
                //}
            }
            else
            {
                Shape shape;
                int n_vertex;
                switch (gs)
                {
                    case gmsh::tri:					    
                        n_vertex = 3;
                        shape = Shape::tri;
                        break;
                    case gmsh::quad:			
                        n_vertex = 4;
                        shape = Shape::quad;
                        break;
                    case gmsh::tet:
                        n_vertex = 4;						
                        shape = Shape::tet;
                        break;
                    case gmsh::hex:
                        n_vertex = 8;						
                        shape = Shape::hex;
                        break;
                    case gmsh::pri:
                        n_vertex = 6;						
                        shape = Shape::pri;
                        break;
                    default:
                        std::cout << "file: " << file_name << std::endl;
                        std::cout << "igs: " << igs << std::endl;
                        std::cout << "invalid n_vertex" << std::endl;
                        std::cout << "igs = " << igs << std::endl;
                        std::cout << "tag = " << tag << std::endl;
                        std::cout << "n_tags = " << n_tags << std::endl;
                        exit(0);
                }

                // vertices.
                std::vector<int> vtx;
                vtx.reserve(n_vertex);

                for (int i=0; i<n_vertex; ++i)
                {
                    in >> tempi;
                    vtx.push_back(tempi);					
                    if (i != 0)
                    {
                        if (vtx[i] == vtx[i-1])
                        {
                            std::cout << "file: " << file_name << std::endl;
                            std::cout << "e: " << e << std::endl;
                            std::cout << "vtx[i]: " << vtx[i] << std::endl;
                            std::cout << "vtx[i-1]: " << vtx[i-1] << std::endl;
                        }
                        assert(vtx[i] != vtx[i-1]);
                    }
                }

                std::vector<MeshPoint> mpts;

                for (int z=0; z<vtx.size(); ++z)
                {
                    //mpts.push_back(mesh.point(Tag(vtx[z])));
                    auto pp = mesh.point(Tag(vtx[z]));
                    assert(pp != nullptr);
                    mpts.push_back(*pp);
                }

                //mesh.add_cell_only(MeshCell(Tag(tag), mesh.tag(), mpts, BouType::interior, shape), n_total);
                mesh.add_cell_only(MeshCell(Tag(tag), mesh.tag(), mpts, phys, shape), n_total);
                if (phys == BouType::wall)
                {
                    if (MeshCell(Tag(tag), mesh.tag(), mpts, phys, shape).face().size() != 1)
                    {
                        std::cout << "gs: " << static_cast<int>(gs) << std::endl;
                        std::cout << "file: " << file_name << std::endl;
                    }
                    assert(MeshCell(Tag(tag), mesh.tag(), mpts, phys, shape).face().size() == 1);
                }

                for (int i=0; i<mesh.cell().back().point().size(); ++i)
                {
                    if (i != 0)
                    {
                        assert(mesh.cell().back().point(i).tag() != mesh.cell().back().point(i-1).tag());
                    }
                }

            }        
        }

        in.close();

        mesh.sort_cells();
        mesh.update_points_from_cell_vertices(-1);
    }

    void read_wall_bou(Mesh& mesh, std::string file_name)
    {
        assert(file_name != "");

        auto aaa = file_name;
        file_name.append("_wall.msh");
        read_mesh_GMSH(mesh, file_name);
        aaa.append("_wall.vtk");
        mesh.print_as_vtk(aaa);
    }

    void read_symmetry_bou(Mesh& mesh, std::string file_name)
    {
        assert(file_name != "");

        file_name.append("_symmetry.msh");
        read_mesh_GMSH(mesh, file_name);
    }

    void read_dirichlet_bou(Mesh& mesh, std::string file_name)
    {
        assert(file_name != "");

        file_name.append("_dirichlet.msh");
        read_mesh_GMSH(mesh, file_name);
    }

    void read_farfield_bou(Mesh& mesh, std::string file_name)
    {
        assert(file_name != "");

        file_name.append("_farfield.msh");
        read_mesh_GMSH(mesh, file_name);
    }

    void read_empty_bou(Mesh& mesh, std::string file_name)
    {
        assert(file_name != "");

        file_name.append("_empty.msh");
        read_mesh_GMSH(mesh, file_name);
    }

    void read_interior_cells(Mesh& mesh, std::string file_name, int rank, bool uniproc)
    {
        assert(file_name != "");
        //if (rank == 0) return;

        /*if (rank > 99)
        {
            file_name.append("_000");
        }
        else if (rank > 9)
        {
            file_name.append("_0000");
        }
        else if (rank != 0)
        {
            file_name.append("_00000");
        }*/
        if (!uniproc)
        {
            file_name.append("_");
            //if (rank != 0)
            {
                file_name.append(std::to_string(rank));
            }
        }
        file_name.append(".msh");
        read_mesh_GMSH(mesh, file_name);
    }

    void read_cells_single(Mesh& mesh, std::string file_name)
    {
        assert(file_name != "");

        read_mesh_GMSH(mesh, file_name);
    }
}
