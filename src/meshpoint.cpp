#include <meshpoint.h>

namespace Tailor
{
    /*size_t MeshPoint::mem() const
    {
        size_t size = 0;

        size += sizeof(bool);
        size += sizeof(int) * 2;
        size += p_.mem();
        size += parent_cell_.mem();

        return size;
    }*/

    bool MeshPoint::erase() const
    {
        return erase_;
    }

    void MeshPoint::mark_to_be_erased()
    {
        erase_ = true;
    }

    bool MeshPoint::operator<(const MeshPoint& other) const
    {
        return tag_ < other.tag();
    }

    bool MeshPoint::operator==(const MeshPoint& other) const
    {
        return tag_ == other.tag();
    }

    void MeshPoint::merge(const MeshPoint& other)
    {
        for (const auto& t: other.parent_cell())
        {
            auto it = std::find(parent_cell_.begin(), parent_cell_.end(), t); 
            if (it == parent_cell_.end())
            {
                parent_cell_.push_back(t);
                //parent_cell_.add(t);
            }
        }
    }

    void MeshPoint::remove_parent_cell(const Tag& ic)
    {
        assert(ic.isvalid());
        parent_cell_.erase(std::remove_if(parent_cell_.begin(), parent_cell_.end(), [&](const auto& a){return a == ic;}), parent_cell_.end());
    }

    const Tag& MeshPoint::parent_mesh() const
    {
        return parent_mesh_;
    }

    void MeshPoint::rotate_point(double angle, int axis, const Vector3& rot_point)
    {
        RotationMatrix rm;

        auto old = p_.r();
        auto z = p_.r() - rot_point;
        auto newz = rm.rotate(angle, axis, z);
        p_.set_r(newz + rot_point - 0.1*sin(angle) + 0.2*cos(angle));

        if (std::isnan(p_.r(0)) || std::isnan(p_.r(1)) || std::isnan(p_.r(2)))
        {
            std::cout << "z: " << z(0) << " " << z(1) << " " << z(2) << std::endl;
            std::cout << "newz: " << newz(0) << " " << newz(1) << " " << newz(2) << std::endl;
            std::cout << "old: " << old(0) << " " << old(1) << " " << old(2) << std::endl;
            std::cout << "angle: " << angle << std::endl;
            std::cout << "axis: " << axis << std::endl;
            std::cout << "rotpoint: " << rot_point(0) << " " << rot_point(1) << " " << rot_point(2) << std::endl;
        }

        assert(!std::isnan(p_.r(0)));
        assert(!std::isnan(p_.r(1)));
        assert(!std::isnan(p_.r(2)));

    }

    void MeshPoint::move_point(const Vector3& v)
    {
        p_.set_r(p_.r() + v);
    }

    //void MeshPoint::set_px(double v)
    //{
        //p_.set_r(v, p_.r(1));
    //}

    void MeshPoint::remove_from_parent_cells()
    {
        parent_cell_.clear();
    }

    //MPI_Datatype mpi_dtype_meshpoint = MPI_DATATYPE_NULL;

    void MeshPoint::set_tag(const Tag& t)
    {
        tag_ = t;
    }
    void MeshPoint::set_parent_mesh(const Tag& ptag)
    {
        parent_mesh_ = ptag;
    }

    const std::vector<Tag>& MeshPoint::parent_cell() const
    {
        return parent_cell_;
    }

    const Tag& MeshPoint::parent_cell(int i) const
    {
        return parent_cell_[i];
    }

    const Tag& MeshPoint::tag() const
    {
        return tag_;
    }

    const Point& MeshPoint::p() const
    {
        return p_;
    }

    void MeshPoint::remove_parent_cells()
    {
        parent_cell_.clear();
    }

    void MeshPoint::add_parent_cell(const Tag& celltag)
    {
        assert(celltag.isvalid());
        auto it = std::find(parent_cell_.begin(), parent_cell_.end(), celltag);
        if (it == parent_cell_.end()) {
            parent_cell_.push_back(celltag);
        }
        //else
        //{
            //assert(celltag.const_addr() != nullptr);
            //it->set_addr(celltag.addr());
        //}

        //auto itt = std::unique(parent_cell_.begin(), parent_cell_.end());
        //assert(itt == parent_cell_.end());
    }

    /*void MeshPoint::add_parent_bface(int i)
      {
      parent_bface.push_back(i);
      }

      void MeshPoint::add_parent_iface(int i)
      {
      parent_iface.push_back(i);
      }*/

    /*void MeshPoint::make_layout()
      {
    // members.
    // 1: int parent_mesh;
    // 1: int partition;
    // 2: Point p;
    // 3: Array<int, 10> parent_cell;
    // 4: Array<int, 10> parent_bface;
    // 5: Array<int, 10> parent_iface;

    int flag;
    MPI_Initialized(&flag);
    if (mpi_dtype_meshpoint == MPI_DATATYPE_NULL && flag == 1)
    {
    int nblock = 5;
    // block count.
    int block_count[nblock];
    block_count[0] = 2;
    block_count[1] = 1;
    block_count[2] = 1;
    block_count[3] = 1;
    block_count[4] = 1;
    // address.
    MPI_Aint addr[nblock];
    MPI_Get_address(&parent_mesh, &addr[0]);
    MPI_Get_address(&p, &addr[1]);
    MPI_Get_address(&parent_cell, &addr[2]);
    MPI_Get_address(&parent_bface, &addr[3]);
    MPI_Get_address(&parent_iface, &addr[4]);
    // offset.
    MPI_Aint offset[nblock];
    offset[0] = 0;
    offset[1] = offset[0] + addr[1] - addr[0];
    offset[2] = offset[1] + addr[2] - addr[1];
    offset[3] = offset[2] + addr[3] - addr[2];
    offset[4] = offset[3] + addr[4] - addr[3];
    // block type.
    MPI_Datatype mdt0;
    MPI_Datatype mdt1;
    MPI_Datatype mdt2;
    parent_cell.make_layout(mdt0, MPI_INT);
    parent_bface.make_layout(mdt1, MPI_INT);
    parent_iface.make_layout(mdt2, MPI_INT);

    MPI_Datatype block_type[nblock];
    block_type[0] = MPI_INT;
    block_type[1] = mpi_dtype_point;
    block_type[2] = mdt0;
    block_type[3] = mdt1;
    block_type[4] = mdt2;
    // define and commit type.
    MPI_Type_create_struct(nblock, block_count, offset, block_type, &mpi_dtype_meshpoint);
    MPI_Type_commit(&mpi_dtype_meshpoint);
    // free MPI datatypes.
    MPI_Type_free(&mdt0);
    MPI_Type_free(&mdt1);
    MPI_Type_free(&mdt2);
    }
    }*/
}
