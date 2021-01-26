/*#include "sp_neighbor.hpp"

int Neighbor::commtag_send() const
{
    return commtag_send_;
}

int Neighbor::commtag_recv() const
{
    return commtag_recv_;
}

void Neighbor::set_commtag_send(int i)
{
    commtag_send_ = i;
}

void Neighbor::set_commtag_recv(int i)
{
    commtag_recv_ = i;
}

int Neighbor::gbi() const
{
    return gbi_;
}

void Neighbor::send_regular(const boost::mpi::communicator& comm, boost::mpi::request& req)
{
    if (rp_send_ != RelativePosition::undefined)
        req = comm.isend(rank_, commtag_send_, cell_);
}

bool Neighbor::send_transit(const boost::mpi::communicator& comm, boost::mpi::request& req)
{
    if (transit_.cell().empty()) return;
    if (transit_.next_dest() == transit_.final_dest()) return;

    bool exact;
    RelativePosition rp = get_relative_position(gbi_, transit_.final_binindex(), nstripe_, binsize_, exact);

    if (exact)
    {
        req = comm.isend(transit_.final_dest(), send_tag_, transit_.cell());
        transit_.clean_cells();
    }
    else
    {
        for (int k=0; k<RM_NNEI; ++k)
        {
            if (rp_recv == nei_[k].rp_recv())
            {
                nei_bin_[k].set_next_transit_dest(nei_bin_[k].rank());
                req = comm.isend(transit_.next_dest(), commtag_send_, transit_.cell());
                break;
            }
        }
    }
}

void Neighbor::recv_regular(const boost::mpi::communicator& comm, boost::mpi::request& req, std::vector<std::vector<MeshCell>>& ac)
{
    assert(static_cast<int>(rp_recv_) >= 0);
    assert(static_cast<int>(rp_recv_) < ac.size());
    if (rp_recv_ != RelativePosition::undefined)
        req = comm.irecv(boost::mpi::any_source, commtag_recv_, ac[static_cast<int>(rp_recv_)]);
}

void Neighbor::recv_transit(const boost::mpi::communicator& comm, boost::mpi::request& req, std::vector<std::vector<MeshCell>>& ac)
{
    req = comm.irecv(boost::mpi::any_source, commtag_recv_, ac[static_cast<int>(rp_recv_)]);
}

void Neighbor::set_gbi(int i)
{
    gbi_ = i;
}

int Neighbor::rank() const
{
    return rank_;
}

const std::vector<MeshCell>& Neighbor::cell() const
{
    return cell_;
}

void Neighbor::add_cell(const MeshCell& c)
{
    cell_.push_back(c);
}

void Neighbor::add_transit_cell(const MeshCell& c)
{
    transit_.add_cell(c);
}

void Transit::add_cell(const MeshCell& c)
{
    cell_.push_back(c);
}

void Neighbor::set_next_transit_dest(unsigned int i)
{
    transit_.set_next_dest(i);
}

void Neighbor::set_final_transit_dest(unsigned int i)
{
    transit_.set_final_dest(i);
}

void Transit::set_next_dest(unsigned int i)
{
    next_dest_ = i;
}

void Transit::set_final_dest(unsigned int i)
{
    final_dest_ = i;
}

size_t Neighbor::size() const
{
    return cell_.size();
}

RelativePosition Neighbor::rp_recv() const
{
    return rp_recv_;
}

RelativePosition Neighbor::rp_send() const
{
    return rp_send_;
}

void Neighbor::set_rp_recv(RelativePosition i)
{
    rp_recv_ = i;
}

void Neighbor::set_rp_send(RelativePosition i)
{
    rp_send_ = i;
}

void Neighbor::set_rank(int rank)
{
    if (this->rank_ == -1)
    {
        this->rank_= rank;
    }
}
int Neighbor::send_to_recv_rp_converter(int rp_send)
{
    if (rp_send > 8 || rp_send < 0)
    {
        std::cout << "invalid input in send_to_recv_rp_converter()" << std::endl;
        exit(-2);
    }

    int recv[9];
    recv[0] = 7;
    recv[1] = 6;
    recv[2] = 5;
    recv[3] = 4;
    recv[4] = 8;
    recv[5] = 3;
    recv[6] = 2;
    recv[7] = 1;
    recv[8] = 0;

    int send[9];
    send[0] = 0;
    send[1] = 1;
    send[2] = 2;
    send[3] = 3;
    send[4] = 8;
    send[5] = 4;
    send[6] = 5;
    send[7] = 6;
    send[8] = 7;

    for (int i=0; i<9; ++i)
    {
        if (rp_send == send[i])
        {
            return recv[i];
        }
    }

    return -1;
}

int Neighbor::recv_to_send_rp_converter(int rp_recv)
{
    if (rp_recv > 8 || rp_recv < 0)
    {
        std::cout << "invalid input in send_to_recv_rp_converter()" << std::endl;
        exit(-2);
    }

    int recv[9];
    recv[0] = 7;
    recv[1] = 6;
    recv[2] = 5;
    recv[3] = 4;
    recv[4] = 8;
    recv[5] = 3;
    recv[6] = 2;
    recv[7] = 1;
    recv[8] = 0;

    int send[9];
    send[0] = 0;
    send[1] = 1;
    send[2] = 2;
    send[3] = 3;
    send[4] = 8;
    send[5] = 4;
    send[6] = 5;
    send[7] = 6;
    send[8] = 7;

    for (int i=0; i<9; ++i)
    {
        if (rp_recv == recv[i])
        {
            return send[i];
        }
    }

    return -1;
}*/
