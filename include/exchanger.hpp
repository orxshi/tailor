#include "exchanger.h"

namespace Tailor
{
    template<typename T>
        Exchanger<T>::Group::Group(int dest_rank, int dest_tag, const T& data): dest_rank_(dest_rank), dest_tag_(dest_tag)
    {
        add(std::move(data));
    }

    template<typename T>
    void Exchanger<T>::Group::add(const T& data)
    {
        data_.push_back(data);
    }
    
    template<typename T>
    void Exchanger<T>::add(int dest_rank, int dest_tag, const T& data)
    {
        auto gr = std::find_if(group_.begin(), group_.end(), [&](const Group& gg){return gg.dest_rank_ == dest_rank;});
        if (gr == group_.end()) {
            group_.push_back(Group(dest_rank, dest_tag, std::move(data)));
        }
        else {
            gr->add(std::move(data));
        }
    }

    template<typename T>
    const std::deque<typename Exchanger<T>::Group>& Exchanger<T>::group() const
    {
        return group_;
    }

    template<typename T>
    std::deque<typename Exchanger<T>::Group>& Exchanger<T>::group()
    {
        return group_;
    }

    template<typename T>
    Exchanger<T>::Exchanger(const boost::mpi::communicator* comm): comm_(comm), nrecv_(0)
    {
    }

    template<typename T>
    void Exchanger<T>::exchange(bool verbose, std::string s, Profiler* profiler, int rank)
    {
        group_.clear();
        prepare(s, profiler);
        send_recv(verbose, rank, s, profiler);
    }

    template<typename T>
    void Exchanger<T>::prepare(std::string s, Profiler* profiler)
    {
        prepare_storage();
        inform_receivers();
    }

    template<typename T>
    void Exchanger<T>::send_recv(bool verbose, int rank, std::string s, Profiler* profiler)
    {
        nrecv_ = 0;

        arrival_.clear();
        
        send(*comm_, verbose, rank, "", profiler);

        while (nrecv_ != total_nrecv_) {
            recv(verbose, rank, "", profiler);
        }

        comm_->barrier();
    }

    template<typename T>
    void Exchanger<T>::inform_receivers()
    {
        global_nrecv_.clear();
        total_nrecv_ = 0;

        std::vector<int> local_nrecv(comm_->size(), 0);

        for (const auto& gr: group_)
        {
            local_nrecv[gr.dest_rank_] += gr.data_.size();
        }

        comm_->barrier();

        for (int i=0; i<comm_->size(); ++i)
        {
            boost::mpi::gather(*comm_, local_nrecv[i], global_nrecv_, i);
        }

        total_nrecv_ = std::accumulate(global_nrecv_.begin(), global_nrecv_.end(), 0);
    }

    template<typename T>
    void Exchanger<T>::recv(bool verbose, int rank, std::string s, Profiler* profiler)
    {
        if (nrecv_ == total_nrecv_) {
            return;
        }

        auto status = comm_->iprobe(boost::mpi::any_source, boost::mpi::any_tag);

        if (!status) {
            return;
        }

        int nincoming = global_nrecv_[status->source()];

        std::vector<T> temp;
        comm_->recv(status->source(), status->tag(), temp);
        nrecv_ += nincoming;
        std::copy(temp.begin(), temp.end(), std::back_inserter(arrival_));
    }

    template<typename T>
    void Exchanger<T>::send(const boost::mpi::communicator& comm, bool verbose, int rank, std::string s, Profiler* profiler)
    {
        for (Group& gr: group_)
        {
            assert(gr.dest_rank_ >= 0);
            assert(!gr.data_.empty());

            boost::mpi::request req = comm.isend(gr.dest_rank_, gr.dest_tag_, gr.data_);

            while (true)
            {
                auto status = req.test();

                if (status)
                {
                    break;
                }

                recv(verbose, rank, s, profiler);
            }
        }
    }

    template<typename T>
        ArrCon<T>& Exchanger<T>::arrival()
    {
        return arrival_;
    }

    template<typename T>
        const ArrCon<T>& Exchanger<T>::arrival() const
    {
        return arrival_;
    }
}
