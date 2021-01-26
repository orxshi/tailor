#ifndef EXCHANGER_H
#define EXCHANGER_H

#include <numeric>
#include <boost/mpi.hpp>
#include "cantor.h"
#include "profiler.h"

namespace Tailor
{
    /*template<typename T>
        class Sender
        {
            public:

                struct Group
                {
                    int dest_rank_;
                    int dest_tag_;
                    std::deque<T> data_;

                    void add(T data);

                    Group(int dest_rank, int dest_tag, T data);
                };

                Sender(int tag);
                void send(const boost::mpi::communicator& comm, std::function<void(bool, int, std::string, Profiler*)> recv_cross, bool verbose=false, int rank=-1, std::string s="", Profiler* profiler=nullptr);
                void add(int dest_rank, int dest_tag, T data);
                int tag() const;
                const std::deque<Group>& group() const;
                std::deque<Group>& group();
                void clear_group();
                void set_tag(int);

            private:

                int tag_;
                std::deque<Group> group_;
        };

    template<typename T> class Exchanger;
        
    template<typename T>
        class Receiver
        {
            public:

                Receiver(int tag);
                int tag() const;
                const ArrCon<T>& arrival() const;
                ArrCon<T>& arrival_p();
                //void copy(const boost::mpi::communicator& comm, std::vector<T>& temp);
                void empty_arrival();
                void set_tag(int);

            private:

                int tag_;
                ArrCon<T> arrival_;
                int pos_;

                friend class Exchanger<T>;
        };*/

    //template<typename T> using SendCon = std::deque<Sender<T>>;
    //template<typename T> using RecvCon = std::deque<Receiver<T>>;

    template<typename T> using ArrCon = std::vector<T>;
    
    template<typename T>
        class Exchanger
        {
            protected:

            struct Group
            {
                int dest_rank_;
                int dest_tag_;
                std::deque<T> data_;

                void add(const T& data);

                Group(int dest_rank, int dest_tag, const T& data);
            };

            public:

                Exchanger(const boost::mpi::communicator*);
                void exchange(bool verbose=false, std::string s="", Profiler* profiler=nullptr, int rank=-1);
                const ArrCon<T>& arrival() const;
                ArrCon<T>& arrival();
                const std::deque<Group>& group() const;

            protected:

                void send_recv(bool verbose=false, int rank=-1, std::string s="", Profiler* profiler=nullptr);
                void add(int dest_rank, int dest_tag, const T& data);
                std::deque<Group>& group();
                ArrCon<T> arrival_;

            private:

                const boost::mpi::communicator* comm_;
                std::vector<int> global_nrecv_;
                int nrecv_;
                std::deque<Group> group_;
                int total_nrecv_;

                void prepare(std::string s="", Profiler* profiler=nullptr);
                void inform_receivers();
                void send(const boost::mpi::communicator& comm, bool verbose=false, int rank=-1, std::string s="", Profiler* profiler=nullptr);
                void recv(bool verbose=false, int rank=-1, std::string s="", Profiler* profiler=nullptr);
                virtual void prepare_storage() = 0;
        };
}

#include "exchanger.hpp"

#endif
