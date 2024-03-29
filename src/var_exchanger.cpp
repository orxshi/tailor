#include "var_exchanger.h"

namespace Tailor
{
    Var::Var(int source_rank, int source_tag, const Vector5& var, int mesh, int cell, const Vector5& dQ): source_rank_(source_rank), source_tag_(source_tag), var_(var), dQ_(dQ)
    {
        mesh_cell_ = std::make_pair(mesh, cell);
    }

    //VarExchanger::VarExchanger(const std::vector<SpatialPartition>* sp, const RecvCon<Request>* other_receiver, const boost::mpi::communicator* comm): comm_(comm), Exchanger(comm), sp_(sp), other_receiver_(other_receiver)
    VarExchanger::VarExchanger(const std::vector<SpatialPartition>* sp, const Exchanger<Request>* other_exc, const boost::mpi::communicator* comm): comm_(comm), Exchanger(comm), sp_(sp), other_exc_(other_exc)
    {
    }

    void VarExchanger::prepare_storage()
    {
        assert(sp_->size() == 1);
        //assert(other_receiver_->size() == sp_->size());

        const SpatialPartition& sp = sp_->front();
        //sender_.push_back(Sender<Var>(sp.tag()()));
        //receiver_.push_back(Receiver<Var>(sp.tag()()));
        //const Receiver<Request>& os = other_receiver_->front();

        //assert(os.tag() == sp.tag()());

        for (const Request& r: other_exc_->arrival())
        {
            auto m = std::find_if(sp.mesh().begin(), sp.mesh().end(), [&](const Mesh& mm){return mm.tag()() == r.mesh_cell_.first;});
            if (m == sp.mesh().end())
            {
                std::cout << comm_->rank() << " sp.mesh().size(): " << sp.mesh().size() << std::endl;
                std::cout << comm_->rank() << " req mesh tag: " << r.mesh_cell_.first << std::endl;
                std::cout << comm_->rank() << " req cell tag: " << r.mesh_cell_.second << std::endl;
            }
            assert(m != sp.mesh().end());
            auto mc = std::find_if(m->cell().begin(), m->cell().end(), [&](const MeshCell& cc){return cc.tag()() == r.mesh_cell_.second;});
            assert(mc != m->cell().end());

            /*if (mc->oga_cell_type() != OGA_cell_type_t::field)
            {
                std::cout << comm_->rank() << " mesh tag: " << m->tag()() << std::endl;
                std::cout << comm_->rank() << " mc tag: " << mc->tag()() << std::endl;
                std::cout << comm_->rank() << " type: " << static_cast<int>(mc->oga_cell_type()) << std::endl;
            }*/
            //assert(mc->oga_cell_type() == OGA_cell_type_t::field);
            assert(!mc->prim().isnan());

            assert(mc->prim(0) > 0.);

            //sender_.front().add(r.source_rank_, r.source_tag_, Var(comm_->rank(), sp.tag()(), mc->prim(), m->tag()(), mc->tag()()));
            add(r.source_rank_, r.source_tag_, Var(comm_->rank(), sp.tag()(), mc->prim(), m->tag()(), mc->tag()(), mc->dQ()));

            //if (mc->parent_mesh()() == 0 && mc->tag()() == 874)
            //{
                //std::cout << "rankkkk: " << comm_->rank() << std::endl;
                //std::cout << "source_rank: " << r.source_rank_ << std::endl;
                //std::cout << "source_tag: " << r.source_tag_ << std::endl;
                //std::cout << "source_mesh: " << r.mesh_cell_.first << std::endl;
                //std::cout << "source_cell: " << r.mesh_cell_.second << std::endl;
            //}
        }
    }

    //void VarExchanger::update(Mesh* mesh, Profiler* profiler, std::string fname)
    void VarExchanger::update(Profiler* profiler, std::string fname)
    {
        if (profiler != nullptr) {profiler->start(fname);}
        //if (mesh != nullptr)
        {
            assert(sp_->size() == 1);
            {
                auto& sp = sp_->front();

                for (auto& gr: group())
                {
                    for (Var& v: gr.data_)
                    {
                        //if (v.mesh_cell_.first != mesh->tag()())
                        //{
                            //continue;
                        //}

                        auto m = std::find_if(sp.mesh().begin(), sp.mesh().end(), [&](const Mesh& mm){return mm.tag()() == v.mesh_cell_.first;});
                        assert(m != sp.mesh().end());
                        //const auto& mc = mesh->cell(Tag(v.mesh_cell_.second));
                        const auto& mc = m->cell(Tag(v.mesh_cell_.second));
                        v.var_ = mc.prim();
                        v.dQ_ = mc.dQ();
                    }
                }
            }
        }
        if (profiler != nullptr) {profiler->stop(fname);}

        send_recv(false, -1, fname, profiler);
    }
}
