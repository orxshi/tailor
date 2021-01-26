#include "sp_exchanger.h"

namespace Tailor
{
    SpExchanger::SpExchanger(const Loadmap* lm, boost::mpi::communicator* comm): comm_(comm), Exchanger(comm), lm_(lm)
    {
    }

    void SpExchanger::prepare_storage()
    {
        const auto& aaa = lm_->aug_bin_tag();

        //sender_.reserve(lm_->rm().bin().size());
        //receiver_.reserve(lm_->rm().bin().size());

        //int ncell = 0;
        //int nwall = 0;
        //int nfar = 0;
        //int ninterog = 0;

        //auto max_bin_tag = lm_->rm().max_bin_tag();
        //sender_.resize(max_bin_tag+1, Sender<SpatialPartition>(-1));
        //receiver_.resize(max_bin_tag+1, Receiver<SpatialPartition>(-1));

        for (int p=0; p<comm_->size(); ++p)
        {
            for (const BinRMTag& t: aaa[p])
            {
                const Bin& bin = lm_->rm().bin(t.bintag());
                
                int btag = bin.tag()();

                //if (btag >= sender_.size())
                //{
                    //std::cout << "btag: " << btag << std::endl;
                    //std::cout << "sender size: " << sender_.size() << std::endl;
                //}
                //assert(btag < sender_.size());
                //assert(btag < receiver_.size());

                //sender_[btag].set_tag(btag);
                //receiver_[btag].set_tag(btag);

                //sender_.push_back(Sender<SpatialPartition>(bin.tag()()));
                //receiver_.push_back(Receiver<SpatialPartition>(bin.tag()()));

                SpatialPartition sp;
                sp.set_aabb(bin.aabb());
                sp.set_tag(bin.tag()); 
                sp.set_nei(bin.nei());

                std::deque<Mesh> mb;
                bin.group_mesh_all(*(lm_->mesh()), mb, comm_->rank());
                //bin.group_mesh_res(*(lm_->mesh()), mb, comm_->rank());

                if (mb.empty()) {
                    continue;
                }

                for (Mesh& m: mb)
                {
                    sp.add_mesh(m);
                    //ncell += m.cell().size();
                    //nwall += m.wall_boundaries().size();
                    //nfar += m.farfield_boundaries().size();
                    //ninterog += m.interog_boundaries().size();
                }


                //sender_.back().add(p, bin.tag()(), std::move(sp));
                //sender_[btag].add(p, btag, std::move(sp));
                add(p, btag, std::move(sp));
            }
        }

        //std::ofstream out;
        //std::string s = "ncll";
        //s.append(".dat");
        //out.open(s, std::ios_base::app);
        //out << comm_->rank() << " " << ncell << " " << nwall << " " << nfar << " " << ninterog << std::endl;
        //out.close();
    }
}
