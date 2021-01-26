#include "profiler.h"

namespace Tailor
{
    void mem_usage(const boost::mpi::communicator* comm, std::string s)
    {
        assert(comm != nullptr);

        std::string fn = "mem-";
        fn.append(std::to_string(comm->rank()));
        fn.append(".dat");
        std::ofstream out;
        out.open(fn, std::ios_base::app);

        struct rusage usage;
        int ret = getrusage(RUSAGE_SELF, &usage);
        out << s << " " << usage.ru_maxrss << std::endl;

        out.close();
    }
    double Profiler::skewness(const std::pair<std::string, int>& func)
    {
        double n = world_->size() - 1;
        double sum = 0.;
        for (int i=1; i<world_->size(); ++i) // not including master.
        {
            sum += elapsed_main_iter_[func.second][i];
        }

        double mean = sum / n;

        double sum_var = 0.;
        double sum_skw = 0.;
        for (int i=1; i<world_->size(); ++i) // not including master.
        {
            sum_var += std::pow((elapsed_main_iter_[func.second][i] - mean), 2);
            sum_skw += std::pow((elapsed_main_iter_[func.second][i] - mean), 3);
        }

        double sdev = std::sqrt(sum_var / (n - 1));
        double skew = (n / ((n-1) * (n-2))) * (sum_skw / std::pow(sdev, 3));

        return skew;
    }

    //Profiler::Profiler(const MPI_Comm& comm, bool with_barrier, const Settings& settings): world_(comm, boost::mpi::comm_attach), with_barrier_(with_barrier)
    //Profiler::Profiler(const MPI_Comm& comm, bool with_barrier): world_(comm, boost::mpi::comm_attach), with_barrier_(with_barrier)
    Profiler::Profiler(const boost::mpi::communicator* comm, bool with_barrier): world_(comm), with_barrier_(with_barrier)
    {
        //for (int j=0; j<elapsed_.size(); ++j)
        //{
        //elapsed_[j].resize(world_.size());
            //elapsed_main_[j].resize(world_.size());
            //elapsed_main_iter_[j].resize(world_.size());
        //}

        //settings_ = std::make_shared<Settings>(settings);

        total_elapsed_.resize(world_->size());
        total_elapsed_iter_.resize(world_->size());

        excluded_ = {
            "add-mesh",
            "first-remap-exchange",
            "first-remap-add-arrival",
            "refine",
            "second-remap-dist",
            "second-remap-send",
            "second-remap-remove",
            "merge-di",
            "donor-info",
            "solver-read",
            "prepare-di",
            "det-undefined-donors",
            "notify-owners",
            "send-ws",
            "recv-ws",
            "recv-req-cells",
            "spc-make-rm",
            "spc-connect-cells"
        };

        included_ = {
            "donor-search"
        };
    }

    const std::vector<std::pair<std::string, int>>& Profiler::func() const
    {
        return func_;
    }

    const size_t Profiler::nprobe() const
    {
        return nprobe_;
    }

    void Profiler::print() const
    {
        //if (world_.rank() != MASTER) return;

        //if (settings_->print_to_file_)
        {
            std::ofstream out;
            out.open("a.csv");
            out << "proc,";
            for (auto it=func_.begin(); it!=func_.end(); ++it)
            {
                out << it->first;
                if (it != std::prev(func_.end()))
                {
                    out << ",";
                }
            }
            out << "\n";
            for (int i=0; i<world_->size(); ++i)
            {
                out << i << ",";
                for (auto it = func_.begin(); it != func_.end(); ++it)
                {
                    out << elapsed_main_[it->second][i];
                    if (it != std::prev(func_.end()))
                        out << ",";
                }
                out << "\n";
            }
            out.close();
        }

        //if (settings_->print_to_file_)
        {
            BOOST_LOG_TRIVIAL(info) << "total elapsed time: " << *std::max_element(total_elapsed_.begin(), total_elapsed_.end());
            auto res = minmax_element(total_elapsed_.begin()+1, total_elapsed_.end());
            double min = *res.first;
            double max = *res.second;
            double dev = (max - min) / max * 100.;
            BOOST_LOG_TRIVIAL(info) << "total elapsed dev: " << dev;
        }

        //if (settings_->print_to_file_)
        {
            std::ofstream out;
            std::string fn = "fn.csv";
            out.open(fn);
            out << "function,elapsed\n";
            for (auto it = func_.begin(); it != func_.end(); ++it)
            {
                out << it->first;
                out << ",";
                double temp = 0.;
                for (int i=0; i<world_->size(); ++i)
                {
                    temp += elapsed_main_[it->second][i];
                }
                out << temp;
                out << "\n";
            }
            out.close();
        }
    }

    std::string Tier::to_string(std::vector<int>::const_iterator begin, std::vector<int>::const_iterator end)
    {
        assert(std::distance(begin, end) != 0);
        std::string out;

        for (auto it = begin; it != end; ++it)
        {
            out.append(std::to_string(*it));
            if (it != std::prev(end)) {
                out.append(".");
            }
        }

        return out;
    }

    std::vector<int> Tier::to_int(std::string s)
    {
        assert(!s.empty());
        std::vector<int> t;

        for (char c: s)
        {
            if (c != '.')
            {
                t.push_back(c - '0');
            }
        }

        return t;
    }

    void Profiler::bstart(std::string func_name)
    {
        world_->barrier();
        start(func_name);
    }

    void Profiler::bstop(std::string func_name)
    {
        world_->barrier();
        stop(func_name);
    }

    void Profiler::start(std::string func_name)
    {
        if (with_barrier_) {
            world_->barrier();
        }

        //std::vector<int> tier_int = Tier::to_int(tier_string);
        //assert(!tier_int.empty());

        //std::cout << "adding tier" << std::endl;
        //add_tier(tier_string);
        //std::cout << "added tier" << std::endl;
                    
        ////tier_tag_address_.insert(std::make_pair(tier_string, ptr));
        //tier_func_.insert(std::make_pair(tier_string, func_name));
            
        int f;
        auto it = std::find_if(func_.begin(), func_.end(), [&](const auto& e){return e.first == func_name;});
        //auto it = std::find_if(func_.begin(), func_.end(), [&](const auto& e){return e.name == func_name;});
        if (it == func_.end())
        {
            nprobe_ = func_.size() + 1;
            elapsed_.push_back(std::vector<double>(world_->size()));
            elapsed_main_.push_back(std::vector<double>(world_->size()));
            elapsed_main_iter_.push_back(std::vector<double>(world_->size()));
            //start_.push_back(time_point());
            start_.push_back(0);
            //stop_.push_back(time_point());
            stop_.push_back(0);
            //func_.insert(std::pair<std::string, int>(func_name, func_.size()));
            func_.push_back(std::pair<std::string, int>(func_name, func_.size()));
            func_rev_.insert(std::pair<int, std::string>(func_.size()-1, func_name));
            //func_.insert(Function(func_name, func_.size(), tier_tag));
            f = func_.size() - 1;
        }
        else
        {
            f = it->second;
        }

        //start_[f] = timer.start_stop();
        //start_[f] = MPI_Wtime();
        start_[f] += MPI_Wtime() - stop_[f];
        assert(elapsed_main_iter_.size() == func_.size());
    }

    void Profiler::add_tier(std::string tier_string)
    {
        auto tier_it = find_if(tier_.begin(), tier_.end(), [&](const Tier& t){return t.tag().front() == tier_string.front();}); 
        std::cout << "tier string: " << tier_string << std::endl;
        if (tier_it == tier_.end())
        {
            std::cout << "tier not existing" << std::endl;
            tier_.push_back(Tier(tier_string));
        }
        else
        {
            std::cout << "tier existing" << std::endl;
            tier_string.erase(0,2);
            std::cout << "taggggg: " << tier_string << std::endl;
            tier_it->add(tier_string);
        }
    }

    Tier::Tier(std::string tier_string)
    {
        tag_ = tier_string;
        std::cout << "tag: " << tier_string << std::endl;

        std::vector<int> tier_int = Tier::to_int(tier_string);
        for (int i: tier_int)
            std::cout << "int: " << i << std::endl;
        assert(tier_int.size() != 0);

        if (tier_int.size() > 1)
        {
            add(tier_int.begin()+1, tier_int.end());
        }
    }

    Tier::Tier(std::vector<int>::const_iterator begin, std::vector<int>::const_iterator end): Tier(to_string(begin, end))
    {
    }

    void Tier::add(std::string tier_string)
    {
        std::vector<int> tier_int = Tier::to_int(tier_string);
        add(tier_int.begin(), tier_int.end());
    }

    void Tier::add(std::vector<int>::const_iterator begin, std::vector<int>::const_iterator end)
    {
        for (auto child: child_)
        {
            assert(child != nullptr);

            //std::cout << "child tag: " << child->tag() << std::endl;
            //std::cout << "to string: " << Tier::to_string(begin, end) << std::endl;
            if (child->tag().front() == Tier::to_string(begin, end).front())
            {
                child->add(begin, end);
                return;
            }
        }

        child_.push_back(std::make_shared<Tier>(begin, end));
    }

    void Profiler::clear_iter()
    {
        elapsed_.clear();
        elapsed_batch_.clear();
        elapsed_main_iter_.clear();
        func_.clear();
        start_.clear();
        stop_.clear();
        /*if (!func_.empty())
        {
            elapsed_.resize(func_.size());
            elapsed_main_iter_.resize(func_.size());
            assert(elapsed_main_iter_.size() == func_.size());
            for (int i=0; i<func_.size(); ++i)
            {
                elapsed_[i].resize(world_.size());
                elapsed_main_iter_[i].resize(world_.size());
            }
        }*/

        total_elapsed_iter_.clear();
        total_elapsed_iter_.resize(world_->size());
        /*if (elapsed_main_iter_.size() != func_.size())
        {
            std::cout << elapsed_main_iter_.size() << std::endl;
            std::cout << func_.size() << std::endl;
        }
        assert(elapsed_main_iter_.size() == func_.size());*/
    }

    void Profiler::print_iter(int iteration)
    {
        //if (settings_->print_to_file_)
        {
            update_batch();

            if (world_->rank() != 0) return;

            std::ofstream out;
            /*out.open("tdev-vs-iter.csv", std::fstream::app);
              if (iteration == 0)
              out << "iter,time\n";
              auto res = minmax_element(total_elapsed_iter_.begin()+1, total_elapsed_iter_.end());
              double min = *res.first;
              double max = *res.second;
              double dev = (max - min) / max * 100.;
              out << iteration << "," << dev << "\n";
              out.close();*/

            std::string fn = "hist-iter-";
            fn.append(std::to_string(iteration));
            fn.append(".csv");
            out.open(fn);
            out << "proc,";
            for (auto it=func_.begin(); it!=func_.end(); ++it)
            {
                out << it->first;
                if (it != std::prev(func_.end()))
                {
                    out << ",";
                }
            }
            out << "\n";
            for (int i=0; i<world_->size(); ++i)
            {
                out << i << ",";
                for (auto it = func_.begin(); it != func_.end(); ++it)
                {
                    out << elapsed_main_iter_[it->second][i];
                    if (it != std::prev(func_.end()))
                        out << ",";
                }
                //out << total_elapsed_iter_[i];
                out << "\n";
            }

            out.close();

            // iter slowest proc time
            //{
                //std::string fn = "slowest-iter.csv";
                //out.open(fn, std::fstream::app);

                //if (iteration == 0) {
                    //out << "iter,time\n";
                //}
                //out << iteration << "," << *std::max_element(std::next(total_elapsed_iter_.begin()), total_elapsed_iter_.end()) << "\n";
                //out.close();
            //}

            // skewness
            //{
            //    std::string fn = "skew-iter-";
            //    fn.append(std::to_string(iteration));
            //    fn.append(".csv");
            //    out.open(fn);

            //    out << "function,time\n";

            //    for (auto it = func_.begin(); it != func_.end(); ++it)
            //    {
            //        out << it->first;
            //        out << ",";

            //        double skew = skewness(*it);

            //        out << skew;
            //        out << "\n";
            //    }

            //    out.close();
            //}
            // impact
            //{
            //    fn = "impact-iter-";
            //    fn.append(std::to_string(iteration));
            //    fn.append(".csv");
            //    out.open(fn);
            //    out << "function,elapsed\n";
            //    std::vector<double> total_elap(func_.size(), 0.);
            //    for (auto it = func_.begin(); it != func_.end(); ++it)
            //    {
            //        for (int i=0; i<world_.size(); ++i)
            //        {
            //            total_elap[std::distance(func_.begin(), it)] += elapsed_main_iter_[it->second][i];
            //        }
            //    }
            //    for (auto it = func_.begin(); it != func_.end(); ++it)
            //    {
            //        out << it->first;
            //        out << ",";
            //        double skew = skewness(*it);
            //        out << total_elap[std::distance(func_.begin(), it)] * skew;
            //        out << "\n";
            //    }
            //    out.close();
            //}
            //{
            //    fn = "fn-iter-";
            //    fn.append(std::to_string(iteration));
            //    fn.append(".csv");
            //    out.open(fn);
            //    out << "function,elapsed\n";
            //    std::vector<double> total_elap(func_.size(), 0.);
            //    for (auto it = func_.begin(); it != func_.end(); ++it)
            //    {
            //        for (int i=0; i<world_.size(); ++i)
            //        {
            //            total_elap[std::distance(func_.begin(), it)] += elapsed_main_iter_[it->second][i];
            //        }
            //    }
            //    for (auto it = func_.begin(); it != func_.end(); ++it)
            //    {
            //        out << it->first;
            //        out << ",";
            //        out << total_elap[std::distance(func_.begin(), it)];
            //        out << "\n";
            //    }
            //    out.close();
            //}
            //{
            //    fn = "fn-sync-elapsed-iter-";
            //    fn.append(std::to_string(iteration));
            //    fn.append(".csv");
            //    out.open(fn);
            //    out << "function,elapsed\n";
            //    std::vector<double> total_elap(func_.size(), 0.);
            //    for (auto it = func_.begin(); it != func_.end(); ++it)
            //    {
            //        out << it->first;
            //        out << ",";
            //        out << elapsed_main_[it->second][0];
            //        out << "\n";
            //    }
            //    out.close();
            //}
        }
    }

    std::string Tier::tag() const
    {
        return tag_;
    }

    bool Profiler::exclude(std::string func_name) const
    {
        auto it = std::find(excluded_.begin(), excluded_.end(), func_name);
        if (it != excluded_.end())
        {
            return true;
        }

        return false;
    }

    bool Profiler::include(std::string func_name) const
    {
        auto it = std::find(included_.begin(), included_.end(), func_name);
        if (it != included_.end())
        {
            return true;
        }

        return false;
    }

    void Profiler::update_batch()
    {
        gather_batch();

        if (world_->rank() != 0) {
            return;
        }

        for (auto it=func_.begin(); it!=func_.end(); ++it)
        //for (int g=0; g<func_.size(); ++g)
        {
            int d = std::distance(func_.begin(), it);
            int f = it->second;
            //int f = g * world_.size();

            for (int i=0; i<world_->size(); ++i)
            {
                int k = d + i * func_.size();
                if (world_->rank() == 13) {
                    std::cout << "k: " << k << std::endl;
                }
                //int k = f + i;
                total_elapsed_[i] += elapsed_batch_[k];
                total_elapsed_iter_[i] += elapsed_batch_[k];
                elapsed_main_[f][i] += elapsed_batch_[k];
                elapsed_main_iter_[f][i] += elapsed_batch_[k];
            }
        }
    }

    void Profiler::update()
    {
        for (auto it=func_.begin(); it!=func_.end(); ++it)
        {
            int f = it->second;

            gather(f);

            assert(total_elapsed_.size() == world_->size());
            assert(total_elapsed_iter_.size() == world_->size());
            assert(elapsed_main_.size() == func_.size());
            assert(elapsed_main_iter_.size() == func_.size());
            if (elapsed_main_.size() <= f)
            {
                std::cout << "f = " << f << std::endl;
                std::cout << "size = " << elapsed_main_.size() << std::endl;
            }
            assert(elapsed_main_.size() > f);
            //if (elapsed_main_iter_.size() <= f)
            //{
            //std::cout << "f = " << f << std::endl;
            //std::cout << "size = " << elapsed_main_iter_.size() << std::endl;
            //}
            assert(elapsed_main_iter_.size() > f);
            assert(total_elapsed_.size() == world_->size());
            assert(total_elapsed_iter_.size() == world_->size());
            //if (world_.rank() == MASTER)
            {
                for (int i=0; i<world_->size(); ++i)
                {
                    total_elapsed_[i] += elapsed_[f][i];
                    //if (!exclude(func_name))
                    //if (include(func_name))
                    //{
                    total_elapsed_iter_[i] += elapsed_[f][i];
                    //}
                    elapsed_main_[f][i] += elapsed_[f][i];
                    elapsed_main_iter_[f][i] += elapsed_[f][i];
                }
            }
        }
    }

    void Profiler::stop(std::string func_name)
    {
        if (with_barrier_) {
            world_->barrier();
        }

        auto it = std::find_if(func_.begin(), func_.end(), [&](const auto& e){return e.first == func_name;});
        assert(it != func_.end());
        int f = it->second;

        //assert(stop_.size() == func_.size());
        assert(stop_.size() > f);
        //stop_[f] = timer.start_stop();
        stop_[f] = MPI_Wtime();
        /*gather(f);
        assert(total_elapsed_.size() == world_.size());
        assert(total_elapsed_iter_.size() == world_.size());
        assert(elapsed_main_.size() == func_.size());
        assert(elapsed_main_iter_.size() == func_.size());
        if (elapsed_main_.size() <= f)
        {
            std::cout << "f = " << f << std::endl;
            std::cout << "size = " << elapsed_main_.size() << std::endl;
        }
        assert(elapsed_main_.size() > f);
        //if (elapsed_main_iter_.size() <= f)
        //{
        //std::cout << "f = " << f << std::endl;
        //std::cout << "size = " << elapsed_main_iter_.size() << std::endl;
        //}
        assert(elapsed_main_iter_.size() > f);
        assert(total_elapsed_.size() == world_.size());
        assert(total_elapsed_iter_.size() == world_.size());
        if (world_.rank() == MASTER)
        {
            for (int i=0; i<world_.size(); ++i)
            {
                total_elapsed_[i] += elapsed_[f][i];
                //if (!exclude(func_name))
                //if (include(func_name))
                //{
                    total_elapsed_iter_[i] += elapsed_[f][i];
                //}
                elapsed_main_[f][i] += elapsed_[f][i];
                elapsed_main_iter_[f][i] += elapsed_[f][i];
            }
        }
        //world_.barrier();*/
    }

    void Profiler::gather_batch()
    {
        std::vector<double> in_values;
        in_values.reserve(func_.size());

        //int tf = -1;
        for (auto it=func_.begin(); it!=func_.end(); ++it)
        {
            int f = it->second;
            //duration elap = stop_[f] - start_[f];
            double elap = stop_[f] - start_[f];
            /*if (world_.rank() == 0)
            {
                if (it->first == "det-best-donor")
                {
                    if (elap.count() > 1.)
                    {
                        std::cout << "best-donor elap for master: " << elap.count() << std::endl;
                    }
                    assert(elap.count() < 1.);
                }
            }
            if (world_.rank() == 0)
            {
                if (it->first == "det-cell-type")
                {
                    tf = in_values.size();
                    if (elap.count() > 1e-5)
                    {
                        std::cout << "det-cell elap for master: " << elap.count() << std::endl;
                    }
                    assert(elap.count() < 1e-5);
                }
            }*/
            //in_values.push_back(elap.count());
            in_values.push_back(elap);
        }

        /*if (world_.rank() == 0)
        {
            if (tf == -1)
            {
                for (auto it=func_.begin(); it!=func_.end(); ++it)
                {
                    std::cout << it->first << std::endl;
                }
            }
            assert(tf != -1);
            if (in_values[tf] > 1e-5)
            {
                std::cout << "inval: " << in_values[tf] << std::endl;
            }
            assert(in_values[tf] < 1e-5);
            std::cout << "tttttttttttttttt: " << tf << std::endl;
        }*/

        boost::mpi::gather(*world_, in_values.data(), in_values.size(), elapsed_batch_, 0);
        //boost::mpi::all_gather(world_, in_values.data(), in_values.size(), elapsed_batch_);
        //assert(std::all_of(elapsed_batch_.begin(), elapsed_batch_.end(), [](double d){return d >= 0.;}));

        /*for (int i=0; i<elapsed_batch_.size(); ++i)
        //for (double d: elapsed_batch_)
        {
            std::cout << i << " " << elapsed_batch_[i] << std::endl;
        }*/
    }

    void Profiler::gather(int f)
    {
        assert(stop_.size() > f);
        assert(start_.size() > f);
        assert(elapsed_.size() > f);
        assert(elapsed_.size() == func_.size());
        //duration elap = stop_[f] - start_[f];
        double elap = stop_[f] - start_[f];
        //boost::mpi::gather(world_, elap.count(), elapsed_[f], MASTER);
        //boost::mpi::gather(world_, elap, elapsed_[f], MASTER);
        boost::mpi::all_gather(*world_, elap, elapsed_[f]);


    }
}
