#include "graph.h"

namespace Tailor
{
    const std::map<int, Graph::Hex*>& Graph::id_to_hex() const
    {
        return id_to_hex_;
    }

    Graph::Vertex& Graph::vertex(int vertexid, int hexid)
    {
        auto it = id_to_hex_.find(hexid);
        assert(it != id_to_hex_.end());
        Hex* hexit = it->second;
        assert(hexit != nullptr);
        assert(hexit->vertex_.size() >= 4);
        for (int i=0; i<nchild_; ++i)
        {
            if (hexit->vertex_[i].id_ == vertexid)
            {
                return hexit->vertex_[i];
            }
        }
        assert(false);
    }

    Graph::Hex::Hex(int id, const std::vector<int>& newvertexid, const std::vector<int>& weight, int nchild): id_(id)
    {
        assert(nchild == 4 || nchild == 8);
        vertex_.resize(nchild);

        vertex_[LLA].parent_hex_ = id_;
        vertex_[LRA].parent_hex_ = id_;
        vertex_[URA].parent_hex_ = id_;
        vertex_[ULA].parent_hex_ = id_;
        if (nchild == 8)
        {
            vertex_[LLF].parent_hex_ = id_;
            vertex_[LRF].parent_hex_ = id_;
            vertex_[URF].parent_hex_ = id_;
            vertex_[ULF].parent_hex_ = id_;
        }

        vertex_[LLA].id_ = newvertexid[0];
        vertex_[LRA].id_ = newvertexid[1];
        vertex_[ULA].id_ = newvertexid[2];
        vertex_[URA].id_ = newvertexid[3];
        if (nchild == 8)
        {
            vertex_[LLF].id_ = newvertexid[4];
            vertex_[LRF].id_ = newvertexid[5];
            vertex_[ULF].id_ = newvertexid[6];
            vertex_[URF].id_ = newvertexid[7];
        }

        vertex_[LLA].pos_ = LLA;
        vertex_[LRA].pos_ = LRA;
        vertex_[URA].pos_ = URA;
        vertex_[ULA].pos_ = ULA;
        if (nchild == 8)
        {
            vertex_[LLF].pos_ = LLF;
            vertex_[LRF].pos_ = LRF;
            vertex_[URF].pos_ = URF;
            vertex_[ULF].pos_ = ULF;
        }

        vertex_[LLA].weight_ = weight[0];
        vertex_[LRA].weight_ = weight[1];
        vertex_[ULA].weight_ = weight[2];
        vertex_[URA].weight_ = weight[3];
        if (nchild == 8)
        {
            vertex_[LLF].weight_ = weight[4];
            vertex_[LRF].weight_ = weight[5];
            vertex_[ULF].weight_ = weight[6];
            vertex_[URF].weight_ = weight[7];
        }
        assert(vertex_.size() >= 4);
        assert(vertex_.size() == 4 || vertex_.size() == 8);
    }

    void Graph::reset_connection(Hex& hex)
    {
        for (int i=0; i<nchild_; ++i)
        {
            hex.vertex_[i].nei_.clear();
        }

        for (int i=0; i<nchild_; ++i)
        {
            if (hex.vertex_[i].hex_) {
                reset_connection(*(hex.vertex_[i].hex_.get()));
            }
        }
    }

    void Graph::reset_connection()
    {
        reset_connection(hex_);
    }

    Graph::Graph(const std::vector<int>& vertex_id, const std::vector<int>& weight, int nchild): current_max_id_(0), hex_(0, vertex_id, weight, nchild), nedge_(0), nchild_(nchild)
    {
        assert(nchild == 4 || nchild == 8);
        // weights are assigned in this order: LLB, LRB, ULB, URB, LLF, LRF, ULF, URF
        id_to_hex_.insert(std::make_pair(0, &hex_));
    }

    void Graph::connect(Hex* hex, int rank)
    {
        if (hex == nullptr) {
            return;
        }

        assert(hex->vertex_.size() >= 4);
        assert(hex->vertex_.size() == 4 || hex->vertex_.size() == 8);

        if (hex->vertex_[LLA].hex_ != nullptr) {
            assert(hex->vertex_[LLA].hex_->vertex_.size() == 4 || hex->vertex_[LLA].hex_->vertex_.size() == 8);
            connect(hex->vertex_[LLA].hex_.get(), rank);
        }
        if (hex->vertex_[LRA].hex_ != nullptr) {
            assert(hex->vertex_[LRA].hex_->vertex_.size() == 4 || hex->vertex_[LRA].hex_->vertex_.size() == 8);
            connect(hex->vertex_[LRA].hex_.get(), rank);
        }
        if (hex->vertex_[ULA].hex_ != nullptr) {
            assert(hex->vertex_[ULA].hex_->vertex_.size() == 4 || hex->vertex_[ULA].hex_->vertex_.size() == 8);
            connect(hex->vertex_[ULA].hex_.get(), rank);
        }
        if (hex->vertex_[URA].hex_ != nullptr) {
            assert(hex->vertex_[URA].hex_->vertex_.size() == 4 || hex->vertex_[URA].hex_->vertex_.size() == 8);
            connect(hex->vertex_[URA].hex_.get(), rank);
        }
        if (nchild_ == 8)
        {
            if (hex->vertex_[LLF].hex_ != nullptr) {
                assert(hex->vertex_[LLF].hex_->vertex_.size() == 4 || hex->vertex_[LLF].hex_->vertex_.size() == 8);
                connect(hex->vertex_[LLF].hex_.get(), rank);
            }
            if (hex->vertex_[LRF].hex_ != nullptr) {
                assert(hex->vertex_[LRF].hex_->vertex_.size() == 4 || hex->vertex_[LRF].hex_->vertex_.size() == 8);
                connect(hex->vertex_[LRF].hex_.get(), rank);
            }
            if (hex->vertex_[ULF].hex_ != nullptr) {
                assert(hex->vertex_[ULF].hex_->vertex_.size() == 4 || hex->vertex_[ULF].hex_->vertex_.size() == 8);
                connect(hex->vertex_[ULF].hex_.get(), rank);
            }
            if (hex->vertex_[URF].hex_ != nullptr) {
                assert(hex->vertex_[URF].hex_->vertex_.size() == 4 || hex->vertex_[URF].hex_->vertex_.size() == 8);
                connect(hex->vertex_[URF].hex_.get(), rank);
            }
        }

        if (hex->vertex_[LLA].hex_ != nullptr) {
            assert(hex->vertex_[LLA].hex_->vertex_.size() == 4 || hex->vertex_[LLA].hex_->vertex_.size() == 8);
        }
        if (hex->vertex_[LRA].hex_ != nullptr) {
            assert(hex->vertex_[LRA].hex_->vertex_.size() == 4 || hex->vertex_[LRA].hex_->vertex_.size() == 8);
        }
        if (hex->vertex_[URA].hex_ != nullptr) {
            assert(hex->vertex_[URA].hex_->vertex_.size() == 4 || hex->vertex_[URA].hex_->vertex_.size() == 8);
        }
        if (hex->vertex_[ULA].hex_ != nullptr) {
            assert(hex->vertex_[ULA].hex_->vertex_.size() == 4 || hex->vertex_[ULA].hex_->vertex_.size() == 8);
        }

        connect(hex->vertex_[LLA], hex->vertex_[LRA], L, rank);
        connect(hex->vertex_[LRA], hex->vertex_[URA], B, rank);
        connect(hex->vertex_[URA], hex->vertex_[ULA], R, rank);
        connect(hex->vertex_[ULA], hex->vertex_[LLA], T, rank);
        if (nchild_ == 8)
        {
            if (hex->vertex_[LLF].hex_ != nullptr) {
                assert(hex->vertex_[LLF].hex_->vertex_.size() == 4 || hex->vertex_[LLF].hex_->vertex_.size() == 8);
            }
            if (hex->vertex_[LRF].hex_ != nullptr) {
                assert(hex->vertex_[LRF].hex_->vertex_.size() == 4 || hex->vertex_[LRF].hex_->vertex_.size() == 8);
            }
            if (hex->vertex_[URF].hex_ != nullptr) {
                assert(hex->vertex_[URF].hex_->vertex_.size() == 4 || hex->vertex_[URF].hex_->vertex_.size() == 8);
            }
            connect(hex->vertex_[LLF], hex->vertex_[LRF], L, rank);
            connect(hex->vertex_[LRF], hex->vertex_[URF], B, rank);
            connect(hex->vertex_[URF], hex->vertex_[ULF], R, rank);
            connect(hex->vertex_[ULF], hex->vertex_[LLF], T, rank);
            connect(hex->vertex_[LLA], hex->vertex_[LLF], A, rank);
            connect(hex->vertex_[LRA], hex->vertex_[LRF], A, rank);
            connect(hex->vertex_[ULA], hex->vertex_[ULF], A, rank);
            connect(hex->vertex_[URA], hex->vertex_[URF], A, rank);
        }
    }

    void Graph::connect(int rank)
    {
        vertex_.clear();
        connect(&hex_, rank);
    }

    void Graph::connect(Vertex& a, Vertex& b, Pos pos, int rank)
    {
        if (a.hex_ == nullptr)
        {
            if (b.hex_ == nullptr)
            {
                a.nei_.push_back(&b);
                b.nei_.push_back(&a);
                ++nedge_;
                //if (a.weight_ != 0.)
                {
                    auto it = std::find_if(vertex_.begin(), vertex_.end(), [&](const Vertex* v){return v->id_ == a.id_;});
                    if (it == vertex_.end()) {
                        vertex_.push_back(&a);
                    }
                }
                //if (b.weight_ != 0.)
                {
                    auto it = std::find_if(vertex_.begin(), vertex_.end(), [&](const Vertex* v){return v->id_ == b.id_;});
                    if (it == vertex_.end()) {
                        vertex_.push_back(&b);
                    }
                }
            }
            else
            {
                assert(b.hex_->vertex_.size() == 4 || b.hex_->vertex_.size() == 8);
                a.nei_.push_back(&b);
                assert(b.hex_ != nullptr);
                if (b.hex_->vertex_.size() < 4)
                {
                    std::cout <<"vertex size:"  << b.hex_->vertex_.size() << std::endl;
                }
                assert(b.hex_->vertex_.size() >= 4);
                if (pos == L)
                {
                    connect(a, b.hex_->vertex_[LLA], pos, rank);
                    connect(a, b.hex_->vertex_[ULA], pos, rank);
                    if (nchild_ == 8)
                    {
                        connect(a, b.hex_->vertex_[LLF], pos, rank);
                        connect(a, b.hex_->vertex_[ULF], pos, rank);
                    }
                }
                else if (pos == B)
                {
                    connect(a, b.hex_->vertex_[LLA], pos, rank);
                    connect(a, b.hex_->vertex_[LRA], pos, rank);
                    if (nchild_ == 8)
                    {
                        connect(a, b.hex_->vertex_[LLF], pos, rank);
                        connect(a, b.hex_->vertex_[LRF], pos, rank);
                    }
                }
                else if (pos == R)
                {
                    assert(b.hex_->vertex_.size() == 4 || b.hex_->vertex_.size() == 8);
                    connect(a, b.hex_->vertex_[LRA], pos, rank);
                    connect(a, b.hex_->vertex_[URA], pos, rank);
                    if (nchild_ == 8)
                    {
                        assert(b.hex_->vertex_.size() == 8);
                        connect(a, b.hex_->vertex_[LRF], pos, rank);
                        connect(a, b.hex_->vertex_[URF], pos, rank);
                    }
                }
                else if (pos == T)
                {
                    connect(a, b.hex_->vertex_[ULA], pos, rank);
                    connect(a, b.hex_->vertex_[URA], pos, rank);
                    if (nchild_ == 8)
                    {
                        connect(a, b.hex_->vertex_[ULF], pos, rank);
                        connect(a, b.hex_->vertex_[URF], pos, rank);
                    }
                }
                else if (pos == A)
                {
                    assert(nchild_ == 8);
                    connect(a, b.hex_->vertex_[LLA], pos, rank);
                    connect(a, b.hex_->vertex_[LRA], pos, rank);
                    connect(a, b.hex_->vertex_[ULA], pos, rank);
                    connect(a, b.hex_->vertex_[URA], pos, rank);
                }
                else if (pos == F)
                {
                    assert(nchild_ == 8);
                    connect(a, b.hex_->vertex_[LLF], pos, rank);
                    connect(a, b.hex_->vertex_[LRF], pos, rank);
                    connect(a, b.hex_->vertex_[ULF], pos, rank);
                    connect(a, b.hex_->vertex_[URF], pos, rank);
                }
            }
        }
        else
        {
            if (a.hex_->vertex_.size() != 4 && a.hex_->vertex_.size() != 8)
            {
                std::cout << rank << " vertex size: " << a.hex_->vertex_.size() << std::endl;
                std::cout << rank << " pos: " << pos << std::endl;
            }
            assert(a.hex_->vertex_.size() == 4 || a.hex_->vertex_.size() == 8);
            if (b.hex_ == nullptr)
            {
                if (pos == L)
                {
                    connect(a.hex_->vertex_[LRA], b, pos, rank);
                    connect(a.hex_->vertex_[URA], b, pos, rank);
                    if (nchild_ == 8)
                    {
                        connect(a.hex_->vertex_[LRF], b, pos, rank);
                        connect(a.hex_->vertex_[URF], b, pos, rank);
                    }
                }
                else if (pos == B)
                {
                    connect(a.hex_->vertex_[ULA], b, pos, rank);
                    connect(a.hex_->vertex_[URA], b, pos, rank);
                    if (nchild_ == 8)
                    {
                        connect(a.hex_->vertex_[ULF], b, pos, rank);
                        connect(a.hex_->vertex_[URF], b, pos, rank);
                    }
                }
                else if (pos == R)
                {
                    assert(a.hex_->vertex_.size() == 4 || a.hex_->vertex_.size() == 8);
                    connect(a.hex_->vertex_[ULA], b, pos, rank);
                    connect(a.hex_->vertex_[LLA], b, pos, rank);
                    if (nchild_ == 8)
                    {
                        assert(a.hex_->vertex_.size() == 8);
                        connect(a.hex_->vertex_[ULF], b, pos, rank);
                        connect(a.hex_->vertex_[LLF], b, pos, rank);
                    }
                }
                else if (pos == T)
                {
                    connect(a.hex_->vertex_[LLA], b, pos, rank);
                    connect(a.hex_->vertex_[LRA], b, pos, rank);
                    if (nchild_ == 8)
                    {
                        connect(a.hex_->vertex_[LLF], b, pos, rank);
                        connect(a.hex_->vertex_[LRF], b, pos, rank);
                    }
                }
                else if (pos == A)
                {
                    assert(nchild_ == 8);
                    connect(a.hex_->vertex_[LLF], b, pos, rank);
                    connect(a.hex_->vertex_[LRF], b, pos, rank);
                    connect(a.hex_->vertex_[ULF], b, pos, rank);
                    connect(a.hex_->vertex_[URF], b, pos, rank);
                }
                else if (pos == F)
                {
                    assert(nchild_ == 8);
                    connect(a.hex_->vertex_[LLA], b, pos, rank);
                    connect(a.hex_->vertex_[LRA], b, pos, rank);
                    connect(a.hex_->vertex_[ULA], b, pos, rank);
                    connect(a.hex_->vertex_[URA], b, pos, rank);
                }
            }
            else
            {
                if (pos == L)
                {
                    assert(a.hex_ != nullptr);
                    assert(b.hex_ != nullptr);
                    assert(a.hex_->vertex_.size() >= 4);
                    connect(a.hex_->vertex_[LRA], b.hex_->vertex_[LLA], pos, rank);
                    connect(a.hex_->vertex_[URA], b.hex_->vertex_[ULA], pos, rank);
                    if (nchild_ == 8)
                    {
                        connect(a.hex_->vertex_[LRF], b.hex_->vertex_[LLF], pos, rank);
                        connect(a.hex_->vertex_[URF], b.hex_->vertex_[ULF], pos, rank);
                    }
                }
                else if (pos == B)
                {
                    assert(a.hex_ != nullptr);
                    assert(b.hex_ != nullptr);
                    assert(a.hex_->vertex_.size() >= 4);
                    assert(b.hex_->vertex_.size() >= 4);
                    connect(a.hex_->vertex_[ULA], b.hex_->vertex_[LLA], pos, rank);
                    connect(a.hex_->vertex_[URA], b.hex_->vertex_[LRA], pos, rank);
                    if (nchild_ == 8)
                    {
                        connect(a.hex_->vertex_[ULF], b.hex_->vertex_[LLF], pos, rank);
                        connect(a.hex_->vertex_[URF], b.hex_->vertex_[LRF], pos, rank);
                    }
                }
                else if (pos == R)
                {
                    assert(a.hex_ != nullptr);
                    assert(b.hex_ != nullptr);
                    assert(a.hex_->vertex_.size() == 4 || a.hex_->vertex_.size() == 8);
                    assert(b.hex_->vertex_.size() == 4 || b.hex_->vertex_.size() == 8);
                    connect(a.hex_->vertex_[LLA], b.hex_->vertex_[LRA], pos, rank);
                    connect(a.hex_->vertex_[ULA], b.hex_->vertex_[URA], pos, rank);
                    if (nchild_ == 8)
                    {
                        connect(a.hex_->vertex_[LLF], b.hex_->vertex_[LRF], pos, rank);
                        connect(a.hex_->vertex_[ULF], b.hex_->vertex_[URF], pos, rank);
                    }
                }
                else if (pos == T)
                {
                    connect(a.hex_->vertex_[LLA], b.hex_->vertex_[ULA], pos, rank);
                    connect(a.hex_->vertex_[LRA], b.hex_->vertex_[URA], pos, rank);
                    if (nchild_ == 8)
                    {
                        connect(a.hex_->vertex_[LLF], b.hex_->vertex_[ULF], pos, rank);
                        connect(a.hex_->vertex_[LRF], b.hex_->vertex_[URF], pos, rank);
                    }
                }
                else if (pos == A)
                {
                    assert(nchild_ == 8);
                    connect(a.hex_->vertex_[LLF], b.hex_->vertex_[LLA], pos, rank);
                    connect(a.hex_->vertex_[LRF], b.hex_->vertex_[LRA], pos, rank);
                    connect(a.hex_->vertex_[ULF], b.hex_->vertex_[ULA], pos, rank);
                    connect(a.hex_->vertex_[URF], b.hex_->vertex_[URA], pos, rank);
                }
                else if (pos == F)
                {
                    assert(nchild_ == 8);
                    connect(a.hex_->vertex_[LLA], b.hex_->vertex_[LLF], pos, rank);
                    connect(a.hex_->vertex_[LRA], b.hex_->vertex_[LRF], pos, rank);
                    connect(a.hex_->vertex_[ULA], b.hex_->vertex_[ULF], pos, rank);
                    connect(a.hex_->vertex_[URA], b.hex_->vertex_[URF], pos, rank);
                }
            }
        }
    }

    void Graph::refine(int vertexid, int hexid, int newhexid, const std::vector<int>& newvertexid, const std::vector<int>& weight)
    {
        reset_connection();
        Vertex& vtx = vertex(vertexid, hexid);
        assert(vtx.id_ == vertexid);
        
        auto it = id_to_hex_.find(hexid);
        assert(it != id_to_hex_.end());
        auto hexit = it->second;

        assert(vtx.hex_ == nullptr);
        vtx.hex_ = std::make_unique<Hex>(Hex(newhexid, newvertexid, weight, nchild_));
        //std::cout << "vtx.quad_->id_: " << vtx.quad_->id_ << std::endl;
        //std::cout << "w0: " << w0 << std::endl;
        //std::cout << "w1: " << w1 << std::endl;
        //std::cout << "w2: " << w2 << std::endl;
        //std::cout << "w3: " << w3 << std::endl;
        //std::cout << "weight0: " << vtx.quad_->vertex_[LL].weight_ << std::endl;
        //std::cout << "weight1: " << vtx.quad_->vertex_[LR].weight_ << std::endl;
        //std::cout << "weight2: " << vtx.quad_->vertex_[UL].weight_ << std::endl;
        //std::cout << "weight2: " << vtx.quad_->vertex_[UR].weight_ << std::endl;
        //id_to_quad_.insert(std::make_pair(newquadid->first+1, vtx.quad_.get()));
        id_to_hex_.insert(std::make_pair(newhexid, vtx.hex_.get()));
    }

    void Graph::print()
    {
        std::cout << vertex_.size() << " " << nedge_ << std::endl;
        for (Vertex* v: vertex_)
        {
            std::cout << v->id_;
            for (Vertex* n: v->nei_)
            {
                std::cout << " ";
                std::cout << n->id_;
            }
            std::cout << std::endl;
        }
    }

    void Graph::print_to_file()
    {
        std::ofstream out;
        std::string fn = "graph.dat";
        out.open(fn);

        for (int i=0; i<group_.size(); ++i)
        {
            out << i << " " << group_[i] << std::endl;;
        }

        out.close();
    }
    
    bool Graph::partition(idx_t nparts, double& dev, std::vector<std::vector<BinRMTag>>& aug_bin_tag, int iter, int rank)
    {
        bool balance = true;

        idx_t nvtxs = vertex_.size();
        idx_t ncon = 1;
        //real_t ubvec[ncon];
        //ubvec[0] = 1.1;
        //ubvec[0] = 1.001;
        //ubvec[0] = 1.01;

        //idx_t ncon = 2; // multi-constraint: 1) load of bin, 2) area of bin

        idx_t xadj[nvtxs+1];
        idx_t adjncy[nedge_ * 2];
        idx_t vwgt[nvtxs * ncon];
        xadj[0] = 0;
        int i=0;
        int j=0;
        for (Vertex* v: vertex_)
        {
            vwgt[i] = v->weight_;
            //vwgt[i+1] = v->area_;
            ++i;
            //i += ncon;
        }
        i=1;
        for (Vertex* v: vertex_)
        {
            int neisize = 0;
            for (Vertex* n: v->nei_)
            {
                assert(n != nullptr);
                auto it = std::find_if(vertex_.begin(), vertex_.end(), [&](const Vertex* v){return v->id_ == n->id_;});
                if (it == vertex_.end()) {
                    continue;
                }
                adjncy[j] = std::distance(vertex_.begin(), it);
                ++j;
                ++neisize;
            }
            xadj[i] = xadj[i-1] + neisize;
            ++i;
        }

        /*std::cout << "nvtxs: " << nvtxs << std::endl;
        std::cout << "ncon: " << ncon << std::endl;
        std::cout << "nparts: " << nparts << std::endl;

        std::cout << "xadj: ";
        for (int i=0; i<nvtxs+1; ++i) {
            std::cout << xadj[i] << " ";
        }

        std::cout << std::endl;

        std::cout << "adjncy: ";
        for (int i=0; i<nedge_*2; ++i) {
            std::cout << adjncy[i] << " ";
        }

        std::cout << std::endl;

        std::cout << "vwgt: ";
        for (int i=0; i<nvtxs; ++i) {
            std::cout << vwgt[i] << " ";
        }

        std::cout << std::endl;*/

        idx_t obj_val;
        idx_t part[nvtxs];


        //std::cout << "nvertices: " << vertex_.size() << std::endl;
        //std::cout << "nparts: " << nparts << std::endl;

        idx_t options[METIS_NOPTIONS];
        METIS_SetDefaultOptions(options);
        options[METIS_OPTION_CONTIG] = 1;
        //options[METIS_OPTION_UFACTOR] = 1;
        //options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;
        //options[METIS_OPTION_NCUTS] = 40;
        //options[METIS_OPTION_NSEPS] = 10;
        //options[METIS_OPTION_DBGLVL] = METIS_DBG_INFO;

        int res = METIS_PartGraphKway(&nvtxs, &ncon, xadj, adjncy, vwgt, NULL, NULL, &nparts, NULL, NULL, options, &obj_val, part);
        std::string resstr;
        if (res == METIS_OK) {
            resstr = "METIS_OK";
        }
        else if (res == METIS_ERROR_INPUT) {
            resstr = "METIS_ERROR_INPUT";
            balance = false;
        }
        else if (res == METIS_ERROR_MEMORY) {
            resstr = "METIS_ERROR_MEMORY";
            balance = false;
        }
        else if (res == METIS_ERROR) {
            resstr = "METIS_ERROR";
            balance = false;
        }
        else
        {
            assert(false);
        }

        //std::ofstream out;
        //if (rank == 0)
        //{
            //std::string fn = "metis";
            //fn.append(std::to_string(iter));
            //fn.append(".dat");
            //out.open(fn);
        //}

        //std::vector<int> group(nparts, 0);
        group_.clear();
        group_ = std::vector<int>(nparts, 0);

        for (int i=0; i<nvtxs; ++i) {
            group_[part[i]] += vertex_[i]->weight_;
        }

        auto result = std::minmax_element(group_.begin(), group_.end());
        double minload = *result.first;
        double maxload = *result.second;
        if (maxload == 0)
        {
            for (int i=0; i<nvtxs; ++i) {
                std::cout << vertex_[i]->parent_hex_ << " " << vertex_[i]->id_ << " " << vertex_[i]->weight_ << " " << part[i] << std::endl;
            }
        }
        assert(maxload != 0);
        dev = ((maxload - minload) / maxload) * 100.;
        //std::cout << "min: " << minload << std::endl;
        //std::cout << "max: " << maxload << std::endl;
        //std::cout << "dev: " << dev << std::endl;
        //std::cout << "balance: " << balance << std::endl;

        //if (rank == 0)
        //{
        //    out << resstr << std::endl;
        //    out << minload << std::endl;
        //    out << maxload << std::endl;
        //    out << dev << std::endl;
        //    for (int i=0; i<nvtxs; ++i) {
        //        out << vertex_[i]->parent_hex_ << " " << vertex_[i]->id_ << " " << vertex_[i]->weight_ << " " << part[i] << std::endl;
        //    }
        //    for (int i=0; i<nparts; ++i) {
        //        out << "group[" << i << "]: " << group[i] << std::endl;
        //    }
        //    out.close();
        //}

        aug_bin_tag.clear();
        aug_bin_tag.resize(nparts);
        for (int i=0; i<nvtxs; ++i)
        {
            //std::cout << "rm: " <<  vertex_[i]->parent_quad_ << std::endl;
            //std::cout << "bt: " <<  vertex_[i]->id_<< std::endl;
            BinRMTag tag(Tag(vertex_[i]->id_), Tag(vertex_[i]->parent_hex_));
            aug_bin_tag[part[i]].push_back(tag);
        }
        assert(!aug_bin_tag.empty());

        /*if (balance)
        {
            std::cout << "aug size: " << aug_bin_tag.size() << std::endl;
            for (int t=0; t<aug_bin_tag.size(); ++t)
            {
                std::cout << "t: " << t << std::endl;
                for (const auto& tt: aug_bin_tag[t])
                {
                    std::cout << "aug[" << t << "].rm: " << tt.rmtag()() << std::endl;
                    std::cout << "aug[" << t << "].bt: " << tt.bintag()() << std::endl;
                }
            }
        }*/

        //std::cout << "done with partitioning" << std::endl;
        return balance;
    }
}
