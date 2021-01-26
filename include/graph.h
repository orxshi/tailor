#ifndef GRAPH_H
#define GRAPH_H

#include <map>
#include <vector>
#include <memory>
#include <set>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "metis.h"
#include "tag.h"

namespace Tailor
{
    class Graph
    {
        struct Hex;

        struct Vertex
        {
            int id_;
            int parent_hex_;
            int pos_;
            int weight_;
            std::unique_ptr<Hex> hex_;
            std::vector<Vertex*> nei_;
        };

        struct Hex
        {
            int id_;
            std::vector<Vertex> vertex_;
            //Vertex vertex_[TAILOR_GRAPH_NCHILD];

            Hex(int id, const std::vector<int>& newvertexid, const std::vector<int>& weight, int nchild);
        };

        public:

            Graph(const std::vector<int>& vertex_id, const std::vector<int>& weight, int nchild);
            void refine(int vertexid, int hexid, int newhexid, const std::vector<int>& newvertexid, const std::vector<int>& weight);
            void connect(int rank);
            void print();
            void print_to_file();
            bool partition(idx_t nparts, double& dev, std::vector<std::vector<BinRMTag>>& aug_bin_tag, int iter, int rank);
            const std::map<int, Hex*>& id_to_hex() const;

        private:

            enum VertexPos
            {
                // B stands for back
                // F stands for front

                LLA=0,
                LRA=1,
                ULA=2,
                URA=3,

                LLF=4,
                LRF=5,
                ULF=6,
                URF=7
            };

            enum Pos
            {
                L, // left
                B, // bottom
                R, // right
                T, // top
                F, // front
                A // back
            };


            std::vector<int> group_;
            int nchild_;
            Hex hex_;
            std::map<int, Hex*> id_to_hex_;
            int nedge_;
            int current_max_id_;
            std::vector<Vertex*> vertex_;

            void connect(Vertex& a, Vertex& b, Pos pos, int rank);
            void connect(Hex* hex, int rank);
            Vertex& vertex(int vertexid, int hexid); // vertexpos such as LL.
            void reset_connection(Hex& hex);
            void reset_connection();
    };
}

#endif
