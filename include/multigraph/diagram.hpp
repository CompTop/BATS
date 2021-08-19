#pragma once

#include <cstddef>
#include <vector>
#include <string>
#include <fstream>
#include <sys/stat.h>

namespace bats {

// diagram class - owns data
template <typename NT, typename ET>
class Diagram {
public:
    struct Edge {
        size_t src;
        size_t targ;

        Edge() {}
        //Edge(const Edge& other) : src(other.src), targ(other.targ) {}
        Edge(size_t s, size_t t) : src(s), targ(t) {}
    };

    std::vector<NT> node;
    std::vector<ET> edata;
    std::vector<Edge> elist;

    Diagram() {}

    Diagram(size_t n, size_t m) {
        node.resize(n);
        edata.resize(m);
        elist.resize(m);
    }

    size_t nnode() const { return node.size(); }
    size_t nedge() const { return edata.size(); }

    inline NT& node_data(size_t i) {return node[i];}
    inline const NT& node_data(size_t i) const {return node[i];}
    inline ET& edge_data(size_t j) {return edata[j];}
    inline const ET& edge_data(size_t j) const {return edata[j];}
    inline size_t edge_source(size_t j) const {return elist[j].src;}
    inline size_t edge_target(size_t j) const {return elist[j].targ;}

    size_t add_node(NT &a) {
        node.emplace_back(a);
        return node.size() - 1;
    }

    size_t add_node(NT &&a) {
        node.emplace_back(a);
        return node.size() - 1;
    }

    // set preallocated node i
    void set_node(size_t i, NT &a) {
        node[i] = a;
    }

    // set preallocated node i
    void set_node(size_t i, NT &&a) {
        node[i] = a;
    }

    size_t add_edge(size_t i, size_t j, ET &data) {
        edata.emplace_back(data);
        elist.emplace_back(Edge(i,j));
        return edata.size() - 1;
    }

    size_t add_edge(size_t i, size_t j, ET &&data) {
        edata.emplace_back(data);
        elist.emplace_back(Edge(i,j));
        return edata.size() - 1;
    }

    // set preallocated edge e
    // void set_edge(size_t i, size_t s, size_t t, ET &data) {
    //     edata[i] = data;
    //     elist[i] = Edge(s, t);
    // }
    void set_edge(size_t i, size_t s, size_t t, const ET &data) {
        edata[i] = data;
        elist[i] = Edge(s, t);
    }

    // set preallocated edge e
    // void set_edge(size_t i, size_t s, size_t t, ET &&data) {
    //     edata[i] = data;
    //     elist[i] = Edge(s, t);
    // }

    void save_metadata(std::string &fname) const {
        std::ofstream file (fname, std::ios::out);
        if (file.is_open()) {
            file << node.size() << '\n'; // number of nodes
            file << elist.size() << '\n'; // number of edges
            for (auto e : elist) {
                // write edge directions
                file << e.src << ',' << e.targ << '\n';
            }
            file.close();
        } else {
            std::cerr << "unable to write metadata at " << fname << std::endl;
        }
    }

    void save(std::string &dname) const {
        mkdir(dname.c_str(), 0777);
        // TODO: write metadata on number of nodes, edges, and edge source/target
        std::string mdname = dname + "/metadata";
        save_metadata(mdname);

        #pragma omp parallel for
        for (size_t i = 0; i < nnode(); i++) {
            std::string fname = dname + "/node" + std::to_string(i);
            node[i].save(fname);
        }
        #pragma omp parallel for
        for (size_t i = 0; i < nedge(); i++) {
            std::string fname = dname + "/edge" + std::to_string(i);
            edata[i].save(fname);
        }
    }
};

} // namespace bats
