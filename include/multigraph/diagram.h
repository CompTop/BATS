#pragma once

#include <cstddef>
#include <vector>

// diagram class - owns data
template <typename NT, typename ET>
class Diagram {
public:
    struct Edge {
        size_t src;
        size_t targ;

        Edge() {}
        Edge(const Edge& other) : src(other.src), targ(other.targ) {}
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
    void set_edge(size_t i, size_t s, size_t t, ET &data) {
        edata[i] = data;
        elist[i] = Edge(s, t);
    }

    // set preallocated edge e
    void set_edge(size_t i, size_t s, size_t t, ET &&data) {
        edata[i] = data;
        elist[i] = Edge(s, t);
    }
};
