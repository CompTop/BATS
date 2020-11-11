#pragma once

#include <cmath>
#include <iostream>
#include <util/sorted.hpp>

namespace bats {

// Struct for holding edge with a value
template <typename TF, typename TI>
struct tedge {
    TF v; // value
    TI s; // source
    TI t; // target

    tedge() {};
    tedge(TF v, TI s, TI t) : v(v), s(s), t(t) {};
};

template <typename TF, typename TI>
inline tedge<TF, TI> make_edge(TI s, TI t, TF v) {
  return tedge<TF, TI> {v, s, t};
}

template <typename TF, typename TI>
inline bool operator<(const tedge<TF, TI> &a, const tedge<TF, TI> &b) {
  return a.v < b.v;
}

// sort edges and value by value
template <typename T>
void sort_edges(std::vector<size_t> &edges, std::vector<T> &v) {
    size_t m = v.size();
    std::vector<tedge<T, size_t>> tmp(m);
    // fill tmp vector;
    for (size_t i = 0; i < m; i++) {
        tmp[i] = tedge(v[i], edges[2*i], edges[2*i + 1]);
    }
    std::sort(
        tmp.begin(),
        tmp.end()
    );
    for (size_t i = 0; i < m; i++) {
        v[i] = tmp[i].v;
        edges[2*i] = tmp[i].s;
        edges[2*i + 1] = tmp[i].t;
    }
}



template <typename TF, typename TI>
std::ostream& operator<<( std::ostream& os, tedge<TF, TI> &x) {
  os << "(" << x.s << ", " << x.t << "): " << x.v;
  return os;
}

// construct edges for rips complex
template <typename T>
void rips_edges(std::vector<T> &x, std::vector<size_t> &edges, std::vector<T> &t) {
    edges.clear();
    size_t nedges = x.size() * (x.size() - 1) / 2;
    edges.reserve(2*nedges);
    t.clear();
    t.reserve(nedges);
    for (size_t i = 0; i < x.size(); i++) {
        for (size_t j = 0; j < i; j++) {
            edges.push_back(j);
            edges.push_back(i);
            t.push_back(std::abs(x[i] - x[j]));
        }
    }
    // TODO: sort by length
    return;
}

// returns all pairs on n vertices
std::vector<size_t> all_pairs(const size_t n) {
    size_t m = ((n-1)*n) / 2;
    std::vector<size_t> edges;
    edges.reserve(m);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < i; j++) {
            edges.push_back(i);
            edges.push_back(j);
        }
    }
    return edges;
}

// euclidean distance in d dimensions
template <typename T, typename TI>
T euclidean(TI a, TI b, size_t d){
    T res = 0;
    for (size_t i = 0; i < d; i++) {
        res += std::pow((*a++ - *b++), 2);
    }
    return std::sqrt(res);
}

// pairwise distances on edges
// todo: take in function
template <typename T>
std::vector<T> pairwise_dist(
    const std::vector<T> &x,
    const std::vector<size_t> &edges,
    const size_t d
) {
    size_t m = (edges.size() / 2);
    std::vector<T> v(m);
    auto vit = v.begin();
    auto eit = edges.cbegin();
    while (vit != v.end()) {
        *vit = euclidean<T>(x.cbegin() + d*(*eit), x.cbegin() + d*(*(eit + 1)), d);
        //std::cout << *vit << std::endl;
        vit++;
        eit+=2;
    }
    return v;
}

} // namespace bats
