#pragma once
/*
Witness complex
*/
#include <vector>
#include "data.h"
#include "neighborhood.h"
#include <util/sorted.h>


template <typename T, typename M>
std::vector<size_t> witness_edges(
    const DataSet<T> &X,
    const DataSet<T> &L,
    const M &dist
) {
    size_t n = L.size();
    std::vector<size_t> edges;
    size_t nedges = (n * (n-1)) / 2;
    edges.reserve(2*nedges);

    // get pairwise distances between L and X
    auto D = dist(L, X);

    // for each column of D, see which edge X[j] witnesses

    return edges;
}

// for weak witness complex, find witnessed edge set for parameterized complex
// template over data type and metric
template <typename T, typename M>
std::vector<size_t> witness_edges(
    const DataSet<T> &X,
    const DataSet<T> &L,
    const M &dist,
    const T rmax
) {

    size_t n = L.size();
    std::vector<size_t> edges;
    size_t nedges = (n * (n-1)) / 2;
    edges.reserve(2*nedges);

    auto nbrs = neighborhoods(L, X, dist, rmax);
    std::vector<size_t> intersection;


    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < i; j++) {

            intersect_sorted(nbrs[i], nbrs[j], intersection);
            if (intersection.size() > 0) {
                edges.push_back(j);
                edges.push_back(i);
            }
        }
    }
    return edges;
}

// for strong witness complex, build list of witnesses for simplices
