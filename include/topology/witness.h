#pragma once
/*
Witness complex
*/
#include <vector>
#include "data.h"
#include "neighborhood.h"
#include <util/sorted.h>

// for weak witness complex, find witnessed edge set
// template over data type and metric
template <typename T, typename M>
std::vector<size_t> weak_witness_edges(
    const DataSet<T> &X,
    const DataSet<T> &L,
    const M &dist,
    const T rmax
) {

    size_t n = L.size();
    std::vector<size_t> edges;
    size_t nedges = (n * (n-1)) / 2;
    edges.reserve(2*nedges);

    nbrs = neighborhoods(L, X, dist, rmax);
    std::vector<size_t> intersection;


    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < i; j++) {

            intersect_sorted(nrbs[i], nbrs[j], intersection);
            if (intersection.size() > 0) {
                edges.push_back(j);
                edges.push_back(i);
            }
        }
    }
    return edges;
}


// for strong witness complex, build list of witnesses for simplices
