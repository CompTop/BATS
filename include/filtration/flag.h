#pragma once
#include <vector>
#include <utility> // make_pair
#include <tuple>
#include <complex/simplicial_complex.h>
#include <filtration/filtration.h>
#include <util/sorted.h>
#include <iostream>
// flag filtrations

// template over filtration type
// TODO: can do unsafe simplex add
template <typename T>
void add_dimension_recursive_flag(
    Filtration<T, SimplicialComplex> &F,
    const std::vector<std::vector<size_t>> &nbrs, // lists of neighbors
    const size_t d, // dimension
    const size_t maxd, // max dimension
    const std::vector<size_t> &iter_idxs,
    std::vector<size_t> &spx_idxs,
    const T &t
) {
    // sorted simplices will end up here
    std::vector<size_t> spx_idxs2(spx_idxs.size() + 1);
    if (d == maxd) {
        // no recursion
        for (auto k : iter_idxs) {
            // append k to spx_idxs, sort
            spx_idxs.push_back(k);
            sort_into(spx_idxs, spx_idxs2);

            // add to F
            F.add_pair(t, spx_idxs2);

            // pop k off spx_idxs
            spx_idxs.pop_back();
        }
    } else { // d < maxd
        // recursion
        std::vector<size_t> iter_idxs2; // indices for recursing on
        iter_idxs2.reserve(iter_idxs.size());
        for (auto k : iter_idxs) {
            // append k to spx_idxs, sort
            spx_idxs.push_back(k);
            sort_into(spx_idxs, spx_idxs2);

            // add to F
            F.add_pair(t, spx_idxs2);

            // recurse
            intersect_sorted_lt(iter_idxs, nbrs[k], k, iter_idxs2);
            add_dimension_recursive_flag(F, nbrs, d+1, maxd, iter_idxs2, spx_idxs2, t);

            // pop k off spx_idxs
            spx_idxs.pop_back();
        }
    }
}


// Flag complex using list of edges
// (edges[2*k], edges[2*k+1]) = (i, j) is an edge
// t - vector of filtration times
// t0 - time for 0-simplices
// n - number of vertices
// maxdim - maximum dimension of simplices
template <typename T>
std::tuple<SimplicialComplex, Filtration<T, SimplicialComplex>> FlagFiltration(
    const std::vector<size_t> &edges,
    const std::vector<T> &t,
    const size_t n, // number of 0-cells
    const size_t maxdim,
    const T t0 = T(0)
) {

    // check that dimensions agree
    size_t m = t.size();
    assert (edges.size() == 2 * m);

    SimplicialComplex X(maxdim);
    Filtration<T, SimplicialComplex> F(X);

    // sets 0-cells
    std::vector<size_t> spx_idxs(1);
    for (size_t k = 0; k < n; k++) {
        spx_idxs[0] = k;
        F.add(t0, spx_idxs); // don't need to try to pair since 0-cells
    }

    std::vector<std::vector<size_t>> nbrs(n);

    spx_idxs.resize(2); // now time to add edges
    std::vector<size_t> iter_idxs;
    iter_idxs.reserve(n); // maximum size

    for (size_t k = 0; k < m; k++) {
        size_t i = edges[2*k];
        size_t j = edges[2*k + 1];
        spx_idxs[0] = i;
        spx_idxs[1] = j;
        F.add_pair(t[k], spx_idxs);

        intersect_sorted(nbrs[i], nbrs[j], iter_idxs);

        nbrs[i].emplace_back(j);
        std::sort(nbrs[i].begin(), nbrs[i].end());
        nbrs[j].emplace_back(i);
        std::sort(nbrs[j].begin(), nbrs[j].end());

        add_dimension_recursive_flag(F, nbrs, 2, maxdim, iter_idxs, spx_idxs, t[k]);
    }

    return std::make_tuple(X, F);
}
