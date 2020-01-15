#pragma once
#include <vector>
#include <set>
#include <utility> // make_pair
#include <algorithm> // sort
#include <tuple>
#include <complex/simplicial_complex.h>
#include <filtration/filtration.h>
#include <util/sorted.h>
#include <iostream>
// flag filtrations

// struct for filtered edges
template <typename T>
struct filtered_edge {
    size_t s; // source
    size_t t; // target
    T r;      // parameter

    filtered_edge() {}
    filtered_edge(const size_t &s, const size_t &t, const T &r) : s(s), t(t), r(r) {}

    // sort by r, then s, then t
    bool operator<(const filtered_edge& other) const {
        return ((r < other.r) ? true :
                (r == other.r && s < other.s) ? true :
                (r == other.r && s == other.s && t < other.t) ? true :
                false
        );
    }
};

// template over filtration type
// TODO: can do unsafe simplex add
template <typename T, typename NT>
void add_dimension_recursive_flag(
    Filtration<T, SimplicialComplex> &F,
    const NT &nbrs, // lists of neighbors
    const size_t d, // dimension
    const size_t maxd, // max dimension
    const std::vector<size_t> &iter_idxs,
    std::vector<size_t> &spx_idxs,
    const T &t
) {
    // sorted simplices will end up here
    std::vector<size_t> spx_idxs2(spx_idxs.size() + 1);
    if (d == maxd) {
        // no recursion - we're adding maximal dimension cells
        for (auto k : iter_idxs) {
            // append k to spx_idxs, sort
            spx_idxs.push_back(k);
            sort_into(spx_idxs, spx_idxs2);

            // add to F
            F.add(t, spx_idxs2);

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
            F.add(t, spx_idxs2);

            // recurse
            intersect_sorted_lt(iter_idxs, nbrs[k], k, iter_idxs2);
            add_dimension_recursive_flag(F, nbrs, d+1, maxd, iter_idxs2, spx_idxs2, t);

            // pop k off spx_idxs
            spx_idxs.pop_back();
        }
    }
}



template <typename T, typename NT>
void add_dimension_recursive_flag_unsafe(
    Filtration<T, SimplicialComplex> &F,
    const NT &nbrs, // lists of neighbors
    const size_t d, // dimension
    const size_t maxd, // max dimension
    const std::vector<size_t> &iter_idxs,
    std::vector<size_t> &spx_idxs,
    const T t
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
            F.add(t, spx_idxs2);

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
            F.add(t, spx_idxs2);

            // recurse
            intersect_sorted_lt(iter_idxs, nbrs[k], k, iter_idxs2);
            add_dimension_recursive_flag_unsafe(F, nbrs, d+1, maxd, iter_idxs2, spx_idxs2, t);

            // pop k off spx_idxs
            spx_idxs.pop_back();
        }
    }
}

// this version keeps track of whether we should attempt to pair
template <typename T, typename NT>
void add_dimension_recursive_flag(
    Filtration<T, SimplicialComplex> &F,
    const NT &nbrs, // lists of neighbors
    const size_t d, // dimension
    const size_t maxd, // max dimension
    const std::vector<size_t> &iter_idxs,
    std::vector<size_t> &spx_idxs,
    const T &t,
    bool face_paired, // was the face paired?
    const cell_ind &fi // face location
) {
    // sorted simplices will end up here
    std::vector<size_t> spx_idxs2(spx_idxs.size() + 1);
    if (d == maxd) {
        // no recursion - we're adding maximal dimension cells
        for (auto k : iter_idxs) {
            if (!face_paired) {
                // append k to spx_idxs, sort
                spx_idxs.push_back(k);
                sort_into(spx_idxs, spx_idxs2);

                // add to F
                cell_ind si = F._add_unsafe_reserve(t, spx_idxs2);

                // pop k off spx_idxs
                spx_idxs.pop_back();

                F._set_pair_unsafe(fi.dim, fi.ind, si.ind);
                face_paired = true;
            } else {
                // face is already paired
                if (!has_intersect_sorted_lt(iter_idxs, nbrs[k], k)) {
                    // homology may be killed by simplex
                    // simplex should be added
                    spx_idxs.push_back(k);
                    sort_into(spx_idxs, spx_idxs2);

                    F._add_unsafe_reserve(t, spx_idxs2);

                    // pop k off spx_idxs
                    spx_idxs.pop_back();
                }
                // else simplex will be paired with simplex one dimension up
                // simplex should not be added since we're not at maximal dimension
                // this will not affect homology
            }
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
            cell_ind si = F._add_unsafe_reserve(t, spx_idxs2);

            // recurse
            intersect_sorted_lt(iter_idxs, nbrs[k], k, iter_idxs2);
            if (!iter_idxs2.empty()) {
                // there will be some recursion
                if (!face_paired) {
                    // we're going to pair with face, so further simplces can't pair with us
                    add_dimension_recursive_flag(F, nbrs, d+1, maxd, iter_idxs2, spx_idxs2, t, true, si);
                } else {
                    // we're not pairing with face, so coface can pair
                    add_dimension_recursive_flag(F, nbrs, d+1, maxd, iter_idxs2, spx_idxs2, t, false, si);
                }
            } // else no cofaces to add

            if (!face_paired) {
                // pair with face first chance
                F._set_pair_unsafe(fi.dim, fi.ind, si.ind);
                face_paired = true;
            }

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
Filtration<T, SimplicialComplex> FlagFiltration(
    std::vector<filtered_edge<T>> &edges,
    const size_t n, // number of 0-cells
    const size_t maxdim,
    const T t0
) {

    std::sort(edges.begin(), edges.end());

    // check that dimensions agree
    size_t m = edges.size();

    // X = SimplicialComplex(maxdim);
    // F = Filtration<T, SimplicialComplex>(X);
    // reset simplicial complex
    Filtration<T, SimplicialComplex> F;

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
        size_t i = edges[k].s;
        size_t j = edges[k].t;
        T t = edges[k].r;
        spx_idxs[0] = i;
        spx_idxs[1] = j;
        auto ret = F.add(t, spx_idxs);
        // std::cout << ret.second << std::endl;
        intersect_sorted(nbrs[i], nbrs[j], iter_idxs);

        if (!iter_idxs.empty()) {
            add_dimension_recursive_flag(F, nbrs, 2, maxdim, iter_idxs, spx_idxs, t);
        }

        // TODO: use std::set for neighbors - insertion is log(n)
        // nbrs[i].emplace(j);
        // nbrs[j].emplace(i);
        // TODO: insertion sort
        nbrs[i].emplace_back(j);
        std::sort(nbrs[i].begin(), nbrs[i].end());
        nbrs[j].emplace_back(i);
        std::sort(nbrs[j].begin(), nbrs[j].end());
    }

    return F;
}


// Flag complex using list of edges
// (edges[2*k], edges[2*k+1]) = (i, j) is an edge
// t - vector of filtration times
// t0 - time for 0-simplices
// n - number of vertices
// maxdim - maximum dimension of simplices
template <typename T>
std::tuple<SimplicialComplex, Filtration<T, SimplicialComplex>> CompleteFlagFiltration(
    const std::vector<size_t> &edges,
    const std::vector<T> &t,
    const size_t n, // number of 0-cells
    const size_t maxdim,
    const T t0
    // SimplicialComplex& X,
    // Filtration<T, SimplicialComplex>& F
) {

    // check that dimensions agree
    size_t m = t.size();
    assert (edges.size() == 2 * m);
    assert (m == (n * (n-1)) / 2);

    // precompute dimensions
    std::vector<size_t> dims = {n, m};
    for (size_t d = 2; d < maxdim+1; d++) {
        // recursive computation of binomial coefficient
        dims.push_back(
            dims[d-1] * (n - d) / (d + 1)
        );
        //std::cout << dims[d] << std::endl;
    }

    // X = SimplicialComplex(maxdim);
    // F = Filtration<T, SimplicialComplex>(X);
    SimplicialComplex X(dims);
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
        F.add(t[k], spx_idxs);

        intersect_sorted(nbrs[i], nbrs[j], iter_idxs);

        // nbrs[i].emplace(j);
        // nbrs[j].emplace(i);
        nbrs[i].emplace_back(j);
        std::sort(nbrs[i].begin(), nbrs[i].end());
        nbrs[j].emplace_back(i);
        std::sort(nbrs[j].begin(), nbrs[j].end());

        add_dimension_recursive_flag_unsafe(F, nbrs, 2, maxdim, iter_idxs, spx_idxs, t[k]);
    }

    return std::make_tuple(X, F);
}
