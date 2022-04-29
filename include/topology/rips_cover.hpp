#pragma once
/*
generate Rips Cover complex
*/

#include <vector>
#include <tuple>
#include <set>
#include "data.hpp"
#include "metric.hpp"
#include "cover.hpp"
#include <complex/simplicial_complex.hpp>
#include <filtration/flag.hpp>

namespace bats {

// template over data type and metric
template <typename T, typename M>
std::vector<filtered_edge<T>> rips_filtration_edges(
    const DataSet<T> &X,
    const bats::Cover &cover,
    const M &dist,
    const T rmax) {
    std::set<filtered_edge<T>> eset;
    std::vector<filtered_edge<T>> edges;

    // loop over each set in cover
    for (auto& nbhd : cover) {
        // loop over each pair of edges in cover
        auto iit = nbhd.cbegin();
        while (iit != nbhd.cend()) {
            auto jit = iit;
            jit++;
            while (jit != nbhd.cend()) {
                T dij = dist(X[*iit], X[*jit]);
                if (dij < rmax) {
                    eset.emplace(filtered_edge<T>(*iit, *jit, dij));
                }
                jit++;
            }
            iit++;
        }
    }
    // put edges into vector
    edges.reserve(eset.size());
    for (auto e : eset) {
        edges.emplace_back(e);
    }
    return edges;
}

// template over data type>
template <typename T, typename M>
Filtration<T, SimplicialComplex> RipsFiltration(
    const DataSet<T> &X,
    const bats::Cover &cover,
    const M &dist,
    T rmax,
    size_t dmax
) {
    size_t n = X.size(); // number of points
    auto edges = rips_filtration_edges(X, cover, dist, rmax);
    return FlagFiltration<SimplicialComplex>(edges, n, dmax, T(0));
}

/**
Strict Cover Filtration
*/
template <typename T, typename CpxT, typename M>
Filtration<T, CpxT> StrictRipsCoverFiltration(
    const DataSet<T> &X,
    const bats::Cover &cover,
    const M &dist,
    T rmax,
    size_t dmax
) {
    size_t n = X.size();
    Filtration<T, CpxT> F(n, dmax);
    std::vector<filtered_edge<T>> edges;

    for (auto& U : cover) {
        edges.clear();
        auto iit = U.cbegin();
        while (iit != U.cend()) {
            auto jit = iit;
            ++jit;
            while (jit != U.cend()) {
                T dij = dist(X[*iit], X[*jit]);
                if (dij < rmax) {
                    edges.emplace_back(filtered_edge<T>(*iit, *jit, dij));
                }
                ++jit;
            }
            ++iit;
        }
        auto FU = FlagFiltration<CpxT>(edges, n, dmax, T(0)); // filtration over this set
        F.union_add(FU);
    }
    return F;
}

} // namespace bats
