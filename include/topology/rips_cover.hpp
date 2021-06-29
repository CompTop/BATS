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
                    eset.emplace(filtered_edge(*iit, *jit, dij));
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

} // namespace bats
