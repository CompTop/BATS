#pragma once

/*
Create Dowker-cover complexes
flag filtration version
*/

#include <vector>
#include <tuple>
#include <set>
#include "data.hpp"
#include "metric.hpp"
#include "cover.hpp"
#include "dowker.hpp"
#include <complex/simplicial_complex.hpp>
#include <filtration/flag.hpp>


namespace bats {

// template over data type and metric
// cover is over the landmark set
template <typename T>
std::vector<filtered_edge<T>> dowker_filtration_edges(
    const Matrix<T> &pdist,
    const bats::Cover &cover,
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
                T dij = dowker_edge_param(pdist, *iit, *jit);
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
template <typename T>
Filtration<T, SimplicialComplex> DowkerFiltration(
    const Matrix<T> &pdist,
    const bats::Cover &cover,
    T rmax,
    size_t dmax
) {
    auto edges = dowker_filtration_edges(pdist, cover, rmax);
    return FlagFiltration<SimplicialComplex>(edges, pdist.nrow(), dmax, T(0));
}

template <typename T, typename M>
Filtration<T, SimplicialComplex> DowkerFiltration(
    const DataSet<T> &L,
    const DataSet<T> &X,
    const M &dist,
    const bats::Cover &cover,
    T rmax,
    size_t dmax
) {
    Matrix<T> pdist = dist(L, X);
    return DowkerFiltration(pdist, cover, rmax, dmax);
}

} // namespace bats
