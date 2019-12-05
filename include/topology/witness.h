#pragma once
/*
Witness complex
*/
#include <vector>
#include "data.h"
#include "neighborhood.h"
#include <util/sorted.h>
#include <complex/simplicial_complex.h>


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
    auto nbrs = neighborhoods(L, X, dist, 2);

    // TODO: finish

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
// template over data type, metric
template <typename T, typename M>
SimplicialComplex WitnessComplex(
    const DataSet<T> &X,
    const DataSet<T> &L,
    const M &dist,
    const size_t dmax
) {


    size_t k = dmax + 1; // number of neighbors to find
    // find the k closest neighbors in
    auto nbrs = neighborhoods(X, L, dist, k);

    // initialize complex
    SimplicialComplex W(dmax);

    // add one vertex to L for every landmark
    for (size_t i = 0; i < L.size(); i++) {
        W.add({i});
    }

    // container for holding simplex
    std::vector<size_t> spx;
    // now we loop over each point in x, and try adding the simplex it witnesses
    for (size_t dim = 1; dim < k; dim++) {
        spx.resize(dim+1);
        // add simplex that X[i] witnesses in dimension dim
        for (size_t i = 0; i < X.size(); i++) {
            for (size_t j = 0; j < dim+1; j++) {
                spx[j] = nbrs[i][j];
            }
            // the add is safe - only succeeds if boundary is present
            W.add(spx);
        }
    }

    return W;
}
