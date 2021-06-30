#pragma once
/*
Witness complex
*/
#include <vector>
#include <algorithm>
#include <limits>
#include "data.hpp"
#include "neighborhood.hpp"
#include "rips.hpp"
#include "nerve.hpp"
#include <util/sorted.hpp>
#include <complex/simplicial_complex.hpp>

namespace bats {

// find next smallest entry of columns of pdist[j]
template <typename T>
std::vector<T>& increment_m(
    const Matrix<T> &pdist,
    std::vector<T> &m
) {
    size_t N = m.size();
    size_t n = pdist.nrow();
    for (size_t j = 0; j < N; j++) {
        T mj = std::numeric_limits<T>::max();
        for (size_t i = 0; i < n; i++) {
            T dij = pdist(i,j);
            mj = (dij > mj) ? mj : (dij > m[j] ? dij : mj);
        }
        // update m[j] to be next smallest element
        m[j] = mj;
    }
    return m;
}

// get nu-th smallest entry in each column
template <typename T>
std::vector<T> get_m(
    const Matrix<T> &pdist,
    size_t nu
) {
    size_t N = pdist.ncol();
    std::vector<T> m(N, T(0));

    while (nu > 0) {
        // update m with next largest entry
        increment_m(pdist, m);
        nu--;
    }

    return m;
}

// get dowker edge parameter for landmarks i and j
// assume that pdist is |L| x |X|
template <typename T>
T dowker_edge_param(
    const Matrix<T> &pdist,
    const size_t i,
    const size_t j
) {
    T Rij = std::numeric_limits<T>::max();
    for (size_t k = 0; k < pdist.ncol(); k++) {
        Rij = std::min(Rij, std::max(pdist(i,k), pdist(j,k)));
    }
    return Rij;
}

// get witness edge parameters
// for witness filtration
template <typename T>
Matrix<T> dowker_edge_param(
    const Matrix<T> &pdist
) {

    auto nL = pdist.nrow();

    // set witness parameters using minmax
    Matrix<T> R(nL, nL);
    for (size_t j = 0; j < nL; j++) {
        for (size_t i = 0; i < j; i++) {
            R(i,j) = R(j,i); // symmetry
        }
        for (size_t i = j+1; i < nL; i++) {
            R(i,j) = dowker_edge_param(pdist, i, j);
        }
    }
    return R;
}

// get witness edge parameters
// for witness filtration
template <typename T, typename M>
Matrix<T> witness_edge_param(
    const DataSet<T> &X,
    const DataSet<T> &L,
    const M &dist,
    const size_t nu
) {

    // pairwise distances
    auto pdist = dist(L, X);

    // get m
    auto m = get_m(pdist, nu);
    auto nX = pdist.ncol();

    // loop over pdist, subtract m[j] from each column, take max with 0
    for (size_t j = 0; j < nX; j++) {
        for (size_t i = 0; i < nX; i++) {
            pdist(i,j) -= m[j];
            pdist(i,j) = (pdist(i,j) > T(0)) ? pdist(i,j) : T(0);
        }
    }

    return dowker_edge_param(pdist);
}


// witness edges
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

// lazy witness complex with parameter rmax
template <typename T, typename M>
SimplicialComplex WitnessComplex(
    const DataSet<T> &X,
    const DataSet<T> &L,
    const M &dist,
    const size_t nu,
    const T rmax,
    const size_t dmax
) {
    auto pdist = witness_edge_param(X, L, dist, nu);
    return RipsComplex(pdist, rmax, dmax);
}

template <typename T, typename M>
inline auto witness_neighborhoods(
    const DataSet<T> &X,
    const DataSet<T> &L,
    const M &dist,
    const size_t nu,
    const T rmax
) {
    auto pdist = witness_edge_param(X, L, dist, nu);
    // get cover
    return neighborhoods(pdist, rmax);
}


// strict witness complex with parameter rmax
// gets
template <typename T, typename M>
SimplicialComplex StrictWitnessComplex(
    const DataSet<T> &X,
    const DataSet<T> &L,
    const M &dist,
    const size_t nu,
    const T rmax,
    const size_t dmax
) {
    // get cover
    auto cover = witness_neighborhoods(X, L, dist, nu, rmax);
    return Nerve(cover, dmax);
}


template <typename T>
std::vector<filtered_edge<T>> flag_filtration_edges(
    const Matrix<T> &pdist,
    const T rmax) {
    size_t n = pdist.ncol();
    std::vector<filtered_edge<T>> edges;
    size_t nedges = (n * (n-1)) / 2;
    edges.reserve(nedges);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < i; j++) {
            T dij = pdist(i,j);
            if (dij < rmax) {
                edges.emplace_back(filtered_edge(i, j, dij));
            }
        }
    }
    // sort edges
    std::sort(edges.begin(), edges.end());
    edges.shrink_to_fit();
    return edges;
}


// template over data type>
template <typename T, typename M>
Filtration<T, SimplicialComplex> WitnessFiltration(
    const DataSet<T> &L,
    const DataSet<T> &X,
    const M &dist,
    T rmax,
    size_t dmax
) {
    size_t n = L.size(); // number of points in landmark set
    auto R = witness_edge_param(X, L, dist, 0); // nu=0
    auto edges = flag_filtration_edges(R, rmax);
    return FlagFiltration<SimplicialComplex>(edges, n, dmax, T(0));
}

// template over data type>
// pdist is nL x nX size
template <typename T>
Filtration<T, SimplicialComplex> DowkerFiltration(
    const Matrix<T> &pdist,
    T rmax,
    size_t dmax
) {
    size_t n = pdist.nrow(); // number of points in landmark set
    auto R = dowker_edge_param(pdist); // nu=0
    auto edges = flag_filtration_edges(R, rmax);
    return FlagFiltration<SimplicialComplex>(edges, n, dmax, T(0));
}


template <typename T, typename M>
Filtration<T, SimplicialComplex> DowkerFiltration(
    const DataSet<T> &L,
    const DataSet<T> &X,
    const M &dist,
    T rmax,
    size_t dmax
) {
    Matrix<T> pdist = dist(L, X);
    return DowkerFiltration(pdist, rmax, dmax);
}

} // namespace bats
