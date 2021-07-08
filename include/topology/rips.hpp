#pragma once
/*
generate Rips complex
*/

#include <vector>
#include <tuple>
#include "data.hpp"
#include "metric.hpp"
#include <complex/simplicial_complex.hpp>
#include "flag.hpp"
#include <filtration/flag.hpp>


namespace bats {

// construct edges for rips complex
// template over data type and metric
template <typename T, typename M>
std::vector<size_t> rips_edges(const DataSet<T> &X, const M &dist, const T rmax) {
    size_t n = X.size();
    std::vector<size_t> edges;
    size_t nedges = (n * (n-1)) / 2;
    edges.reserve(2*nedges);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < i; j++) {
            T dij = dist(X[i], X[j]);
            if (dij <= rmax) {
                edges.push_back(j);
                edges.push_back(i);
            }
        }
    }
    return edges;
}

// template over data type and metric
template <typename T, typename M>
std::vector<filtered_edge<T>> rips_filtration_edges(
    const DataSet<T> &X, const M &dist, const T rmax) {
    size_t n = X.size();
    std::vector<filtered_edge<T>> edges;
    size_t nedges = (n * (n-1)) / 2;
    edges.reserve(nedges);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < i; j++) {
            T dij = dist(X[i], X[j]);
            if (dij <= rmax) {
                edges.emplace_back(filtered_edge(i, j, dij));
            }
        }
    }
    // sort edges
    std::sort(edges.begin(), edges.end());
    edges.shrink_to_fit();
    return edges;
}

// work on pairwise distance matrix
template <typename T>
std::vector<filtered_edge<T>> rips_filtration_edges(
    const Matrix<T> &pdist,
    const T rmax)
{
    size_t n = pdist.ncol();
    std::vector<filtered_edge<T>> edges;
    size_t nedges = (n * (n-1)) / 2;
    edges.reserve(nedges);
    for (size_t j = 0; j < n; j++) {
        for (size_t i = 0; i < j; i++) {
            if (pdist(i,j) <= rmax) {
                edges.emplace_back(filtered_edge(i, j, pdist(i,j)));
            }
        }
    }
    return edges;
}

// work on pairwise distance matrix
template <typename T>
std::vector<size_t> rips_edges(const Matrix<T> &pdist, const T rmax) {
    size_t n = pdist.ncol();
    std::vector<size_t> edges;
    size_t nedges = (n * (n-1)) / 2;
    edges.reserve(nedges);
    for (size_t j = 0; j < n; j++) {
        for (size_t i = 0; i < j; i++) {
            if (pdist(i,j) <= rmax) {
                edges.push_back(i);
                edges.push_back(j);
            }
        }
    }
    return edges;
}

/*
Generate Rips complex from data
x - point cloud
dimension - dimension of point cloud
rmax - maximum radius for Rips complex
dmax - maximum simplex dimension
*/

// template over data type>
template <typename CpxT, typename T, typename M>
CpxT RipsComplex(
    const DataSet<T> &X,
    const M &dist,
    T rmax,
    size_t dmax
) {
    size_t n = X.size(); // number of points
    auto redges = rips_edges(X, dist, rmax);
    return FlagComplex<CpxT>(redges, n, dmax);
}

template <typename CpxT, typename T>
CpxT RipsComplex(
    const Matrix<T> &pdist,
    T rmax,
    size_t dmax
) {
    size_t n = pdist.ncol(); // number of points
    auto redges = rips_edges(pdist, rmax);
    return FlagComplex<CpxT>(redges, n, dmax);
}

// template over data type>
template <typename CpxT, typename T, typename M>
Filtration<T, CpxT> RipsFiltration(
    const DataSet<T> &X,
    const M &dist,
    T rmax,
    size_t dmax
) {
    size_t n = X.size(); // number of points
    auto edges = rips_filtration_edges(X, dist, rmax);
    return FlagFiltration<CpxT>(edges, n, dmax, T(0));
}

template <typename CpxT, typename T>
Filtration<T, CpxT> RipsFiltration(
    const Matrix<T> &pdist,
    T rmax,
    size_t dmax
) {
    size_t n = pdist.ncol(); // number of points
    auto redges = rips_filtration_edges(pdist, rmax);
    return FlagFiltration<CpxT>(redges, n, dmax, T(0));
}


// template over data type>
template <typename CpxT, typename T, typename M>
auto RipsFiltration_extension(
    const DataSet<T> &X,
    const M &dist,
    T rmax,
    size_t dmax
) {
    size_t n = X.size(); // number of points
    auto edges = rips_filtration_edges(X, dist, rmax);
    return FlagFiltration_extension<CpxT>(edges, n, dmax, T(0));
}

template <typename CpxT, typename T>
auto RipsFiltration_extension(
    const Matrix<T> &pdist,
    T rmax,
    size_t dmax
) {
    size_t n = pdist.ncol(); // number of points
    auto redges = rips_filtration_edges(pdist, rmax);
    return FlagFiltration_extension<CpxT>(redges, n, dmax, T(0));
}

} // namespace bats
