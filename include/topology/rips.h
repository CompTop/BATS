#pragma once
/*
generate Rips complex
*/

#include <vector>
#include "data.h"
#include "metric.h"
#include <complex/simplicial_complex.h>

// construct edges for rips complex
// template over data type and metric
template <typename T, typename M>
std::vector<size_t> rips_edges(const Matrix<T> &X, const M &dist, const T rmax) {
    size_t n = X.ncol();
    std::vector<size_t> edges;
    size_t nedges = (n * (n-1)) / 2;
    edges.reserve(2*nedges);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < i; j++) {
            T dij = dist(X[i], X[j]);
            if (dij < rmax) {
                edges.push_back(j);
                edges.push_back(i);
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
template <typename T, typename M>
SimplicialComplex RipsComplex(
    const DataSet<T> &X,
    const M &dist,
    T rmax,
    size_t dmax
) {
    size_t n = X.size(); // number of points
    auto redges = rips_edges(X.data, dist, rmax);
    return FlagComplex(redges, n, dmax);
}
