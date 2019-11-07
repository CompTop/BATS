#pragma once
/*
generate Rips complex
*/

#include <vector>
#include "metric.h"
#include <complex/simplicial_complex.h>

// construct edges for rips complex
template <typename T>
std::vector<size_t> rips_edges(const std::vector<T> &x, const size_t d, const T rmax) {
    std::vector<size_t> edges;
    size_t nedges = x.size() * (x.size() - d) / (d * d * 2);
    edges.reserve(2*nedges);
    for (size_t i = 0; i < (x.size() / d); i++) {
        for (size_t j = 0; j < i; j++) {
            T dij = euclidean<T>(x.cbegin() + d*i, x.cbegin() + d*j, d);
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
template <typename T>
SimplicialComplex RipsComplex(
    std::vector<T> x,
    size_t dimension,
    T rmax,
    size_t dmax
) {
    size_t n = x.size() / dimension; // number of points
    auto redges = rips_edges(x, dimension, rmax);
    return FlagComplex(redges, n, dmax);
}
