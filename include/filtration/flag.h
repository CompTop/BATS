#pragma once
#include <vector>
#include <filtration/filtration.h>
#include <complex/simplicial_complex.h>
// flag filtrations

// Flag complex using list of edges
// (edges[2*k], edges[2*k+1]) = (i, j) is an edge
// t - vector of filtration times
// t0 - time for 0-simplices
// n - number of vertices
// maxdim - maximum dimension of simplices
template <typename T>
Filtration<SimplicialComplex, T> FlagFiltration(std::vector<size_t> edges, std::vector<T> t, T t0, size_t n, size_t maxdim) {

  Filtration<SimplicialComplex, T> F = Filtration<SimplicialComplex, T>(maxdim);

  // sets 0-cells
  for (size_t k = 0; k < n; k++) {
    F.add_unsafe({k}, t0);
  }

  std::vector<std::vector<size_t>> nbrs(n);

  std::vector<size_t> spx_idxs(2);
  std::vector<size_t> iter_idxs;
  iter_idxs.reserve(n); // maximum size

  size_t m = edges.size() / 2;
  for (size_t k = 0; k < m; k++) {
    size_t i = edges[2*k];
    size_t j = edges[2*k + 1];
    spx_idxs[0] = i;
    spx_idxs[1] = j;
    F.add_unsafe(spx_idxs, t[k]);

    intersect_sorted(nbrs[i], nbrs[j], iter_idxs);

    nbrs[i].emplace_back(j);
    std::sort(nbrs[i].begin(), nbrs[i].end());
    nbrs[j].emplace_back(i);
    std::sort(nbrs[j].begin(), nbrs[j].end());

    F.add_dimension_recursive_flag_unsafe(nbrs, 2, maxdim, iter_idxs, spx_idxs, t[k]);
  }

  return F;
}
