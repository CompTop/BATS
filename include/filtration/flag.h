#pragma once
#include <vector>
#include <filtraiton/filtration.h>
#include <complex/simplicial_complex.h>
// flag filtrations

// Flag complex using list of edges
// (edges[2*k], edges[2*k+1]) = (i, j) is an edge
// n - number of vertices
// maxdim - maximum dimension of simplices
Filtration<SimplicialComplex, float> FlagComplex(
  std::vector<size_t> edges,
  std::vector<float> vals;
  size_t n,
  size_t maxdim
) {
  Filtration<SimplicialComplex, float> X(maxdim);
  // sets 0-cells
  X.set_ncells0(n);

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
    X.add_unsafe(spx_idxs);

    intersect_sorted(nbrs[i], nbrs[j], iter_idxs);

    nbrs[i].emplace_back(j);
    std::sort(nbrs[i].begin(), nbrs[i].end());
    nbrs[j].emplace_back(i);
    std::sort(nbrs[j].begin(), nbrs[j].end());

    X.add_dimension_recursive_flag_unsafe(nbrs, 2, maxdim, iter_idxs, spx_idxs);
  }

  return X;
}
