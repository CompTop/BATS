#pragma once

#include <vector>
#include <algorithm>
#include <numeric>
#include <complex/abstract_complex.h>
#include <complex/simplicial_complex.h>
#include <linalg/sparse_vector.h>
#include <linalg/col_matrix.h>
#include <util/sorted.h>

// template over
//  TC - complex type
//  TF - filtration type
template <class TC, typename TF>
class Filtration
{
private:
  TC cpx;
  std::vector<std::vector<TF>> val; // filtration value for each cell
  std::vector<std::vector<size_t>> filt_perm; // filtration perumation for each dimension
  std::vector<std::vector<size_t>> inv_perm; // inverse permutation for each dimension
public:

  // empty initializer
  Filtration() : cpx(TC()) {}

  // initialize with maxdim
  Filtration(size_t maxdim) : cpx(TC(maxdim)) {
    val = std::vector<std::vector<TF>>(maxdim+1);
  }

  // initialize with only complex
  Filtration(TC cpx) : cpx(cpx) {}

  Filtration(TC cpx,
    std::vector<std::vector<TF>> val) : cpx(cpx), val(val) {}

  bool add_unsafe(std::vector<size_t>, const TF t);

  bool add(std::vector<size_t>, const TF t);

  void add_dimension_recursive_flag_unsafe(std::vector<std::vector<size_t>> &nbrs,
    size_t d, size_t maxd,
    std::vector<size_t> &iter_idxs,
    std::vector<size_t> &spx_idxs,
    const TF t
  );

  // assume filtration value has been set for each cell.  Find sort permutaion in each dimension
  void sort() {
    filt_perm.resize(cpx.maxdim() + 1);
    inv_perm.resize(cpx.maxdim() + 1);
    // loop over dimension
    for (size_t dim = 0; dim < cpx.maxdim() + 1; dim++) {
      // std::cout << "sorting dim = " << dim << std::endl;
      size_t ncells_dim = cpx.ncells(dim);
      // std::cout << "ncells dim = " << ncells_dim << std::endl;
      filt_perm[dim] = std::vector<size_t>(ncells_dim);
      // std::cout << "created vector" << std::endl;
      std::iota(filt_perm[dim].begin(), filt_perm[dim].end(), 0);
      // std::cout << "filled iota" << std::endl;
      std::sort(filt_perm[dim].begin(), filt_perm[dim].end(),
        [&](const int& a, const int& b) {
          return (val[dim][a] < val[dim][b]);
        }
      );
      // std::cout << "sorted" << std::endl;
      // now filt_perm[dim] holds the sort permutation on cells for the filtration
      inv_perm[dim] = std::vector<size_t>(ncells_dim);
      for (size_t i = 0; i < ncells_dim; i++) {
        inv_perm[dim][filt_perm[dim][i]] = i;
      }
      // std::cout << "inv sorted" << std::endl;
      // now inv_perm[dim] holds the inverse permutation
    }
    return;
  }

  void print(size_t dim) {
    for (auto i : filt_perm[dim]) {
      std::cout << val[dim][i] << " | ";
      cpx.print_cell(dim, i);
    }
  }

  void print() {
    for (size_t dim = 0; dim < cpx.maxdim() + 1; dim ++) {
      std::cout << "dim " << dim << " : " << cpx.ncells(dim) << " cells" << std::endl;
      print(dim);
    }
  }

  // create boundary in dimension dim
  // template over column type
  template <class TVec>
  ColumnMatrix<TVec> boundary(size_t dim) {
    // first get column matrix of boundary
    ColumnMatrix<TVec> B = cpx.template boundary<TVec>(dim);
    // // TODO: permute rows and columns
    B.permute(filt_perm[dim-1], filt_perm[dim]);
    return B;
  }

};

// add cell c to complex at value t
// Special method fo Simplicial Complexes
// template <class TC, typename TF> class Filtration; // primary template
// template <typename TF>
template <>
bool Filtration<SimplicialComplex, float>::add(std::vector<size_t> c, const float t) {
  bool added = cpx.add(c);
  if (added) {
    size_t dim = c.size() - 1;
    while (val.size() < dim + 1) {
      val.emplace_back(std::vector<float>());
    }
    val[dim].emplace_back(t);
  }
  return added;
}

template <>
bool Filtration<SimplicialComplex, float>::add_unsafe(std::vector<size_t> c, const float t) {
  bool added = cpx.add_unsafe(c);
  if (added) {
    size_t dim = c.size() - 1;
    val[dim].emplace_back(t);
  }
  return added;
}

template<>
void Filtration<SimplicialComplex, float>::add_dimension_recursive_flag_unsafe(
  std::vector<std::vector<size_t>> &nbrs,
  size_t d, size_t maxd,
  std::vector<size_t> &iter_idxs,
  std::vector<size_t> &spx_idxs,
  const float t
) {
  cpx.add_dimension_recursive_flag_unsafe(nbrs, d, maxd, iter_idxs, spx_idxs);
  // add filtration time until val arrays have same number of elements as complex
  for (size_t dim = d; dim < maxd + 1; dim++) {
    while (val[dim].size() < cpx.ncells(dim)) {
      val[dim].emplace_back(t);
    }
  }
}

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
