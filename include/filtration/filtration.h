#pragma once

#include <vector>
#include <algorithm>
#include <numeric>
#include <complex/abstract_complex.h>
#include <complex/simplicial_complex.h>
#include <linalg/sparse_vector.h>
#include <linalg/col_matrix.h>

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

  // initialize with only complex
  Filtration(TC cpx) : cpx(cpx) {}

  Filtration(TC cpx,
    std::vector<std::vector<TF>> val) : cpx(cpx), val(val) {}

  bool add(std::vector<size_t>, const TF t);

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
      val.push_back(std::vector<float>());
    }
    val[dim].push_back(t);
  }
  return added;
}

// construct clique complex on n vertices
