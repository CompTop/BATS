#pragma once

#include "abstract_complex.h"
#include <linalg/sparse_vector.h>
#include <linalg/col_matrix.h>
#include <vector>
#include <cstddef>
#include <iostream>
#include <map>
#include <limits>
#include <algorithm>

// no index
#define NO_IND std::numeric_limits<size_t>::max();

/*
simplicial complex implementation
*/
class SimplicialComplex : public AbstractComplex
{
private:
  // hold number of 0-cells
  size_t ncells0;
  // spx_list[k][i] holds i'th simplex in dimension k+1
  std::vector<std::vector<std::vector<size_t>>> spx_list;
  // map to find simplices when forming boundary
  std::vector<std::map<std::vector<size_t>, size_t>> spx_to_idx;
  // TODO: if forming full clique complex
  // can give all simplices unique location for lookup

  // returns index of simplex
  size_t find_idx(std::vector<size_t> s) {
    size_t dim = s.size() - 1;
    if (dim == 0) {
      return s[0] < ncells0 ? s[0] : NO_IND;
    } else {
      if (spx_to_idx[dim-1].count(s)) {
        return spx_to_idx[dim-1][s];
      }
    }
    return NO_IND;
  }

public:

  // default constructor
  SimplicialComplex() : ncells0(0) {}

  // constructor that initializes to set dimension
  SimplicialComplex(size_t maxdim) : ncells0(0), spx_list(maxdim), spx_to_idx(maxdim) {}

  // adds a simplex to the complex without doing any checks
  bool add_unsafe(const std::vector<size_t> &s) {
    size_t dim = s.size() - 1;
    if (dim == 0) {
      ncells0 = std::max(s[0] + 1, ncells0);
    } else {
      // add simplex to appropriate dimension
      spx_list[dim-1].push_back(s);
      // set reverse map
      spx_to_idx[dim-1][s] = spx_list[dim-1].size() - 1;
    }
    return true;
  }

  // add simplex to complex with appropriate checks
  bool add(std::vector<size_t> s) {
    size_t dim = s.size() - 1;
    if (dim == 0){
      ncells0 = std::max(s[0] + 1, ncells0);
      return true;
    }

    // ensure simplex is sorted
    std::sort(s.begin(), s.end());

    // add dimensions if necessary
    while (spx_list.size() < dim) {
      spx_list.push_back(std::vector<std::vector<size_t>>());
      spx_to_idx.push_back(std::map<std::vector<size_t>, size_t>());
    }

    // check if simplex is already in complex
    if (spx_to_idx[dim-1].count(s) > 0) {
      // simplex is already in complex
      return false;
    }

    return add_unsafe(s);
  }

  // return 0-skeleton of cell
  std::vector<size_t> skeleton0(size_t dim, size_t i) {
    if (dim == 0) {
      return {i};
    }
    return spx_list[dim-1][i];
  }



  size_t ncells(size_t dim) {
    return dim == 0 ? ncells0 : spx_list[dim - 1].size();
  }

  size_t maxdim() {
    return spx_list.size();
  }

  void print_dims() {
    std::cout << "maxdim = " << maxdim() << std::endl;
    for (size_t dim = 0; dim < maxdim() + 1; dim ++) {
      std::cout << "dim " << dim << " : " << ncells(dim) << std::endl;
    }
  }

  void print_cell(size_t dim, size_t i) {
    if (dim == 0) {
      std::cout << i << std::endl;
      return;
    }
    for (size_t j = 0; j < dim; j++) {
      std::cout << spx_list[dim - 1][i][j] << ",";
    }
    std::cout << spx_list[dim - 1][i][dim] << std::endl;
  }

  void print(size_t dim) {
    if (dim == 0) {
      for (size_t i = 0; i < ncells0; i++) {
        print_cell(0, i);
      }
    } else {
      for (size_t i = 0; i < ncells(dim); i++) {
        print_cell(dim, i);
      }
    }
  }

  void print() {
    std::cout << "SimplicialComplex of dimension " << maxdim() << std::endl;
    for (size_t dim = 0; dim < maxdim() + 1; dim ++) {
      std::cout << "dim " << dim << " : " << ncells(dim) << " cells" << std::endl;
      print(dim);
    }
  }


  // get boundary of simplex i in dimension dim
  template <typename TV>
  SparseVector<size_t, TV> boundary(size_t dim, size_t i) {
    //assert(dim > 0);
    std::vector<size_t> bdr_ind;
    std::vector<int> bdr_val;
    std::vector<size_t> face;
    int coeff = -1;
    // loop over faces in lexicographical order
    for (size_t k = 0; k < dim+1; k++) {
      size_t k2 = dim-k; // index to skip
      face.clear();
      for (size_t j = 0; j < k2; j++) {
        face.push_back(spx_list[dim-1][i][j]);
      }
      for (size_t j = k2+1; j < dim+1; j++) {
        face.push_back(spx_list[dim-1][i][j]);
      }

      bdr_ind.push_back(find_idx(face));
      bdr_val.push_back(coeff);
      coeff = -coeff;
    }
    return SparseVector<size_t, TV>(bdr_ind, bdr_val);
  }

  // template over column type
  template <class TVec>
  ColumnMatrix<TVec> boundary(size_t dim) {
    //assert(dim > 0);
    std::vector<TVec> col;
    std::vector<size_t> bdr_ind;
    std::vector<int> bdr_val;
    std::vector<size_t> face;
    for (size_t i = 0; i < ncells(dim); i++) {
      bdr_ind.clear();
      bdr_val.clear();
      int coeff = -1;
      // loop over faces in lexicographical order
      for (size_t k = 0; k < dim+1; k++) {
        size_t k2 = dim-k; // index to skip
        face.clear();
        for (size_t j = 0; j < k2; j++) {
          face.push_back(spx_list[dim-1][i][j]);
        }
        for (size_t j = k2+1; j < dim+1; j++) {
          face.push_back(spx_list[dim-1][i][j]);
        }

        bdr_ind.push_back(find_idx(face));
        bdr_val.push_back(coeff);
        coeff = -coeff;
      }
      col.push_back(TVec(bdr_ind, bdr_val));
    }
    return ColumnMatrix<TVec>(col);
  }



};
