#pragma once

#include <iostream>
#include <algorithm>

// template over type of values
template <typename TI, typename TV>
class SparseVector
{
private:
  // store index-value pairs
  std::vector<std::pair<TI, TV>> indval;
  // // non-zero indices
  // std::vector<TI> ind;
  // // non-zero values
  // std::vector<TV> val;
public:

  SparseVector(std::vector<TI> ind, std::vector<TV> val) {
    size_t nz = ind.size();
    indval = std::vector<std::pair<TI, TV>>(nz);
    for (size_t i = 0; i < nz; i++) {
      indval[i] = std::pair<TI, TV>(ind[i], val[i]);
    }
  }
  // get index and set index

  // nnz
  inline size_t nnz() {
    return indval.size();
  }

  // sort in-place
  void sort() {
    // sort in-place
    std::sort(indval.begin(), indval.end());
  }

  // permute in-place
  void permute(const std::vector<size_t>  &perm) {
    for (size_t i = 0; i < nnz(); i++) {
      indval[i].first = perm[indval[i].first];
    }
    sort();
  }

  // axpy - in place
  // scal - in place

  // add, subtract, multiply by scalar

  void print() {
    for (size_t i = 0; i < nnz(); i++) {
      std::cout << indval[i].first << " : " << indval[i].second << std::endl;
    }
  }
};

// TODO: sparse F2 vector implementation
template <typename TI>
class SparseF2Vector
{
private:
  // non-zero indices
  std::vector<TI> ind;
public:


  // get index and set index
};
