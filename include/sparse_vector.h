#pragma once

#include <iostream>

// template over type of values
template <typename TI, typename TV>
class SparseVector
{
private:
  // non-zero indices
  std::vector<TI> ind;
  // non-zero values
  std::vector<TV> val;
public:

  SparseVector(std::vector<TI> ind, std::vector<TV> val) : ind(ind), val(val) {}
  // get index and set index

  // sort in-place
  void sort() {

  }

  // axpy - in place
  // scal - in place

  // add, subtract, multiply by scalar
  void print() {
    for (size_t i = 0; i < ind.size(); i++) {
      std::cout << ind[i] << " : " << val[i] << std::endl;
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
