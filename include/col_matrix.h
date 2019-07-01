#pragma once

#include <vector>
#include <cstddef>
#include "abstract_matrix.h"

// standard list of list implementation for sparse matrix

// template over column type
template <class TC>
class ColumnMatrix : public AbstractMatrix
{
private:
  std::vector<TC> col;
public:

  ColumnMatrix(std::vector<TC> col) : col(col) {}

  // addition, substraction, scalar multiplication

  // gemm

  // triangular solve

  // schur complement friend

  // dense matrix
};
