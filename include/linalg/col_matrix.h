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

  inline size_t width() {
    return col.size();
  }

  TC & operator[](size_t index)  {
    return col[index];
  }

  // permutations: permute, permute_rows, permute_cols

  // permute columns in-place
  // TODO: evaluate if this is best method.
  void permute_cols(const std::vector<size_t> &colperm) {
    size_t ncols = col.size();
    // record swapped columns
    bool visited[ncols] = {0};
    for (size_t j = 0; j < ncols; j++)   {
       size_t next_j = j;
       while(!visited[next_j] && !visited[colperm[next_j]]) {
         std::swap(col[next_j], col[colperm[next_j]]);
         visited[next_j] = true;
         next_j = colperm[next_j];
       }
    }
  }

  // permute rows in-place
  void permute_rows(const std::vector<size_t> &rowperm) {
    // TODO: this is trivially parallelizable
    for (size_t i = 0; i < col.size(); i++) {
      col[i].permute(rowperm);
    }
  }

  // permute both rows and columns
  void permute(const std::vector<size_t> &rowperm,
               const std::vector<size_t> &colperm) {
    permute_cols(colperm);
    permute_rows(rowperm);
  }


  // addition, substraction, scalar multiplication

  // gemm

  // triangular solve

  // schur complement friend

  void print() {
    std::cout << "transpose: " << std::endl;
    for (size_t i = 0; i < col.size(); i++) {
      std::cout << i << " : ";
      col[i].print_row();
    }
  }

  // dense matrix
};

template <class TC>
ColumnMatrix<TC> identity(size_t n) {
  std::vector<TC> col(n);
  for (size_t j = 0; j < n; j++) {
    col[j] = TC(j);
  }
  return ColumnMatrix<TC>(col);
}
