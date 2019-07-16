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
  //
  // TC& constcol const (size_t index) {
  //   return col[index];
  // }

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


// return y = A*x
template <class TC>
TC gemv(ColumnMatrix<TC> &A, const TC &x) {
  TC y;  // zero initializer
  // loop over nonzero indices of v
  for (auto p = x.nzbegin(); p < x.nzend(); ++p) {
    y.axpy(p->second, A[p->first]); // y <- x[j]*A[j]
  }
  return y;
}

// return y = U \ x
// solves x = U * y
// Assumes U is upper triangular, with unit diagonal
// pseudo-code
// for i = n:-1:1
//    y[i] = y[i] - U[i,j]x[j]
template <class TC>
TC ut_solve(ColumnMatrix<TC> &U, const TC &x) {
  //std::cout << "entering solve" << std::endl;
  TC y(x);
  auto yi = y.nzend() - 1;
  //std::cout << yi->first << " : " << yi->second << std::endl;
  while (yi >= y.nzbegin() && yi < y.nzend() && yi->first > 0) {
    //std::cout << yi->first << " : " << yi->second << std::endl;
    auto i = yi->first;
    //auto yp = y.nzpair(yi);
    // y[yi->first] /= U[i][i]
    // but U is unit triangular
    y.axpy(-(yi->second), U[i], 1); //y.axpy(-yp.second, U[i][:i-1])
    // find next nonzero index
    yi = y.find_last_nz(i - 1);
    //y.print_row();
  }
  return y;
}

// solve U \ A
template <class TC>
ColumnMatrix<TC> ut_solve(ColumnMatrix<TC> &U, ColumnMatrix<TC> &A) {
  //std::cout << "entering solve" << std::endl;
  std::vector<TC> col;
  for (size_t j = 0; j < A.width(); j++) {
    std::cout << j << std::endl;
    A[j].print_row();
    std::cout << "\nU" << std::endl;
    U.print();
    auto Uinvj = ut_solve(U, A[j]);
    std::cout << "sol:" << std::endl;
    Uinvj.print_row();
    col.push_back(Uinvj);
  }
  return ColumnMatrix<TC>(col);
}
