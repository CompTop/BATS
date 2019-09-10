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
    size_t m = 0; // number of rows
    size_t n = 0; // number of columns
    std::vector<TC> col;
public:

    // default constructor
    ColumnMatrix() {}

    // construct empty matrix.  same as zero matrix
    ColumnMatrix(size_t m, size_t n) : m(m), n(n) {
        col.resize(n, TC());
    }

    ColumnMatrix(std::vector<TC> col) : col(col) {
        n = col.size();
    }

    ColumnMatrix(size_t m, size_t n, std::vector<TC> col) : m(m), n(n), col(col) {}

    // constructor from CSCMatrix over the integers
    ColumnMatrix(const CSCMatrix<int, size_t> &A) : m(A.nrow()), n(A.ncol()) {
        col.reserve(n);
        auto colptr = A.get_colptr();
        auto rowind = A.get_rowind();
        auto val = A.get_rowval();
        auto rptr = rowind.cbegin();
        auto vptr = val.cbegin();
        for (size_t j = 0; j < n; j++) {
            // use iterator constructor
            col.emplace_back(TC(
                rptr + colptr[j],
                vptr + colptr[j],
                colptr[j+1] - colptr[j]
            ));
        }
    }

    inline size_t nrow() { return m; }
    inline size_t ncol() { return n; }

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
        std::vector<bool> visited(ncols, false);
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

    // gemv
    TC gemv(const TC &x) {
        TC y;  // zero initializer
        // loop over nonzero indices of x
        for (size_t i = 0; i < x.nnz(); i++) {
            y.axpy(x.nzval(i), col[x.nzind(i)]); // y <- x[j]*A[j]
        }
        return y;
    }

    // gemm C = self * B
    ColumnMatrix operator*(const ColumnMatrix &B) {
        ColumnMatrix C(m, B.n);
        for (size_t j = 0; j < B.n; j++) {
            C.col[j] = gemv(B.col[j]);
        }
        return C;
    }


    ColumnMatrix operator+(const ColumnMatrix &B) {
        ColumnMatrix C(m, n);
        for (size_t j = 0; j < n; j++) {
            C.col[j] = col[j] + B.col[j];
        }
        return C;
    }

    ColumnMatrix operator-(const ColumnMatrix &B) {
        ColumnMatrix C(m, n);
        for (size_t j = 0; j < n; j++) {
            C.col[j] = col[j] - B.col[j];
        }
        return C;
    }

    // triangular solve

    // schur complement friend

    void print() {
        std::cout << m << " x " << n << " matrix. transpose: " << std::endl;
        for (size_t i = 0; i < col.size(); i++) {
            std::cout << i << " : ";
            col[i].print_row();
        }
    }

    // expose as static member function?
    static ColumnMatrix identity(size_t n) {
        std::vector<TC> col(n);
        for (size_t j = 0; j < n; j++) {
            col[j] = TC(j);
        }
        return ColumnMatrix(n, n, col);
    }
    // dense matrix
};



// matvec

// matmul

// return y = A*x
template <class TC>
TC gemv(ColumnMatrix<TC> &A, const TC &x) {
    return A.gemv(x);
}

// return y = U \ x
// solves x = U * y
// Assumes U is upper triangular, with unit diagonal
// pseudo-code
// for i = n:-1:1
//    y[i] = y[i] - U[i,j]x[j]
// this is not generic over TC!
template <class TC>
TC ut_solve(ColumnMatrix<TC> &U, const TC &x) {
    //std::cout << "entering solve" << std::endl;
    TC y(x);
    auto yi = y.nzend() - 1;
    //std::cout << yi->first << " : " << yi->second << std::endl;
    while (yi >= y.nzbegin() && yi < y.nzend()) {
        size_t i = yi - y.nzbegin();
        size_t iind = y.nzind(i);
        if (iind == 0) {
            break;
        }
        y.axpy(-y.nzval(i), U[iind], 1); //y.axpy(-yp.second, U[i][:i-1])
        // find next nonzero index
        yi = y.find_last_nz(iind - 1);
    }
    return y;
}

// solve U \ A
template <class TC>
ColumnMatrix<TC> ut_solve(ColumnMatrix<TC> &U, ColumnMatrix<TC> &A) {
    //std::cout << "entering solve" << std::endl;
    size_t m = A.nrow();
    size_t n = A.ncol();
    std::vector<TC> col;
    for (size_t j = 0; j < A.width(); j++) {

        auto Uinvj = ut_solve(U, A[j]);

        col.push_back(Uinvj);
    }
    return ColumnMatrix<TC>(col);
}
