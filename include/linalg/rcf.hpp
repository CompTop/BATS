#pragma once
/*
Compute the rational canonical form of a matrix over a field

Naive algorithm:
c.f. Dummit & Foote section 12.2
*/
#include <vector>
#include "sparse_vector.hpp"
#include "col_matrix.hpp"
#include "polynomial.hpp"
#include "pid.hpp"

// generate basis for k[T] module
// from R factor of Smith factorization of A
template <typename T>
ColumnMatrix<SparseVector<T>> generating_basis(
    const ColumnMatrix<SparseVector<UnivariatePolynomial<T>>>& R,
    const ColumnMatrix<SparseVector<T>> &A
) {
    using VT = SparseVector<T>;
    auto m = R.nrow();
    std::vector<VT> col(m);
    for (size_t j = 0; j < m; j++) {
        auto p = R[j].nzbegin();
        while (p != R[j].nzend()) {
            col[j].axpy(T(1), p->val(A, VT(p->ind))); // Cj = R[i,j](A) * ei
            p++;
        }
    }
    return ColumnMatrix(m, m, col);
}
// Sanity check: if smith form S(j,j) = 1,
// then basis vector j should be 0.

// generate whole basis for RCF
// S is smith form of xI - A
// B is generating basis for k[A]
template <typename T>
ColumnMatrix<SparseVector<T>> RCF_basis(
    const ColumnMatrix<SparseVector<UnivariatePolynomial<T>>>& S,
    const ColumnMatrix<SparseVector<T>> &A,
    const ColumnMatrix<SparseVector<T>> &B
) {
    using VT = SparseVector<T>;

    auto m = S.nrow();
    std::vector<VT> col;
    col.reserve(m);
    // loop over columns of R
    for (size_t j = 0; j < m; j++) {
        if (S(j,j).degree() > 0) {
            auto v = B[j];
            col.emplace_back(v);
            for (size_t k = 0; k < S(j,j).degree(); k++) {
                v = A * v;
                col.emplace_back(v);
            }
        }

    }
    return ColumnMatrix(m, m, col);
}
