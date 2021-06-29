#pragma once
// schur complement
#include <iostream>
#include "col_matrix.hpp"

// schur complement over a matrix type
// S = D - C A^{-1} B
// assume column-dominant - e.g. ColumnMatrix or dense column major
// assume A is upper triangular
// TODO: dispatch based on triangle type of A
// template <typename MatT>
// inline MatT ut_schur(MatT &A, MatT &B, MatT &C, MatT &D) {
//     return D - C * ut_solve(A, B);
// }

// S = D - C * A \ B
// uses operations defined for CSCMatrix
template <typename T>
T schur(T& A, T& B, T& C, T& D) {
    T R, S;
    std::cout << "trilu" << std::endl;
    trilu(A, B, S);   // S <- A \ B
    std::cout << "gemm" << std::endl;
    gemm(C, S, R);    // R <- C * S = C * A \ B
    std::cout << "sum" << std::endl;
    sum(-1, R, D, S); // S <- D - R = D - C * A \ B
    return S;
}

// column matrix specialization
template <typename CT>
ColumnMatrix<CT> schur(
    ColumnMatrix<CT> &A,
    ColumnMatrix<CT> &B,
    ColumnMatrix<CT> &C,
    ColumnMatrix<CT> &D
) {
    using MT = ColumnMatrix<CT>;

    // A.print_size();
    // B.print_size();
    // C.print_size();
    // D.print_size();
    // std::cout << "solve" << std::endl;
    MT S = solve_U(A, B); // S <- A \ B
    MT R = C * S;         // R <- C * S = C * A \ B
    return D - R;         // D - R = D - C * A \ B
}

/*
Schur complement with respect to A[i0:i1, j0:j1]
A is modified to produce AS
      A                 AS
[A11 A12 A13]     [S11  0  S13]
[A21 A22 A23] ->  [ 0  A22  0 ]
[A31 A32 A33]     [S31  0  S33]
where A22 = A.block(i0,i1, j0,j1)
returns L1 U1, AS, L2 U2
so that A = L1 U1 AS L2 U2

// assume A22 is invertible
*/
template <typename TV>
auto schur(
    const ColumnMatrix<TV>& A,
    const size_t i0,
    const size_t i1,
    const size_t j0,
    const size_t j1
) {
    // TODO: check that A22 is square

    auto m = A.nrow();
    auto n = A.ncol();

    // obtain block sub-matrices
    auto A11 = A.block(0, i0, 0, j0);
    auto A12 = A.block(0, i0, j0, j1);
    auto A13 = A.block(0, i0, j1, n);

    auto A21 = A.block(i0, i1, 0, j0);
    auto A22 = A.block(i0, i1, j0,j1);
    auto A23 = A.block(i0, i1, j1, n);

    auto A31 = A.block(i1, m, 0, j0);
    auto A32 = A.block(i1, m, j0, j1);
    auto A33 = A.block(i1, m, j1, n);

    auto A22inv = inv(A22);


    // A32 * A22^{-1}
    auto L32 =  A32 * A22inv;
    // A12 * A22^{-1}
    auto U12 = A12 * A22inv;
    // A22^{-1} * A21
    auto L21 = A22inv * A21;
    // A22^{-1} * A23
    auto U23 = A22inv * A23;

    // S11 = A11 - A12 * A22^{-1} * A21
    auto S11 = A11 - A12 * L21;
    // S31 = A31 - A32 * A22^{-1} * A21
    auto S31 = A31 - A32 * L21;
    // S13 = A13 - A12 * A22^{-1} * A23
    auto S13 = A13 - A12 * U23;
    // S33 = A33 - A32 * A22^{-1} * A23
    auto S33 = A33 - A32 * U23;

    // now construct the matrices for return
    auto L1 = ColumnMatrix<TV>::identity(m);
    // set 32 block to L32
    L1.set_block(i1, m, j0, j1, L32);

    auto U1 = ColumnMatrix<TV>::identity(m);
    // set 12 block to U12
    U1.set_block(0, i0, j0, j1, U12);

    auto L2 = ColumnMatrix<TV>::identity(n);
    // set 21 block to L21
    L2.set_block(i0, i1, 0, j0, L21);

    auto U2 = ColumnMatrix<TV>::identity(n);
    // set 23 block to U23
    U2.set_block(i0, i1, j1, n, U23);

    ColumnMatrix<TV> AS(m,n);
    AS.set_block(0, i0, 0, j0, S11);
    AS.set_block(0, i0, j1, n, S13);
    AS.set_block(i0, i1, j0,j1, A22);
    AS.set_block(i1, m, 0, j0, S31);
    AS.set_blcok(i1, m, j1, n, S33);

    return std::tuple(L1, U1, AS, L2, U2);
}
