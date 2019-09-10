#pragma once
// schur complement
#include <iostream>


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
