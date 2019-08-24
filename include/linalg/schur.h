#pragma once
// schur complement

// schur complement over a matrix type
// S = D - C A^{-1} B
// assume column-dominant - e.g. ColumnMatrix or dense column major
// assume A is upper triangular
// TODO: dispatch based on triangle type of A
template <typename MatT>
inline MatT ut_schur(MatT &A, MatT &B, MatT &C, MatT &D) {
    return D - C * ut_solve(A, B);
}
