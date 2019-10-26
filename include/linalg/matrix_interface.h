
#pragma once

// parent of all matrices of all implementations
struct MAT{};

template<typename Impl> struct A : MAT {}; // arbitrary matrix
template<typename Impl> struct T : A<Impl> {}; // arbitrary triangular matrix
template<typename Impl> struct L : T<Impl> {}; // lower triangular
template<typename Impl> struct U : T<Impl> {}; // upper triangular
template<typename Impl> struct E : A<Impl> {}; // arbitrary echelon matrix
template<typename Impl> struct D : E<Impl>, L<Impl>, U<Impl> {}; // diagonal
template<typename Impl> struct EL : E<Impl>, L<Impl> {}; // echelon L matrix
template<typename Impl> struct EU : E<Impl>, U<Impl> {}; // echelon U matrix
template<typename Impl> struct P : E<Impl> {}; // permutation matrix


MAT matmul(MAT,MAT);

std::tuple<MAT,MAT,MAT,MAT> LEUP_Fact(MAT);

MAT commute_el(MAT,MAT);

MAT l_solve(MAT,MAT);
