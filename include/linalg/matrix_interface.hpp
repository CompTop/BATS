
#pragma once

// parent of all matrices of all implementations
struct MAT{};

template<typename Impl> struct A : MAT {}; // arbitrary matrix
template<typename Impl> struct T : A<Impl> {}; // arbitrary triangular matrix
template<typename Impl> struct L : T<Impl> {}; // lower triangular
template<typename Impl> struct U : T<Impl> {}; // upper triangular
template<typename Impl> struct E : A<Impl> {}; // arbitrary echelon matrix
template<typename Impl> struct D : E<Impl>, L<Impl>, U<Impl> {}; // diagonal
template<typename Impl> struct EL : E<Impl>, L<Impl> {}; // echelon L pivot matrix
template<typename Impl> struct EU : E<Impl>, U<Impl> {}; // echelon U pivot matrix
template<typename Impl> struct ELH : E<Impl>, L<Impl> {}; // echelon L hat pivot matrix
template<typename Impl> struct EUH : E<Impl>, U<Impl> {}; // echelon U hat pivot matrix
template<typename Impl> struct P : E<Impl> {}; // permutation matrix


MAT matmul(MAT,MAT);

std::tuple<MAT,MAT,MAT,MAT> LEUP_Fact(MAT);
std::tuple<MAT,MAT,MAT,MAT> PLEU_Fact(MAT);
std::tuple<MAT,MAT,MAT,MAT> UELP_Fact(MAT);
std::tuple<MAT,MAT,MAT,MAT> PUEL_Fact(MAT);

MAT commute(MAT,MAT);

MAT apply_inverse(MAT,MAT);
