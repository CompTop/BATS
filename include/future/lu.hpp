#pragma once

#include <utility>

namespace bats {
namespace future {

/*
function to perform computation in-place

L is a unit-lower triangular matrix.

yit and xit are random access iterators.

yit is over-written with residual from forward substitution

xit is over-written to solve
Ly = x

returns true if system is consistent.  Otherwise, returns false
*/
template <typename MT, typename I1, typename I2>
bool l_solve(const MT &L, I1& y, I2& x) {

	// forward substitution
	for (size_t j = 0; j < L.ncol(); j++) {
		auto coeff = y[j];
		if (coeff != 0) {
			x[j] = coeff;
			for (size_t i = j; i < L.nrow(); i++) {
				y[i] -= coeff * L(i,j);
			}
		} else {
			x[j] = 0;
		}
	}

	return true;
}

// find residual in-place
template <typename MT, typename I1>
bool l_residual(const MT &L, I1& y) {

	// forward substitution
	for (size_t j = 0; j < L.ncol(); j++) {
		auto coeff = y[j];
		if (coeff != 0) {
			for (size_t i = j; i < L.nrow(); i++) {
				y[i] -= coeff * L(i,j);
			}
		}
	}

	return true;
}
// find residual in-place
template <typename MT, typename I1>
bool l_residual(const MT &&L, I1&& y) {

	// forward substitution
	for (size_t j = 0; j < L.ncol(); j++) {
		auto coeff = y[j];
		if (coeff != 0) {
			for (size_t i = j; i < L.nrow(); i++) {
				y[i] -= coeff * L(i,j);
			}
		}
	}

	return true;
}

// solve P L x = y -> L x = P^T y
// over-write y with residual
template <typename MT1, typename MT2, typename I1, typename I2>
bool l_solve(const MT1 &Pt, const MT2 &L, I1& y, I2& x) {
	// Pt = P^T
	Pt.matvec_inplace(y);
	return l_solve(L, y, x);
}



/*
LU Factorization struct

A = PLUQ
*/

template <typename MT>
struct LUFact {
	using val_type = typename MT::value_type;

	MT P; // row pivoting
	MT L; // unit lower triangular
	MT U; // upper triangular
	MT Q; // column pivoting

	LUFact(MT &&P, MT &&L, MT &&U, MT &&Q) : P(P), L(L), U(U), Q(Q) {}

	MT prod() const {return (P * L) * (U * Q);}

	void print_info() const {
		std::cout << "[" << this << "] : " << \
        " LUFact" << std::endl;
	}
	void print() const {
		print_info();
		P.print();
		L.print();
		U.print();
		Q.print();
	}

};

// find pivot in A[k:, k:]
// we'll just search for a non-zero
template <typename MT>
std::pair<size_t, size_t> find_pivot_complete(const MT &A, size_t k) {

	// loop in row-major order
	for (size_t i = k; i < A.nrow(); i++) {
		for (size_t j = k; j < A.ncol(); j++) {
			if (A(i,j) != 0) {return std::make_pair(i,j);}
		}
	}
	return std::make_pair(k,k);
}

// return unit lower triangular matrix part of A
template <typename MT>
MT unit_lower(const MT &A) {
	auto [m, n] = A.dims();
	// initialize to lower triangular
	MT L = MT::identity(m);
	// loop in row major order
	for (size_t i = 1; i < m; i++) {
		for (size_t j = 0; j < i; j++ ) {
			L(i,j) = A(i,j);
		}
	}
	return L;
}

// return upper triangular matrix part of A
template <typename MT>
MT upper(const MT &A) {

	using T = typename MT::value_type;

	auto [m, n] = A.dims();

	MT U(m, n, T(0));
	for (size_t i = 0; i < m; i++) {
		for (size_t j = i; j < n; j++ ) {
			U(i,j) = A(i,j);
		}
	}
	return U;

}

// complete pivoting
// GvL alg. 3.4.3
template <typename MT>
LUFact<MT> LU(const MT &A0) {

	MT A = A0;
	auto [m, n] = A.dims();

	MT P = MT::identity(m);
	MT Q = MT::identity(n);
	// maintain invariant P^T A Q^T

	for (size_t k = 0; k < n-1; k++) {
		auto [mu, lam] = find_pivot_complete(A, k);
		// do pivoting
		A.swap_rows(k, mu);
		P.swap_columns(k, mu); // record row pivots in P
		A.swap_columns(k, lam);
		Q.swap_rows(k,lam); // record column pivots in Q

		// TODO: record pivots in P, Q
		if (A(k,k) != 0) {
			for (size_t i = k+1; i < m; i++) {
				A(i,k) /= A(k,k);
			}
			for (size_t i = k+1; i < m; i++) {
				for (size_t j = k+1; j < n; j++) {
					A(i, j) -= A(i,k) * A(k,j);
				}
			}
		}
	}

	return LUFact(P.transpose(), unit_lower(A), upper(A), Q.transpose());
}

}
}
