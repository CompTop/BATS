#pragma once

#include <map>
#include <unordered_map>
#include <set>
#include <linalg/sparse_vector.h>
#include <linalg/col_matrix.h>

#define p2c_type std::unordered_map<size_t, size_t>

// perform reduction algorithm on a column matrix in-place
template <class TVec>
p2c_type reduce_matrix(ColumnMatrix<TVec> &M) {

	p2c_type pivot_to_col;

	// loop over columns
	for (size_t j = 0; j < M.ncol(); j++) {
		while(M[j].nnz() > 0) {
			// std::cout << j << " : ";
			// M[j].print_row();
			// piv is index-value nzpair
			auto piv = M[j].lastnz();
			if (pivot_to_col.count(piv.ind) > 0) {
				// eliminate pivot
				size_t k = pivot_to_col[piv.ind];
				auto a = piv.val / M[k].lastnz().val;
				M[j].axpy(-a, M[k]);
			} else {
				// new pivot
				pivot_to_col[piv.ind] = j;
				break;
			}
		}
	}
	return pivot_to_col;
}

// perform reduction algorithm on a column matrix in-place
// apply change of basis to U
// invariant M * U = R
template <class TVec>
p2c_type reduce_matrix(ColumnMatrix<TVec> &M, ColumnMatrix<TVec> &U) {

	p2c_type pivot_to_col;

	// loop over columns
	for (size_t j = 0; j < M.ncol(); j++) {
		while(M[j].nnz() > 0) {
			// std::cout << j << " : ";
			// M[j].print_row();
			// piv is index-value nzpair
			auto piv = M[j].lastnz();
			if (pivot_to_col.count(piv.ind) > 0) {
				// eliminate pivot
				size_t k = pivot_to_col[piv.ind];
				auto a = piv.val / M[k].lastnz().val;
				M[j].axpy(-a, M[k]);
				U[j].axpy(-a, U[k]); // update change of basis
			} else {
				// new pivot
				pivot_to_col[piv.ind] = j;
				break;
			}
		}
	}
	return pivot_to_col;
}

// extract indices for homology-revealing basis
// Rk - reduced boundary matrix
// p2ck1 - p2c returned by reduction of B{k+1}
// dimk - number of cells in dimension k
template <typename MT>
std::set<size_t> extract_basis_indices(MT &Rk, p2c_type &p2ck1) {

	std::set<size_t> I;
	size_t dimk = Rk.ncol();

	for (size_t j = 0; j < dimk; j++) {
		if (Rk[j].nnz() == 0) {
			// column is zero
			if (p2ck1.count(j) == 0) {
				// column is not a pivot, so add
				I.emplace(j); // TODO: can use emplace_hint with pointer to back
			}
		}
	}

	return I;
}

// extract indices for homology-revealing basis
// Rk - reduced boundary matrix
// no p2c for dimension above
// dimk - number of cells in dimension k
template <typename MT>
std::set<size_t> extract_basis_indices(MT &Rk) {

	std::set<size_t> I;
	size_t dimk = Rk.ncol();

	for (size_t j = 0; j < dimk; j++) {
		if (Rk[j].nnz() == 0) {
			// column is zero
			// column is not a pivot, so add
			I.emplace(j); // TODO: can use emplace_hint with pointer to back
		}
	}

	return I;
}
