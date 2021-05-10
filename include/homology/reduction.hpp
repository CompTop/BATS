#pragma once

#include <map>
#include <unordered_map>
#include <set>
#include <linalg/sparse_vector.hpp>
#include <linalg/col_matrix.hpp>
// #include "tsl/hopscotch_map.h"
// #include "tsl/robin_map.h"
// #include "tsl/sparse_map.h"

// template <typename T>
// struct vector_to_dict {
// 	std::vector<T> vals;
//
// 	template <typename ...Args>
// 	vector_to_dict(Args (&...args)) : vals(args...) {}
//
// 	vector_to_dict(size_t n) : vals(n, bats::NO_IND) {}
//
// 	inline T& operator=(size_t i) {return vals[i];}
// 	inline const T& operator=(size_t i) const {return vals[i];}
// 	inline size_t count(size_t i) const {return vals[i] != bats::NO_IND;}
// };

// #define p2c_type std::unordered_map<size_t, size_t>
// #define p2c_type tsl::hopscotch_map<size_t, size_t>
// #define p2c_type tsl::robin_map<size_t, size_t>
// #define p2c_type tsl::sparse_map<size_t, size_t>
#define p2c_type std::vector<size_t>


namespace bats {

// flags for selecting optimizations
struct no_optimization_flag {};
struct clearing_flag {};
struct compression_flag {};

// flags for selecting algorithms
struct standard_reduction_flag {};
struct extra_reduction_flag {};


// perform reduction algorithm on a column matrix in-place
template <class TVec>
p2c_type reduce_matrix_standard(ColumnMatrix<TVec> &M) {

	// p2c_type pivot_to_col;
	p2c_type pivot_to_col(M.nrow(), bats::NO_IND);
	// create a temporary vector for use in axpys
	typename TVec::tmp_type tmp;

	// loop over columns
	for (size_t j = 0; j < M.ncol(); j++) {
		while(M[j].nnz() > 0) {
			// std::cout << j << " : ";
			// M[j].print_row();
			// piv is index-value nzpair
			auto piv = M[j].lastnz();
			// if (pivot_to_col.count(piv.ind) > 0) {
			if (pivot_to_col[piv.ind] != bats::NO_IND) {
				// eliminate pivot
				size_t k = pivot_to_col[piv.ind];
				auto a = piv.val / M[k].lastnz().val;
				M[j].axpy(-a, M[k], tmp);
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
// do extra reduction to eliminate nonzeros even after pivot has been found
template <class TVec>
p2c_type reduce_matrix_extra(ColumnMatrix<TVec> &M) {

	// p2c_type pivot_to_col;
	p2c_type pivot_to_col(M.nrow(), bats::NO_IND);
	// create a temporary vector for use in axpys
	typename TVec::tmp_type tmp;

	// loop over columns
	for (size_t j = 0; j < M.ncol(); j++) {
		bool found_pivot = false;
		size_t end_offset = 1;
		auto piv = M[j].nzend() - end_offset; // nonzero location
		while(piv - M[j].nzbegin() > 0) {
			// while the nonzero we are looking at is in the vector
			// piv is index-value nzpair
			// auto piv = M[j].lastnz();
			// if (pivot_to_col.count(piv->ind) > 0) {
			if (pivot_to_col[piv->ind] != bats::NO_IND) {
				// eliminate pivot
				size_t k = pivot_to_col[piv->ind];
				auto a = piv->val / M[k].lastnz().val;
				M[j].axpy(-a, M[k], tmp);
				piv = M[j].nzend() - end_offset; // next nonzero location
			} else if (!found_pivot) {
				// new pivot
				pivot_to_col[piv->ind] = j;
				found_pivot = true;
				end_offset++;
				piv--;
			} else {
				// we skip zeroing out this entry
				// need to increment offset
				end_offset++;
				piv--;
			}
		}
	}
	return pivot_to_col;
}

// default behavior
template <class TVec>
inline p2c_type reduce_matrix(ColumnMatrix<TVec> &M) {
	return reduce_matrix_standard(M);
}
// behavior with standard flag
template <class TVec>
inline p2c_type reduce_matrix(ColumnMatrix<TVec> &M, bats::standard_reduction_flag) {
	return reduce_matrix_standard(M);
}
// behavior with extra flag
template <class TVec>
inline p2c_type reduce_matrix(ColumnMatrix<TVec> &M, bats::extra_reduction_flag) {
	return reduce_matrix_extra(M);
}

// perform reduction algorithm on a column matrix in-place
// apply change of basis to U
// invariant M * U = R
template <class TVec>
p2c_type reduce_matrix_standard(ColumnMatrix<TVec> &M, ColumnMatrix<TVec> &U) {

	if (M.ncol() != U.ncol()) {throw std::runtime_error("Number of columns are not the same!");}

	// p2c_type pivot_to_col;
	p2c_type pivot_to_col(M.nrow(), bats::NO_IND);
	// create a temporary vector for use in axpys
	typename TVec::tmp_type tmp;

	// loop over columns
	for (size_t j = 0; j < M.ncol(); j++) {
		while(M[j].nnz() > 0) {
			// std::cout << j << " : ";
			// M[j].print_row();
			// piv is index-value nzpair
			auto piv = M[j].lastnz();
			// if (pivot_to_col.count(piv.ind) > 0) {
			if (pivot_to_col[piv.ind] != bats::NO_IND) {
				// eliminate pivot
				size_t k = pivot_to_col[piv.ind];
				auto a = piv.val / M[k].lastnz().val;
				M[j].axpy(-a, M[k], tmp);
				U[j].axpy(-a, U[k], tmp); // update change of basis
			} else {
				// new pivot
				pivot_to_col[piv.ind] = j;
				break;
			}
		}
	}
	return pivot_to_col;
}

template <class TVec>
p2c_type reduce_matrix_extra(ColumnMatrix<TVec> &M, ColumnMatrix<TVec> &U) {

	if (M.ncol() != U.ncol()) {throw std::runtime_error("Number of columns are not the same!");}

	// p2c_type pivot_to_col;
	p2c_type pivot_to_col(M.nrow(), bats::NO_IND);
	// create a temporary vector for use in axpys
	typename TVec::tmp_type tmp;

	// loop over columns
	for (size_t j = 0; j < M.ncol(); j++) {
		bool found_pivot = false;
		size_t end_offset = 1;
		auto piv = M[j].nzend() - end_offset; // nonzero location
		while(piv - M[j].nzbegin() > 0) {
			// while the nonzero we are looking at is in the vector
			// piv is index-value nzpair
			// auto piv = M[j].lastnz();
			// if (pivot_to_col.count(piv->ind) > 0) {
			if (pivot_to_col[piv->ind] != bats::NO_IND) {
				// eliminate pivot
				size_t k = pivot_to_col[piv->ind];
				auto a = piv->val / M[k].lastnz().val;
				M[j].axpy(-a, M[k], tmp);
				U[j].axpy(-a, U[k], tmp); // update change of basis
				piv = M[j].nzend() - end_offset; // next nonzero location
			} else if (!found_pivot) {
				// new pivot
				pivot_to_col[piv->ind] = j;
				found_pivot = true;
				end_offset++;
				piv--;
			} else {
				// we skip zeroing out this entry
				// need to increment offset
				end_offset++;
				piv--;
			}
		}
	}
	return pivot_to_col;
}

// default behavior
template <class TVec>
inline p2c_type reduce_matrix(ColumnMatrix<TVec> &M, ColumnMatrix<TVec> &U) {
	return reduce_matrix_standard(M, U);
}
// behavior with standard flag
template <class TVec>
inline p2c_type reduce_matrix(ColumnMatrix<TVec> &M, ColumnMatrix<TVec> &U, bats::standard_reduction_flag) {
	return reduce_matrix_standard(M, U);
}
// behavior with extra flag
template <class TVec>
inline p2c_type reduce_matrix(ColumnMatrix<TVec> &M, ColumnMatrix<TVec> &U, bats::extra_reduction_flag) {
	return reduce_matrix_extra(M, U);
}

// get clearing indices from pivots
std::vector<size_t> get_clearing_inds(const p2c_type &p2c) {
	std::vector<size_t> inds;
	inds.reserve(p2c.size());
	size_t i = 0;
	auto it = p2c.begin();
	while (it != p2c.end()) {
		// put keys into indices - these are the clearing indices for
		// one dimension down.
		// inds.emplace_back(it->first);
		if (*it != bats::NO_IND) { inds.emplace_back(i); }
		it++;
		i++;
	}
	std::sort(inds.begin(), inds.end());
	return inds;
}

// perform reduction algorithm on a column matrix in-place
// use clearing optimization
// assume clear_inds is in sorted order
// also known as the twist algorithm
// get clearing indices from one dimension above
template <class TVec, typename flag>
p2c_type reduce_matrix_clearing(
	ColumnMatrix<TVec> &M,
	const std::vector<size_t> &clear_inds,
	flag
) {

	// std::cout << "clearing cols" << std::endl;
	// just zero-out columns ahead of time
	for (auto cind : clear_inds) {
		// we hit an index we want to clear
		M[cind].clear(); // just zero out column without doing work
	}
	// std::cout << "reducing matrix" << std::endl;

	// now run standard reduction algorithm with flag
	return reduce_matrix(M, flag());
}

// get compression indices from reduced matrix R
// compression indices are column indices
// where column has not been zeroed out.
template <class TVec>
std::vector<bool> get_compression_inds(const ColumnMatrix<TVec> &R) {
	std::vector<bool> inds(R.ncol());
	for (size_t j = 0; j < R.ncol(); j++) {
		inds[j] = R[j].nnz() > 0;
	}
	return inds;
}

// reduce matrix using compression
// comp_inds hold indices of rows which will never be pivots
// we zero-out these rows to reduce number of operations
template <class TVec, typename flag>
p2c_type reduce_matrix_compression(
	ColumnMatrix<TVec> &M,
	const std::vector<bool> &comp_inds,
	flag
) {
	// first zero-out rows designated by comp_inds
	M.clear_rows(comp_inds);

	// then run standard reduction algorithm
	return reduce_matrix(M, flag());
}


template <class TVec, typename flag>
p2c_type reduce_matrix_compression(
	ColumnMatrix<TVec> &M,
	ColumnMatrix<TVec> &U,
	const std::vector<bool> &comp_inds,
	flag
) {
	// first zero-out rows designated by comp_inds
	M.clear_rows(comp_inds);

	// we can still record basis.
	return reduce_matrix(M, U, flag());
}
// TODO: add flags



// extract indices for homology-revealing basis
// Rk - reduced boundary matrix
// p2ck1 - p2c returned by reduction of B{k+1}
// dimk - number of cells in dimension k
template <typename MT>
std::vector<size_t> extract_basis_indices(const MT& Rk, const p2c_type &p2ck1) {

	std::vector<size_t> I;
	size_t dimk = Rk.ncol();

	for (size_t j = 0; j < dimk; j++) {
		if (Rk[j].nnz() == 0) {
			// column is zero
			// if (p2ck1.count(j) == 0) {
			if (p2ck1[j] == bats::NO_IND) {
				// column is not a pivot, so add
				I.emplace_back(j); // TODO: can use emplace_hint with pointer to back
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
std::vector<size_t> extract_basis_indices(const MT &Rk) {

	std::vector<size_t> I;
	size_t dimk = Rk.ncol();
	// I.reserve(dimk);

	for (size_t j = 0; j < dimk; j++) {
		if (Rk[j].nnz() == 0) {
			// column is zero
			// column is not a pivot, so add
			I.emplace_back(j); // TODO: can use emplace_hint with pointer to back
		}
	}

	return I;
}

} // namespace bats
