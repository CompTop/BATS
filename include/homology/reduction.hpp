#pragma once

#include <map>
#include <unordered_map>
#include <set>
#include <utility>
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
struct apparent_pairs_flag {};
struct no_apparent_pairs_flag {};

// flags for selecting algorithms
struct standard_reduction_flag {};
struct extra_reduction_flag {};
struct sparse_reduction_flag {};


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
	// size_t ct = 0;
	for (size_t j = 0; j < M.ncol(); j++) {
		while(M[j].nnz() > 0) {
			// ++ct;
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
	// std::cout << "# iterations = " << ct << std::endl;
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
			} else if (end_offset == 1) {
				// new pivot
				pivot_to_col[piv->ind] = j;
				++end_offset;
				--piv;
			} else {
				// we skip zeroing out this entry
				// need to increment offset
				++end_offset;
				--piv;
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

/**
greedily introduce sparsity into column j of M
using columns k < j
has the effect of reducing column j if not already reduced.

@param M matrix
@param j column to reduce
@param pivot_to_col maps pivots to columns
@param coeff preallocated map to use when sparsifying.
@param tmp preallocated to use with axpy
*/
template <typename TVec, typename F>
void reduce_column_sparsify(
	ColumnMatrix<TVec>& M,
	const size_t j,
	p2c_type& pivot_to_col,
	std::map<F, size_t>& coeff,
	typename TVec::tmp_type& tmp
) {
	size_t end_offset = 1; // start looking at last nonzero
	auto piv = M[j].nzend() - end_offset; // nonzero location
	while (piv - M[j].nzbegin() > 0) {
		if (pivot_to_col[piv->ind] == j) {
			// do nothing, because this is already a pivot
			++end_offset;
			--piv;
		} else if (pivot_to_col[piv->ind] < j && end_offset == 1) {
			// we simply eliminate this non-zero
			size_t k = pivot_to_col[piv->ind];
			auto a = piv->val / M[k].lastnz().val;
			M[j].axpy(-a, M[k], tmp);
			piv = M[j].nzend() - end_offset; // next nonzero location
		} else if (pivot_to_col[piv->ind] < j) {
			// determine if we should eliminate this non-zero
			size_t i = piv->ind;
			size_t k = pivot_to_col[i];
			// determine largest number of non-zeros we might eliminate
			coeff.clear();
			M[j].coeff_intersection(M[k], coeff); // calculate intersections
			F a(0);
			size_t ct = 0;
			for (auto it : coeff) {
				if (it.second > ct) {
					ct = it.second;
					a = it.first;
				}
			}
			// determine whether or not we would introduce more non-zeros than we eliminate
			if (ct > (M[k].nnz() / 2)) {
				M[j].axpy(-a, M[k], tmp);
			}
			piv = M[j].nzend() - end_offset; // next nonzero location
			if (piv->ind == i) {++end_offset; --piv;} // if we didn't eliminate this nonzero, continue
		} else if (end_offset == 1) {
			// we can set a new pivot
			pivot_to_col[piv->ind] = j;
			++end_offset;
			--piv;
		} else {
			// we skip zeroing out this entry
			// need to increment offset
			++end_offset;
			--piv;
		}
	}
	return;
}

/**
greedily introduce sparsity into columns j of U and R
using columns k < j
assumes columns k <= j are already reduced

objective to greedily minimize is nnz(U[j]) + nnz(R[j])

@param R reduced matrix
@param U change of basis matrix
@param j column
@param coeff preallocated map to use when sparsifying.
@param tmp preallocated to use with axpy
*/
template <typename TVec, typename F>
void sparsify_basis(
	ColumnMatrix<TVec>& R,
	ColumnMatrix<TVec>& U,
	const size_t j,
	std::map<F, size_t>& coeff,
	typename TVec::tmp_type& tmp
) {
	size_t end_offset = 2; // entry above diagonal in U[j]
	auto piv = U[j].nzend() - end_offset; // nonzero location
	while (piv - U[j].nzbegin() > 0) {
		size_t k = piv->ind;
		// std::cout << k << std::endl;
		if (R[k].nnz() == 0 || R[k].lastnz().ind < j) {
			// can potentially modify using this entry
			coeff.clear();
			U[j].coeff_intersection(U[k], coeff); // calculate intersections
			R[j].coeff_intersection(R[k], coeff); // add counts for matrix R
			F a(0);
			size_t ct = 0;
			for (auto it : coeff) {
				if (it.second > ct) {
					ct = it.second;
					a = it.first;
				}
			}
			// determine whether or not we would introduce more non-zeros than we eliminate
			if (ct > ((U[k].nnz() + R[k].nnz()) / 2)) {
				U[j].axpy(-a, U[k], tmp);
				R[j].axpy(-a, R[k], tmp);
			}
			piv = U[j].nzend() - end_offset; // next nonzero location
			if (piv->ind == k) {++end_offset; --piv;} // if we didn't eliminate this nonzero, continue
		} else {
			++end_offset;
			--piv;
		}
	}
	return;
}

/**
greedily introduce sparsity into columns of U and R
assumes R is already reduced

objective to greedily minimize is nnz(U[j]) + nnz(R[j])

@param R reduced matrix
@param U change of basis matrix
*/
template <typename TVec>
void sparsify_basis(
	ColumnMatrix<TVec>& R,
	ColumnMatrix<TVec>& U
) {

	if (R.ncol() != U.ncol()) {throw std::runtime_error("Number of columns are not the same!");}

	using F = typename TVec::val_type;
	std::map<F, size_t> coeff;
	typename TVec::tmp_type tmp;

	for (size_t j = 0; j < R.ncol(); ++j) {
		sparsify_basis(R, U, j, coeff, tmp);
	}

	return;
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



/**
Update change of basis matrix to not be as dense by
removing lower grade cycles from U.

Let j1 < j2, and R[j1] = 0
Then, we can add a linear combination of U[j1] to U[j2]
without changing the matrix invatiant B*U = R

@param R reduced matrix
@param U change of basis matrix
*/
template <class TVec>
void remove_extra_cycles(
	const ColumnMatrix<TVec> &R,
	ColumnMatrix<TVec> &U
) {

	if (R.ncol() != U.ncol()) {throw std::runtime_error("Number of columns are not the same!");}

	// create a temporary vector for use in axpys
	typename TVec::tmp_type tmp;

	size_t n = U.ncol(); // number of columns
	for (size_t j = 0; j < n; ++j) {
		size_t end_offset = 2; // entry above diagonal
		auto piv = U[j].nzend() - end_offset; // nonzero location
		while (piv - U[j].nzbegin() > 0) {
			size_t k = piv->ind;
			if (R[k].nnz() == 0) {
				// we can eliminate entry
				auto a = piv->val / U[k].lastnz().val;
				U[j].axpy(-a, U[k], tmp); // update change of basis
				// R[j].axpy(R[k]) does nothing because R[k] = 0
				piv = U[j].nzend() - end_offset; // next nonzero location
			} else {
				// we skip zeroing out this entry
				// need to increment offset
				end_offset++;
				piv--;
			}
		}
	}

}

// TODO: write a function on SparseVector which
// will calculate the coefficient to use to greedily reduce the nnz
// by adding another vector (may be zero)

// might also write a function that will try to minimize
// nnz(U) + nnz(R)
// since we can add copies of columns of R with smaller pivots as well

namespace detail {

/**
Calculate

@return c coefficient for a += c*b
@return d difference in number of nonzeros of a using this coefficient
*/
template <typename TVec>
auto pivot_coeff(const TVec& a, const TVec& b) {

}

} // namespace detail



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
