#pragma once
#include <unordered_map>
#include <vector>
#include <utility>
#include <linalg/sparse_vector.hpp>
#include <linalg/col_matrix.hpp>
#include "reduction.hpp"

#include <omp.h> // openMP header

namespace bats {

struct divide_conquer_flag {};


/**
reduce a single column of the matrix M

@param M matrix
@param U change of basis matrix
@param j column to reduce
@param p2c map from pivots to columns
@param tmp preallocated for faster axpys

assumes that p2c only contains columns < j
*/
template <typename TVec>
void reduce_column_standard(
	ColumnMatrix<TVec>& M,
	ColumnMatrix<TVec>& U,
	const size_t j,
	std::unordered_map<size_t, size_t>& p2c,
	typename TVec::tmp_type& tmp
) {
	while(M[j].nnz() > 0) {
		// piv is index-value nzpair
		auto piv = M[j].lastnz();

		if (p2c.count(piv.ind) > 0) {
			// eliminate pivot
			size_t k = p2c[piv.ind];
			auto a = piv.val / M[k].lastnz().val;
			M[j].axpy(-a, M[k], tmp);
			U[j].axpy(-a, U[k], tmp); // update change of basis
		} else {
			// new pivot
			p2c[piv.ind] = j;
			break;
		}
	}
}

/**
reduce a block of columns sequentially

@param M matrix to be reduced
@param U change of basis matrix
@param j0 start of column range
@param j1 upper bound of column range
@return p2c map from pivots to columns for this block
@return nzcol vector of nonzero columns after reduction
*/
template <typename TVec>
std::tuple<std::unordered_map<size_t,size_t>, std::vector<size_t>> reduce_block_sequential(
	ColumnMatrix<TVec>& M,
	ColumnMatrix<TVec>& U,
	const size_t j0,
	const size_t j1
) {
	std::vector<size_t> nzcol;
	std::unordered_map<size_t, size_t> p2c;
	typename TVec::tmp_type tmp; // preallocate

	for (size_t j = j0; j < j1; ++j) {
		reduce_column_standard(M, U, j, p2c, tmp);
		if (M[j].nnz() > 0) { nzcol.emplace_back(j); }
	}

	return std::make_tuple(p2c, nzcol);
}



/**
reduce a block of columns via divide and conquer

@param M matrix to be reduced
@param U change of basis matrix
@param j0 start of column range
@param j1 upper bound of column range
@param max_block_size maximum block size for sequential base case
@return p2c map from pivots to columns for this block
@return nzcol vector of nonzero columns after reduction
*/
template <typename TVec>
std::tuple<std::unordered_map<size_t,size_t>, std::vector<size_t>> reduce_block_dq(
	ColumnMatrix<TVec>& M,
	ColumnMatrix<TVec>& U,
	const size_t j0,
	const size_t j1,
	const size_t max_block_size
) {
	// determine if we should run base case
	if (j1 - j0 < max_block_size) {
		return reduce_block_sequential(M, U, j0, j1);
	}

	// else we divide and conquer
	size_t j2 = (j0 + j1) >> 1; // midpoint
	std::vector<size_t> nzcol1, nzcol2;
	std::unordered_map<size_t, size_t> p2c;

	// left half
	// TODO: OMP task
	// std::cout << "dq: " << j0 << " -- " << j2 << std::endl;
	std::tie(p2c, nzcol1) = reduce_block_dq(M, U, j0, j2, max_block_size);
	// right half
	// TODO: OMP task
	// std::cout << "dq: " << j2 << " -- " << j1 << std::endl;
	std::tie(std::ignore, nzcol2) = reduce_block_dq(M, U, j2, j1, max_block_size);

	// now we combine
	// can simply iterate over the right block adding to the data for the
	// left block
	typename TVec::tmp_type tmp; // preallocate
	for (auto j : nzcol2) {
		reduce_column_standard(M, U, j, p2c, tmp);
		if (M[j].nnz() > 0) { nzcol1.emplace_back(j); }
	}

	return std::make_tuple(p2c, nzcol1);
}


// behavior with divide conquer flag
template <class TVec>
inline p2c_type reduce_matrix(
	ColumnMatrix<TVec> &M,
	ColumnMatrix<TVec>& U,
	divide_conquer_flag
) {
	size_t n = M.ncol();

	size_t max_block_size = 1024;

	std::unordered_map<size_t, size_t> p2c;
	std::tie(p2c, std::ignore) = reduce_block_dq(M, U, 0, n, max_block_size);

	// fill in return vector
	p2c_type pivot_to_col(M.nrow(), bats::NO_IND);
	for (auto& it : p2c) {
		pivot_to_col[it.first] = it.second;
	}
	return pivot_to_col;
}

/**
reduce a single column of the matrix M

@param M matrix
@param j column to reduce
@param p2c map from pivots to columns
@param tmp preallocated for faster axpys

assumes that p2c only contains columns < j
*/
template <typename TVec>
void reduce_column_standard(
	ColumnMatrix<TVec>& M,
	const size_t j,
	std::unordered_map<size_t, size_t>& p2c,
	typename TVec::tmp_type& tmp
) {
	while(M[j].nnz() > 0) {
		// piv is index-value nzpair
		auto piv = M[j].lastnz();

		if (p2c.count(piv.ind) > 0) {
			// eliminate pivot
			size_t k = p2c[piv.ind];
			auto a = piv.val / M[k].lastnz().val;
			M[j].axpy(-a, M[k], tmp);
		} else {
			// new pivot
			p2c[piv.ind] = j;
			break;
		}
	}
}

/**
reduce a single column of the matrix M
continue eliminating entries after finding pivot

@param M matrix
@param j column to reduce
@param p2c map from pivots to columns
@param tmp preallocated for faster axpys

assumes that p2c only contains columns < j
*/
template <typename TVec>
void reduce_column_extra(
	ColumnMatrix<TVec>& M,
	const size_t j,
	std::unordered_map<size_t, size_t>& p2c,
	typename TVec::tmp_type& tmp
) {
	size_t end_offset = 1;
	auto piv = M[j].nzend() - end_offset; // nonzero location
	while(piv - M[j].nzbegin() > 0) {
		// piv is index-value nzpair

		if (p2c.count(piv->ind) > 0) {
			// eliminate pivot
			size_t k = p2c[piv->ind];
			auto a = piv->val / M[k].lastnz().val;
			M[j].axpy(-a, M[k], tmp);
			piv = M[j].nzend() - end_offset; // next nonzero location
		} else if (end_offset == 1){
			// new pivot
			p2c[piv->ind] = j;
			++end_offset;
			--piv;
		} else {
			++end_offset;
			--piv;
		}
	}
}

/**
reduce a block of columns sequentially

@param M matrix to be reduced
@param j0 start of column range
@param j1 upper bound of column range
@return p2c map from pivots to columns for this block
@return nzcol vector of nonzero columns after reduction
*/
template <typename TVec>
std::tuple<std::unordered_map<size_t,size_t>, std::vector<size_t>> reduce_block_sequential(
	ColumnMatrix<TVec>& M,
	const size_t j0,
	const size_t j1
) {
	std::vector<size_t> nzcol;
	std::unordered_map<size_t, size_t> p2c;
	typename TVec::tmp_type tmp; // preallocate

	for (size_t j = j0; j < j1; ++j) {
		reduce_column_standard(M, j, p2c, tmp);
		// reduce_column_extra(M, j, p2c, tmp);
		if (M[j].nnz() > 0) { nzcol.emplace_back(j); }
	}

	return std::make_tuple(p2c, nzcol);
}

/**
reduce a block of columns via divide and conquer

@param M matrix to be reduced
@param j0 start of column range
@param j1 upper bound of column range
@param max_block_size maximum block size for sequential base case
@return p2c map from pivots to columns for this block
@return nzcol vector of nonzero columns after reduction
*/
template <typename TVec>
std::tuple<std::unordered_map<size_t,size_t>, std::vector<size_t>> reduce_block_dq(
	ColumnMatrix<TVec>& M,
	const size_t j0,
	const size_t j1,
	const size_t max_block_size,
	const size_t level,
	const size_t max_depth
) {
	// std::cout << j0 << " -- " << j1 << "\t level=" << level << "/" << max_depth << std::endl;
	// determine if we should run base case
	if (j1 - j0 < max_block_size || level > max_depth) {
		return reduce_block_sequential(M, j0, j1);
	}

	// else we divide and conquer
	size_t j2 = (j0 + j1) >> 1; // midpoint
	std::vector<size_t> nzcol1, nzcol2;
	std::unordered_map<size_t, size_t> p2c, p2c_ignore;

	// left half
	// std::cout << "dq: " << j0 << " -- " << j2 << std::endl;
	#pragma omp task shared(M) shared(p2c, nzcol1) firstprivate(j0, j2, max_block_size, level, max_depth)
	std::tie(p2c, nzcol1) = reduce_block_dq(M, j0, j2, max_block_size, level+1, max_depth);
	// right half
	// std::cout << "dq: " << j2 << " -- " << j1 << std::endl;
	#pragma omp task shared(M) shared(nzcol2, p2c_ignore) firstprivate(j2, j1, max_block_size, level, max_depth)
	std::tie(p2c_ignore, nzcol2) = reduce_block_dq(M, j2, j1, max_block_size, level+1, max_depth);
	#pragma omp taskwait

	// now we combine
	// can simply iterate over the right block adding to the data for the
	// left block
	typename TVec::tmp_type tmp; // preallocate
	for (auto j : nzcol2) {
		reduce_column_standard(M, j, p2c, tmp);
		if (M[j].nnz() > 0) { nzcol1.emplace_back(j); }
	}

	return std::make_tuple(p2c, nzcol1);
}

// behavior with divide conquer flag
template <class TVec>
inline p2c_type reduce_matrix(
	ColumnMatrix<TVec> &M,
	divide_conquer_flag
) {
	size_t n = M.ncol();

	size_t max_block_size = 256;
	size_t max_depth = 2;

	std::unordered_map<size_t, size_t> p2c;
	std::vector<size_t> nz_ignore;
	// begin parallel region
    #pragma omp parallel default(none) shared(p2c, nz_ignore, M, n, max_block_size, max_depth)
    {
        #pragma omp single nowait
        {
			std::tie(p2c, nz_ignore) = reduce_block_dq(M, 0, n, max_block_size, 0, max_depth);
		}
	}

	// fill in return vector
	p2c_type pivot_to_col(M.nrow(), bats::NO_IND);
	for (auto& it : p2c) {
		pivot_to_col[it.first] = it.second;
	}
	return pivot_to_col;
}

/**
do an initial parallel sweep to zero out columns as possible

@param M matrix to reduce
@param block_size size of column blocks to be processed in parallel
*/
template <typename TVec>
void partial_reduce_parallel(
	ColumnMatrix<TVec>& M,
	const size_t block_size
) {
	size_t n = M.ncol();

	#pragma omp parallel for
	for (size_t j = 0; j < n-block_size; j+=block_size) {
		reduce_block_sequential(M, j, j+block_size);
	}
	return;
}
/**
do an initial parallel sweep to zero out columns as possible

@param M matrix to reduce
@param U change of basis
@param block_size size of column blocks to be processed in parallel
*/
template <typename TVec>
void partial_reduce_parallel(
	ColumnMatrix<TVec>& M,
	ColumnMatrix<TVec>& U,
	const size_t block_size
) {
	size_t n = M.ncol();

	#pragma omp parallel for
	for (size_t j = 0; j < n-block_size; j+=block_size) {
		reduce_block_sequential(M, U, j, j+block_size);
	}
	return;
}

} // namespace bats
