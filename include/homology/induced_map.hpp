#pragma once
/*
compute induced maps on homology
*/
#include <iostream>
#include "reduction.hpp"
#include "basis.hpp"
#include "cohom_basis.hpp"
#include <chain/chain_map.hpp>

#include <linalg/sparse_vector.hpp>
#include <linalg/col_matrix.hpp>

namespace bats {


// get induced map on homology for dimension k
template <class TVec>
ColumnMatrix<TVec> induced_map(
	const ChainMap<ColumnMatrix<TVec>> &F,
	const ReducedChainComplex<ColumnMatrix<TVec>> &C,
	const ReducedChainComplex<ColumnMatrix<TVec>> &D,
	size_t k
) {

	// Check dimensions
	// F: C -> D
	if (F[k].ncol() != C.U[k].nrow() || F[k].nrow() != D.U[k].ncol()) {
		throw std::runtime_error("ChainMap dimensions do not agree with basis dimensions!");
	}

	std::vector<TVec> col;
	// iterate over homology generators in C
	for (auto it = C.I[k].cbegin(); it != C.I[k].cend(); it++) {
		// put image in homology revealing basis in D
		auto y = u_solve(D.U[k], F[k] * C.U[k][*it]);
		// find preferred representative for homology class
		D.find_preferred_representative(y, k);
		// extract indices
		col.emplace_back(y[D.I[k]]);
	}
	return ColumnMatrix<TVec>(D.I[k].size(), C.I[k].size(), col);
}

/**
obtain induced map on cohomology for dimension k

We assume that F is a chain map, so it is dualized before computing the map

@param F ChainMap
@param C ReducedCoChainComplex
@param D ReducedCoChainComplex
@parak k dimension
*/
template <class TVec>
ColumnMatrix<TVec> induced_map(
	const ChainMap<ColumnMatrix<TVec>> &F,
	const ReducedCochainComplex<ColumnMatrix<TVec>> &C,
	const ReducedCochainComplex<ColumnMatrix<TVec>> &D,
	size_t k
) {

	// Check dimensions
	// F*: D -> C
	if (F[k].ncol() != C.R[k].nrow() || F[k].nrow() != D.R[k].nrow()) {
		throw std::runtime_error("ChainMap dimensions do not agree with basis dimensions!");
	}

	std::vector<TVec> col;
	// iterate over homology generators in D
	// F*: D -> C
	for (auto it = D.I[k].cbegin(); it != D.I[k].cend(); ++it) {
		// put image in homology revealing basis in D
		auto y = u_solve(C.U[k+1], D.U[k+1][*it] * F[k]); // F[k].T * D.U[k]
		// find preferred representative for homology class
		C.find_preferred_representative(y, k);
		// extract indices
		col.emplace_back(y[C.I[k]]);
	}
	return ColumnMatrix<TVec>(C.I[k].size(), D.I[k].size(), col);
}

/**
obtain induced map on cohomology for dimension k

We assume that F is a chain map, so it is dualized before computing the map

@param F ChainMap
@param C ReducedCoChainComplex
@param D ReducedCoChainComplex
@parak k dimension
*/
template <class TVec>
ColumnMatrix<TVec> induced_map(
	const DGLinearMap<ColumnMatrix<TVec>> &F,
	const ReducedDGVectorSpace<ColumnMatrix<TVec>> &C,
	const ReducedDGVectorSpace<ColumnMatrix<TVec>> &D,
	size_t k
) {

	if (C.degree != D.degree) {
		throw std::runtime_error("Degree of DGVectorSpaces must agree!");
	}
	int degree = C.degree;
	int dg_offset = (degree == -1) ? 0 : 1;

	// check dimensions
	if (F[k].ncol() != C.R[k+dg_offset].ncol() || F[k].nrow() != D.R[k+dg_offset].ncol()) {
		throw std::runtime_error("ChainMap dimensions do not agree with basis dimensions!");
	}

	std::vector<TVec> col;
	col.reserve(C.hdim(k));

	for (auto it = C.I[k+dg_offset].cbegin(); it != C.I[k+dg_offset].cend(); ++it) {
		// put image in homology revealing basis in D
		auto y = u_solve(D.U[k+dg_offset], F[k] * C.U[k+dg_offset][*it]);
		// find preferred representative for homology class
		D.find_preferred_representative(y, k);
		// extract indices
		col.emplace_back(y[D.I[k+dg_offset]]);
	}
	return ColumnMatrix<TVec>(D.hdim(k), C.hdim(k), col);

}

} // namespace bats
