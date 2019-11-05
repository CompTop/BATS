#pragma once
/*
compute induced maps on homology
*/
#include <iostream>
#include "reduction.h"
#include "basis.h"
#include <chain/chain_map.h>

#include <linalg/sparse_vector.h>
#include <linalg/col_matrix.h>


// get induced map on homology for dimension k
template <class TVec>
ColumnMatrix<TVec> induced_map(
	const ChainMap<ColumnMatrix<TVec>> &F,
	const ReducedChainComplex<ColumnMatrix<TVec>> &C,
	const ReducedChainComplex<ColumnMatrix<TVec>> &D,
	size_t k
) {

	std::vector<TVec> col;
	ColumnMatrix<TVec> bdry = u_solve(D.U[k], D.R[k+1]); // TODO: if k is max dimension
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
