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

template <class TVec>
void find_preferred_representative(
	TVec &y,
	const ColumnMatrix<TVec> &bdry,
	const p2c_type &p2c
) {
	size_t j = bdry.nrow();
	auto yit = y.upper_bound(j); // find iterator past last nonzero
	while (yit != y.nzbegin()) {
		--yit;
		j = yit->ind;
		// find if j is a pivot of bdry
		if (p2c.count(j) > 0) {
			auto i = p2c.at(j);
			auto Ui_it = bdry[i].lastnz();
			auto a = (yit->val) / Ui_it.val;
			y.axpy(-a, bdry[i]);

			// get next nonzero
			yit = y.upper_bound(j-1);
		}
		// else j is a preferred representative
		// so we do nothing
	}
}

// get induced map on homology for dimension k
template <class TVec>
ColumnMatrix<TVec> induced_map_maxdim(
	const ChainMap<ColumnMatrix<TVec>> &F,
	const ReducedChainComplex<ColumnMatrix<TVec>> &C,
	const ReducedChainComplex<ColumnMatrix<TVec>> &D,
	size_t k
) {

	std::vector<TVec> col;
	// iterate over homology generators in C
	for (auto it = C.I[k].cbegin(); it != C.I[k].cend(); it++) {
		auto y = F[k] * C.U[k][*it];
		// extract indices
		// every representative is preferred because no boundary k+1
		col.emplace_back(y[D.I[k]]);
	}
	return ColumnMatrix<TVec>(D.I[k].size(), C.I[k].size(), col);
}


// get induced map on homology for dimension k
template <class TVec>
ColumnMatrix<TVec> induced_map(
	const ChainMap<ColumnMatrix<TVec>> &F,
	const ReducedChainComplex<ColumnMatrix<TVec>> &C,
	const ReducedChainComplex<ColumnMatrix<TVec>> &D,
	size_t k
) {
	if ( k == D.maxdim() ) {
		return induced_map_maxdim(F, C, D, k);
	}

	std::vector<TVec> col;
	ColumnMatrix<TVec> bdry = u_solve(D.U[k], D.R[k+1]); // TODO: if k is max dimension
	// iterate over homology generators in C
	for (auto it = C.I[k].cbegin(); it != C.I[k].cend(); it++) {
		auto y = F[k] * C.U[k][*it];
		find_preferred_representative(y, bdry, D.p2c[k+1]);
		// extract indices
		col.emplace_back(y[D.I[k]]);
	}
	return ColumnMatrix<TVec>(D.I[k].size(), C.I[k].size(), col);
}
