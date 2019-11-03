#pragma once
/*
compute homology-revealing bases for a chain complex
*/
#include <vector>
#include <set>
#include <chain/chain_complex.h>
#include "reduction.h"

template <typename MT>
struct ReducedChainComplex {

	std::vector<MT> U; // basis matrices
	std::vector<MT> R; // reduced matrix
	std::vector<std::set<size_t>> I;
	std::vector<p2c_type> p2c;

	// compute reduced chain complex from chain complex
	ReducedChainComplex(ChainComplex<MT> &C) {
		size_t dmax = C.maxdim() + 1;
		U.resize(dmax);
		R.resize(dmax);
		I.resize(dmax);
		p2c.resize(dmax);

		// TODO: can parallelize this
		for (size_t k = 0; k < dmax; k++) {
			U[k] = MT::identity(C.dim[k]);
			R[k] = C.boundary[k];
			p2c[k] = reduce_matrix(R[k], U[k]);
		}
		// TODO: can parallelize this
		for (size_t k = 0; k < dmax-1; k++) {
			I[k] = extract_basis_indices(R[k], p2c[k+1]);
		}
		I[dmax-1] = extract_basis_indices(R[dmax-1]);
	}

	// size of homology vector space in dimension k
	inline size_t hdim(size_t k) const { return I[k].size(); }
	inline size_t maxdim() const { return R.size() - 1; }
};
