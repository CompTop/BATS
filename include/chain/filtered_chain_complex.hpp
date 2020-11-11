#pragma once

#include "chain_complex.h"
#include <filtration/filtration.h>

// template over filtration and matrix type
template <typename FT, typename MT>
struct FilteredChainComplex {
	std::vector<std::vector<FT>> val; // stored in original order
	ChainComplex<MT> C; // stored in permutation order
	std::vector<std::vector<size_t>> iperm; // inverse permutation from permutaiton order to original order

	FilteredChainComplex() {}

	template <typename CpxT>
	FilteredChainComplex(const Filtration<FT, CpxT> &F) : val(F.vals()), C(F.complex()) {
		// step 1: compute permutation to put val in order

		// step 2: put ChainComplex C in permutation order

		// step 3: store inverse perumutation to map back to original order

	}

	inline size_t dim(const size_t k) { return C.dim[k]; }

	inline const ChainComplex<MT>& complex() const {return C;}
	inline const std::vector<std::vector<FT>>& vals() const { return val; }

	// update filtration
	void update_filtration(const std::vector<std::vector<T>> newval) {
		// step 1: determine permutation order for newval

		// step 2: determine update to old permutation

		// step 3: apply permutation updates to ChainComplex C

	}

};

// defualt return
template <typename FT, typename T, typename CpxT>
inline auto __FilteredChainComplex(const Filtration<FT, CpxT> &F, T) {
	using VT = SparseVector<T, size_t>;
	using MT = ColumnMatrix<VT>;

	return FilteredChainComplex<FT, MT>(F);
}
