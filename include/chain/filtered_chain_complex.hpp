#pragma once

#include "chain_complex.hpp"
#include <filtration/filtration.hpp>

namespace bats {

// template over filtration and matrix type
template <typename FT, typename MT>
struct FilteredChainComplex {
	std::vector<std::vector<FT>> val; // stored in original order
	ChainComplex<MT> C; // stored in permutation order
	std::vector<std::vector<size_t>> perm; // map from column order in chain complex to val

	FilteredChainComplex() {}

	template <typename CpxT>
	FilteredChainComplex(const Filtration<FT, CpxT> &F) : val(F.vals()), C(F.complex()) {
		// step 1: compute permutation to put val in order
		auto perms = filtration_sortperm(val);

		// step 2: put ChainComplex C in permutation order
		C.ipermute_basis(perms);

		// step 3: store inverse perumutation to map back to original order
		perm = perms;
	}

	inline size_t dim(const size_t k) { return C.dim(k); }

	inline const ChainComplex<MT>& complex() const {return C;}
	inline const std::vector<std::vector<FT>>& vals() const { return val; }

	// update filtration
	void update_filtration(const std::vector<std::vector<FT>> newval) {
		// step 1: determine permutation order for newval
		auto perms = filtration_sortperm(newval);

		// get inverse permutation for current permutation
		auto iperm = filtration_iperm(perm);

		// step 2: determine update to old permutation
		// iperm[k] will hold the updated permutation temporarily
		std::vector<FT> tmp; // temporary vector
		for (size_t k = 0; k < iperm.size(); k++) {
			bats::util::apply_perm(iperm[k].data(), tmp, perms[k]);
		}

		// step 3: apply permutation updates to ChainComplex C
		C.permute_basis(iperm);

		// step 4: store new permutation
		perm = perms;

		// store values as well
		val = newval;

	}

	template <typename I>
	void update_filtration_general(const I & updating_information){
		// step 1: update ChainComplex C by new filtered filtration information
		C.update_basis_general(updating_information);
		
		// step 2: store new permutation
		perm = updating_information.F_Y_perms;

		// store values as well
		val = updating_information.F_Y_vals;
	}

};

// defualt return
template <typename FT, typename T, typename CpxT>
inline auto __FilteredChainComplex(const Filtration<FT, CpxT> &F, T) {
	using VT = SparseVector<T, size_t>;
	using MT = ColumnMatrix<VT>;

	return FilteredChainComplex<FT, MT>(F);
}

// defualt return
template <typename FT, typename T, typename CpxT>
inline auto Chain(const Filtration<FT, CpxT> &F, T) {return __FilteredChainComplex(F, T());}


} // namespace bats
