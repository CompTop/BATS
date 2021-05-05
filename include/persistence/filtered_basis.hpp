#pragma once

#include <vector>
#include <limits>
#include <chain/filtered_chain_complex.hpp>
#include <homology/basis.hpp>
#include <util/common.hpp>
#include "barcode.hpp"

namespace bats {

template <typename T, typename MT>
struct ReducedFilteredChainComplex {

	ReducedChainComplex<MT> RC; // reduced in permutation order
	std::vector<std::vector<T>> val; // stored in original order
	std::vector<std::vector<size_t>> iperm; // inverse permutation from permutaiton order to original order

	ReducedFilteredChainComplex() {}

	// variadic template for passing arguments
	template <typename... Args>
	ReducedFilteredChainComplex(const FilteredChainComplex<T, MT>& C, Args ...args) :
		RC(C.complex(), args...),
		val(C.val),
		iperm(C.iperm) {}

	inline size_t maxdim() const { return RC.maxdim(); }
	inline size_t dim(const size_t k) const {return RC.dim(k);}
	inline size_t hdim(const size_t k) const {return RC.hdim(k);}

	// persistence pairs in dimension k
	std::vector<PersistencePair<T>> persistence_pairs(const size_t k) {
		std::vector<PersistencePair<T>> pairs;
		for (size_t i =0; i < dim(k); i++) {
			if (RC.R[k][i].nnz() == 0) {
				// homology generated
				if (k == maxdim() || RC.p2c[k+1][i] == bats::NO_IND)  {
					// infinite bar
					pairs.emplace_back(
						PersistencePair(k, i, bats::NO_IND,
							val[k][iperm[k][i]], std::numeric_limits<T>::infinity()
							// val[k][i], std::numeric_limits<T>::infinity()
						)
					);
				} else {
					size_t j = RC.p2c[k+1][i];
					// finite bar
					pairs.emplace_back(
						PersistencePair(k, i, j,
							val[k][iperm[k][i]], val[k+1][iperm[k+1][j]]
							// val[k][i], val[k+1][j]
						)
					);
				}
			}

		}
		return pairs;
	}

	// return representative for a pair
	inline auto representative(const PersistencePair<T>& p) {
		return RC.U[p.dim][p.birth_ind];
	}

	// barcode w/out critical inds
	std::vector<T> barcode(const size_t k);

	// critical cells for barcode in dimension k
	std::vector<size_t> critical_cells(const size_t k);

	// update filtration
	void update_filtration(const std::vector<std::vector<T>> newval) {
		// step 1: determine permutation order for newval
		auto perms = filtration_sortperm(newval);

		// step 2: determine update to old permutation
		// iperm[k] will hold the updated permutation temporarily
		std::vector<T> tmp; // temporary vector
		for (size_t k = 0; k < iperm.size(); k++) {
			bats::util::apply_perm(iperm[k].data(), tmp, perms[k]);
		}

		// step 3: apply permutation updates to ReducedChainComplex RC
		RC.permute_basis(iperm);


		// step 4: store new inverse permutation
		iperm = filtration_iperm(perms);

		// store values as well
		val = newval;

	}

	// get subcomplex
	ReducedFilteredChainComplex get_subcomplex() const;

	void print_summary() const {
		std::cout << "ReducedFilteredChainComplex with " << maxdim() << " dimensions:" << std::endl;
		for (size_t k = 0; k < maxdim() + 1; k++) {
			std::cout << "\tdim " << k << ": " << dim(k)
			<< ", betti_" << k << ": " << hdim(k) << "\n";
		}
	}

};

// defualt return
template <typename FT, typename T, typename CpxT, typename... Args>
inline auto __ReducedFilteredChainComplex(const Filtration<FT, CpxT> &F, T, Args ...args) {
	using VT = SparseVector<T, size_t>;
	using MT = ColumnMatrix<VT>;

	return ReducedFilteredChainComplex(FilteredChainComplex<FT, MT>(F), args...);
}

template <typename FT, typename T, typename CpxT, typename... Args>
inline auto Reduce(const Filtration<FT, CpxT> &F, T, Args ...args) {
	return __ReducedFilteredChainComplex(F, T(), args...);
}

template <typename T, typename MT, typename... Args>
inline auto Reduce(const FilteredChainComplex<T, MT>& C, Args ...args) {
	return ReducedFilteredChainComplex(C, args...);
}


} // namespace bats
