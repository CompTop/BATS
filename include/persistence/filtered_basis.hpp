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

	std::vector<std::vector<T>> val;
	ReducedChainComplex<MT> RC;

	ReducedFilteredChainComplex() {}

	// variadic template for passing arguments
	template <typename... Args>
	ReducedFilteredChainComplex(const FilteredChainComplex<T, MT>& C, Args... args) :
		val(C.vals()),
		RC(C.complex(), args...) {}

	inline size_t maxdim() const { return RC.maxdim(); }
	inline size_t dim(const size_t k) const {return RC.dim(k);}

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
							val[k][i], std::numeric_limits<T>::infinity()
						)
					);
				} else {
					size_t j = RC.p2c[k+1][i];
					// finite bar
					pairs.emplace_back(
						PersistencePair(k, i, j,
							val[k][i], val[k+1][j]
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

		// step 2: determine update to old permutation

		// step 3: apply permutation updates to ReducedChainComplex RC

		// step 4: update reduction in RC

	}

	// get subcomplex
	ReducedFilteredChainComplex get_subcomplex() const;

};

// defualt return
template <typename FT, typename T, typename CpxT, typename... Args>
inline auto __ReducedFilteredChainComplex(const Filtration<FT, CpxT> &F, T, Args... args) {
	using VT = SparseVector<T, size_t>;
	using MT = ColumnMatrix<VT>;

	return ReducedFilteredChainComplex(FilteredChainComplex<FT, MT>(F), args...);
}

} // namespace bats
