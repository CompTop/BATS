#pragma once

#include <vector>
#include <limits>
#include <dgvs/filtered_dgvs.hpp>
#include <homology/dgbasis.hpp>
#include <util/common.hpp>
#include "barcode.hpp"
#include <filtration/update_information.hpp>

namespace bats {

template <typename T, typename MT>
struct ReducedFilteredDGVectorSpace {

    ReducedDGVectorSpace<MT> RC; // reduced in permutation order
	std::vector<std::vector<T>> val; // stored in original order
	std::vector<std::vector<size_t>> perm; //from permutaiton order to original order

	ReducedFilteredDGVectorSpace() {}

	// variadic template for passing arguments
	template <typename... Args>
	ReducedFilteredDGVectorSpace(const FilteredDGVectorSpace<T, MT>& C, Args ...args) :
		RC(C.complex(), args...),
		val(C.val),
		perm(C.perm) {}

	inline size_t maxdim() const { return RC.maxdim(); }
	inline size_t dim(const size_t k) const {return RC.dim(k);}
	inline size_t hdim(const size_t k) const {return RC.hdim(k);}

	/**
	persistence pairs in dimension k

	@param k	homology dimension
	@param permuted	set to true to return critical indices permuted by filtration parameter
	set to false to return with indices in original order.  Default: false
	*/
	std::vector<PersistencePair<T>> persistence_pairs(
		const size_t k
	)const {
		std::vector<PersistencePair<T>> pairs;
		if (RC.degree == -1) {
			pairs.reserve(RC.R[k].ncol());
			for (size_t i =0; i < dim(k); i++) {
				if (RC.R[k][i].nnz() == 0) {
					// homology generated
					if (k == maxdim() || RC.p2c[k+1][i] == bats::NO_IND)  {
						// infinite bar
						pairs.emplace_back(
							PersistencePair<T>(k, perm[k][i], bats::NO_IND,
								val[k][perm[k][i]], std::numeric_limits<T>::infinity()
							)
						);
					} else {
						size_t j = RC.p2c[k+1][i];
						// finite bar
						pairs.emplace_back(
							PersistencePair<T>(k, perm[k][i], perm[k+1][j],
								val[k][perm[k][i]], val[k+1][perm[k+1][j]]
							)
						);
					}
				}
			}
		} else { // degree == +1
			pairs.reserve(RC.R[k+1].ncol());
			for (size_t i =0; i < dim(k); i++) {
				if (RC.R[k+1][i].nnz() == 0) {
					// homology generated
					if (RC.p2c[k][i] == bats::NO_IND)  {
						// infinite bar
						pairs.emplace_back(
							PersistencePair<T>(k, perm[k+1][i], bats::NO_IND,
								val[k+1][perm[k+1][i]], std::numeric_limits<T>::infinity()
							)
						);
					} else {
						size_t j = RC.p2c[k][i];
						// finite bar
						pairs.emplace_back(
							PersistencePair<T>(k, perm[k+1][i], perm[k][j],
								val[k+1][perm[k+1][i]], val[k][perm[k][j]]
							)
						);
					}
				}
			}
		}

		return pairs;
	}

};

} // namespace bats
