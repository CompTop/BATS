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
						// infinite bar.
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
			/*
			TODO: this is actually *relative* relative barcode
			Should maybe have persistence_pairs_relative

			See "Dualitites in persistent (co)homology" DeSilva et al 2011
			Section 3.4
			*/
			for (size_t i = 0; i < dim(k); ++i) {
				if (RC.R[k+1][i].nnz() > 0) {
					// kills bar
					size_t j = RC.R[k+1][i].lastnz().ind;
					pairs.emplace_back(
						PersistencePair<T>(
							k,
							perm[k][i], // birth in filtration is cohomology death
							perm[k+1][j], // death in filtration is cohomology birth
							val[k][perm[k][i]], // birth in filtration is cohomology death
							val[k+1][perm[k+1][j]] // death in filtration is cohomology birth
						)
					);
				} else {
					// check to see if infinte bar
					size_t j = RC.p2c[k][i];
					if (j == bats::NO_IND) {
						pairs.emplace_back(
							PersistencePair<T>(
								k,
								perm[k][i],
								bats::NO_IND,
								val[k][perm[k][i]],
								std::numeric_limits<T>::infinity()
							)
						);
					}
				}
			}
			// pairs.reserve(RC.R[k+2].ncol());
			// for (size_t i =0; i < dim(k+1); i++) {
			// 	if (RC.R[k+2][i].nnz() == 0) {
			// 		// homology generated
			// 		size_t j = RC.p2c[k+1][i];
			// 		if (j == bats::NO_IND)  {
			// 			// infinite bar.  Born at column index RC.R[k+1][i]
			// 			pairs.emplace_back(
			// 				PersistencePair<T>(
			// 					k,
			// 					perm[k+1][i],
			// 					bats::NO_IND,
			// 					val[k+1][perm[k+1][i]],
			// 					std::numeric_limits<T>::infinity()
			// 				)
			// 			);
			// 		} else {
			// 			// finite bar.  Dies at column index RC[k][j]
			// 			pairs.emplace_back(
			// 				PersistencePair<T>(
			// 					k,
			// 					perm[k][j], // birth in filtration is cohomology death
			// 					perm[k+1][i], // death in filtration is cohomology birth
			// 					val[k][perm[k][j]], // birth in filtration is cohomology death
			// 					val[k+1][perm[k+1][i]] // death in filtration is cohomology birth
			// 				)
			// 			);
			// 		}
			// 	}
			// }
		}

		return pairs;
	}

	// update filtration
	void update_filtration(const std::vector<std::vector<T>> newval) {
		// step 1: determine permutation order for newval
		auto perms = filtration_sortperm(newval);

		// if degree is +1, reverse
		if (RC.degree == +1) {
			// for cohomology reverse rows and columns
			for (auto& pi : perms) {
				std::reverse(pi.begin(), pi.end());
			}
		}

		// get inverse permutation for current permutation
		auto iperm = filtration_iperm(perm);

		// step 2: determine update to old permutation
		// iperm[k] will hold the updated permutation temporarily
		std::vector<size_t> tmp; // temporary vector
		for (size_t k = 0; k < iperm.size(); k++) {
			bats::util::apply_perm(iperm[k].data(), tmp, perms[k]);
		}

		// step 3: apply permutation updates to ReducedChainComplex RC
		RC.permute_basis(iperm);

		// step 4: store new permutation
		perm = perms;

		// store values as well
		val = newval;

	}

};

} // namespace bats
