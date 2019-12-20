#pragma once

#include <vector>
#include <limits>
#include <chain/filtered_chain_complex.h>
#include <homology/basis.h>
#include <util/common.h>
#include "barcode.h"


template <typename T, typename MT>
struct ReducedFilteredChainComplex {

	std::vector<std::vector<T>> val;
	ReducedChainComplex<MT> RC;

	ReducedFilteredChainComplex() {}

	ReducedFilteredChainComplex(const FilteredChainComplex<T, MT>& C) :
		val(C.vals()),
		RC(C.complex()) {}

	inline size_t maxdim() const { return RC.maxdim(); }
	inline size_t dim(const size_t k) const {return RC.dim[k];}

	// persistence pairs in dimension k
	std::vector<PersistencePair<T>> persistence_pairs(const size_t k) {
		std::vector<PersistencePair<T>> pairs;
		for (size_t i =0; i < dim(k); i++) {
			if (RC.R[k][i].nnz() == 0) {
				// homology generated
				if (k == maxdim() || RC.p2c[k+1].count(i) == 0)  {
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

	// barcode w/out critical inds
	std::vector<T> barcode(const size_t k);

	// critical cells for barcode in dimension k
	std::vector<size_t> critical_cells(const size_t k);

};
