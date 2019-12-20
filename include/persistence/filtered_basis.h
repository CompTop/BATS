#pragma once

#include <vector>
#include <chain/filtered_chain_complex.h>
#include <homology/basis.h>
#include "barcode.h"


template <typename T, typename MT>
struct ReducedFilteredChainComplex {

	std::vector<std::vector<T>> val;
	ReducedChainComplex<MT> RC;

	ReducedFilteredChainComplex() {}

	ReducedFilteredChainComplex(const FilteredChainComplex<T, MT>& C) :
		val(C.vals()),
		RC(C.complex()) {}

	// persistence pairs in dimension k
	std::vector<PersistencePair<T>> persistence_pairs(const size_t k);

	// barcode w/out critical inds
	std::vector<T> barcode(const size_t k);

	// critical cells for barcode in dimension k
	std::vector<size_t> critical_cells(const size_t k);

};
