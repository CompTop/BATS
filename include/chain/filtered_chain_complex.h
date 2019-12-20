#pragma once

#include "chain_complex.h"
#include <filtration/filtration.h>

// template over filtration and matrix type
template <typename FT, typename MT>
struct FilteredChainComplex {
	std::vector<std::vector<FT>> val;
	ChainComplex<MT> C;

	FilteredChainComplex() {}

	template <typename CpxT>
	FilteredChainComplex(const Filtration<FT, CpxT> &F) : val(F.vals()), C(F.complex()) {}

	inline const ChainComplex<MT>& complex() const {return C;}
	inline const std::vector<std::vector<FT>>& vals() const { return val; }

};
