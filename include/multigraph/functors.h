#pragma once
/*
Functors from one type of diagram to another
*/
#include "diagram.h"
#include <chain/chain_complex.h>
#include <chain/chain_map.h>
#include <homology/basis.h>
#include <homology/induced_map.h>

// ChainComplex functor
// template over matrix type, diagram type
template <typename TM, typename DT>
Diagram<ChainComplex<TM>, ChainMap<TM>> Chain(const DT &D) {
	size_t n = D.nnode();
	size_t m = D.nedge();
	// Diagram of chain complexes and chain maps
	Diagram<ChainComplex<TM>, ChainMap<TM>> CD(n, m);

	// apply chain functor to nodes
	#pragma omp parallel for
	for (size_t i = 0; i < n; i++) {
		CD.set_node(i, ChainComplex<TM>(D.node[i]));
	}

	// apply chain functor to edges
	#pragma omp parallel for
	for (size_t j = 0; j < m; j++) {
		auto s = D.elist[j].src;
	    auto t = D.elist[j].targ;
	    CD.set_edge(j, s, t, ChainMap<TM>(D.edata[j]));
	}

	return CD;
}


// Homology functor for dimension k
// template over matrix type
template <typename TM>
Diagram<ReducedChainComplex<TM>, TM> Hom(
	const Diagram<ChainComplex<TM>, ChainMap<TM>> &D,
	size_t k
) {

	size_t n = D.nnode();
	size_t m = D.nedge();
	// Diagram of chain complexes and chain maps
	Diagram<ReducedChainComplex<TM>, TM> HD(n, m);

	// apply hom functor to nodes
	#pragma omp parallel for
	for (size_t i = 0; i < n; i++) {
		HD.set_node(i, ReducedChainComplex<TM>(D.node[i]));
	}

	// apply hom functor to edges
	#pragma omp parallel for
	for (size_t j = 0; j < m; j++) {
		auto s = D.elist[j].src;
		auto t = D.elist[j].targ;
		HD.set_edge(j, s, t,
			induced_map(D.edata[j], HD.node[s], HD.node[t], k)
		);
	}

	return HD;

}
