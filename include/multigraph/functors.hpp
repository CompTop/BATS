#pragma once
/*
Functors from one type of diagram to another
*/
#include "diagram.hpp"
#include <topology/data.hpp>
#include <topology/cover.hpp>
#include <topology/rips.hpp>
#include <topology/nerve.hpp>
#include <complex/simplicial_complex.hpp>
#include <complex/simplicial_map.hpp>
#include <chain/chain_complex.hpp>
#include <chain/chain_map.hpp>
#include <homology/basis.hpp>
#include <homology/induced_map.hpp>

#include <stdexcept>

Diagram<SimplicialComplex, CellularMap> Nerve(
	const Diagram<bats::Cover, std::vector<size_t>> &D,
	const size_t dmax // maximum simplex dimension
) {
	size_t n = D.nnode();
	size_t m = D.nedge();
	// Diagram of simplicial complexes and maps
	Diagram<SimplicialComplex, CellularMap> TD(n, m);

	// apply functor to nodes
	#pragma omp parallel for
	for (size_t i = 0; i < n; i++) {
		TD.set_node(i, Nerve(D.node[i], dmax));
	}

	// apply functor to edges
	#pragma omp parallel for
	for (size_t i = 0; i < m; i++) {
		auto s = D.elist[i].src;
		auto t = D.elist[i].targ;
		TD.set_edge(i, s, t,
						SimplicialMap(TD.node[s] , TD.node[t], D.edata[i])
					   );
	}

	return TD;


}


// Create diagram of Rips complexes from subsets
template <typename T, typename M>
Diagram<SimplicialComplex, CellularMap> Rips(
	const Diagram<std::set<size_t>, std::vector<size_t>> &D,
	const DataSet<T> &X,
	const M &dist, // distance
	const T rmax, // maximum radius
	const size_t dmax // maximum simplex dimension
) {

	size_t n = D.nnode();
	size_t m = D.nedge();
	// Diagram of simplicial complexes and maps
	Diagram<SimplicialComplex, CellularMap> TD(n, m);

	// apply functor to nodes
	#pragma omp parallel for
	for (size_t i = 0; i < n; i++) {
	    auto XI = get_subset(X, D.node[i]);
	    TD.set_node(i, RipsComplex(XI, dist, rmax, dmax));
	}

	// apply functor to edges
	#pragma omp parallel for
	for (size_t i = 0; i < m; i++) {
	    auto s = D.elist[i].src;
	    auto t = D.elist[i].targ;
	    TD.set_edge(i, s, t,
	                    SimplicialMap(TD.node[s] , TD.node[t], D.edata[i])
	                   );
	}

	return TD;
}


// Create diagram of Rips complexes from subsets
// uses a different rips parameter for each node
template <typename T, typename M>
Diagram<SimplicialComplex, CellularMap> Rips(
	const Diagram<std::set<size_t>, std::vector<size_t>> &D,
	const DataSet<T> &X,
	const M &dist, // distance
	const std::vector<T> &rmax, // maximum radius for each node
	const size_t dmax // maximum simplex dimension
) {

	size_t n = D.nnode();
	size_t m = D.nedge();
	// Diagram of simplicial complexes and maps
	Diagram<SimplicialComplex, CellularMap> TD(n, m);

	// apply functor to nodes
	#pragma omp parallel for
	for (size_t i = 0; i < n; i++) {
	    auto XI = get_subset(X, D.node[i]);
	    TD.set_node(i, RipsComplex(XI, dist, rmax[i], dmax));
	}

	// apply functor to edges
	#pragma omp parallel for
	for (size_t i = 0; i < m; i++) {
	    auto s = D.elist[i].src;
	    auto t = D.elist[i].targ;
		if (rmax[t] < rmax[s]) {
			throw std::range_error("Rips parameter must be non-decreasing from source to target.");
		}
	    TD.set_edge(i, s, t,
	                    SimplicialMap(TD.node[s] , TD.node[t], D.edata[i])
	                   );
	}

	return TD;
}


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

// easy chain functor
template <typename DT, typename T>
inline auto __Chain(const DT &D, T) {
	using VT = SparseVector<T, size_t>;
	using MT = ColumnMatrix<VT>;

	return Chain<MT, DT>(D);
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
