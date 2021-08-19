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
#include <dgvs/dgvs.hpp>
#include <dgvs/dgmap.hpp>
#include <homology/dgbasis.hpp>

#include <stdexcept>

namespace bats {

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
	    TD.set_node(i, RipsComplex<SimplicialComplex>(XI, dist, rmax, dmax));
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
	    TD.set_node(i, RipsComplex<SimplicialComplex>(XI, dist, rmax[i], dmax));
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
Diagram<ChainComplex<TM>, ChainMap<TM>> ChainFunctor(const DT &D) {
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

template <typename TF, typename DT>
inline auto ChainFunctor(const DT &D, TF) {
	using TM = ColumnMatrix<SparseVector<TF>>;
	return ChainFunctor<TM>(D);
}

// easy chain functor
template <typename DT, typename T>
inline auto __ChainFunctor(const DT &D, T) {
	using VT = SparseVector<T, size_t>;
	using MT = ColumnMatrix<VT>;

	return ChainFunctor<MT, DT>(D);
}

template <typename CpxT, typename T>
inline auto Chain(const Diagram<CpxT, CellularMap>& D, T) {
	using VT = SparseVector<T, size_t>;
	using MT = ColumnMatrix<VT>;

	return ChainFunctor<MT, Diagram<CpxT, CellularMap>>(D);
}


/**
Homology functor for dimension k
template over matrix type
*/
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

/**
Homology functor in all dimensions
template over matrix type

@param D diagram of ChainComplexes and ChainMaps
@param topd (optional, default: false) if true will compute top dimensional homology.

when topd is true, a k-dimensional Chain complex will be assumed to be 0 in dimension k+1
and H_k will be computed.

Assumes that all chain complexes have same dimension.
*/
template <typename TM>
Diagram<ReducedChainComplex<TM>, std::vector<TM>> Hom(
	const Diagram<ChainComplex<TM>, ChainMap<TM>> &D,
	bool topd=false
) {

	size_t n = D.nnode();
	size_t m = D.nedge();
	// Diagram of chain complexes and chain maps
	Diagram<ReducedChainComplex<TM>, std::vector<TM>> HD(n, m);
	size_t maxdim = D.node_data(0).maxdim();
	maxdim = topd ? maxdim + 1 : maxdim; // compute 1-dimension higher if topd is true
	if (maxdim == 0) { throw std::runtime_error("No homology to compute!");}

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
		std::vector<TM> F(maxdim);
		for (size_t k = 0; k < maxdim; ++k) {
			F[k] = induced_map(D.edata[j], HD.node[s], HD.node[t], k);
		}
		HD.set_edge(j, s, t, F);
	}

	return HD;
}

// Create diagram of Rips complexes from subsets
// uses a different rips parameter for each node
template <typename T, typename M, typename FT>
auto RipsHom(
	const Diagram<std::set<size_t>, std::vector<size_t>> &D,
	const DataSet<T> &X,
	const M &dist, // distance
	const std::vector<T> &rmax, // maximum radius for each node
	const size_t hdim, // homology dimension
	FT // field
) {

	using VT = SparseVector<FT, size_t>;
	using MT = ColumnMatrix<VT>;
	using CpxT = SimplicialComplex;

	size_t n = D.nnode();
	size_t m = D.nedge();
	// Diagram of dimensions and induced maps
	Diagram<size_t, MT> TD(n, m);

	// apply functor to edges
	// note this will do every node computation twice
	// but the parallelization is trivial
	#pragma omp parallel for
	for (size_t i = 0; i < m; i++) {
	    auto s = D.elist[i].src;
	    auto t = D.elist[i].targ;
		if (rmax[t] < rmax[s]) {
			throw std::range_error("Rips parameter must be non-decreasing from source to target.");
		}
		// handle nodes
		auto Is = get_subset(X, D.node[s]);
		auto It = get_subset(X, D.node[t]);
		auto Xs = RipsComplex<CpxT>(Is, dist, rmax[s], hdim+1);
		auto Xt = RipsComplex<CpxT>(It, dist, rmax[t], hdim+1);
		// skip chain construction in memory
		auto Rs = Reduce(Xs, FT());
		auto Rt = Reduce(Xt, FT());

		// handle edge
		auto F = SimplicialMap(Xs , Xt, D.edata[i]);
		auto CF = ChainMap<MT>(F);
		auto HF = induced_map(CF, Rs, Rt, hdim);

	    TD.set_edge(i, s, t, HF);
	}

	// extract dimensions for nodes
	// we do this by iterating over edges
	// do this sequentially to avoid race conditions
	for (size_t i = 0; i < m; i++) {
		auto s = D.elist[i].src;
	    auto t = D.elist[i].targ;
		auto ds = TD.edata[i].ncol();
		auto dt = TD.edata[i].nrow();
		TD.set_node(s, ds);
		TD.set_node(t, dt);
	}

	return TD;
}


/**
Functor from topological category to category of
differential graded vector spaces

Chain functor is degree = -1 (default)
Cochain functor is degree = +1.  This is contravariant (reverses arrows)
*/
template <typename TM, typename DT>
Diagram<DGVectorSpace<TM>, DGLinearMap<TM>> DGLinearFunctor(
	const DT &D,
	int degree=-1
) {
	size_t n = D.nnode();
	size_t m = D.nedge();
	// Diagram of chain complexes and chain maps
	Diagram<DGVectorSpace<TM>, DGLinearMap<TM>> DGD(n, m);

	// apply chain functor to nodes
	#pragma omp parallel for
	for (size_t i = 0; i < n; i++) {
		DGD.set_node(i, DGVectorSpace<TM>(D.node[i], degree));
	}

	// apply chain functor to edges
	#pragma omp parallel for
	for (size_t j = 0; j < m; j++) {
		auto s = (degree == -1) ? D.elist[j].src : D.elist[j].targ;
	    auto t = (degree == -1) ? D.elist[j].targ : D.elist[j].src;
	    DGD.set_edge(j, s, t, DGLinearMap<TM>(D.edata[j], degree));
	}

	return DGD;
}

/**
Homology functor for dimension k
template over matrix type
*/
template <typename TM>
Diagram<ReducedDGVectorSpace<TM>, TM> Hom(
	const Diagram<DGVectorSpace<TM>, DGLinearMap<TM>> &D,
	size_t k
) {

	size_t n = D.nnode();
	size_t m = D.nedge();
	// Diagram of chain complexes and chain maps
	Diagram<ReducedDGVectorSpace<TM>, TM> HD(n, m);

	// apply hom functor to nodes
	#pragma omp parallel for
	for (size_t i = 0; i < n; i++) {
		HD.set_node(i, ReducedDGVectorSpace<TM>(D.node[i]));
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

} // namespace bats
