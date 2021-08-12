#pragma once
/*
helpers for constructions in "Zigzag Zoology: Rips Zigzags for Homology Inference"
by Oudot and Sheehy (2015) [OS15]
*/

#include "data.hpp"
#include "inclusion.hpp" 				// linear_subset_union_diagram
#include "landmark.hpp"				// greedy_landmarks_hausdorff
#include <multigraph/diagram.hpp>
#include <multigraph/functors.hpp> 	// Rips functor
#include <tuple>
#include <vector>
#include <iostream>

namespace bats {

/**
determine set of indices to use for dM-ZZ construction

@param hds hausdorff distance of each point to point set.
	assume in decreasing order.
@param rho rips multiplier - should be >10 [OS15 thm 4.6]
@param theta depends on dimension of data.  In [0.5, 1/sqrt(2))[OS15 eq. 1]
	  - if unsure, set to 0.5 for finer discretization
*/
template <typename T>
std::vector<size_t> get_dM_ZZ_inds(
	const std::vector<T> &hds,
	const T rho,
	const T theta
) {
	T eta = (3 * rho + 20) / (10 * theta * rho + 20); // [OS15 eq. 12]

	// to return
	std::vector<size_t> inds2;
	std::vector<T> hds2;

	inds2.emplace_back(0); 		// emplace first point
	hds2.emplace_back(hds[0]); 	// Hausdorff distance of first point to full set

	for (size_t i = 1; i < hds.size(); i++) {
		if ( !(hds[i] > eta * hds2.back()) ) {
			inds2.emplace_back(i);
			hds2.emplace_back(hds[i]);
		}
	}

	return inds2; //make_tuple(inds2, hds2);
}


template <typename T, typename M>
std::tuple<Diagram<std::set<size_t>, std::vector<size_t>>, std::vector<T>> DiscreteMorozovZigzagSets(
	const DataSet<T> &D,
	const M &dist,
	T rho
) {
	size_t i0 = approx_center(D, dist);

	auto [ inds, hds ] = greedy_landmarks_hausdorff(D, dist, i0);

	T theta = std::sqrt(T(D.dim()) / (2 * T(D.dim()) + 1)); // [OS15 eq. 1]

	auto inds2 = get_dM_ZZ_inds(hds, rho, theta);
	size_t n = inds2.size();

	// construct sets for Rips
	std::set<size_t> s;
	std::vector<T> rmax(2*n - 1); // stores rips parameters

	Diagram<std::set<size_t>, std::vector<size_t>> SetDgm(2*n - 1, 2*n - 2);


	auto indit  = inds.cbegin();
	auto indit2 = inds2.cbegin();

	// std::cout << "constructing sets" << std::endl;
	// loop over each set size
	size_t j = 0; // node index
	while (indit2 != inds2.cend()) {
		// std::cout << j << std::endl;

		// fill in sets up to level specified by inds2[i]
		while (s.size() <= *indit2) {
			s.emplace(*indit++);
		}

		// std::cout << j << " : " << s.size() << "," << *indit2 << std::endl;

		if (j > 0) {
			SetDgm.set_node(j, s); // set node
			rmax[j] = hds[*(indit2 -1)] * rho; // rips level
			j++;
		}

		SetDgm.set_node(j, s); // set node
		rmax[j] = hds[*indit2++] * rho; // rips level
		j++;

	}

	// std::cout << "constructing edges" << std::endl;
	// set edges to be inclusion maps
	// #pragma omp parallel for
	for (size_t i = 0; i < n - 1; i++) {
		// ->
		j = 2*i; // edge_index
		size_t src = 2*i;
		size_t targ = 2*i + 1;

		SetDgm.set_edge(j, src, targ, vertex_inclusion_map(SetDgm.node[src], SetDgm.node[targ]));

		// <-
		j = 2*i+1; // update edge index
		src = 2*i + 2; // update source
		SetDgm.set_edge(j, src, targ, vertex_inclusion_map(SetDgm.node[src], SetDgm.node[targ]));
	}

	return make_tuple(SetDgm, rmax);
}

/**
Construct diagram of sets with a vector of Rips parameters

Creates this for Oscillating Rips Zigzag construction

Ref Oudot-Sheehy '15

@param D data set
@param dist metric
@param rho multiplier
@param eta multiplier.  Must have eta <= rho
*/
template <typename T, typename M>
std::tuple<Diagram<std::set<size_t>, std::vector<size_t>>, std::vector<T>> OscillatingRipsZigzagSets(
	const DataSet<T> &D,
	const M &dist,
	T rho,
	T eta
) {
	if (rho < eta) {
		throw std::runtime_error("Must have eta <= rho.");
	}

	size_t i0 = approx_center(D, dist);

	auto [ inds, hds ] = greedy_landmarks_hausdorff(D, dist, i0);

	size_t n = inds.size();

	// construct sets for Rips
	std::set<size_t> s;
	std::vector<T> rmax(2*n - 1); // stores rips parameters

	Diagram<std::set<size_t>, std::vector<size_t>> SetDgm(2*n - 1, 2*n - 2);

	auto indit  = inds.cbegin();
	auto hdit = hds.cbegin();

	// std::cout << "constructing sets" << std::endl;
	// loop over each set size
	T eps1; // prev eps
	size_t j = 0; // node index
	while (indit != inds.cend()) {
		// std::cout << *indit << ','<< j << std::endl;

		// put in next element of set
		s.emplace(*indit++);
		T eps = *hdit++; // Hausdorff distance to full set

		if (j > 0) {
			SetDgm.set_node(j, s);
			rmax[j] = rho * eps1;
			j++;
		}

		SetDgm.set_node(j, s); // set node
		rmax[j] = eta * eps;
		j++;

		eps1 = eps;

	}

	// std::cout << "constructing edges" << std::endl;
	// set edges to be inclusion maps
	// #pragma omp parallel for
	for (size_t i = 0; i < n - 1; i++) {
		// ->
		j = 2*i; // edge_index
		size_t src = 2*i;
		size_t targ = 2*i + 1;
		// std::cout << j << ":" << src << "->" << targ << std::endl;

		SetDgm.set_edge(j, src, targ, vertex_inclusion_map(SetDgm.node[src], SetDgm.node[targ]));

		// <-
		j = 2*i+1; // update edge index
		src = 2*i + 2; // update source
		// std::cout << j << ":" << src << "->" << targ << std::endl;
		SetDgm.set_edge(j, src, targ, vertex_inclusion_map(SetDgm.node[src], SetDgm.node[targ]));
	}

	return make_tuple(SetDgm, rmax);
}


// rho - rips multiplier - should be >10 [OS15 thm 4.6]
// returns diagram of Rips complexes, rips parameter for each complex.
template <typename T, typename M>
std::tuple<Diagram<SimplicialComplex, CellularMap>, std::vector<T>> DiscreteMorozovZigzag(
	const DataSet<T> &D,
	const M &dist,
	T rho,
	size_t dmax
) {

	auto [SetDgm, rmax] = DiscreteMorozovZigzagSets(D, dist, rho);

	// std::cout << "applying Rips" << std::endl;
	// now apply Rips functor
	auto RipsDgm = Rips(SetDgm, D, dist, rmax, dmax);

	// std::cout << "constructing tuple" << std::endl;
	return make_tuple(RipsDgm, rmax);
}

// TODO: add magic function that computes barcode for you

} // namespace bats
