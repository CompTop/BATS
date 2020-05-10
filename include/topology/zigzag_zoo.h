#pragma once
/*
helpers for constructions in "Zigzag Zoology: Rips Zigzags for Homology Inference"
by Oudot and Sheehy (2015) [OS15]
*/

#include "data.h"
#include "inclusion.h" 				// linear_subset_union_diagram
#include "landmark.h"				// greedy_landmarks_hausdorff
#include <multigraph/diagram.h>
#include <multigraph/functors.h> 	// Rips functor
#include <tuple>
#include <vector>
#include <iostream>

/*
determine set of indices to use for dM-ZZ construction
hds - hausdorff distance of each point to point set.
	assume in decreasing order.
rho - rips multiplier - should be >10 [OS15 thm 4.6]
theta - depends on dimension of data.  In [0.5, 1/sqrt(2))[OS15 eq. 1]
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


// rho - rips multiplier - should be >10 [OS15 thm 4.6]
// returns diagram of Rips complexes, rips parameter for each complex.
template <typename T, typename M>
std::tuple<Diagram<SimplicialComplex, CellularMap>, std::vector<T>> DiscreteMorozovZigzag(
	const DataSet<T> &D,
	const M &dist,
	T rho,
	size_t dmax
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
		size_t j = 2*i+1; // update edge index
		src = 2*i + 2; // update source
		SetDgm.set_edge(j, src, targ, vertex_inclusion_map(SetDgm.node[src], SetDgm.node[targ]));
	}

	// std::cout << "applying Rips" << std::endl;
	// now apply Rips functor
	auto RipsDgm = Rips(SetDgm, D, dist, rmax, dmax);

	// std::cout << "constructing tuple" << std::endl;
	return make_tuple(RipsDgm, rmax);
}

// TODO: add magic function that computes barcode for you
