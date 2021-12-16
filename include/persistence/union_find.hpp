#include <iostream>
#include <vector>
#include <utility>
#include <limits>
#include <algorithm>
#include <map>
#include <cstddef>
#include "barcode.hpp"

namespace bats {

/**
Find the parent of node i.
Performs path compression.
*/
int find_parent(std::vector<size_t> &parent, size_t i) {
	if (i != parent[i]) {
		parent[i] = find_parent(parent, parent[i]);
	}
	return parent[i];
}


/**

Compute the 0-dimensional barcode using the Union-Find algorithm

returns vector of PersistencePairs containing 0-dimensional barcode
*/
template <typename T, typename MT>
std::vector<PersistencePair<T>> union_find_pairs(
	const FilteredChainComplex<T, MT> &F
) {

	// first get boundary in dimension 1
	// assumed to be sorted in filtration order
	auto& B = F.complex()[1];
	// vals are not stored in filtraiton order
	const auto& val = F.val;
	// perm is used for mapping
	const auto& perm = F.perm;

	size_t N = B.nrow();
	size_t M = B.ncol();



	// initialize parent vector
	std::vector<size_t> parent(N);
	std::iota(parent.begin(), parent.end(), 0);

	std::vector<PersistencePair<T>> pairs;
	pairs.reserve(N);
	//first fill birth times, death times
	for (size_t k = 0; k < N; ++k) {
		pairs.emplace_back(PersistencePair(
			0,
			k,
			bats::NO_IND,
			val[0][perm[0][k]],
			std::numeric_limits<T>::infinity()
		));
	}

	// loop over columns of boundary matrix
	size_t nfinite = 0; // number of finite bars
	for (size_t k = 0; k < M; ++k) {
		// first, get simplices on boundary

		// This is not safe if we do not know we're working with 1-dimensional boundary
		auto nzit = B[k].nzbegin();
		size_t i = nzit->ind; ++nzit;
		size_t j = nzit->ind;

		size_t pi = find_parent(parent, i);
		size_t pj = find_parent(parent, j);

		// if parents are same, H1 created.  Nothing to do here.
		if (pi == pj) { continue; }

		// get birth times in filtration order
		// float bi = pairs[pi].birth;
		// float bj = pairs[pj].birth;
		if (pi < pj) { // bi < bj || (bi == bj && pi < pj)
			// should be ok because boundary is in filtration order
			// component with j merges into component with i
			// kill bar at pj
			pairs[pj].death_ind = k;
			pairs[pj].death = val[1][perm[1][k]];

			// merge components
			parent[pj] = pi;
			parent[j] = pi;
		} else {
			// component with i merges into component with j
			// kill bar at pi
			pairs[pi].death_ind = k;
			pairs[pi].death = val[1][perm[1][k]];

			// merge components
			parent[pi] = pj;
			parent[i] = pj;
		}
		nfinite++;
		// if we've found spanning tree, then break
		if (nfinite == N - 1) { break; }




	}

	return pairs;
}

} // namespace bats
