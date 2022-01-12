#include <iostream>
#include <vector>
#include <utility>
#include <limits>
#include <algorithm>
#include <map>
#include <cstddef>
#include <math.h>
#include "barcode.hpp"

namespace bats {


/**
Find parent of node i
Performs path compression
Implementation that avoids recursive function calls

Modified from https://github.com/stat37411/tda/blob/main/include/union_find.hpp
*/
size_t find_parent(
    std::vector<size_t> &parent,
    size_t i
) {
    size_t root = i;

    // first find cluster representative, i.e. root of tree
    while (parent[root] != root) {
        root = parent[root];
    }

    // now do path compression
    while (parent[i] != root) {
        auto p = parent[i];
        parent[i] = root;
        i = p;
    }

    return root;
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

/**
compute k-choose-2 = k * (k-1) / 2
*/
template <typename T>
inline T k_choose_2(T k) { return (k * (k-1)) >> 1; }

/**
returns next smallest k such that k-choose-2 is >= targ
*/
inline size_t k_choose_2_inv(size_t targ){
	return (size_t) ceil((1 + (size_t) sqrt((double)(1 + 8*targ))) / 2);
}

/**
compute the squareform indices of binomial index k
inputs:
	k, index between 0 and N-choose-2
	N size of set
returns:
	(i,j) where i < j

 k = N-choose-2 - (N-i)-choose-2 + (j-i-1)
 (N-i)-choose-2 = N-choose-2 - k - (j-i-1)
 (N-i)-choose-2 >= N-choose-2 - k
*/
template <typename T>
inline std::pair<T, T> binom_to_inds(
	T k,
	T N
) {
	T N2 = k_choose_2(N);
	T targ = N2 - k;
	T ub = k_choose_2_inv(targ);
	T ub2 = k_choose_2(ub);
	T i = N - ub;
	T j = (k + ub2 + i + 1) - N2;
	return std::make_pair(i,j);
}
/**

Compute the 0-dimensional Rips barcode using the Union-Find algorithm

Inputs:
inds: sort indices of n-choose-2 pairwise distances
vals: values of n-choose-2 pairwise distances (sorted by inds)

outputs:
0-dimensional persistence pairs

returns vector of PersistencePairs containing 0-dimensional barcode
*/
template <typename T>
std::vector<PersistencePair<T>> rips_union_find_pairs(
	const std::vector<size_t> &inds,
	const std::vector<T>& vals
) {

	size_t M = inds.size(); // number of edges
	// calculate number of nodes.
	// N-choose-2 = M = N (N - 1) // 2
	// N^2 - N - 2M = 0
	// N = [1 + sqrt(1 + 8*M)] // 2 by quadratic formula
	size_t N = (1 + (size_t) sqrt((double)(1 + 8*M))) / 2;

	// initialize output
	std::vector<PersistencePair<T>> pairs;
	pairs.reserve(N);

	//first fill birth times, death times
	for (size_t k = 0; k < N; ++k) {
		pairs.emplace_back(PersistencePair(
			0,
			k,
			bats::NO_IND,
			0.0, // standard birth at 0.0
			std::numeric_limits<T>::infinity()
		));
	}

	// initialize parent vector
	std::vector<size_t> parent(N);
	std::iota(parent.begin(), parent.end(), 0);

	// loop over columns of boundary matrix
	size_t nfinite = 0; // number of finite bars
	for (size_t k = 0; k < inds.size(); ++k) {

		// calulate i, j from edge index ei
		auto ei = inds[k];
		auto [i, j] = binom_to_inds(ei, N);

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
			pairs[pj].death = vals[ei];

			// merge components
			parent[pj] = pi;
			parent[j] = pi;
		} else {
			// component with i merges into component with j
			// kill bar at pi
			pairs[pi].death_ind = k;
			pairs[pi].death = vals[ei];

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
