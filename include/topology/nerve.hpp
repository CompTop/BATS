#pragma once
/*
Nerve of a cover
*/

#include <vector>
#include <set>
#include <util/sorted.hpp>
#include <complex/simplicial_complex.hpp>
#include "cover.hpp"

namespace bats {

void add_dimension_recursive_nerve(
	SimplicialComplex &N,
	std::vector<size_t> spx,
	const bats::Cover &cover,
	std::vector<size_t> intersection,
	const size_t dmax
) {
	size_t dim = spx.size() - 1;
	if (dim >= dmax) { return; }

	// finds maximum k
	size_t kmax = *std::min_element(spx.begin(), spx.end());
	std::vector<size_t> intersectionk;
	std::vector<size_t> spxk;

	for (size_t k = 0; k < kmax; k++) {
		bats::util::intersect_sorted(cover[k], intersection, intersectionk);
		if (intersectionk.size() > 0) {
			// add simplex
			spx.emplace_back(k);
			bats::util::sort_into(spx, spxk);
			// add to nerve
			N.add(spxk);
			if (dmax > dim+1) {
				add_dimension_recursive_nerve(N, spxk, cover, intersectionk, dmax);
			}
			spx.pop_back(); // remove k
		}
	}

}

/*
Build nerve of cover up to dmax skeleton
*/
SimplicialComplex Nerve(
	const bats::Cover &cover,
	const size_t dmax
) {

	// initialize complex
    SimplicialComplex N(dmax);

	std::vector<size_t> intersection;
	std::vector<size_t> spx;

	for (size_t i = 0; i < cover.size(); i++) {
		// add 0-simplex i
		N.add({i});

		for (size_t j = 0; j < i; j++) {
			// check to add edge between i and j
			bats::util::intersect_sorted(cover[i], cover[j], intersection);
			if (intersection.size() > 0) {
				// add edge
				spx = {j, i};
				N.add(spx);
				// add higher order simplices recursively, but only for elements < j
				if (dmax > 1) {
					add_dimension_recursive_nerve(N, spx, cover, intersection, dmax);
				}
			}
		}
	}
	return N;
}

} // namespace bats
