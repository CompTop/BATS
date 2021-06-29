#pragma once
/*
utilities to compute density on point clouds
*/

#include "data.hpp"
#include "metric.hpp"
#include "neighborhood.hpp"

namespace bats {

// get distance to kth nearest neighbor
template <typename T, typename M>
std::vector<T> kdist(
	const DataSet<T> &X,
	const M& dist,
	const size_t k
) {
	std::vector<T> dk(X.size());

	auto pdist = dist(X, X);

	std::vector<T> dj(X.size());

	for (size_t j = 0; j < X.size(); j++) {
		// copy column into vector
		for (size_t i = 0; i < X.size(); i++) {
			dj[i] = pdist(i,j);
		}
		// find kth largest element
		std::nth_element(dj.begin(), dj.begin() + k, dj.end());
		dk[j] = dj[k];
	}

	return dk;
}

} // namespace bats
