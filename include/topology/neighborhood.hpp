#pragma once
/*
points for neighborhoods
*/
#include <vector>
#include <util/permutation.hpp>
#include "data.hpp"

namespace bats {

/*
get neighbors of point x in DataSet X within radius r in metric dist
*/
template <typename T, typename M>
std::vector<size_t> neighborhood(
	const VectorView<T> &x,
	const DataSet<T> &X,
	const M &dist,
	const T r
) {

	std::vector<size_t> ind;

	for (size_t i = 0; i < X.size(); i++) {
		if (dist(x, X[i]) < r) {
			ind.emplace_back(i);
		}
	}

	return ind;
}

/*
get k nearest neighbors of point x in DataSet X
*/
template<typename T, typename M>
std::vector<size_t> neighborhood(
	const VectorView<T> &x,
	const DataSet<T> &X,
	const M &dist,
	const size_t k
) {
	auto dxX = dist(x, X);
	return bats::util::firstk(dxX.cbegin(), dxX.cend(), k);
}

/*
get neighborhoods of points in X in DataSet Y within radius r in metric dist
*/
template <typename T, typename M>
std::vector<std::vector<size_t>> neighborhoods(
	const DataSet<T> &X,
	const DataSet<T> &Y,
	const M &dist,
	const T r
) {

	std::vector<std::vector<size_t>> nbrs(X.size());

	#pragma omp parallel for
	for (size_t i = 0; i < X.size(); i++) {
		nbrs.emplace_back(neighborhood(X[i], Y, dist, r));
	}

	return nbrs;
}

/*
get k nearest neighborhoods of points in X in DataSet Y
*/
template <typename T, typename M>
std::vector<std::vector<size_t>> neighborhoods(
	const DataSet<T> &X,
	const DataSet<T> &Y,
	const M &dist,
	const size_t k
) {

	std::vector<std::vector<size_t>> nbrs(X.size());

	#pragma omp parallel for
	for (size_t i = 0; i < X.size(); i++) {
		// k nearest neighbors of X[i] in L
		nbrs[i] = neighborhood(X[i], Y, dist, k);
	}

	return nbrs;
}


/*
put columns of pdist in neighborhoods of rows of pdist
each column is a point in X
each row is a point in L
assign each column to neighborhood of L if the distance to L is <= eps
*/
template <typename T>
std::vector<std::set<size_t>> neighborhoods(
    const Matrix<T> &pdist,
    const T eps = T(0)
) {

	std::vector<std::set<size_t>> nbrs(pdist.nrow());

	for (size_t j = 0; j < pdist.ncol(); j++) {
		for (size_t i = 0; i < pdist.nrow(); i++) {
			if (pdist(i,j) <= eps) {
				nbrs[i].emplace_hint(nbrs[i].end(), j);
			}
		}
	}
	return nbrs;
}

} // namespace bats
