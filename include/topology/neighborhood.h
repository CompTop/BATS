#pragma once
/*
points for neighborhoods
*/
#include <vector>
#include "data.h"

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
get neighborhoods of points in L in DataSet X within radius r in metric dist
*/
template <typename T, typename M>
std::vector<std::vector<size_t>> neighborhoods(
	const DataSet<T> &L,
	const DataSet<T> &X,
	const M &dist,
	const T r
) {

	std::vector<std::vector<size_t>> nbrs;
	nbrs.reserve(L.size());

	for (size_t i = 0; i < L.size(); i++) {
		nbrs.emplace_back(neighborhood(L[i], X, dist, r));
	}

	return nbrs;
}
