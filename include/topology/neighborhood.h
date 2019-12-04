#pragma once
/*
points for neighborhoods
*/
#include <vector>
#include <util/permutation.h>
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
	return bats::firstk(dxX.cbegin(), dxX.cend(), k);
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

/*
get k nearest neighborhoods of points in L in DataSet X
*/
template <typename T, typename M>
std::vector<std::vector<size_t>> neighborhoods(
	const DataSet<T> &L,
	const DataSet<T> &X,
	const M &dist,
	const size_t k
) {

	std::vector<std::vector<size_t>> nbrs;
	nbrs.reserve(L.size());

	for (size_t i = 0; i < L.size(); i++) {
		nbrs.emplace_back(neighborhood(L[i], X, dist, k));
	}

	return nbrs;
}
