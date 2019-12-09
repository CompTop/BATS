#pragma once
/*
utility for landmarking pointclouds
*/
#include <random>
#include <vector>
#include <set>
#include <algorithm>

#include "data.h"
#include "metric.h"
#include <util/set.h>

// random landmarking with k landmarks
// assume k < n
template <typename T>
DataSet<T> random_landmarks(
	const DataSet<T> &D,
	const size_t k
 ) {
	size_t n = D.size();

	auto inds = random_subset(n, k);

	return D[inds];
 }

// greedy landmarking
// D is dataset
// k is number of landmarks
// dist is Metric struct
// i0 is seed point for landmark set
// template over data type T and metric M
template <typename T, typename M>
DataSet<T> greedy_landmarks(
	const DataSet<T> &D,
	const size_t k,
	const M &dist,
	const size_t i0=0
) {
	std::set<size_t> inds;
	inds.insert(i0);

	std::vector<T> d = dist(D[i0], D);

	while (inds.size() < k) {
		// get furthest point from landmark set
		auto it = std::max_element(d.cbegin(), d.cend());
		size_t i = it - d.cbegin();

		// insert point
		inds.insert(i);

		// get all distances
		std::vector<T> di = dist(D[i], D);
		// update distances from set
		for (size_t k = 0; k < d.size(); k++) {
			d[k] = (d[k] < di[k]) ? d[k] : di[k];
		}
	}

	return D[inds];
}
