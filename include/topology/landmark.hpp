#pragma once
/*
utility for landmarking pointclouds
*/
#include <random>
#include <vector>
#include <set>
#include <algorithm>
#include <tuple>
#include <iterator>

#include "data.hpp"
#include "metric.hpp"
#include <util/set.hpp>

namespace bats {

// random landmarking with k landmarks
// assume k < n
template <typename T>
DataSet<T> random_landmarks(
	const DataSet<T> &D,
	const size_t k
 ) {
	size_t n = D.size();

	auto inds = bats::util::random_subset(n, k);

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
		size_t i = std::distance(d.cbegin(), it);

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

// hausdorff landmarking
// D is dataset
// k is number of landmarks
// dist is Metric struct
// i0 is seed point for landmark set
// template over data type T and metric M
template <typename T, typename M>
DataSet<T> hausdorff_landmarks(
	const DataSet<T> &D,
	const T r,
	const M &dist,
	const size_t i0=0
) {

	std::set<size_t> inds;
	inds.insert(i0);

	std::vector<T> d = dist(D[i0], D);
	auto it = std::max_element(d.cbegin(), d.cend());
	T dDL = *it;

	while (dDL > r) {
		// get furthest point from landmark set
		// assume we've already found max_element with it
		size_t i = std::distance(d.cbegin(), it);

		// insert point
		inds.insert(i);

		// get all distances from point i
		std::vector<T> di = dist(D[i], D);
		// update distances from set
		for (size_t k = 0; k < d.size(); k++) {
			d[k] = (d[k] < di[k]) ? d[k] : di[k];
		}

		// update distance to rest of data
		it = std::max_element(d.cbegin(), d.cend());
		dDL = *it;
	}

	return D[inds];
}


// greedy landmarking
// D is dataset
// k is number of landmarks
// dist is Metric struct
// i0 is seed point for landmark set
// template over data type T and metric M
// return indices of landmarks in order, hausdorff distance to rest of set
template <typename T, typename M>
std::tuple<std::vector<size_t>, std::vector<T>> greedy_landmarks_hausdorff(
	const DataSet<T> &D,
	const M &dist,
	const size_t i0=0
) {
	std::vector<size_t> inds; // index order
	std::vector<T>      hds; // Hausdorff distance to full data set of indices up that point
	inds.emplace_back(i0);

	std::vector<T> d = dist(D[i0], D);

	while (inds.size() < D.size()) {
		// get furthest point from landmark set
		auto it = std::max_element(d.cbegin(), d.cend());
		hds.emplace_back(*it); // Hausdorff distance of set up to that point.
		size_t i = std::distance(d.cbegin(), it);

		// insert new point
		inds.emplace_back(i);

		// get all distances
		std::vector<T> di = dist(D[i], D);
		// update distances from set
		for (size_t k = 0; k < d.size(); k++) {
			d[k] = (d[k] < di[k]) ? d[k] : di[k];
		}
	}

	hds.emplace_back(T(0));

	return std::make_tuple(inds, hds);
}

// greedy landmarking
// D is dataset
// k is number of landmarks
// dist is Metric struct
// i0 is seed point for landmark set
// template over data type T and metric M
// return indices of landmarks in order, hausdorff distance to rest of set
template <typename T>
std::tuple<std::vector<size_t>, std::vector<T>> greedy_landmarks_hausdorff(
	const Matrix<T> &pdist,
	const size_t i0=0
) {
	std::vector<size_t> inds; // index order
	std::vector<T>      hds; // Hausdorff distance to full data set of indices up that point
	inds.emplace_back(i0);

	std::vector<T> d(pdist.ncol());
	auto dit = d.begin();
	for (size_t i = 0; i < pdist.ncol(); i++) {
		*dit++ = pdist(i0,i);
	}

	while (inds.size() < pdist.ncol()) {
		// get furthest point from landmark set
		auto it = std::max_element(d.cbegin(), d.cend());
		hds.emplace_back(*it); // Hausdorff distance of set up to that point.
		size_t i = std::distance(d.cbegin(), it);

		// insert new point
		inds.emplace_back(i);

		// get all distances
		auto di = pdist[i]; // distance from point i
		// update distances from set
		for (size_t k = 0; k < d.size(); k++) {
			d[k] = (d[k] < di[k]) ? d[k] : di[k];
		}
	}

	hds.emplace_back(T(0));

	return std::make_tuple(inds, hds);
}


// find approximate center index by doing k steps of max-min landmarking
// then finding point that minimizes maximum distance to all landmarks
// if k = 0, use , k=D.dim()+1 for triangulation.
// i0 is first landmark point (default 0)
template <typename T, typename M>
size_t approx_center(
	const DataSet<T> &D,
	const M &dist,
	size_t k=0,
	size_t i0=0
) {
	// if k = 0, set k = D.dim + 1;
	k = (k== 0) ? D.dim() + 1 : k;

	std::vector<T> maxd = dist(D[i0], D); // minimum maximum distance to any landmark
	std::vector<T> mind = dist(D[i0], D); // maximum minimum distance to any landmark

	for (size_t j = 0; j < k; j++) {

		// find next landmark
		auto it = std::max_element(maxd.cbegin(), maxd.cend());
		size_t i = std::distance(maxd.cbegin(), it);

		// distances to landmark i
		std::vector<T> di = dist(D[i], D);

		// update max and min distances
		for (size_t jj = 0; jj < di.size(); jj++) {
			maxd[jj] = (maxd[jj] < di[jj]) ? maxd[jj] : di[jj];
			mind[jj] = (mind[jj] > di[jj]) ? mind[jj] : di[jj];
		}

	}

	return std::distance(
		mind.cbegin(),
		std::min_element(mind.cbegin(), mind.cend())
	);
}

} // namespace bats
