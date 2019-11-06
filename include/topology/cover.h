#pragma once
/*
utilities for creating covers
*/

#include <cstddef>
#include <vector>
#include <set>
#include <algorithm>
#include <cmath>


// helper function for assigning lower set
template <typename T>
int assign_set_lower(const T v, const T min_val, const double bin_width, const size_t nsets) {
	double offset = (v - min_val) / bin_width;
	int bin = (int) floor(offset) * 2;
	if (!(bin < nsets)) {
		return -1;
	}
	return bin;
}

// helper function for assigning upper set
template <typename T>
int assign_set_upper(const T v, const T min_val, const double bin_width, const size_t nsets) {
	double offset = (v - min_val) / bin_width - 0.5;
	int bin = 1 + (int) floor(offset) * 2;
	if (!(bin > 0) || !(bin < nsets)) {
		return -1;
	}
	return bin;
}

/*
Create a cover of a space from projection onto the real line
INPUTS:
	x - vector of coordinates of points
	nsets - number of open sets to produce
OUTPUT:
	vector of sets, each set holds indices of x
	each subsequent set shares half the width
*/
template <typename T>
std::vector<std::set<size_t>> uniform_interval_cover(
	const std::vector<T> &x,
	const size_t nsets
) {
	auto [min, max] = std::minmax_element(x.cbegin(), x.cend());
	T min_val = *min;
	T max_val = *max;
	T length = max_val - min_val;

	double bin_width = (double) length * 2.0 / (nsets);

	// start a cover with nsets
	std::vector<std::set<size_t>> cover(nsets);

	for (size_t i = 0; i < x.size(); i++) {
		int bl = assign_set_lower(x[i], min_val, bin_width, nsets);
		int bu = assign_set_upper(x[i], min_val, bin_width, nsets);
		if (bl != -1) {
			cover[bl].emplace(i);
		}
		if (bu != -1) {
			cover[bu].emplace(i);
		}
	}
	return cover;
}
