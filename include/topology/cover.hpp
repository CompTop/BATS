#pragma once
/*
utilities for creating covers
*/

#include <cstddef>
#include <vector>
#include <set>
#include <algorithm>
#include <cmath>
#include <tuple>

#include <util/set.hpp>
#include <multigraph/diagram.hpp>
#include "data.hpp"
#include "inclusion.hpp"
#include "landmark.hpp"
#include "neighborhood.hpp"

namespace bats {

using Cover = std::vector<std::set<size_t>>;




// helper function for assigning lower set
template <typename T>
int assign_set_lower(const T v, const T min_val, const double bin_width, const size_t nsets) {
	double offset = (v - min_val) / bin_width;
	size_t bin = (size_t) floor(offset) * 2;
	if (!(bin < nsets)) {
		return -1;
	}
	return bin;
}

// helper function for assigning upper set
template <typename T>
int assign_set_upper(const T v, const T min_val, const double bin_width, const size_t nsets) {
	double offset = (v - min_val) / bin_width - 0.5;
	size_t bin = 1 + (size_t) floor(offset) * 2;
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
bats::Cover uniform_interval_cover(
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




// project data x in d dimensions onto coordinate i
template <typename T>
std::vector<T> coordinate_projection(const DataSet<T> &X, const size_t i) {
    std::vector<T> p;
    p.reserve(X.size());
    for (size_t j = 0; j < X.size(); j++) {
		// std::cout << j << ": " << X(j, i) << std::endl;
		p.emplace_back(X(j, i));
    }
    return p;
}

// get subset of data
template <typename T>
DataSet<T> get_subset(
    const DataSet<T> &X,
    const std::set<size_t> &ind
) {
	size_t d = X.dim();
    Matrix<T> XS(ind.size(), d);
	size_t i = 0;
    for (auto it = ind.cbegin(); it != ind.cend(); it++) {
        for (size_t j = 0; j < d; j++) {
            XS(i, j) = X(*it, j);
        }
		i++;
    }
    return DataSet(XS);
}

/*
Assign points in data sets to k nearest landmarks
*/
// template over data type, distance
template <typename T, typename M>
bats::Cover landmark_cover(
	const DataSet<T> &X,
	const DataSet<T> &L,
	const M &dist,
	size_t k
) {
	auto ns = neighborhoods(X, L, dist, k);
	// ns[i] contains k closest points in L to X[i]
	bats::Cover cover(L.size());
	for (size_t i = 0; i < ns.size(); i++) {
		for (auto j : ns[i]) {
			cover[j].emplace(i);
		}
	}
	return cover;
}

// get cover that assigns to closest landmark
// also any other landmark that is within eps distance of closest landmark
// template over data type, distance
template <typename T, typename M>
bats::Cover landmark_eps_cover(
	const DataSet<T> &X,
	const DataSet<T> &L,
	const M &dist,
	T eps
) {
	// initialize cover
	bats::Cover cover(L.size());
	// loop over points in X
	for (size_t i = 0; i < X.size(); i++) {
		auto di = dist(X[i], L); // distance from X[i] to points in L
		auto it = std::min_element(di.cbegin(), di.cend());
		for (size_t j = 0; j < L.size(); j++) {
			if (di[j] < *it + eps) {
				// add i to any landmark set that is within eps of closest landmark
				// this should include closest point
				cover[j].emplace(i);
			}
		}
	}
	return cover;
}


/*
produce bivariate cover of cover1, cover2
i.e. restriction to diagonal in cover1 x cover2
also return maps to cover1 and cover2 (extends to SimplicialMap on Nerves)
*/
auto bivariate_cover(
	const std::vector<std::set<size_t>> &cover1,
	const std::vector<std::set<size_t>> &cover2
) {
	std::vector<std::set<size_t>> cover12;
	std::vector<size_t> f1;
	std::vector<size_t> f2;

	for (size_t i = 0; i < cover1.size(); i++) {
		for (size_t j = 0; j < cover2.size(); j++) {
			std::set<size_t> intersection;
			bats::util::intersect_sorted(
				cover1[i],
				cover2[j],
				intersection
			);
			if (intersection.size() > 0) {
				// cover1[i] and cover2[j] have nonempy intersection
				cover12.emplace_back(intersection);
				f1.emplace_back(i);
				f2.emplace_back(j);
			}
		}
	}

	return std::make_tuple(cover12, f1, f2);
}

} // namespace bats
