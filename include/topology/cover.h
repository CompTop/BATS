#pragma once
/*
utilities for creating covers
*/

#include <cstddef>
#include <vector>
#include <set>
#include <algorithm>
#include <cmath>

#include <multigraph/diagram.h>
#include "inclusion.h"

// return union of two sets
// template over containter types
template <typename CT1, typename CT2>
std::set<size_t> set_union(const CT1 &s1, const CT2 &s2) {
    std::set<size_t> s;
    for (auto x : s1) {
        s.insert(x);
    }
    for (auto x : s2) {
        s.insert(x);
    }
    return s;
}

// intersection of two sets
// template over container type
template <typename C1, typename C2>
std::set<size_t> set_intersection(const C1 &a, const C2 &b) {
	std::set<size_t> s;

	auto is = s.end();
    auto ia = a.cbegin();
    auto ib = b.cbegin();
    while (ia != a.cend() && ib != b.cend()) {
        if (*ia < *ib) {
            ++ia;
        } else if (*ib < *ia) {
            ++ib;
        } else {
            s.emplace_hint(is, *ia);
			is = s.end();
            ++ia;
            ++ib;
        }
	}
	return s;
}

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

Diagram<std::set<size_t>, std::vector<size_t>> linear_cover_union_diagram(
	std::vector<std::set<size_t>> &cover
) {
	Diagram<std::set<size_t>, std::vector<size_t>> SetDgm;
	for (size_t i = 0; i < cover.size(); i++) {
	    auto i1 = SetDgm.add_node(cover[i]);
	    if (i > 0) {
	        // add backward arrow
	        SetDgm.add_edge(i1, i1-1, vertex_inclusion_map(SetDgm.node[i1], SetDgm.node[i1-1]));
	    }
	    if (i == (cover.size() - 1)) { break; }
	    auto i2 = SetDgm.add_node(set_union(cover[i], cover[i+1]));
	    // map from i1 to i2
	    SetDgm.add_edge(i1, i2, vertex_inclusion_map(SetDgm.node[i1], SetDgm.node[i2]));
	}

	return SetDgm;
}

Diagram<std::set<size_t>, std::vector<size_t>> linear_cover_intersection_diagram(
	std::vector<std::set<size_t>> &cover
) {
	Diagram<std::set<size_t>, std::vector<size_t>> SetDgm;
	for (size_t i = 0; i < cover.size(); i++) {
	    auto i1 = SetDgm.add_node(cover[i]);
	    if (i > 0) {
	        // add arrow from inclusion
	        SetDgm.add_edge(i1-1, i, vertex_inclusion_map(SetDgm.node[i1-1], SetDgm.node[i1]));
	    }
	    if (i == (cover.size() - 1)) { break; }
	    auto i2 = SetDgm.add_node(set_intersection(cover[i], cover[i+1]));
	    // map from i2 to i1
	    SetDgm.add_edge(i2, i1, vertex_inclusion_map(SetDgm.node[i2], SetDgm.node[i1]));
	}

	return SetDgm;
}





// project data x in d dimensions onto coordinate i
template <typename T>
std::vector<T> coordinate_projection(const std::vector<T> &x, const size_t d, const size_t i) {
    std::vector<T> p;
    p.reserve(x.size() / d);
    for (auto it = x.cbegin(); it != x.cend(); it += d) {
        p.emplace_back(*(it + i));
    }
    return p;
}

// get subset of data
template <typename T>
std::vector<T> get_subset(
    const std::vector<T> &x,
    const size_t d,
    const std::set<size_t> &ind
) {
    std::vector<T> xs;
    xs.reserve(ind.size() * d);
    for (auto it = ind.cbegin(); it != ind.cend(); it++) {
        for (size_t i = 0; i < d; i++) {
            xs.emplace_back(x[(*it) * d + i]);
        }
    }
    return xs;
}
