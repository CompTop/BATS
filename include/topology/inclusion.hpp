#pragma once

#include <vector>
#include <util/set.hpp>

namespace bats {


// vertex inclusion map of set s into set t
// template over container type
template <typename T1, typename T2>
std::vector<size_t> vertex_inclusion_map(const T1 &s, const T2 &t) {
    std::vector<size_t> i;
    i.reserve(s.size());
    size_t ct = 0;
    auto it = t.cbegin();
    for (auto x : s) {
        while ( *it != x) {
            it++;
            ct++;
        }
        i.emplace_back(ct);
    }
    return i;
}

Diagram<std::set<size_t>, std::vector<size_t>> linear_subset_union_diagram(
	std::vector<std::set<size_t>> &subsets
) {
	Diagram<std::set<size_t>, std::vector<size_t>> SetDgm;
	for (size_t i = 0; i < subsets.size(); i++) {
	    auto i1 = SetDgm.add_node(subsets[i]);
	    if (i > 0) {
	        // add backward arrow
	        SetDgm.add_edge(i1, i1-1, vertex_inclusion_map(SetDgm.node[i1], SetDgm.node[i1-1]));
	    }
	    if (i == (subsets.size() - 1)) { break; }
	    auto i2 = SetDgm.add_node(bats::util::set_union(subsets[i], subsets[i+1]));
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
	    auto i2 = SetDgm.add_node(bats::util::set_intersection(cover[i], cover[i+1]));
	    // map from i2 to i1
	    SetDgm.add_edge(i2, i1, vertex_inclusion_map(SetDgm.node[i2], SetDgm.node[i1]));
	}

	return SetDgm;
}

} // namespace bats
