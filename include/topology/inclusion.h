#pragma once

#include <vector>

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
