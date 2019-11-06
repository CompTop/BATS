#pragma once

#include <cstddef>
#include <cmath>
#include <vector>

template <typename TI>
auto norm(TI start, const TI end) {
    using T = typename std::iterator_traits<TI>::value_type;
    T res = 0;
    while (start != end) {
        res += (*start) * (*start);
        ++start;
    }
    return std::sqrt(res);
}

// euclidean distance in d dimensions
template <typename T, typename TI>
T euclidean(TI a, TI b, size_t d){
    T res = 0;
    for (size_t i = 0; i < d; i++) {
        res += std::pow((*a++ - *b++), 2);
    }
    return std::sqrt(res);
}

// construct edges for rips complex
template <typename T>
std::vector<size_t> rips_edges(const std::vector<T> &x, const size_t d, const T rmax) {
    std::vector<size_t> edges;
    size_t nedges = x.size() * (x.size() - d) / (d * d * 2);
    edges.reserve(2*nedges);
    for (size_t i = 0; i < (x.size() / d); i++) {
        for (size_t j = 0; j < i; j++) {
            T dij = euclidean<T>(x.cbegin() + d*i, x.cbegin() + d*j, d);
            if (dij < rmax) {
                edges.push_back(j);
                edges.push_back(i);
            }
        }
    }
    return edges;
}
