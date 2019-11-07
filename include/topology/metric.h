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
