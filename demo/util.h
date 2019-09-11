# pragma once
// utilities for generating examples

#include <cmath>
#include <vector>
#include <random>
#include <iterator>
#include <iostream>

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

/*
    generate d-dimensional sphere with n samples
*/
template <typename T>
std::vector<T> gen_sphere(
    const size_t d,
    const size_t n
) {
    std::default_random_engine generator;
    std::normal_distribution<T> distribution(0.0,1.0);

    // fill with gaussian random numbers
    std::vector<T> x(d*n);
    auto it = x.begin();
    while (it < x.end()) {
        *it = distribution(generator);
        ++it;
    }

    // normalize
    for (size_t i = 0; i < n; i++) {
        auto vst = x.begin() + (i*d);
        auto vend = x.begin() + (i*d + d);
        T vnorm = norm(vst, vend);
        while (vst < vend) {
            *vst /= vnorm;
            ++vst;
        }
    }

    return x;

}

// print bars with length > minlen
template <typename T>
void print_barcodes(
    std::vector<std::vector<T>> bars,
    T minlen=T(0)
) {
    for (size_t d = 0; d < bars.size(); d++) {
        std::cout << "dimension " << d << std::endl;
        auto it = bars[d].cbegin();
        while (it != bars[d].cend()) {
            T b = *it++;
            T d = *it++;
            T len = d - b;
            if (len > minlen) {
                std::cout << "   (" << b << ',' << d << ')' << std::endl;
            }
        }
    }
    return;

}
