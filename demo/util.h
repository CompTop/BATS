# pragma once
// utilities for generating examples

#include <cmath>
#include <vector>
#include <random>
#include <iterator>

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
