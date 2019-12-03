#pragma once
/*
various metrics for use with geometric constructions
TODO: add pairwise function
*/

#include <cstddef>
#include <cmath>
#include <vector>
#include "data.h"

struct AbstractMetric {};

// metrics are implemented as structs
struct Euclidean : AbstractMetric {

    template <typename T>
    T operator() (const VectorView<T> &x, const VectorView<T> &y) const {
        T n = T(0);
        for (size_t i = 0; i < x.size(); i++) {
            T diff = x(i) - y(i);
            n += diff * diff;
        }
        return std::sqrt(n);
    }
};

struct L1Dist : AbstractMetric {

    template <typename T>
    T operator() (const VectorView<T> &x, const VectorView<T> &y) const {
        T n = T(0);
        for (size_t i = 0; i < x.size(); i++) {
            n += std::abs(x(i) - y(i));
        }
        return n;
    }
};

struct LInfDist : AbstractMetric {

    template <typename T>
    T operator() (const VectorView<T> &x, const VectorView<T> &y) const {
        T n = T(0);
        for (size_t i = 0; i < x.size(); i++) {
            T diff = std::abs(x(i) - y(i));
            n = (n > diff) ? n : diff;
        }
        return n;
    }
};
