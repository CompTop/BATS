#pragma once
/*
various metrics for use with geometric constructions
TODO: add pairwise function
*/

#include <cstddef>
#include <cmath>
#include <vector>
#include "data.h"

// CRTP over derived class D
template <class D>
struct AbstractMetric {

    template <typename T>
    inline T dist(const VectorView<T> &x, const VectorView<T> &y) const {
        return static_cast<const D*>(this)->dist(x, y);
    }

    // distance between two points
    template <typename T>
    inline T operator() (const VectorView<T> &x, const VectorView<T> &y) const {
        return static_cast<const D*>(this)->dist(x, y);
    };

    // all distances from point x
    template <typename T>
    std::vector<T> operator()(const VectorView<T> &x, const DataSet<T> &X) const {
        std::vector<T> dists(X.size());
        for (int i = 0; i < X.size(); i++) {
            dists[i] = static_cast<const D*>(this)->dist(x, X[i]);
        }
        return dists;
    }

    // all pair distances
    template <typename T>
    Matrix<T> operator()(const DataSet<T> &X, const DataSet<T> &Y) const {
        Matrix<T> dists(X.size(), Y.size());
        for (size_t i = 0; i < X.size(); i++) {
            for (size_t j = 0; j < Y.size(); j++) {
                dists(i,j) = static_cast<const D*>(this)->dist(X[i], Y[j]);
            }
        }
        return dists;
    }
};


// metrics are implemented as structs
struct Euclidean : AbstractMetric<Euclidean> {

    template <typename T>
    T dist (const VectorView<T> &x, const VectorView<T> &y) const {
        T n = T(0);
        for (size_t i = 0; i < x.size(); i++) {
            T diff = x(i) - y(i);
            n += diff * diff;
        }
        return std::sqrt(n);
    }

};

struct L1Dist : AbstractMetric<L1Dist> {

    template <typename T>
    T dist (const VectorView<T> &x, const VectorView<T> &y) const {
        T n = T(0);
        for (size_t i = 0; i < x.size(); i++) {
            n += std::abs(x(i) - y(i));
        }
        return n;
    }
};

struct LInfDist : AbstractMetric<LInfDist> {

    template <typename T>
    T dist (const VectorView<T> &x, const VectorView<T> &y) const {
        T n = T(0);
        for (size_t i = 0; i < x.size(); i++) {
            T diff = std::abs(x(i) - y(i));
            n = (n > diff) ? n : diff;
        }
        return n;
    }
};
