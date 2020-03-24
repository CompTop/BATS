#pragma once
/*
various metrics for use with geometric constructions
TODO: add pairwise function
*/
#pragma once

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
        for (size_t i = 0; i < X.size(); i++) {
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

// standard Euclidean distance
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

// L1 metric
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

// L-infinity metric
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


// Angular metric
struct AngleDist : AbstractMetric<AngleDist> {

    template <typename T>
    T dist (const VectorView<T> &x, const VectorView<T> &y) const {
        T nx = T(0);
        T ny = T(0);
        T ip = T(0);
        for (size_t i = 0; i < x.size(); i++) {
            nx += (x(i) * x(i));
            ip += (x(i) * y(i));
            ny += (y(i) * y(i));
        }
        T cxy = ip / (std::sqrt(nx) * std::sqrt(ny)); // cosine of angle
        return std::acos(cxy > 1.0 ? 1.0 : cxy);
    }
};

// metric on RP based on angular distance between representatives
struct RPAngleDist : AbstractMetric<RPAngleDist> {

    template <typename T>
    T dist (const VectorView<T> &x, const VectorView<T> &y) const {
        T nx = T(0);
        T ny = T(0);
        T ip = T(0);
        for (size_t i = 0; i < x.size(); i++) {
            nx += x(i) * x(i);
            ip += x(i) * y(i);
            ny += y(i) * y(i);
        }
        ip = std::abs(ip); // positive inner produce minimizes distance between RP representatives
        T cxy = ip / (std::sqrt(nx) * std::sqrt(ny)); // cosine of angle
        return std::acos(cxy);
    }

};
