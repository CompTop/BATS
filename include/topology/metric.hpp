#pragma once
/*
various metrics for use with geometric constructions
TODO: add pairwise function
*/
#pragma once

#include <cstddef>
#include <cmath>
#include <vector>
#include <limits>
#include "data.hpp"

namespace bats {

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

    // all pair distances
    template <typename T>
    Matrix<T> operator()(const DataSet<T> &X) const {
        Matrix<T> dists(X.size(), X.size());
        for (size_t i = 0; i < X.size(); i++) {
            for (size_t j = 0; j < X.size(); j++) {
                dists(i,j) = static_cast<const D*>(this)->dist(X[i], X[j]);
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


// Cosine metric 1 - cos(x,y) = 1 - x^T y / ||x||*||y||
struct CosineDist : AbstractMetric<CosineDist> {

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
        return T(1) - cxy;
    }
};

// Cosine metric 1 - cos(x,y) = 1 - x^T y / ||x||*||y||
struct RPCosineDist : AbstractMetric<CosineDist> {

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
        return T(1) - std::abs(cxy);
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
        return std::acos(cxy > T(1) ? T(1) : cxy);
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

/**
compute the enclosing radius from matrix of pairwise distances

this is the minimum of the largest entry of each row
*/
template <typename T>
T enclosing_radius(const Matrix<T>& D) {
    T r_enc = std::numeric_limits<T>::max(); // smallest largest entry
    for (size_t i = 0; i < D.nrow(); ++i) {
        T r_row = T(0); // largest entry in row
        for (size_t j = 0; j < D.ncol(); ++j) {
            r_row = (r_row < D(i,j)) ? D(i,j) : r_row;
        }
        r_enc = (r_enc < r_row) ? r_enc : r_row;
    }

    return r_enc;
}

} // namespace bats
