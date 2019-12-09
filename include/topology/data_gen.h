#pragma once
/*
utilities for generating data
*/

#include <cmath>
#include <vector>
#include <random>
#include <iterator>
#include <iostream>

#include "data.h"



// add normal_noise
// operates in-place on X
template <typename T>
Matrix<T>& add_normal_noise(
    Matrix<T> &X,
    const T mu=T(0),
    const T sigma=T(1)
) {
    std::default_random_engine generator;
    std::normal_distribution distribution(mu,sigma);

    #pragma omp simd
    for (size_t i = 0; i < X.size(); i++) {
        X(i) += distribution(generator);
    }
    return X;
}

template <typename T>
inline DataSet<T>& add_normal_noise(
    DataSet<T> &X,
    const T mu=T(0),
    const T sigma=T(1)
) {
    add_normal_noise(X.data, mu, sigma);
    return X;
}


template <typename T>
Matrix<T>& add_uniform_noise(Matrix<T> &X,
    const T lb = T(0),
    const T ub = T(1)
) {
    std::default_random_engine generator;
    std::uniform_real_distribution distribution(lb, ub);

    #pragma omp simd
    for (size_t i = 0; i < X.size(); i++) {
        X(i) += distribution(generator);
    }
    return X;
}

template <typename T>
inline DataSet<T>& add_uniform_noise(
    DataSet<T> &X,
    const T lb = T(0),
    const T ub = T(1)
) {
    return add_uniform_noise(X.data, lb, ub);
}


// form product space X x Y
template <typename T>
DataSet<T> product_space(
    const DataSet<T> &X,
    const DataSet<T> &Y
) {
    size_t dx = X.dim();
    size_t nx = X.size();
    size_t dy = Y.dim();
    size_t ny = Y.size();

    size_t nz = nx * ny;
    size_t dz = dx + dy;
    //std::cout << 'x' << nx << 'y' << ny << 'z' <<  nz << 'd' << dz << std::endl;
    Matrix<T> Z(dz, nz);
    //z.reserve(nz * dz);
    size_t k = 0;
    for (size_t i = 0; i < nx; i++) {
        for (size_t j = 0; j < ny; j++) {
            // put x[i] in first coordinate block
            for (size_t di = 0; di < dx; di++) {
                Z(k++) = X(di, i);
            }
            // put y[j] in second coordinate block
            for (size_t dj = 0; dj < dy; dj++) {
                Z(k++) = Y(dj, j);
            }
        }
    }
    return DataSet(Z);
}


/*
    generate (d-1)-dimensional sphere with n samples
*/
template <typename T>
DataSet<T> sample_sphere(
    const size_t d,
    const size_t n
) {

    // fill with gaussian random numbers
    Matrix<T> X(d, n);
    add_normal_noise(X);

    // normalize
    for (size_t j = 0; j < n; j++) {
        T vnorm = norm(X[j]);
        X[j] /= vnorm;
    }

    return DataSet(X);

}

// sample n points uniformly at random from d dimensional cube
template <typename T>
DataSet<T> sample_cube(
    const size_t d,
    const size_t n
) {
    Matrix<T> X(d, n);
    add_uniform_noise(X);
    return DataSet(X);
}

// n equispaced points between min and max
template <typename T>
DataSet<T> interval(
    const T min,
    const T max,
    const size_t n
) {
    Matrix<T> x(1,n);
    T stride = (max - min) / (n-1);
    x(0) = min;
    for (size_t i = 1; i < n; i++) {
        x(i) = x(i-1) + stride;
    }
    return DataSet(x);
}

// circle embedded in 2d on n points
template <typename T>
DataSet<T> circle(
    const T rad,
    const size_t n
) {
    Matrix<T> x(2, n);
    T dtheta = 2 * M_PI / n;
    T theta = T(0);
    for (size_t i = 0; i < n; i++) {
        x(0, i) = rad * std::cos(theta);
        x(1, i) = rad * std::sin(theta);
        theta += dtheta;
    }
    return DataSet(x);
}

/*
    generate 3-d cylinder that is I x S^1
    total number of points is n_len x n_cir
*/
DataSet<double> gen_cylinder(
    const size_t n_len,
    const size_t n_cir
) {
    auto I = interval(0.0, 1.0, n_len);
    auto S = circle(1.0, n_cir);

    return product_space(I, S);
}
