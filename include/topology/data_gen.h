#pragma once
/*
utilities for generating data
*/

#include <cmath>
#include <vector>
#include <random>
#include <iterator>
#include <iostream>

#include "metric.h"

template <typename T>
void add_normal_noise(std::vector<T> &x,
    const T mu=T(0),
    const T sigma=T(1)
) {
    std::default_random_engine generator;
    std::normal_distribution<T> distribution(mu,sigma);

    #pragma omp simd
    for (size_t i = 0; i < x.size(); i++) {
        x[i] += distribution(generator);
    }
}

// form product space X x Y
template <typename T>
std::vector<T> product_space(
    const std::vector<T> &x,
    size_t dx,
    const std::vector<T> &y,
    size_t dy
) {
    size_t nx = x.size()/dx;
    size_t ny = y.size()/dy;
    size_t nz = nx * ny;
    size_t dz = dx + dy;
    //std::cout << 'x' << nx << 'y' << ny << 'z' <<  nz << 'd' << dz << std::endl;
    std::vector<T> z;
    z.reserve(nz * dz);
    for (size_t i = 0; i < nx*dx; i+=dx) {
        for (size_t j = 0; j < ny*dy; j+=dy) {
            // put x[i] in first coordinate block
            for (size_t ii = i; ii < i + dx; ii++) {
                z.emplace_back(x[ii]);
            }
            // put y[j] in second coordinate block
            for (size_t jj = j; jj < j + dy; jj++) {
                z.emplace_back(y[jj]);
            }
        }
    }
    return z;
}


/*
    generate (d-1)-dimensional sphere with n samples
*/
template <typename T>
std::vector<T> sample_sphere(
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

// sample n points uniformly at random from d dimensional cube
template <typename T>
std::vector<T> sample_cube(
    const size_t d,
    const size_t n
) {
    std::default_random_engine generator;
    std::uniform_real_distribution<T> distribution(0.0,1.0);

    // fill with gaussian random numbers
    std::vector<T> x(d*n);
    auto it = x.begin();
    while (it < x.end()) {
        *it = distribution(generator);
        ++it;
    }

    return x;

}

// n equispaced points between min and max
template <typename T>
std::vector<T> interval(
    const T min,
    const T max,
    const size_t n
) {
    std::vector<T> x(n);
    T stride = (max - min) / (n-1);
    x[0] = min;
    for (size_t i = 1; i < n; i++) {
        x[i] = x[i-1] + stride;
    }
    return x;
}

// circle embedded in 2d on n points
template <typename T>
std::vector<T> circle(
    const T rad,
    const size_t n
) {
    std::vector<T> x(2*n);
    T dtheta = 2 * M_PI / n;
    T theta = T(0);
    for (size_t i = 0; i < n; i++) {
        x[2*i] = rad * std::cos(theta);
        x[2*i + 1] = rad * std::sin(theta);
        theta += dtheta;
    }
    return x;
}

/*
    generate 3-d cylinder that is I x S^1
    total number of points is n_len x n_cir
*/
std::vector<double> gen_cylinder(
    const size_t n_len,
    const size_t n_cir
) {
    auto I = interval(0.0, 1.0, n_len);
    auto S = circle(1.0, n_cir);

    return product_space(I, 1, S, 2);
}
