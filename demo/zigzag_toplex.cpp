#include <bats.hpp>
#include <vector>
#include <iostream>
#include <chrono>

#include <cmath>

/*
Generate a smooth cube function for image levelset

output will be stored in column-major order

This generates same smooth function as in Morse zigzag paper.

Function is on cube [-2,2]^3
*/
template <typename T=double>
auto smooth_cube_fn(
    const size_t n
) {
    std::vector<std::vector<std::vector<T>>> f0(n);
    // initialize
    for (size_t i = 0; i < n; ++i) {
        f0[i] = std::vector<std::vector<T>>(n);
        for (size_t j = 0; j < n; ++j) {
            f0[i][j] = std::vector<T>(n);
        }
    }


    for (size_t i = 0; i < n; ++i) {
        T x = i * T(4) / (n - 1) - T(2);
        T sx = std::sin(x);
        T s2x = std::sin(2*x);
        T s3x = std::sin(3*x);
        T cx = std::cos(x);
        T c2x = std::cos(2*x);
        T c3x = std::cos(3*x);
        for (size_t j = 0; j < n; ++j) {
            T y = j * T(4) / (n - 1) - T(2);
            T sy = std::sin(y);
            T s2y = std::sin(2*y);
            T s3y = std::sin(3*y);
            T cy = std::cos(y);
            T c2y = std::cos(2*y);
            T c3y = std::cos(3*y);
            for (size_t k = 0; k < n; ++k) {
                T z = k * T(4) / (n - 1) - T(2);
                T sz = std::sin(z);
                T s2z = std::sin(2*z);
                T s3z = std::sin(3*z);
                T cz = std::cos(z);
                T c2z = std::cos(2*z);
                T c3z = std::cos(3*z);
                f0[i][j][k] =
                    T(1) * sx * s2y * s3z +
                    T(2) * s2x * sy * s3z +
                    T(3) * s3x * s2y * sz +
                    T(4) * sx * s3y * s2z +
                    T(5) * s2x * s3y * sz +
                    T(6) * s3x * sy * s2z +
                    T(1) * c3x * cy * c2z +
                    T(2) * c2x * cy * c3z +
                    T(3) * cx * c2y * c3z +
                    T(4) * c3x * c2y * cz +
                    T(5) * c2x * c3y * cz +
                    T(6) * cx * c3y * c2z;
            }
        }
    }
    return f0;
}

int main() {

    size_t n = 32;
    auto img = smooth_cube_fn(n);
    auto X = bats::zigzag_toplex(img);

    auto X0 = X.levelset(-1.0,1.0);
    X0.print_summary();

    auto D = bats::zigzag_levelsets(X, 0.1, -15.0, 15.0);


    return EXIT_SUCCESS;

}
