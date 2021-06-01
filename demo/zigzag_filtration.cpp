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
std::vector<T> smooth_cube_fn(
    const size_t n
) {
    std::vector<T> f0(n*n*n);

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
                f0[k + (j + i*n)*n] =
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

template <typename T>
auto gen_cube_zigzag(
    const size_t n,
    const T eps
) {
    auto f0 = smooth_cube_fn(n);
    bats::CubicalComplex X = bats::CubicalComplex::generate_cube(n);

    return bats::extend_zigzag_filtration(f0, X, eps, n);
}

auto gen_rips_cylinder(
    const size_t n_len, // number of points in interval
    const size_t n_cir, // number of points around circle
    const double r, // Rips parameter
    const double eps // levelset radius
) {
    auto X = bats::gen_cylinder(n_len, n_cir, 0.3);
    auto p = bats::coordinate_projection(X, 0);

    auto R = bats::RipsComplex<bats::SimplicialComplex>(
        X, bats::Euclidean(), r, 2
    );
    R.print_summary();

    return bats::extend_zigzag_filtration(p, R, eps);
}

auto gen_rips_line(
    const size_t n_len, // number of points in interval
    const double r, // Rips parameter
    const double eps // levelset radius
) {
    auto X = bats::interval(0.0, 1.0, n_len);
    auto p = bats::coordinate_projection(X, 0);

    auto R = bats::RipsComplex<bats::SimplicialComplex>(
        X, bats::Euclidean(), r, 2
    );
    R.print_summary();

    return bats::extend_zigzag_filtration(p, R, eps);
}

// determine Lipschitz constant
template <typename T>
T lipschitz_constant(
    const std::vector<T>& f0,
    const size_t n
) {
    T lc = std::numeric_limits<T>::min();
    for (size_t i = 0; i < n-1; ++i) {
        for (size_t j = 0; j < n-1; ++j) {
            for (size_t k = 0; k < n-1; ++k) {
                T cval = bats::detail::cube_val(f0, i, j, k, n);
                T cvali =bats::detail::cube_val(f0, i+1, j, k, n);
                T cvalj =bats::detail::cube_val(f0, i, j+1, k, n);
                T cvalk =bats::detail::cube_val(f0, i, j, k+1, n);
                lc = std::max(lc, std::abs(cvali - cval));
                lc = std::max(lc, std::abs(cvalj - cval));
                lc = std::max(lc, std::abs(cvalk - cval));
            }
        }
    }
    return lc;
}

int main() {

    bats::ZigzagFiltration<bats::SimplicialComplex> F;

    std::vector<size_t> spx;
    // create a cycle that persists for a while
    spx = {0,1}; F.add_recursive(0.0, 10.0, spx);
    spx = {0,2}; F.add_recursive(0.0, 10.0, spx);
    spx = {1,2}; F.add_recursive(0.0, 10.0, spx);

    {
        F.complex().print_summary();

        using F2 = ModP<int, 2>;
        auto ps = bats::barcode(F, F2(),
            bats::no_optimization_flag(),
            bats::standard_reduction_flag()
        );

        for (auto& pk : ps) {
            for (auto p : pk) {
                if (p.length() > 0)
                    std::cout << p.str() << std::endl;
            }
        }
    }

    // now block cycle for some period of time
    std::cout <<"\nadding block:\n";
    spx = {0,1,2}; F.add(2.0, 4.0, spx);

    {
        F.complex().print_summary();

        using F2 = ModP<int, 2>;
        auto ps = bats::barcode(F, F2(),
            bats::no_optimization_flag(),
            bats::standard_reduction_flag()
        );

        for (auto& pk : ps) {
            for (auto p : pk) {
                if (p.length() > 0)
                    std::cout << p << std::endl;
            }
        }
    }

    std::cout << "\n\nRips on line:\n";
    {

        using F2 = ModP<int, 2>;
        double r = 0.3;
        double eps = 0.15;
        auto start = std::chrono::steady_clock::now();
        // auto F = gen_rips_cylinder(10, 20,  r, eps);
        auto F = gen_rips_line(10, r, eps);
        auto end = std::chrono::steady_clock::now();
        std::cout << "\nSetup: "
            << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
            << "\u03BCs" << std::endl;
        F.complex().print_summary();

        auto R = bats::Reduce(F.complex(), F2());
        for (size_t k = 0; k < R.maxdim()+1; ++k) {
            std::cout <<"betti " << k << ": " << R.hdim(k) << std::endl;
        }


        start = std::chrono::steady_clock::now();
        using F2 = ModP<int, 2>;
        auto ps = bats::barcode(F, F2(),
            bats::no_optimization_flag(),
            bats::standard_reduction_flag()
        );
        end = std::chrono::steady_clock::now();
        std::cout << "\nCompute barcode: "
            << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
            << "\u03BCs" << std::endl;

        for (auto& pk : ps) {
            for (auto p : pk) {
                if (p.length() > 0)
                    std::cout << p.str() << std::endl;
            }
        }
    }

    std::cout << "\n\nRips on cylinder:\n";
    {

        using F2 = ModP<int, 2>;
        double r = 0.3;
        double eps = 0.15;
        auto start = std::chrono::steady_clock::now();
        auto F = gen_rips_cylinder(10, 10,  r, eps);
        // auto F = gen_rips_line(10, r, eps);
        auto end = std::chrono::steady_clock::now();
        std::cout << "\nSetup: "
            << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
            << "\u03BCs" << std::endl;
        F.complex().print_summary();

        auto R = bats::Reduce(F.complex(), F2());
        for (size_t k = 0; k < R.maxdim()+1; ++k) {
            std::cout <<"betti " << k << ": " << R.hdim(k) << std::endl;
        }


        start = std::chrono::steady_clock::now();
        using F2 = ModP<int, 2>;
        auto ps = bats::barcode(F, F2(),
            bats::no_optimization_flag(),
            bats::standard_reduction_flag()
        );
        end = std::chrono::steady_clock::now();
        std::cout << "\nCompute barcode: "
            << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
            << "\u03BCs" << std::endl;

        for (auto& pk : ps) {
            for (auto p : pk) {
                if (p.length() > 0)
                    std::cout << p.str() << std::endl;
            }
        }
    }

    // TODO: create Rips complex on samples
    // extend function on vertices

    // // test for generating cubes
    // // for (auto n : {3, 5, 9, 17, 33, 65, 129}) {
    // {
    //     size_t n = 65;
    //     double eps = 8.0;
    //     auto start = std::chrono::steady_clock::now();
    //     auto scon = start;
    //     bats::CubicalComplex X = bats::CubicalComplex::generate_cube(n);
    //     auto end = std::chrono::steady_clock::now();
    //     std::cout << "\nBuild cube on " << n << "^3 vertices: "
    //         << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
    //         << "\u03BCs" << std::endl;
    //     X.print_summary();
    //
    //     start = std::chrono::steady_clock::now();
    //     auto f0 = smooth_cube_fn(n);
    //     end = std::chrono::steady_clock::now();
    //     std::cout << "\nBuild smooth cube function: "
    //         << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
    //         << "\u03BCs" << std::endl;
    //
    //     std::cout << "Lipschitz constant: " << lipschitz_constant(f0, n) << std::endl;
    //     auto it1 = std::min_element(f0.begin(), f0.end());
    //     std::cout << "minimum: " << *it1 << std::endl;
    //     auto it2 = std::max_element(f0.begin(), f0.end());
    //     std::cout << "maximum: " << *it2 << std::endl;
    //
    //
    //     start = std::chrono::steady_clock::now();
    //     auto val = bats::extend_levelset(f0, X, eps, n);
    //     end = std::chrono::steady_clock::now();
    //     std::cout << "\nExtend to zigzag filtration: "
    //         << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
    //         << "\u03BCs" << std::endl;
    //
    //     start = std::chrono::steady_clock::now();
    //     auto XF = bats::ZigzagFiltration(X, val);
    //     end = std::chrono::steady_clock::now();
    //     std::cout << "\nPut in ZigzagFiltration: "
    //         << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
    //         << "\u03BCs" << std::endl;
    //     auto econ = end;
    //     std::cout << "\nTotal setup: "
    //         << std::chrono::duration_cast<std::chrono::microseconds>(econ - scon).count()
    //         << "\u03BCs" << std::endl;
    //
    //
    //     start = std::chrono::steady_clock::now();
    //     using F2 = ModP<int, 2>;
    //     auto ps = bats::barcode(XF, F2(),
    //         bats::no_optimization_flag(),
    //         bats::standard_reduction_flag(),
    //         bats::no_apparent_pairs_flag()
    //     );
    //     end = std::chrono::steady_clock::now();
    //     std::cout << "\nCompute barcode: "
    //         << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
    //         << "\u03BCs" << std::endl;
    // }
    //
    // {
    //     size_t n = 8;
    //     double eps = 8.0;
    //     auto start = std::chrono::steady_clock::now();
    //     auto F = gen_cube_zigzag(n, eps);
    //     auto end = std::chrono::steady_clock::now();
    //     std::cout << "\nSetup: "
    //         << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
    //         << "\u03BCs" << std::endl;
    //     F.complex().print_summary();
    //
    //
    //     start = std::chrono::steady_clock::now();
    //     using F2 = ModP<int, 2>;
    //     auto ps = bats::barcode(F, F2(),
    //         bats::no_optimization_flag(),
    //         bats::standard_reduction_flag()
    //     );
    //     end = std::chrono::steady_clock::now();
    //     std::cout << "\nCompute barcode: "
    //         << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
    //         << "\u03BCs" << std::endl;
    // }




    return EXIT_SUCCESS;
}
