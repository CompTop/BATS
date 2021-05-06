#include <bats.hpp>
#include <vector>
#include <iostream>
#include <math.h>
#include <chrono>
#include <random>

using FT = ModP<int, 2>;

// using CpxT = bats::LightSimplicialComplex<size_t, std::unordered_map<size_t, size_t>>;
using CpxT = bats::SimplicialComplex;

inline size_t _get_idx(size_t i, size_t j, size_t n) {return j + n * i;}

CpxT freudenthal_2d(size_t m, size_t n) {
    CpxT X(m*n, 2);

    for (size_t i = 0; i < m-1; i++) {
        for (size_t j = 0; j < n-1; j++) {
            auto k1 = _get_idx(i,j,n);
            auto k2 = _get_idx(i+1, j, n);
            auto k3 = _get_idx(i, j+1, n);
            auto k4 = _get_idx(i+1, j+1, n);
            X.add_recursive({k1, k2, k4});
            X.add_recursive({k1, k3, k4});
        }
    }

    return X;
}

std::vector<double> get_fn(size_t m, size_t n) {
    std::vector<double> f(m*n);
    for (size_t i = 0; i < m-1; i++) {
        for (size_t j = 0; j < n-1; j++) {
            auto k = _get_idx(i,j,n);
            f[k] = cos(i) + sin(j);
        }
    }
    return f;
}

std::vector<double> get_fn_noisy(size_t m, size_t n, double sigma) {
    auto f = get_fn(m, n);
    std::uniform_real_distribution<double> unif(0, sigma);
    std::default_random_engine re;
    for (auto& fi : f) {
        fi += unif(re);
    }
    return f;
}


int main() {

    /*
    First we build a simplicial complex which can be used to extend filtrations
    */
    size_t m = 256;
    size_t n = 256;
    auto start = std::chrono::steady_clock::now();
    auto  X = freudenthal_2d(m, n);
    auto end = std::chrono::steady_clock::now();
    std::cout << "Construction of Complex: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
        << "ms" << std::endl;
    X.print_summary();

    /*
    Now let's build a filtration and reduce
    */
    auto f = get_fn(m, n);
    start = std::chrono::steady_clock::now();
    auto vals = lower_star_filtration(X, f);
    auto F = bats::Filtration(X, vals);
    auto C = bats::Chain(F, FT());
    auto R = bats::Reduce(C);
    end = std::chrono::steady_clock::now();
    std::cout << "\nBuild and reduce filtration: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
        << "ms" << std::endl;
    R.print_summary();

    /*
    Now, let's change the filtration value
    */
    auto k = _get_idx(m/2, n/2, n);
    // std::cout << "\n" << f[k] << " -> " << 1-f[k] << " at " << k << "\n";
    // f[k] = 1-f[k]; // flip this function value
    f = get_fn_noisy(m, n, 0.05);
    /*
    If we want to change the filtration we have a variety of options.
    First is to simply build a new filtration
    Then, we reduce the complex again
    */
    {
        start = std::chrono::steady_clock::now();
        auto vals = lower_star_filtration(X, f);
        end = std::chrono::steady_clock::now();
        auto t0 = end - start;
        start = std::chrono::steady_clock::now();
        auto F2 = bats::Filtration(X, vals);
        end = std::chrono::steady_clock::now();
        auto t1 = end - start;
        start = std::chrono::steady_clock::now();
        auto C2 = bats::Chain(F2, FT());
        end = std::chrono::steady_clock::now();
        auto t2 = end - start;
        start = std::chrono::steady_clock::now();
        auto R2 = bats::Reduce(C2);
        end = std::chrono::steady_clock::now();
        auto t3 = end - start;
        std::cout << "\nBuild and reduce filtration again: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t0 + t1 + t2 + t3).count()
            << "ms" << std::endl;
        std::cout << "\tFiltration extension: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t0).count()
            << "ms" << std::endl;
        std::cout << "\tFiltration creation: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t1).count()
            << "ms" << std::endl;
        std::cout << "\tFiltered Chain complex creation: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t2).count()
            << "ms" << std::endl;
        std::cout << "\treduction: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t3).count()
            << "ms" << std::endl;
        R2.print_summary();
    }
    /*
    Another option is to update the filtered chain complex
    */
    {
        start = std::chrono::steady_clock::now();
        vals = lower_star_filtration(X, f);
        end = std::chrono::steady_clock::now();
        auto t0 = end - start;
        start = std::chrono::steady_clock::now();
        C.update_filtration(vals);
        end = std::chrono::steady_clock::now();
        auto t1 = end - start;
        start = std::chrono::steady_clock::now();
        auto R3 = bats::Reduce(C);
        end = std::chrono::steady_clock::now();
        auto t2 = end - start;
        std::cout << "\nUpdate filtered chain complex then reduce: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t0 + t1 + t2).count()
            << "ms" << std::endl;
        std::cout << "\tFiltration extension: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t0).count()
            << "ms" << std::endl;
        std::cout << "\tpermute chain complex: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t1).count()
            << "ms" << std::endl;
        std::cout << "\treduction: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t2).count()
            << "ms" << std::endl;
        R3.print_summary();
    }
    /*
    The final option is to update the reduced filtered chain complex
    */
    {
        start = std::chrono::steady_clock::now();
        vals = lower_star_filtration(X, f);
        end = std::chrono::steady_clock::now();
        auto t0 = end - start;
        start = std::chrono::steady_clock::now();
        R.update_filtration(vals);
        end = std::chrono::steady_clock::now();
        auto t1 = end - start;
        std::cout << "\nUpdate reduced filtered chain complex: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t1 + t0).count()
            << "ms" << std::endl;
        std::cout << "\tFiltration extension: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t0).count()
            << "ms" << std::endl;
        std::cout << "\tupdate reduced chain complex: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t1).count()
            << "ms" << std::endl;
        R.print_summary();
    }
    // Let's now flip the filtration value back to get more detail on update
    std::cout << "\n" << f[k] << " -> " << 1-f[k] << " at " << k << "\n";
    f[k] = 1-f[k]; // flip this function value
    {
        start = std::chrono::steady_clock::now();
        vals = lower_star_filtration(X, f);
        end = std::chrono::steady_clock::now();
        auto t0 = end - start;
        std::cout << "\nTime to extend filtration: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t0).count()
            << "ms" << std::endl;

        start = std::chrono::steady_clock::now();
        R.update_filtration(vals);
        end = std::chrono::steady_clock::now();
        auto t1 = end-start;
        std::cout << "\nTime to update reduced filtered chain complex: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t1).count()
            << "ms" << std::endl;
    }


    // /*
    // First we build a simplicial complex which can be used to extend filtrations
    // */
    // CpxT X(3,2);
    //
    // std::vector<size_t> s;
    //
    // s = {0}; X.add(s);
    // s = {1}; X.add(s);
    // s = {2}; X.add(s);
    // s = {0,1}; X.add(s);
    // s = {1,2}; X.add(s);
    // s = {0,2}; X.add(s);
    // s = {0,1,2}; X.add(s);
    //
    // X.print_summary();
    //
    // /*
    // Now let's exend a filtration
    // */
    // std::vector<double> f0 = {0.0, 0.1, 0.2};
    // // lower star filtration
    // std::function<double(const std::vector<size_t>&)> filtfn = [&f0](
    //     const std::vector<size_t>& s
    // ) -> double {
    //     return f0[*std::max_element(s.begin(), s.end(), [&f0](size_t i, size_t j) {return f0[i] < f0[j];})];
    // };
    // auto vals = extend_filtration(X, filtfn);
    // for (auto& valsk: vals) {
    //     for (auto& v: valsk) {
    //         std::cout << v << ", ";
    //     }
    //     std::cout << "\n";
    // }
    //
    // /*
    // The above is wrapped by the lower_star_filtration function
    // */
    // vals = lower_star_filtration(X, f0);
    // for (auto& valsk: vals) {
    //     for (auto& v: valsk) {
    //         std::cout << v << ", ";
    //     }
    //     std::cout << "\n";
    // }
    //
    /*
    Now let's build a filtration and reduce
    */
    // auto F = bats::Filtration(X, vals);
    // auto C = bats::Chain(F, FT());
    // auto R = bats::Reduce(C);
    // R.print_summary();
    //
    // for (auto& p: R.persistence_pairs(0)) {
    //     std::cout << p.str() << std::endl;
    // }
    //
    // /*
    // If we want to change the filtration we have a variety of options.
    // First is to update the filtration on the FilteredChainComplex
    // Then, we reduce the complex again
    // */
    // // update filtration
    // f0 = {1.1, 1.0, 1.2};
    // vals = lower_star_filtration(X, f0);
    // C.update_filtration(vals);
    // R = bats::Reduce(C);
    // R.print_summary();
    //
    // for (auto& p: R.persistence_pairs(0)) {
    //     std::cout << p.str() << std::endl;
    // }
    //
    // /*
    // The second option is to update the ReducedFilteredChainComplex
    // In some situations, this may be faster
    // */
    // // update filtration again
    // f0 = {0.0, 0.1, 0.2};
    // vals = lower_star_filtration(X, f0);
    // R.update_filtration(vals);
    // R.print_summary();
    //
    // for (auto& p: R.persistence_pairs(0)) {
    //     std::cout << p.str() << std::endl;
    // }

    /*
    The final option would be to just call Chain() again to build
    the updated chain complex.
    */

    return EXIT_SUCCESS;
}
