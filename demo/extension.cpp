#include <bats.hpp>
#include <vector>
#include <iostream>
#include <math.h>
#include <chrono>
#include <random>

using FT = ModP<int, 2>;

// using CpxT = bats::LightSimplicialComplex<size_t, std::unordered_map<size_t, size_t>>;
using CpxT = bats::SimplicialComplex;

std::vector<double> get_fn(size_t m, size_t n) {
    std::vector<double> f(m*n);
    for (size_t i = 0; i < m-1; i++) {
        for (size_t j = 0; j < n-1; j++) {
            auto k = bats::rowmajor::get_idx(i,j,n);
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
    auto  X = bats::Freudenthal<CpxT>(m, n);
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
    auto [vals, cinds] = lower_star_filtration(X, f);
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
    auto k = bats::rowmajor::get_idx(m/2, n/2, n);
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
        auto [vals, cinds] = lower_star_filtration(X, f);
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
        auto [vals, cinds] = lower_star_filtration(X, f);
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
        auto [vals, cinds] = lower_star_filtration(X, f);
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
        auto [vals, cinds] = lower_star_filtration(X, f);
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


    return EXIT_SUCCESS;
}
