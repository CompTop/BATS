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
            f[k] = cos(4.0*3.1415926*i) + sin(4.0*3.1415926*j);
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
    size_t m = 100;
    size_t n = 100;
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
    auto vipair = lower_star_filtration(X, f);
    auto vals = std::get<0>(vipair);

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
    f = get_fn_noisy(m, n, 0.05);
    vipair = lower_star_filtration(X, f);
    vals = std::get<0>(vipair);
    // Compute Persistent Homology
    auto F_Y = bats::Filtration(X, vals);
    auto C_Y = bats::Chain(F_Y, FT());
    auto R_Y = bats::Reduce(C_Y);
    /*
    There are two options one is used for Filtration with complex that does not change
    */

    
    {
        C = bats::Chain(F, FT());
        R = bats::Reduce(C);
        R.update_filtration(vals);

        if(test_reduce_result(R_Y, R)){
            std::cout << "By compare two RFCC, this method success!!"  << std::endl;
        }
    }

    {
        C = bats::Chain(F, FT());
        R = bats::Reduce(C);

        auto UI = Update_info(F, F_Y);
        //the complex in a lower star filtration is not sorted by their filtration values
        UI.filtered_info(C.perm); 

        R.update_filtration_general(UI);

        if(test_reduce_result(R_Y, R)){
            std::cout << "By compare two RFCC, this method success!!"  << std::endl;
        }
    }


    return EXIT_SUCCESS;
}