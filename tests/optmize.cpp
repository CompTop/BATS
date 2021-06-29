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
    auto [vals, imap] = lower_star_filtration(X, f);
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
    // auto k = _get_idx(m/2, n/2, n);
    // std::cout << "\n" << f[k] << " -> " << 1-f[k] << " at " << k << "\n";
    // f[k] = 1-f[k]; // flip this function value
    f = get_fn_noisy(m, n, 0.05);
    /*
    The final option is to update the reduced filtered chain complex
    */
    double lr = 0.1;
    for (size_t i = 0; i< 6; i++){
        std::cout << "\niter "<< i << std::endl;
        start = std::chrono::steady_clock::now();
        auto [vals, imap] = lower_star_filtration(X, f);
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
        auto nnz_of_R = R.get_nnz_R();
        auto nnz_of_U = R.get_nnz_U();
        std::cout << "nnz of R: ";
        print_1D_vectors(nnz_of_R);
        std::cout << "nnz of U: ";
        print_1D_vectors(nnz_of_U);

        if(i == 3){
            R.sparsify_basis();
            std::cout << "\ni = 3, after sparsify_basis()" << std::endl;
            std::cout << "nnz of R: ";
            print_1D_vectors(R.get_nnz_R());
            std::cout << "nnz of U: ";
            print_1D_vectors(R.get_nnz_U());
        }

        auto ps = R.persistence_pairs(1);
        // std::cout << "\nmaximum of imap[1] ";
        // std::cout << *std::max_element(imap[1].begin(), imap[1].end()) << std::endl;
        // std::cout << "\nmaximum of imap[2] ";
        // std::cout << *std::max_element(imap[2].begin(), imap[2].end()) << std::endl;
        // std::cout << "\nf.size() = "<< f.size() << std::endl;

        for(auto p: ps){
            size_t d = p.dim;
            size_t bi = imap[d][p.birth_ind]; // maps birth_ind to pixel where it appeared
            size_t di = imap[d+1][p.death_ind]; // maps death_ind to pixel where it appeared
            f[bi] = f[bi] - lr;
            f[di] = f[di] + lr;
        }
    }

    return EXIT_SUCCESS;
}
