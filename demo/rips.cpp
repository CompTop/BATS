#include <iostream>
#include <filtration/flag.h>
#include <filtration/rips.h>
#include <linalg/field.h>
#include <linalg/sparse_vector.h>
#include <linalg/set_vector.h>
#include <linalg/col_matrix.h>
#include <persistence/homology_reduction.h>
#include <vector>
#include <random>
#include <chrono>
#include "util.h"

#define TE tedge<float, size_t>
// filtration type
#define TF double

// field type
#define FT ModP<int, 3>
#define CT SparseVector<FT>
// #define CT SetVector<FT>
#define MT ColumnMatrix<CT>

int main() {

    size_t d = 2; // dimension of sphere + 1
    size_t n = 100;
    // size_t n = 4;
    size_t maxdim = d;

    std::cout << "creating pairwise distances" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<TF> x = gen_sphere<TF>(d, n);
    // for (auto xi : x) {
    //     std::cout << xi << ',';
    // }
    // std::cout << std::endl;
    // std::vector<TF> x = {1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0};
    std::vector<size_t> edges = all_pairs(n);
    std::vector<TF> t = pairwise_dist(x, edges, d);
    // for (size_t k = 0; k < t.size(); k++) {
    //     std::cout << edges[d*k] << ',' << edges[d*k + 1] << ',' <<  t[k] << std::endl;
    // }
    // std::cout << std::endl;
    sort_edges(edges, t);
    // for (size_t k = 0; k < t.size(); k++) {
    //     std::cout << edges[d*k] << ',' << edges[d*k + 1] << ',' <<  t[k] << std::endl;
    // }
    // std::cout << std::endl;
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "microseconds: " << duration.count() << std::endl;
    std::cout << "n edges =" << t.size() << std::endl;

    // std::cout << "entering Flag construction" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    // SimplicialComplex X;
    // Filtration<float, SimplicialComplex> F;
    auto [X, F] = FlagFiltration(edges, t, n, maxdim, TF(0));
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "microseconds: " << duration.count() << std::endl;

  // std::cout << "exiting Flag construction" << std::endl;

    for (size_t k = 0; k < X.maxdim() + 1; k++) {
        std::cout << "cells in dim " << k << " = " << X.ncells(k) << std::endl;
    }
    std::cout << F.maxdim() << std::endl;

    std::cout << "forming boundary 1" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    auto B1 = X.boundary_csc(1);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "microseconds: " << duration.count() << std::endl;
    std::cout << B1.nrow() << ',' << B1.ncol() << std::endl;
    // B1.print();

    std::cout << "Copying to ColumnMatrix" << std::endl;
    MT M1(B1);
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "microseconds: " << duration.count() << std::endl;
    std::cout << M1.nrow() << ',' << M1.ncol() << std::endl;
    // M1.print();

    std::cout << "Running Reduction Alg" << std::endl;
    auto p2c1 = reduce_matrix(M1);
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "microseconds: " << duration.count() << std::endl;
    std::cout << p2c1.size() << std::endl;
    // for (auto& [k, v] : p2c1) {
    //     std::cout << k << ':' << v << ',' << std::endl;
    // }
    // M1.print();

    std::cout << "forming boundary 2" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    auto B2 = X.boundary_csc(2);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "microseconds: " << duration.count() << std::endl;
    std::cout << B2.nrow() << ',' << B2.ncol() << std::endl;
    // B2.print();

    std::cout << "Copying to ColumnMatrix" << std::endl;
    MT M2(B2);
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "microseconds: " << duration.count() << std::endl;
    std::cout << M2.nrow() << ',' << M2.ncol() << std::endl;
    // M2.print();

    std::cout << "Running Reduction Alg" << std::endl;
    auto p2c2 = reduce_matrix(M2);
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "microseconds: " << duration.count() << std::endl;
    std::cout << p2c2.size() << std::endl;
    // for (auto& [k, v] : p2c2) {
    //     std::cout << k << ':' << v << ',' << std::endl;
    // }

    return 0;
}
