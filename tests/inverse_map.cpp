#include <bats.hpp>
#include <util/io.hpp>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include <chrono>
#include <random>

using FT = ModP<int, 2>;

using CpxT = bats::LightSimplicialComplex<size_t, std::unordered_map<size_t, size_t>>;
// using CpxT = bats::SimplicialComplex;


int main() {

    
    /*
    First we build a simplicial complex which can be used to extend filtrations
    */
    CpxT X(3,2);
    
    std::vector<size_t> s;
    
    s = {2}; X.add(s);
    s = {0}; X.add(s);
    s = {1}; X.add(s);
    s = {0,1}; X.add(s);
    s = {1,2}; X.add(s);
    s = {0,2}; X.add(s);
    s = {0,1,2}; X.add(s);
    
    X.print_summary();
    
    /*
    Now let's exend a filtration
    */
    std::vector<double> f0 = {0.0, 0.1, 0.2};
    // lower star filtration
    
    auto [vals, inds] = lower_star_filtration(X, f0);
    for (auto& valsk: vals) {
        for (auto& v: valsk) {
            std::cout << v << ", ";
        }
        std::cout << "\n";
    }
    std::cout << "\ninds is" << std::endl;
    print_2D_vectors(inds);

    // inverse map of Rips Filtration
    auto Y = bats::sample_sphere<double>(3, 7);
	auto dist = bats::Euclidean(); // metric

    auto [F_Y, inds_Y] = bats::RipsFiltration_extension<CpxT>(Y, dist, 1.5, size_t(2));
    // auto [F_Y, inds_Y] = RipsFiltration_extension<CpxT>(Y, dist, 1.5, size_t(2));
    std::cout << "\nprint inds_Y" << std::endl;
    print_2D_vectors(inds_Y);
    print_filtration_info(F_Y);
}