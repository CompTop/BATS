#include <bats.hpp>
#include <vector>
#include <iostream>
#include <chrono>

// TODO: make this a static method in cubical_complex
// generates cube on n^3 vertices.
// should have n > 1
// bats::CubicalComplex generate_cube(
//     const size_t n
// ) {
//     if (n < 2) {throw std::runtime_error("Cube must have at least 2 verices in each dimension.");}
//     bats::CubicalComplex X(3);
//     for (size_t i = 0; i < n-1; ++i) {
//         for (size_t j = 0; j < n-1; ++j) {
//             for (size_t k = 0; k < n-1; ++k) {
//                 X.add_recursive({i, i+1, j, j+1, k, k+1});
//             }
//         }
//     }
//     return X;
// }

int main() {

    bats::RightFiltration<bats::SimplicialComplex> F;

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
            bats::extra_reduction_flag()
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
            bats::extra_reduction_flag()
        );

        for (auto& pk : ps) {
            for (auto p : pk) {
                if (p.length() > 0)
                    std::cout << p << std::endl;
            }
        }
    }

    // test for generating cubes
    for (auto n : {3, 5, 9, 17, 33, 65, 129}) {
        auto start = std::chrono::steady_clock::now();
        bats::CubicalComplex X = bats::CubicalComplex::generate_cube(n);
        auto end = std::chrono::steady_clock::now();
        std::cout << "\nBuild cube on " << n << "^3 vertices: "
            << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
            << "\u03BCs" << std::endl;
        X.print_summary();
    }




    return EXIT_SUCCESS;
}
