#include <bats.hpp>
#include <vector>
#include <iostream>

int main() {

    bats::LightSimplicialComplex<size_t> X(100,2);

    std::vector<size_t> s;

    s = {0}; std::cout << X.simplex_key(s) << "\n"; X.add(s);
    s = {1}; std::cout << X.simplex_key(s) << "\n"; X.add(s);
    s = {2}; std::cout << X.simplex_key(s) << "\n"; X.add(s);
    s = {0,1}; std::cout << X.simplex_key(s) << "\n"; X.add(s);
    s = {1,2}; std::cout << X.simplex_key(s) << "\n"; X.add(s);
    s = {0,2}; std::cout << X.simplex_key(s) << "\n"; X.add(s);
    s = {0,1,2}; std::cout << X.simplex_key(s) << "\n"; X.add(s);
    s = {1,2,3}; std::cout << X.simplex_key(s) << "\n"; X.add_recursive(s);


    {
        auto w = X.key_to_simplex(1, 1);
        std::cout << "{"; for (auto i : w) {std::cout << i << ","; } std::cout << "}\n";
    }


    std::cout << X.max_vertex(0, 0) << ", "
              << X.max_vertex(1, 0) << ", "
              << X.max_vertex(2, 0) << ", "
              << X.max_vertex(0, 1) << ", "
              << X.max_vertex(1, 1) << ", "
              << X.max_vertex(2, 1) << ", "
              << X.max_vertex(0, 2) << std::endl;



    X.print_summary();

    auto B1 = X.boundary_csc(1);
    B1.print();

    auto B2 = X.boundary_csc(2);
    B2.print();

    for (const auto& spx : X.get_simplices(1)) {
        std::cout << "{"; for (auto i : spx) {std::cout << i << ","; } std::cout << "}\n";
    }

    auto [ind, val] = X.boundary(2,0);
    auto b2 = SparseVector(ind, val);
    b2.print();

    return EXIT_SUCCESS;
}
