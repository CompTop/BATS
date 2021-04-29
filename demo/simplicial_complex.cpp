#include <bats.hpp>
#include <vector>
#include <iostream>

int main() {

    bats::SimplicialComplex X, Y;

    X.add_recursive({0,1,2});
    X.print_cells();
    Y.add_recursive({1,2,3});
    Y.print_cells();

    auto I = intersection(X, Y);
    I.print_cells();

    auto U = simplicial_union(X, Y);
    U.print_cells();

    auto B = mayer_vietoris_boundary(X, U, I, 0);
    B.print();

    auto [ind, val] = X.boundary(1, 1);
    SparseVector(ind, val).print();

    X.boundary_csc(1).print();

    return 0;
}
