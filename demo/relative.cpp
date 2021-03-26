#include <iostream>
#include <vector>
#include <chrono>

#include <bats.hpp>


#define FT ModP<int, 3>
#define VT SparseVector<FT>
#define MT ColumnMatrix<VT>

int main(int argc, char* argv[]) {

    bats::SimplicialComplex X, A;

    X.add({0});
    X.add({1});
    X.add({2});
    X.add({0,1});
    X.add({0,2});
    X.add({1,2});
    X.add({0,1,2});

    X.print_summary();
    A.add({0});
    A.add({1});
    A.add({2});
    A.add({0,1});
    A.add({0,2});
    A.add({1,2});
    A.print_summary();

    auto inds = X.get_indices(A);
    std::vector<std::vector<size_t>> cinds;
    for (size_t d = 0; d < X.maxdim()+1; d++) {
        auto& ind = inds[d];
        cinds.emplace_back(bats::util::sorted_complement(ind, X.ncells(d)));
        std::cout << "indicies: ";
        for (auto i : ind) {
            std::cout << i << ',';
        }
        std::cout << std::endl;
        std::cout << "complement: ";
        for (auto i : bats::util::sorted_complement(ind, X.ncells(d))) {
            std::cout << i << ',';
        }
        std::cout << std::endl;
    }

    auto B = X.boundary_csc(1);
    B.print();

    auto R = B.submatrix(cinds[0], cinds[1]);
    R.print();

    auto C = __ChainComplex(X, A, FT());
    C[0].print();
    C[1].print();

    auto RC = bats::ReducedChainComplex(C);
    for (size_t d = 0; d <= RC.maxdim(); d++) {
        std::cout << "betti_" << d << ": " << RC.betti(d) << std::endl;
    }


    return 0;
}
