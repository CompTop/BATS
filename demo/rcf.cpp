#include <bats.hpp>
#include <vector>
#include <iostream>

// compute rational canonical form of a matrix A
// A = X * R * X^{-1}
// returns tuple(X, R, Xinv)
template <typename T>
auto rational_canonical_form(
    const ColumnMatrix<SparseVector<T>> &A
) {

    std::cout << "forming characteristic matrix" << std::endl;
    auto pA = characteristic_matrix(A);
    pA.print();
    std::cout << "computing Smith factorization" << std::endl;
    auto F = smith_factorization(pA);
    F.S.print();

    std::cout << "computing Smith rows..." << std::endl;
    auto [R, S] = smith_rows(pA);
    std::cout << "computing generating basis..." << std::endl;
    auto B = generating_basis(R, A);

    std::cout << "computing RCF basis..." << std::endl;
    auto X = RCF_basis(S, A, B);
    X.print();
    std::cout << "inverting..." << std::endl;
    auto Xi = inv(X);

    std::cout << "computing product..." << std::endl;
    auto RCF = Xi * A * X;
    std::cout << "returning..." << std::endl;
    return std::tuple(X, RCF, Xi);
}



int main() {
    using T = ModP<int, 3>;

    UnivariatePolynomial p{T(0), T(0), T(1), T(1), T(1)};

    auto A = p.companion_matrix();
    A.print();

    auto [X, R, Xi] = rational_canonical_form(A);
    //
    (Xi * A * X).print();

    return 0;
}
