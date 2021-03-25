#include <bats.hpp>
#include <vector>
#include <iostream>

#define T ModP<int, 3>

int main() {

    UnivariatePolynomial p{T(1), T(2), T(2), T(1)};

    auto A = p.companion_matrix();
    A.print();

    auto pA = characteristic_matrix(A);
    auto F = smith_factorization(pA);

    std::cout << "R:  "; F.R.print();
    std::cout << "S:  "; F.S.print();
    std::cout << "C:  "; F.C.print();
    F.prod().print();

    pA.print();

    auto [R, S] = smith_rows(pA);
    R.print();

    auto B = generating_basis(R, A);
    B.print();

    auto B2 = RCF_basis(S, A, B);
    B2.print();

    auto B2i = inv(B2);
    B2i.print();

    (B2i * A * B2).print();

    return 0;
}
