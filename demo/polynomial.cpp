#include <bats.hpp>
#include <vector>
#include <iostream>

#define T ModP<int, 3>

int main() {

    UnivariatePolynomial p{T(2), T(1)};
    std::cout << p << std::endl;

    UnivariatePolynomial p2(T(1));
    std::cout << p2 << std::endl;

    auto p3 = -p;
    std::cout << p3 << std::endl;

    std::cout << (p2 + p) << std::endl;
    std::cout << (p2 - p) << std::endl;

    std::cout << (p * p) << std::endl;

    std::cout << UnivariatePolynomial<T>::monomial(2) << std::endl;
    std::cout << UnivariatePolynomial<T>::monomial(2, T(2)) << std::endl;

    auto pp = p * p;
    std::cout << "poly: " << p << std::endl;
    std::cout << "square: " << pp << std::endl;
    auto [q, r] = pp.divrem(p);
    std::cout << "q: " << q << std::endl;
    std::cout << "r: " << r << std::endl;

    std::cout << "gcd: " << pp.gcd(p) << std::endl;

    auto [x,y,g,ax,by] = extended_euclidean(pp, p);

    std::cout << "\nextended euclidean" << std::endl;
    std::cout << x << std::endl;
    std::cout << y << std::endl;
    std::cout << g << std::endl;

    std::cout << "\ntesting companion matrices/Smith form" << std::endl;
    p = UnivariatePolynomial{T(1), T(2), T(2), T(1)};
    std::cout << p << std::endl;

    auto A = p.companion_matrix();
    A.print();

    auto pA = characteristic_matrix(A);
    pA.print();

    smith_normal_form(pA);
    pA.print();

    pA = characteristic_matrix(A);
    auto F = smith_factorization(pA);

    std::cout << "R:  "; F.R.print();
    std::cout << "S:  "; F.S.print();
    std::cout << "C:  "; F.C.print();
    F.prod().print();

    pA.print();

    return 0;
}
