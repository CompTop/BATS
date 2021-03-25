#pragma once
#include <tuple>
#include "col_matrix.hpp"
#include <iostream>
#include "polynomial.hpp"


int sgn(const int &a) {
    return (a < 0) ? -1 : (a > 0) ? 1 : 0;
}

/*
compute extended euclidean algorithm.
given a and b, compute x and y so that
ax + by = gcd(a,b)
returns x, y, gcd(a,b), a/gcd(a,b), b/gcd(a,b)
*/
template <typename PID>
std::tuple<PID, PID, PID, PID, PID> extended_euclidean(const PID& a, const PID& b) {
    // TODO: check if a or b is 0.
    // if (a == 0) {
    //
    //     return std::tie(old_s, old_t, old_r, t, s); // x, y, gcd, a/gcd, b/gcd
    // } else if (b == 0) {
    //
    // }

    PID s(0), old_s(1);
    PID t(1), old_t(0);
    PID r(b), old_r(a);

    int i = 0;
    while (r != 0) {
        i++;
        PID quotient = old_r / r; // division without remainder
        PID tmp = old_r;
        old_r = r; // remember remainder
        r = tmp - quotient * r; // new remainder

        tmp = old_s;
        old_s = s;
        s = tmp - quotient * s;

        tmp = old_t;
        old_t = t;
        t = tmp - quotient * t;
    }

    // correct sign of division by gcd
    if ((i & 1)) {t = -t;}
    if (!(i & 1)) {s = -s;}

    return std::tie(old_s, old_t, old_r, t, s); // x, y, gcd, a/gcd, b/gcd
}

/*
Compute the Smith normal form of a matrix A
Do this in-place
note that row operations are generally inefficent
*/
template <typename TC>
void smith_normal_form(ColumnMatrix<TC> &A) {
    using T = typename TC::val_type;

    size_t rank = 0;
    // loop over columns
    for (size_t j = 0; j < A.ncol(); j++) {

        // STEP 1: put pivot in top-left corner
        if (A[j].nnz() == 0) {
            // swap columns so there is a non-zero
            for (size_t j2 = j+1; j2 < A.ncol(); j2++) {
                if (A[j2].nnz() > 0) {
                    A.swap_cols(j, j2);
                    break;
                }
            }
            // if A[j] is still zero, we're done
            if (A[j].nnz() == 0) break;
        }
        rank++; // increment the rank

        // now we have at least one non-zero in column
        // pivot non-zero to pivot location
        auto pind = A[j].nzbegin()->ind;
        if (pind != j) {
            //swap rows to put in in pivot location
            A.swap_rows(j, pind);
        }


        // STEP 2: improve pivot
        // and eliminate entries in row and column
        bool continue_reduction = true;
        while (continue_reduction) {
            // first, we loop over columns to the right
            auto p = A[j].nzbegin();
            for (size_t j2 = j+1; j2 < A.ncol(); j2++) {
                if (A[j2].nnz() > 0 && A[j2].firstnz().ind == j) {
                    auto [x, y, g, w, z] = extended_euclidean(p->val, A[j2].firstnz().val);
                    A.mix_cols(j, j2, x, -z, y, w);
                }
            }

            // now we loop over rows
            continue_reduction = false;
            p = A[j].nzbegin();
            auto rp = A[j].nzend() - 1;
            while (rp != p) {
                auto [x, y, g, w, z] = extended_euclidean(p->val, rp->val);
                if (!continue_reduction && p->val == g) {
                    A[j].mix_rows(p->ind, rp->ind, x, y, -z, w);
                } else {
                    // we're going to mix a row with the pivot row.
                    continue_reduction = true;
                    A.mix_rows(p->ind, rp->ind, x, y, -z, w);
                }
                rp = A[j].nzend() - 1;
            }


            // if continue_reduction is true, we mixed a row with the pivot row
            // so need to continue while-loop
        }

    }

    // finalize by ensuring divisibility condition
    // these loops take O(rank^2)
    for (size_t j = 0; j < rank; j++) {
        // check that j divides everything after it
        for (size_t j2 = j+1; j2 < rank; j2++) {
            auto [x, y, g, w, z] = extended_euclidean(A(j,j), A(j2,j2));
            if (A(j,j) != g) {
                // we need to put gcd in A(j,j)
                A.mix_cols(j, j2, x, -z, y, w);
                auto yz = y*z;
                A.mix_rows(j, j2, T(1), T(1), -yz, T(1) - yz);
            }
        }
        // if working over polynomials, make monic
        if constexpr (is_UnivariatePolynomial<T>::value) {
            // std::cout << "UnivariatePolynomial detected!" << std::endl;
            auto coeff = A(j,j).leading_coeff().inv();
            A[j].scale_inplace(coeff); // only one non-zero in row, so just scale column
            // R[j].scale_inplace(coeff.inv()); // corresponding inverse op
        }
    }

    return;
}

// factorization for Smith normal form
// A = R * S * C
// S is smith normal form of A
template <class TC>
struct SmithFact {
    // struct to hold matrices in factorization
    ColumnMatrix<TC> R;
    ColumnMatrix<TC> S;
    ColumnMatrix<TC> C;

    inline ColumnMatrix<TC> prod() const {
        return R * S * C;
    }


};

template <class TC>
SmithFact<TC> smith_factorization(const ColumnMatrix<TC> &A) {
    // Smith factorization of matrix A

    using MatT = ColumnMatrix<TC>;
    using T = typename TC::val_type;

    size_t m = A.nrow();
    size_t n = A.ncol();

    MatT R = MatT::identity(m); // row operations
    MatT S = A;
    MatT C = MatT::identity(n); // column operations

    size_t rank = 0;
    // loop over columns
    for (size_t j = 0; j < S.ncol(); j++) {

        // STEP 1: put pivot in top-left corner
        if (S[j].nnz() == 0) {
            // swap columns so there is a non-zero
            for (size_t j2 = j+1; j2 < S.ncol(); j2++) {
                if (S[j2].nnz() > 0) {
                    S.swap_cols(j, j2);
                    C.swap_rows(j, j2);
                    break;
                }
            }
            // if A[j] is still zero, we're done
            if (S[j].nnz() == 0) break;
        }
        rank++; // increment the rank

        // now we have at least one non-zero in column
        // pivot non-zero to pivot location
        auto pind = S[j].nzbegin()->ind;
        if (pind != j) {
            //swap rows to put in in pivot location
            S.swap_rows(j, pind);
            R.swap_cols(j, pind);
        }


        // STEP 2: improve pivot
        // and eliminate entries in row and column
        bool continue_reduction = true;
        while (continue_reduction) {
            // first, we loop over columns to the right
            auto p = S[j].nzbegin();
            for (size_t j2 = j+1; j2 < S.ncol(); j2++) {
                if (S[j2].nnz() > 0 && S[j2].firstnz().ind == j) {
                    auto [x, y, g, w, z] = extended_euclidean(p->val, S[j2].firstnz().val);
                    S.mix_cols(j, j2, x, -z, y, w);
                    C.mix_rows(j, j2, w, z, -y, x);
                }
            }

            // now we loop over rows
            continue_reduction = false;
            auto rp = S[j].nzend() - 1;
            p = S[j].nzbegin();
            while (rp != p) {
                auto [x, y, g, w, z] = extended_euclidean(p->val, rp->val);
                if (!continue_reduction && p->val == g) {
                    S[j].mix_rows(p->ind, rp->ind, x, y, -z, w);
                    R.mix_cols(p->ind, rp->ind, w, -y, z, x);
                } else {
                    // we're going to mix a row with the pivot row.
                    continue_reduction = true;
                    S.mix_rows(p->ind, rp->ind, x, y, -z, w);
                    R.mix_cols(p->ind, rp->ind, w, -y, z, x);
                }
                rp = S[j].nzend() - 1;
            }
            // if continue_reduction is true, we mixed a row with the pivot row
            // so need to continue while-loop
        }

    }

    // finalize by ensuring divisibility condition
    // these loops take O(rank^2)
    for (size_t j = 0; j < rank; j++) {
        // check that j divides everything after it
        for (size_t j2 = j+1; j2 < rank; j2++) {
            auto [x, y, g, w, z] = extended_euclidean(S(j,j), S(j2,j2));
            if (S(j,j) != g) {
                // we need to put gcd in A(j,j)
                S.mix_cols(j, j2, x, -z, y, w);
                C.mix_rows(j, j2, w, z, -y, x);

                auto yz = y*z;
                S.mix_rows(j, j2, T(1), T(1), -yz, T(1) - yz);
                R.mix_cols(j, j2, T(1) - yz, -T(1), yz, T(1));
            }
        }
        // if working over polynomials, make monic
        if constexpr (is_UnivariatePolynomial<T>::value) {
            // std::cout << "UnivariatePolynomial detected!" << std::endl;
            auto coeff = S(j,j).leading_coeff().inv();
            S[j].scale_inplace(coeff); // only one non-zero in row, so just scale column
            R[j].scale_inplace(coeff.inv()); // corresponding inverse op
        }
    }

    return SmithFact<TC>{R, S, C};
}

// smith normal form of matrix A
// only keep track of inverse row operations
// useful for RCF
template <class TC>
std::tuple<ColumnMatrix<TC>, ColumnMatrix<TC>> smith_rows(const ColumnMatrix<TC> &A) {

    using MatT = ColumnMatrix<TC>;
    using T = typename TC::val_type;

    size_t m = A.nrow();

    MatT R = MatT::identity(m); // row operations
    MatT S = A;

    size_t rank = 0;
    // loop over columns
    for (size_t j = 0; j < S.ncol(); j++) {

        // STEP 1: put pivot in top-left corner
        if (S[j].nnz() == 0) {
            // swap columns so there is a non-zero
            for (size_t j2 = j+1; j2 < S.ncol(); j2++) {
                if (S[j2].nnz() > 0) {
                    S.swap_cols(j, j2);
                    break;
                }
            }
            // if A[j] is still zero, we're done
            if (S[j].nnz() == 0) break;
        }
        rank++; // increment the rank

        // now we have at least one non-zero in column
        // pivot non-zero to pivot location
        auto pind = S[j].nzbegin()->ind;
        if (pind != j) {
            //swap rows to put in in pivot location
            S.swap_rows(j, pind);
            R.swap_cols(j, pind);
        }


        // STEP 2: improve pivot
        // and eliminate entries in row and column
        bool continue_reduction = true;
        while (continue_reduction) {
            // first, we loop over columns to the right
            auto p = S[j].nzbegin();
            for (size_t j2 = j+1; j2 < S.ncol(); j2++) {
                if (S[j2].nnz() > 0 && S[j2].firstnz().ind == j) {
                    auto [x, y, g, w, z] = extended_euclidean(p->val, S[j2].firstnz().val);
                    S.mix_cols(j, j2, x, -z, y, w);
                }
            }

            // now we loop over rows
            continue_reduction = false;
            auto rp = S[j].nzend() - 1;
            p = S[j].nzbegin();
            while (rp != p) {
                auto [x, y, g, w, z] = extended_euclidean(p->val, rp->val);
                if (!continue_reduction && p->val == g) {
                    S[j].mix_rows(p->ind, rp->ind, x, y, -z, w);
                    R.mix_cols(p->ind, rp->ind, w, -y, z, x);
                } else {
                    // we're going to mix a row with the pivot row.
                    continue_reduction = true;
                    S.mix_rows(p->ind, rp->ind, x, y, -z, w);
                    R.mix_cols(p->ind, rp->ind, w, -y, z, x);
                }
                rp = S[j].nzend() - 1;
            }
            // if continue_reduction is true, we mixed a row with the pivot row
            // so need to continue while-loop
        }

    }

    // finalize by ensuring divisibility condition
    // these loops take O(rank^2)
    for (size_t j = 0; j < rank; j++) {
        // check that j divides everything after it
        for (size_t j2 = j+1; j2 < rank; j2++) {
            auto [x, y, g, w, z] = extended_euclidean(S(j,j), S(j2,j2));
            if (S(j,j) != g) {
                // we need to put gcd in A(j,j)
                S.mix_cols(j, j2, x, -z, y, w);

                auto yz = y*z;
                S.mix_rows(j, j2, T(1), T(1), -yz, T(1) - yz);
                R.mix_cols(j, j2, T(1) - yz, -T(1), yz, T(1));
            }
        }
        // if working over polynomials, make monic
        if constexpr (is_UnivariatePolynomial<T>::value) {
            // std::cout << "UnivariatePolynomial detected!" << std::endl;
            auto coeff = S(j,j).leading_coeff().inv();
            S[j].scale_inplace(coeff); // only one non-zero in row, so just scale column
            R[j].scale_inplace(coeff.inv()); // corresponding inverse op
        }
    }

    return std::tie(R, S);
}
