#pragma once

#include <cstddef>
#include <cassert>
#include <vector>
#include "abstract_matrix.hpp"
#include <util/sorted.hpp>

#include <iostream>
#include <iomanip>

// template over value type and index type
template <typename TV, typename TI=size_t>
class CSCMatrix : public AbstractMatrix
{
private:
    std::vector<TI> colptr;
    std::vector<TI> rowind;
    std::vector<TV> val;

    // dimensions
    size_t m; // number of rows
    size_t n; // number of columns

    // sort column j so row indices are in increasing order
    void sort_col(
        size_t j,
        std::vector<size_t> &perm,
        std::vector<TI> &tmp1,
        std::vector<TV> &tmp2
    ) {
        // get sortperm on row indices
        bats::util::fill_sortperm(
            rowind.cbegin() + colptr[j],
            rowind.cbegin() + colptr[j+1],
            perm
        );
        // apply permutation to row indices
        bats::util::apply_perm(
            rowind.data() + colptr[j],
            tmp1,
            perm
        );
        // apply permutation to row values
        bats::util::apply_perm(
            val.data() + colptr[j],
            tmp2,
            perm
        );
    }

    void sort_col(size_t j) {
        std::vector<size_t> perm;
        std::vector<TI> tmp1;
        std::vector<TV> tmp2;
        sort_col(j, perm, tmp1, tmp2);
    }

    // sort all columns
    void sort() {
        std::vector<size_t> perm;
        std::vector<TI> tmp1;
        std::vector<TV> tmp2;
        for (size_t j = 0; j < n; j++) {
            sort_col(j, perm, tmp1, tmp2);
        }
    }

public:

    CSCMatrix() : m(0), n(0) {};

    CSCMatrix(
        size_t m,
        size_t n,
        const std::vector<TI> &colptr,
        const std::vector<TI> &rowind,
        const std::vector<TV> &val
    ) : colptr(colptr), rowind(rowind), val(val), m(m), n(n) {
        sort();
    }

    CSCMatrix(
        const std::vector<TI> &colptr,
        const std::vector<TI> &rowind,
        const std::vector<TV> &val
    ) : colptr(colptr), rowind(rowind), val(val), n(colptr.size()-1) {
        sort();
    }

    TV getval(size_t i, size_t j) const {
        for (size_t k = colptr[j]; k < colptr[j+1]; k++) {
            if (rowind[k] == i) { return val[k]; }
        }
        return TV(0);
    }

    // access to data structures
    // std::vector<TI> colptr;
    // std::vector<TI> rowind;
    // std::vector<TV> val;
    //
    // // dimensions
    // size_t m; // number of rows
    // size_t n; // number of columns
    const std::vector<TI>& get_colptr() const { return colptr; }
    const std::vector<TI>& get_rowind() const { return rowind; }
    const std::vector<TV>& get_val() const {return val; }
    size_t nrow() const { return m; }
    size_t ncol() const { return n; }

    void print_size() const {
        std::cout << "[" << this << "] : " << m << " x " << n <<\
        " CSCMatrix" << std::endl;
    }

    void print(
        size_t rowmin,
        size_t rowmax,
        size_t colmin,
        size_t colmax
    ) const {
        // loop over rows
        for (size_t i = rowmin; i < rowmax; i++) {
            for (size_t j = colmin; j < colmax; j++) {
                std::cout << std::setw(3) << getval(i, j) << " ";
            }
            std::cout << std::endl;
        }
        return;
    }

    void print() const {
        print_size();
        print(0, m, 0, n);
        return;
    }

    // return submatrix indexed by rows and columns
    CSCMatrix submatrix(
        const std::vector<size_t> &rind,
        const std::vector<size_t> &cind
    ) const {
        CSCMatrix A;
        A.m = rind.size();
        A.n = cind.size();
        A.colptr.resize(A.n + 1);
        A.colptr[0] = 0;
        A.rowind.clear();
        A.val.clear();

        auto prow = bats::util::partial_perm(rind, nrow());

        size_t j = 0;
        // loop in column permutation order
        for ( size_t Mj : cind) {
            // loop over nzs in the jth column of M
            for (size_t Mp = colptr[Mj]; Mp < colptr[Mj+1]; Mp++) {
                if (prow[rowind[Mp]] != bats::NO_IND) {
                    // this item makes it to the block
                    A.rowind.emplace_back(prow[rowind[Mp]]);
                    A.val.emplace_back(val[Mp]);
                }
            }
            j++;
            A.colptr[j] = A.rowind.size();
        }
        return A;

    }

    // obtain block A from M from columns cind and
    // row permutation prow
    friend void block_select(
        const CSCMatrix &M,
        const std::vector<size_t> &cind,
        const std::vector<size_t> &prow,
        const size_t m,
        CSCMatrix &A
    ) {
        size_t n = cind.size(); // number of columns
        A.colptr.resize(n+1);
        A.colptr[0] = 0;
        // reserve enough space for 1 nonzero in each column.  This can always be expanded
        A.rowind.clear();
        A.rowind.reserve(m);
        A.val.clear();
        A.val.reserve(m);
        A.m = m;
        A.n = n;

        size_t j = 0;
        // loop in column permutation order
        for ( size_t Mj : cind) {
            // loop over nzs in the jth column of M
            for (size_t Mp = M.colptr[j]; Mp < M.colptr[j+1]; Mp++) {
                if (prow[M.rowind[Mp]] != bats::NO_IND) {
                    // this item makes it to the block
                    A.rowind.emplace_back(prow[M.rowind[Mp]]);
                    A.val.emplace_back(M.val[Mp]);
                }
            }
            j++;
            A.colptr[j] = A.rowind.size();
        }
    }

    /*
    operates on multiple row blocks
    WARNING: assumes that row blocks have unique indices
    */
    template <size_t N>
    friend void block_select(
        const CSCMatrix &M,
        const std::vector<size_t> &cind,
        const std::vector<size_t>* (&&prow)[N],
        const size_t (&&m)[N],
        CSCMatrix* (&&A)[N]
    ) {
        size_t n = cind.size(); // number of columns
        for (size_t i = 0; i < N; i++) {
            (*A[i]).colptr.resize(n+1);
            (*A[i]).colptr[0] = 0;
            // reserve enough space for 1 nonzero in each column.  This can always be expanded
            (*A[i]).rowind.clear();
            (*A[i]).rowind.reserve(m[i]);
            (*A[i]).val.clear();
            (*A[i]).val.reserve(m[i]);
            (*A[i]).m = m[i];
            (*A[i]).n = n;
        }

        size_t j = 0;
        // loop in column permutation order
        for ( size_t Mj : cind) {
            // loop over nzs in the jth column of M
            for (size_t Mp = M.colptr[Mj]; Mp < M.colptr[Mj+1]; Mp++) {
                for (size_t i = 0; i < N; i++) {
                    if ((*prow[i])[M.rowind[Mp]] != bats::NO_IND) {
                        // this item makes it to the block
                        (*A[i]).rowind.emplace_back((*prow[i])[M.rowind[Mp]]);
                        (*A[i]).val.emplace_back(M.val[Mp]);
                        break;
                    }
                }
            }
            j++;
            for (size_t i = 0; i < N; i++) {
                (*A[i]).colptr[j] = (*A[i]).rowind.size();
            }
        }
        for (size_t i = 0; i < N; i++) {
            (*A[i]).sort();
        }
    }

    // gemm implementation
    // TV can be arbitrary ring
    // C = A * B
    // C is over-written
    friend void gemm(const CSCMatrix &A, const CSCMatrix &B, CSCMatrix &C) {

        assert (A.n == B.m);
        size_t m = A.m;
        size_t n = B.n;

        C.m = m;
        C.n = n;
        C.colptr.resize(n+1);
        C.colptr[0] = 0;
        C.rowind.clear();
        C.val.clear();

        // loop over columns of C
        // C[:,j] = A[:,k] * B[k, j]
        // allocations for temp
        std::vector<size_t> perm;
        std::vector<TI> tmp1;
        std::vector<TV> tmp2;
        for (size_t j = 0; j < n; j++) {
            // loop over nonzeros in jth column of B
            for (size_t pBj = B.colptr[j]; pBj < B.colptr[j+1]; pBj++) {
                TI k = B.rowind[pBj];
                TV v = B.val[pBj];
                // push onto jth column of C
                for (size_t pAk = A.colptr[k]; pAk < A.colptr[k+1]; pAk++) {
                    C.rowind.emplace_back(A.rowind[pAk]);
                    C.val.emplace_back(v * A.val[pAk]);
                }
            }
            // sort & sum reduce column C[j]
            bats::util::sort_sum_reduce(C.rowind, C.val, C.colptr[j], perm, tmp1, tmp2);
            C.colptr[j+1] = C.rowind.size();
        }

    }

    CSCMatrix operator*(const CSCMatrix &other) const {
        CSCMatrix C;
        gemm(*this, other, C);
        return C;
    }

    // gemm implementation
    // TV can be arbitrary ring
    // C = A + B
    // C is over-written
    friend void sum(const CSCMatrix &A, const CSCMatrix &B, CSCMatrix &C) {

        assert (A.n == B.n);
        assert (A.m == B.m);
        size_t n = B.n;
        size_t m = A.m;

        C.m = m;
        C.n = n;
        C.colptr.resize(n+1);
        C.colptr[0] = 0;
        C.rowind.clear();
        C.val.clear();

        // loop over columns of C
        // C[:,j] = A[:,k] * B[k, j]
        // allocations for temp
        std::vector<size_t> perm;
        std::vector<TI> tmp1;
        std::vector<TV> tmp2;
        for (size_t j = 0; j < n; j++) {
            // loop over noneros in jth column of B
            for (size_t pBj = B.colptr[j]; pBj < B.colptr[j+1]; pBj++) {
                C.rowind.emplace_back(B.rowind[pBj]);
                C.val.emplace_back(B.val[pBj]);
            }
            for (size_t pAj = A.colptr[j]; pAj < A.colptr[j+1]; pAj++) {
                C.rowind.emplace_back(A.rowind[pAj]);
                C.val.emplace_back(A.val[pAj]);
            }

            // sort & sum reduce column C[j]
            bats::util::sort_sum_reduce(C.rowind, C.val, C.colptr[j], perm, tmp1, tmp2);
            C.colptr[j+1] = C.rowind.size();
        }

    }

    // c = alpha * A + B
    friend void sum(const TV& alpha, const CSCMatrix &A, const CSCMatrix &B, CSCMatrix &C) {

        assert (A.n == B.n);
        assert (A.m == B.m);
        size_t n = B.n;
        size_t m = A.m;

        C.m = m;
        C.n = n;
        C.colptr.resize(n+1);
        C.colptr[0] = 0;
        C.rowind.clear();
        C.val.clear();

        // loop over columns of C
        // C[:,j] = alpha * A[:,j] + B[:,j]
        // allocations for temp
        std::vector<size_t> perm;
        std::vector<TI> tmp1;
        std::vector<TV> tmp2;
        for (size_t j = 0; j < n; j++) {
            // loop over noneros in jth column of B
            for (size_t pBj = B.colptr[j]; pBj < B.colptr[j+1]; pBj++) {
                C.rowind.emplace_back(B.rowind[pBj]);
                C.val.emplace_back(B.val[pBj]);
            }
            for (size_t pAj = A.colptr[j]; pAj < A.colptr[j+1]; pAj++) {
                C.rowind.emplace_back(A.rowind[pAj]);
                C.val.emplace_back(alpha * A.val[pAj]);
            }

            // sort & sum reduce column C[j]
            bats::util::sort_sum_reduce(C.rowind, C.val, C.colptr[j], perm, tmp1, tmp2);
            C.colptr[j+1] = C.rowind.size();
        }

    }

    // upper triangular solve implementation
    // TV can be arbitrary ring, assuming diagonal has units
    // C = A \ B
    // C is over-written
    // WARNING: assumes A is upper-triangular, with units on diagonal
    // TODO: this may be faster with a symbolic factorization first
    // or using temporary arrays for storage as with sorted sparse vectors
    friend void trilu(const CSCMatrix &A, const CSCMatrix &B, CSCMatrix &C) {

        assert (A.n == B.m);
        assert (A.n == A.m); // square = invertible
        size_t m = A.m;
        size_t n = B.n;

        C.m = m;
        C.n = n;
        C.colptr.resize(n+1);
        C.colptr[0] = 0;
        C.rowind.clear();
        C.val.clear();

        // loop over columns of C
        // C[:,j] = A \ B[:,j]
        // allocations for temp
        std::vector<size_t> perm;
        std::vector<TI> tmp1;
        std::vector<TV> tmp2;
        for (size_t j = 0; j < n; j++) {
            // y = U \ x:
            // y = x
            // for i = m:-1:1:
            //      y[i] = y[i] / U[i,i]
            //      y[:i-1] += y[i]*U[:i-1,i]
            // first copy B[:,j] into C[:,j]
            size_t offsetCj = C.colptr[j];
            for (size_t pBj = B.colptr[j]; pBj < B.colptr[j+1]; pBj++) {
                C.rowind.emplace_back(B.rowind[pBj]);
                C.val.emplace_back(B.val[pBj]);
            }

            size_t lastnz = find_sorted_lt(C.rowind.cbegin() + offsetCj, C.rowind.cend(), m);

            while (lastnz != bats::NO_IND) {

                TI i = *(C.rowind.cbegin() + offsetCj + lastnz);
                TV yi = *(C.val.cbegin() + offsetCj + lastnz);
                *(C.val.begin() + offsetCj + lastnz) /= A.val[A.colptr[i+1] - 1]; // y[i] = y[i] / U[i,i]
                // y[:i-1] += y[i]*U[:i-1, i]
                size_t offsetAi = A.colptr[i];

                auto Cr = C.rowind.begin() + offsetCj;
                auto Cv = C.val.begin() + offsetCj;
                auto Ar = A.rowind.cbegin() + offsetAi;
                auto Av = A.val.cbegin() + offsetAi;
                while (*Cr < i || *Ar < i) {

                    if (*Cr == *Ar) {
                        *Cv -= *Av * yi;
                        ++Cr; ++Cv; ++Ar; ++Av;
                    } else if (*Cr < *Ar) {
                        // no update
                        ++Cr; ++Cv;
                    } else { // *Cr > *Ar, and we need to insert
                        Cr = C.rowind.insert(Cr, *Ar) + 1;
                        Cv = C.val.insert(Cv, -(*Av)*yi) + 1;
                        ++Ar; ++Av;
                    }
                }

                lastnz = find_sorted_lt(C.rowind.cbegin() + offsetCj, C.rowind.cend(), i);
            }
            C.colptr[j+1] = C.rowind.size();
        }
    }

};
