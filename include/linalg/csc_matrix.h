#pragma once

#include <cstddef>
#include <vector>
#include "abstract_matrix.h"
#include <util/sorted.h>

#include <iostream>

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
        fill_sortperm(
            rowind.cbegin() + colptr[j],
            rowind.cbegin() + colptr[j+1],
            perm
        );
        // apply permutation to row indices
        apply_perm(
            rowind.data() + colptr[j],
            tmp1,
            perm
        );
        // apply permutation to row values
        apply_perm(
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
        std::vector<TI> &colptr,
        std::vector<TI> &rowind,
        std::vector<TV> &val
    ) : colptr(colptr), rowind(rowind), val(val), m(m), n(n) {
        sort();
    }

    CSCMatrix(
        std::vector<TI> &colptr,
        std::vector<TI> &rowind,
        std::vector<TV> &val
    ) : colptr(colptr), rowind(rowind), val(val), n(colptr.size()-1) {
        sort();
    }

    TV getval(size_t i, size_t j) const {
        for (size_t k = colptr[j]; k < colptr[j+1]; k++) {
            if (rowind[k] == i) { return val[k]; }
        }
        return TV(0);
    }


    void print() {
        std::cout << "[" << this << "] : " << m << " x " << n <<\
        " CSCMatrix" << std::endl;
        // loop over rows
        for (size_t i = 0; i < m; i++) {
            for (size_t j = 0; j < n; j++) {
                std::cout << std::setw(3) << getval(i, j) << " ";
            }
            std::cout << std::endl;
        }
        return;
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
            for (size_t Mp = M.colptr[j]; Mp < M.colptr[j+1]; Mp++) {
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
    }

};
