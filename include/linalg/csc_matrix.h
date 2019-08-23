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

    CSCMatrix(
        size_t m,
        size_t n,
        std::vector<TI> &colptr,
        std::vector<TI> &rowind,
        std::vector<TV> &val
    ) : colptr(colptr), rowind(rowind), val(val), m(m), n(n) {
        //sort();
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

};
