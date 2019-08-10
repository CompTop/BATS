#pragma once

#include <cstddef>
#include <vector>
#include "abstract_matrix.h"

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

public:

    CSCMatrix(
        size_t m,
        size_t n,
        std::vector<TI> &colptr,
        std::vector<TI> &rowind,
        std::vector<TV> &val
    ) : colptr(colptr), rowind(rowind), val(val), m(m), n(n) {}

    CSCMatrix(
        std::vector<TI> &colptr,
        std::vector<TI> &rowind,
        std::vector<TV> &val
    ) : colptr(colptr), rowind(rowind), val(val), n(colptr.size()-1) {}

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
