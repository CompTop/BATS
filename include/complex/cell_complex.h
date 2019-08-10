#pragma once

#include <cstddef>
#include <limits>
#include <vector>
#include <algorithm>
#include <iterator>
#include "abstract_complex.h"
#include <linalg/csc_matrix.h>


class CellComplex : public AbstractComplex
{
private:
    std::vector<size_t> _ncells;

    // boundaries are stored in CSC format
    // only store boundaries > 0
    std::vector<std::vector<size_t>> ptr; // pointer to cell start
    std::vector<std::vector<size_t>> bdr; // boundaries of cells
    std::vector<std::vector<int>> coeff; // coefficients of boundaries

    // reserve for maxdim dimensional cells
    void reserve(size_t maxdim) {
        if ( _ncells.size() < maxdim + 1 ) { _ncells.resize(maxdim+1, 0); }
        if ( ptr.size() < maxdim ) { ptr.resize(maxdim, std::vector<size_t>(1, 0)); }
        if ( bdr.size() < maxdim ) { bdr.resize(maxdim); }
        if ( coeff.size() < maxdim ) { coeff.resize(maxdim); }
        return;
    }
    // reserve for k cells in dimension dim
    // assume dim > 0
    void reserve(size_t dim, size_t k) {
        // first reserve for dimension
        reserve(dim);
        dim--;
        if ( ptr[dim].size() < k ) { ptr[dim].reserve(k); }
        if ( bdr[dim].size() < k ) { bdr[dim].reserve(k); }
        if ( coeff[dim].size() < k ) { coeff[dim].reserve(k); }
        return;
    }

    // return the identifier for the added cell
    size_t add_unsafe(
        const std::vector<size_t> &b,
        const std::vector<int> &c,
        size_t k
    ) {
        // copy boundary
        std::copy(
            b.cbegin(),
            b.cend(),
            std::back_inserter(bdr[k-1])
        );
        // copy coefficients
        std::copy(
            c.cbegin(),
            c.cend(),
            std::back_inserter(coeff[k-1])
        );
        // update pointer
        ptr[k-1].emplace_back(bdr[k-1].size());
        return _ncells[k]++;
    }

    // add vertex, return index
    inline size_t _add_vertex() {
        return _ncells[0]++;
    }

    // add k vertices to complex, return total number at end
    inline size_t _add_vertices(size_t k) {
        _ncells[0] += k;
        return _ncells[0];
    }

    size_t add_safe(
        const std::vector<size_t> &b,
        const std::vector<int> &c,
        size_t k
    ) {
        if (k == 0) { return add_vertex(); }
        // first check that b and c are same length
        assert(b.size() == c.size());
        // now check that we have reserved memory for dimension k
        reserve(k);

        // we're clear to add
        return add_unsafe(b, c, k);
    }


    // get CSC integer matrix boundary restricted to indices
public:

    CellComplex() { reserve(0); }

    CellComplex(size_t maxdim) { reserve(maxdim); }

    // maximum dimension of cell
    inline size_t maxdim() const { return _ncells.size() - 1; }
    // number of cells in dimension k
    inline size_t ncells(size_t k) const { return _ncells[k]; }

    // add cell with boundary b, coeffeicents c to dimension k
    inline size_t add(
        const std::vector<size_t> &b,
        const std::vector<int> &c,
        size_t k
    ) {
        return add_safe(b, c, k);
    }

    inline size_t add_vertex() { return _add_vertex(); }
    inline size_t add_vertices(size_t k) { return _add_vertices(k); }

    // get CSC integer matrix boundary in dimension dim
    inline CSCMatrix<int, size_t> boundary_csc(size_t dim) {
        return CSCMatrix<int, size_t>(
            _ncells[dim-1],
            _ncells[dim],
            ptr[dim-1],
            bdr[dim-1],
            coeff[dim-1]
        );
    }

    // create boundary for dimension k
    // template over matrix type
    // can specialize for F2
    void boundary(size_t k);

    // create boundary for dimension k restricted to indices
    // template over matrix type
    // can specialize for F2
    void boundary(
        size_t k,
        std::vector<size_t> row,
        std::vector<size_t> col
    );

};
