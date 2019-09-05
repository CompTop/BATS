#pragma once

#include <cstddef>
#include <limits>
#include <vector>
#include <algorithm>
#include <iterator>

#include <util/common.h>
#include "abstract_complex.h"
#include <linalg/sparse_vector.h>
#include <linalg/col_matrix.h>
#include <linalg/csc_matrix.h>
#include <correspondence/correspondence.h>
#include <util/sorted.h>
#include <util/trie.h>
#include <iostream>
#include <map>
#include <filtration/filtration.h>
#include <morse/pairing.h>



#include <util/simplex.h> // for hash function


/*
simplicial complex implementation
TODO: if forming full clique complex
can give all simplices unique location for lookup
*/
class SimplicialComplex
{
private:
    // maybe make template parameter for container?
    // using spx_map = std::unordered_map<std::vector<size_t>, size_t, SimplexHasher>;
    // using spx_map = std::map<std::vector<size_t>, size_t>;
    using spx_map = SparseTrie<size_t, size_t>;

    // container for simplices of a fixed dimension
    std::vector<std::vector<size_t>> faces;
    std::vector<std::vector<int>> coeff;

    // keeps track of how many cells are in given dimension
    std::vector<size_t> _ncells;
    // keeps track of reserved capacity
    std::vector<size_t> _reserved;

    // map to find simplices when forming boundary
    spx_map spx_to_idx;


    // returns index of simplex
    // TODO: find a way to do this with a single traversal of spx_to_idx
    size_t find_idx(std::vector<size_t> &s) {
        if (spx_to_idx.count(s)) {
            return spx_to_idx[s];
        }
        return bats::NO_IND;
    }

    // reserve for maxdim dimensional simplices
    void reserve(size_t maxdim) {
        if ( _ncells.size() < maxdim + 1 ) { _ncells.resize(maxdim+1, 0); }
        if ( _reserved.size() < maxdim + 1 ) { _reserved.resize(maxdim+1, 0); }
        if ( faces.size() < maxdim ) { faces.resize(maxdim); }
        if ( coeff.size() < maxdim ) { coeff.resize(maxdim); }
        return;
    }

    // assume dim > 0
    void reserve(size_t dim, size_t k) {
        reserve(dim);
        _reserved[dim] = std::max(_reserved[dim], k);
        if (dim == 0) { return; }
        if ( faces[dim-1].capacity() < k * (dim + 1) ) {
          faces[dim-1].reserve(k * (dim + 1));
          coeff[dim-1].reserve(k * (dim + 1));
        }
        return;
    }

    // // add vertex, return index
    // inline size_t _add_vertex() {
    //     return _ncells[0]++;
    // }
    //
    // // add k vertices to complex, return total number at end
    // inline size_t _add_vertices(size_t k) {
    //     _ncells[0] += k;
    //     return _ncells[0];
    // }

    // adds a simplex to the complex without doing any checks
    // assume dim > 0
    cell_ind _add_unsafe(std::vector<size_t> &s) {
        size_t dim = s.size() - 1;

        // get index of simplex
        size_t ind = _ncells[dim]++;
        // add to map
        spx_to_idx.emplace(s, ind);

        // determine faces
        if (dim > 0){
            int c = -1;
            // loop over faces in lexicographical order
            std::vector<size_t> face;
            //face.reserve(dim);
            size_t spx_len = dim + 1;
            for (size_t k = 0; k < spx_len; k++) {
                size_t k2 = dim-k; // index to skip
                face.clear();
                for (size_t j = 0; j < k2; j++) {
                    face.emplace_back(s[j]);
                }
                for (size_t j = k2+1; j < spx_len; j++) {
                    face.emplace_back(s[j]);
                }

                faces[dim-1].emplace_back(find_idx(face));
                coeff[dim-1].emplace_back(c);
                c = -c;
            }
        }

        return cell_ind(dim, ind);
    }

    // add simplex to complex with appropriate checks
    cell_ind add_safe(std::vector<size_t> &s) {
        size_t dim = s.size() - 1;
        // check that we have reserved memory for dimension
        reserve(dim);

        // ensure simplex is sorted
        std::sort(s.begin(), s.end());
        // check if simplex is already in complex
        if (spx_to_idx.count(s) > 0) {
            // simplex is already in complex
            return cell_ind(dim, spx_to_idx[s]);
        }
        // we're clear to add
        return _add_unsafe(s);
    }


public:

    // default constructor
    SimplicialComplex() { reserve(0); }

    // constructor that initializes to set dimension
    SimplicialComplex(size_t maxdim) { reserve(maxdim); }

    // constructor that reserves space for given number of simplices in each dimension
    SimplicialComplex(const std::vector<size_t> &dim) {
        for (size_t d = 0; d < dim.size(); d++) {
            reserve(d, dim[d]);
        }
    }

    // maximum dimension of cell
    inline size_t maxdim() const { return _ncells.size() - 1; }
    // number of cells in dimension k
    inline size_t ncells(const size_t k) const { return _ncells[k]; }
    // total number of cells
    size_t ncells() const {
      size_t ct = 0;
      for (size_t k = 0; k < maxdim() + 1; k++) {
        ct += ncells(k);
      }
      return ct;
    }

    // inline cell_ind add_unsafe(
    //     std::vector<size_t> &s
    // ) { return add_unsafe(s); }

    // add simplex to complex
    inline cell_ind add(
        std::vector<size_t> &s
    ) { return add_safe(s); }

    // inline size_t add_vertex() { return _add_vertex(); }
    // inline size_t add_vertices(size_t k) { return _add_vertices(k); }

    inline auto faces_begin(const size_t dim, const size_t i) const {
        return faces[dim-1].cbegin() + ((dim + 1) * i);
    }
    inline auto faces_begin(const cell_ind &ci) const { return faces_begin(ci.dim, ci.ind); }
    inline auto faces_end(const size_t dim, const size_t i) const {
        return faces[dim-1].cbegin() + ((dim + 1) * (i+1));
    }
    inline auto faces_end(const cell_ind &ci) const { return faces_end(ci.dim, ci.ind); }

    // get CSC integer matrix boundary in dimension dim
    CSCMatrix<int, size_t> boundary_csc(const size_t dim) const {
        size_t m = ncells(dim-1);
        size_t n = ncells(dim);
        // create colptr
        std::vector<size_t> colptr(n+1);
        auto vbeg = colptr.begin();
        auto vend = colptr.end();
        size_t val = 0;
        while (vbeg != vend) {
            *vbeg++ = val;
            val += (dim + 1);
        }

        // create rowind
        std::vector<size_t> rowind;
        rowind.reserve((dim+1)*n);

        // create val
        std::vector<int> ival;
        ival.reserve((dim+1)*n);

        std::copy(
            faces[dim-1].cbegin(),
            faces[dim-1].cend(),
            std::back_inserter(rowind)
        );
        std::copy(
            coeff[dim-1].cbegin(),
            coeff[dim-1].cend(),
            std::back_inserter(ival)
        );

        return CSCMatrix<int, size_t>(
            m,
            n,
            colptr,
            rowind,
            ival
        );
    }

    friend class MorsePairing<SimplicialComplex>;

    ~SimplicialComplex() {
        //std::cout << "in simplicial complex destructor" << std::endl;
    }

};
