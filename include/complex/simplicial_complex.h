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

    // container that holds simplices
    // TODO: can combine this with Trie pointer
    std::vector<std::vector<size_t>> spx;

    // container for simplices of a fixed dimension
    std::vector<std::vector<size_t>> faces;
    std::vector<std::vector<int>> coeff;



    // keeps track of how many cells are in given dimension
    std::vector<size_t> _ncells;
    // keeps track of reserved capacity
    std::vector<size_t> _reserved;
	// for iterating over faces
    std::vector<size_t> __face;

    // map to find simplices when forming boundary
    spx_map spx_to_idx;


    // reserve for maxdim dimensional simplices
    void reserve(size_t maxdim) {
        if ( _ncells.size() < maxdim + 1 ) { _ncells.resize(maxdim+1, 0); }
        if ( _reserved.size() < maxdim + 1 ) { _reserved.resize(maxdim+1, 0); }
        if ( spx.size() < maxdim ) { spx.resize(maxdim + 1); }
        if ( faces.size() < maxdim ) { faces.resize(maxdim); }
        if ( coeff.size() < maxdim ) { coeff.resize(maxdim); }
        return;
    }

    // assume dim > 0
    void reserve(size_t dim, size_t k) {
        reserve(dim);
        _reserved[dim] = std::max(_reserved[dim], k);
        if ( spx[dim].capacity() < k * (dim + 1) ) {
            spx[dim].reserve(k * (dim + 1));
        }
        if (dim == 0) { return; }
        if ( faces[dim-1].capacity() < k * (dim + 1) ) {
            faces[dim-1].reserve(k * (dim + 1));
            coeff[dim-1].reserve(k * (dim + 1));
        }
        return;
    }

    // adds a simplex to the complex without doing any checks
    // assume dim > 0
    cell_ind _add_unsafe(std::vector<size_t> &s) {
        size_t dim = s.size() - 1;

        // get index of simplex
        size_t ind = _ncells[dim]++;
        // add to map
        spx_to_idx.emplace(s, ind);
        for (auto v : s) {
            spx[dim].emplace_back(v);
        }

        // determine faces
        if (dim > 0){
            int c = -1;
            // loop over faces in lexicographical order
            size_t spx_len = dim + 1;
            for (size_t k = 0; k < spx_len; k++) {

                size_t k2 = dim-k; // index to skip
                __face.clear();
                for (size_t j = 0; j < k2; j++) {
                    __face.emplace_back(s[j]);
                }
                for (size_t j = k2+1; j < spx_len; j++) {
                    __face.emplace_back(s[j]);
                }

                faces[dim-1].emplace_back(find_idx(__face));
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

    // inline spx_map& trie() { return spx_to_idx; }

    // returns index of simplex
    // TODO: find a way to do this with a single traversal of spx_to_idx
    size_t find_idx(const std::vector<size_t> &s) {
        return spx_to_idx.get(s, bats::NO_IND);
    }

    size_t find_idx(const std::vector<size_t> &s) const {
        return spx_to_idx.get(s, bats::NO_IND);
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
    inline cell_ind add(
        std::vector<size_t> &&s
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

    // get simplex from index
    // get simplex i in dimension dim
    inline auto simplex_begin(const size_t dim, const size_t i) const {
        return spx[dim].cbegin() + ((dim + 1) * i);
    }
    inline auto simplex_end(const size_t dim, const size_t i) const {
        return spx[dim].cbegin() + ((dim + 1) * (i+1));
    }

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


// template over filtration type
// TODO: can do unsafe simplex add
template <typename NT>
void add_dimension_recursive_flag(
    SimplicialComplex &X,
    const NT &nbrs, // lists of neighbors
    const size_t d, // dimension
    const size_t maxd, // max dimension
    const std::vector<size_t> &iter_idxs,
    std::vector<size_t> &spx_idxs
) {
    // sorted simplices will end up here
    std::vector<size_t> spx_idxs2(spx_idxs.size() + 1);
    if (d == maxd) {
        // no recursion - we're adding maximal dimension cells
        for (auto k : iter_idxs) {
            // append k to spx_idxs, sort
            spx_idxs.push_back(k);
            sort_into(spx_idxs, spx_idxs2);

            // add to X
            X.add(spx_idxs2);

            // pop k off spx_idxs
            spx_idxs.pop_back();
        }
    } else if (d < maxd) { // d < maxd
        // recursion
        std::vector<size_t> iter_idxs2; // indices for recursing on
        iter_idxs2.reserve(iter_idxs.size());
        for (auto k : iter_idxs) {
            // append k to spx_idxs, sort
            spx_idxs.push_back(k);
            sort_into(spx_idxs, spx_idxs2);

            // add to X
            X.add(spx_idxs2);

            // recurse
            intersect_sorted_lt(iter_idxs, nbrs[k], k, iter_idxs2);
            add_dimension_recursive_flag(X, nbrs, d+1, maxd, iter_idxs2, spx_idxs2);

            // pop k off spx_idxs
            spx_idxs.pop_back();
        }
    }
    // else do nothing
}


// Flag complex using list of edges
// (edges[2*k], edges[2*k+1]) = (i, j) is an edge
// n - number of vertices
// maxdim - maximum dimension of simplices
SimplicialComplex FlagComplex(
    const std::vector<size_t> &edges,
    const size_t n, // number of 0-cells
    const size_t maxdim
) {

    // check that dimensions agree
    size_t m = edges.size() / 2;
    if (!(edges.size() == 2 * m)) {
        throw std::logic_error("edge vector must have length multiple of 2!");
    }

    // X = SimplicialComplex(maxdim);
    // F = Filtration<T, SimplicialComplex>(X);
    // reset simplicial complex
    SimplicialComplex X(maxdim);

    // sets 0-cells
    std::vector<size_t> spx_idxs(1);
    for (size_t k = 0; k < n; k++) {
        spx_idxs[0] = k;
        X.add(spx_idxs);
    }

    std::vector<std::vector<size_t>> nbrs(n);

    spx_idxs.resize(2); // now time to add edges
    std::vector<size_t> iter_idxs;
    iter_idxs.reserve(n); // maximum size

    for (size_t k = 0; k < m; k++) {
        size_t i = edges[2*k];
        size_t j = edges[2*k + 1];
        spx_idxs[0] = i;
        spx_idxs[1] = j;
        X.add(spx_idxs);
        // std::cout << ret.second << std::endl;
        intersect_sorted(nbrs[i], nbrs[j], iter_idxs);

        if (!iter_idxs.empty()) {
            add_dimension_recursive_flag(X, nbrs, 2, maxdim, iter_idxs, spx_idxs);
        }

        // TODO: use std::set for neighbors - insertion is log(n)
        // nbrs[i].emplace(j);
        // nbrs[j].emplace(i);
        // TODO: insertion sort
        nbrs[i].emplace_back(j);
        std::sort(nbrs[i].begin(), nbrs[i].end());
        nbrs[j].emplace_back(i);
        std::sort(nbrs[j].begin(), nbrs[j].end());
    }

    return X;
}
