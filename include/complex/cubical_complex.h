#pragma once

/*
	Cubical complex implementation.

	Similar to simplical complex implementation,
	but computation of faces, boundary is different
*/

#include <cstddef>
#include <limits>
#include <vector>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <stdexcept>

#include "abstract_complex.h"
#include <util/common.h>
#include <linalg/sparse_vector.h>
#include <linalg/col_matrix.h>
#include <linalg/csc_matrix.h>
#include <util/sorted.h>
#include <util/trie.h>
#include <filtration/filtration.h>
#include <morse/pairing.h>

#include <util/simplex.h> // for hash function

/*
	cubical complex of dimension d
	cubes are vectors of length 2*d, storing degeneracies
*/
class CubicalComplex
{
private:
    // maybe make template parameter for container?
    // using spx_map = std::unordered_map<std::vector<size_t>, size_t, SimplexHasher>;
    // using spx_map = std::map<std::vector<size_t>, size_t>;
    using spx_map = SparseTrie<size_t, size_t>;
	// used for looking up index of a cube

    // container that holds simplices
	// this is useful for getting order of cubes
    // TODO: can combine this with Trie pointer
    std::vector<std::vector<size_t>> spx;

    // container for cubes of a fixed dimension
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

	size_t __maxdim;


    // reserve for maxdim dimensional simplices
    void reserve(size_t maxdim) {
        if ( _ncells.size() < maxdim + 1 ) { _ncells.resize(maxdim+1, 0); }
        if ( _reserved.size() < maxdim + 1 ) { _reserved.resize(maxdim+1, 0); }
        if ( spx.size() < maxdim + 1 ) { spx.resize(maxdim + 1); }
        if ( faces.size() < maxdim ) { faces.resize(maxdim); }
        if ( coeff.size() < maxdim ) { coeff.resize(maxdim); }
		__maxdim = maxdim;
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

	// check that s is a cube that can be added to complex
	bool is_valid_cube(const std::vector<size_t> &s) {
		if (s.size() != 2*__maxdim) { return false; }
		// check that each interval is either degenerate or unit length
		for (size_t i = 0; i < 2*__maxdim; i+=2) {
			if (s[i+1] < s[i] || s[i+1] > s[i] + 1) { return false; }
		}
		return true;
	}

	// get dimension of cube by removing degeneracies
	size_t cube_dim(const std::vector<size_t> &s) {
		size_t dim = 0;
		for (size_t i = 0; i < 2*__maxdim; i+=2) {
			dim += s[i+1] - s[i]; // assume valid cube
		}
		return dim;
	}

    // adds a cube to the complex without doing any checks
    // assume dim > 0
    cell_ind _add_unsafe(const std::vector<size_t> &s) {
        size_t dim = cube_dim(s);
		__face.resize(s.size());

        // determine faces
        if (dim > 0){
            // first we'll add faces
            size_t old_size = faces[dim-1].size();

			// copy cube into __face
			std::copy(s.cbegin(), s.cend(), __face.begin());

            int c = 1; // coefficient
			// loop over possible faces
			for (size_t k = 0; k < 2*__maxdim; k+=2) {
				if (s[k] == s[k+1]) { continue; } // degenerate
				// else, we have two faces
				__face[k] = s[k+1]; // face # 1
				size_t ind = find_idx(__face);
				if (ind == bats::NO_IND) {
                    // remove boundary placed so far
                    faces[dim-1].resize(old_size);
                    coeff[dim-1].resize(old_size);
                    return cell_ind(dim, bats::NO_IND);
                }
				faces[dim-1].emplace_back(ind);
                coeff[dim-1].emplace_back(c);

				__face[k] = s[k]; __face[k+1] = s[k]; // face # 2
				ind = find_idx(__face);
				if (ind == bats::NO_IND) {
                    // remove boundary placed so far
                    faces[dim-1].resize(old_size);
                    coeff[dim-1].resize(old_size);
                    return cell_ind(dim, bats::NO_IND);
                }
				faces[dim-1].emplace_back(ind);
                coeff[dim-1].emplace_back(-c);

				__face[k] = s[k]; __face[k+1] = s[k+1]; // return face to original state
				c = -c; // alternate coefficient for next dimension
			}
        }

        // if we succeeded (faces were present), we add the cube
        // get index of cube
        size_t ind = _ncells[dim]++;
        // add to map
        spx_to_idx.emplace(s, ind);
        for (auto v : s) {
            spx[dim].emplace_back(v);
        }

        return cell_ind(dim, ind);
    }

    // add simplex to complex with appropriate checks
    cell_ind add_safe(const std::vector<size_t> &s) {
		if (!is_valid_cube(s)) {throw std::runtime_error("Not a valid cube for dimension of complex!");}

        // check if cube is already in complex
        if (spx_to_idx.count(s) > 0) {
            // simplex is already in complex
            return cell_ind(cube_dim(s), spx_to_idx[s]);
        }
        // we're clear to add
        return _add_unsafe(s);
    }

	// adds a cube to the complex without doing any checks
    // assume dim > 0
    cell_ind _add_unsafe_toplex(const std::vector<size_t> &s) {
        size_t dim = cube_dim(s);
		std::vector<size_t> cface(s.size()); // current face

        // determine faces
        if (dim > 0){
            // first we'll add faces

			// copy cube into cface
			std::copy(s.cbegin(), s.cend(), cface.begin());

            int c = 1; // coefficient
			// loop over possible faces
			for (size_t k = 0; k < 2*__maxdim; k+=2) {
				if (s[k] == s[k+1]) { continue; } // degenerate
				// else, we have two faces
				cface[k] = s[k+1]; // face # 1
				size_t ind = find_idx(cface);
				if (ind == bats::NO_IND) {
                    auto ret = _add_unsafe_toplex(cface);
					ind = ret.ind;
                }
				faces[dim-1].emplace_back(ind);
                coeff[dim-1].emplace_back(c);

				cface[k] = s[k]; cface[k+1] = s[k]; // face # 2
				ind = find_idx(cface);
				if (ind == bats::NO_IND) {
                    auto ret = _add_unsafe_toplex(cface);
					ind = ret.ind;
                }
				faces[dim-1].emplace_back(ind);
                coeff[dim-1].emplace_back(-c);

				cface[k] = s[k]; cface[k+1] = s[k+1]; // return face to original state
				c = -c; // alternate coefficient for next dimension
			}
        }

        // if we succeeded (faces were present), we add the cube
        // get index of cube
        size_t ind = _ncells[dim]++;
        // add to map
        spx_to_idx.emplace(s, ind);
        for (auto v : s) {
            spx[dim].emplace_back(v);
        }

        return cell_ind(dim, ind);
    }

	// add simplex to complex with appropriate checks
    cell_ind add_safe_toplex(const std::vector<size_t> &s) {
		if (!is_valid_cube(s)) {throw std::runtime_error("Not a valid cube for dimension of complex!");}

        // check if cube is already in complex
        if (spx_to_idx.count(s) > 0) {
            // simplex is already in complex
            return cell_ind(cube_dim(s), spx_to_idx[s]);
        }
        // we're clear to add
        return _add_unsafe_toplex(s);
    }

public:

    // default constructor
    CubicalComplex() { reserve(0); }

    // constructor that initializes to set dimension
    CubicalComplex(size_t maxdim) { reserve(maxdim); }


    // returns index of simplex
    // TODO: find a way to do this with a single traversal of spx_to_idx
    size_t find_idx(const std::vector<size_t> &s) {
        return spx_to_idx.get(s, bats::NO_IND);
    }

    size_t find_idx(const std::vector<size_t> &s) const {
        return spx_to_idx.get(s, bats::NO_IND);
    }

    // maximum dimension of cell
    inline size_t maxdim() const { return __maxdim; }
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

	// set dimension of complex
	void set_dimension(size_t maxdim) {
		if (ncells() > 0) {
			throw std::runtime_error("can not set dimension - cubes already added!");
			return;
		}
		reserve(maxdim);
	}

    void print_summary() {
        std::cout << "CubicalComplex, maxdim = " << maxdim() << std::endl;
        for (size_t k = 0; k < maxdim() + 1; k++) {
            std::cout << "\tdim " << k << " : " << ncells(k) << " cells" << std::endl;
        }
        std::cout << ncells() << " total" << std::endl;
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

	inline cell_ind add_toplex(
		const std::vector<size_t> &s
	) { return add_safe_toplex(s); }
	inline cell_ind add_toplex(
		const std::vector<size_t> &&s
	) { return add_safe_toplex(s); }

    // inline size_t add_vertex() { return _add_vertex(); }
    // inline size_t add_vertices(size_t k) { return _add_vertices(k); }

    inline auto faces_begin(const size_t dim, const size_t i) const {
        return faces[dim-1].cbegin() + (2*dim * i);
    }
    inline auto faces_begin(const cell_ind &ci) const { return faces_begin(ci.dim, ci.ind); }
    inline auto faces_end(const size_t dim, const size_t i) const {
        return faces[dim-1].cbegin() + (2*dim * (i+1));
    }
    inline auto faces_end(const cell_ind &ci) const { return faces_end(ci.dim, ci.ind); }

    // get simplex from index
    // get simplex i in dimension dim
    inline auto simplex_begin(const size_t dim, const size_t i) const {
        return spx[dim].cbegin() + (2*__maxdim * i);
    }
    inline auto simplex_end(const size_t dim, const size_t i) const {
        return spx[dim].cbegin() + ((2*__maxdim) * (i+1));
    }

	// return simplex i in dimension dim
	std::vector<size_t> get_cube(size_t dim, size_t i) const {
		std::vector<size_t> s;
		s.reserve(dim+1);
		for (auto it = simplex_begin(dim, i); it != simplex_end(dim, i); it++) {
			s.emplace_back(*it);
		}
		return s;
	}

	// return simplices in dimension dim
	std::vector<std::vector<size_t>> get_cubess(const size_t dim) const {
		std::vector<std::vector<size_t>> simplices;
		std::vector<size_t> s;
		for (auto i : spx[dim]) {
			s.emplace_back(i);
			if (s.size() == dim + 1) {
				simplices.emplace_back(s);
				s.clear();
			}
		}
		return simplices;
	}

	// return all simplices
	std::vector<std::vector<size_t>> get_cubes() const {
		std::vector<std::vector<size_t>> simplices;
		std::vector<size_t> s;
		for (size_t dim = 0; dim < maxdim() + 1; dim++){
			for (auto i : spx[dim]) {
				s.emplace_back(i);
				if (s.size() == dim + 1) {
					simplices.emplace_back(s);
					s.clear();
				}
			}
		}
		return simplices;
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

    void save(std::string &fname) const {
        std::ofstream file (fname, std::ios::out);
        if (file.is_open()) {
            file << "SimplicialComplex\n";
            for (size_t dim = 0; dim < maxdim() + 1; dim++) {
                for (size_t i =0; i < ncells(dim); i++) {
                    write_simplex(file, simplex_begin(dim, i), simplex_end(dim, i));
                }
            }
            file.close();
        } else {
            std::cerr << "unable to write simplicial complex to " << fname << std::endl;
        }
    }

};
