#pragma once

#include <cstddef>
#include <limits>
#include <vector>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <tuple>

#include "abstract_complex.hpp"
#include <util/common.hpp>
#include <linalg/sparse_vector.hpp>
#include <linalg/col_matrix.hpp>
#include <linalg/csc_matrix.hpp>
#include <util/sorted.hpp>
#include <util/trie.hpp>
#include <filtration/filtration.hpp>
#include <morse/pairing.hpp>

#include <util/simplex.hpp> // for hash function

namespace bats {

/**
@brief A simplicial complex based using a trie data structure

A class which can be used to hold simplicial complexes on large or exapanding vertex sets.
For a lighter-weight option, see LightSimplicialComplex

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
	std::vector<size_t> __perm; // for sortperms
	std::vector<int> __tmpc; // temporary coeff


    // map to find simplices when forming boundary
    spx_map spx_to_idx;


    // reserve for maxdim dimensional simplices
    void reserve(size_t maxdim) {
        if ( _ncells.size() < maxdim + 1 ) { _ncells.resize(maxdim+1, 0); }
        if ( _reserved.size() < maxdim + 1 ) { _reserved.resize(maxdim+1, 0); }
        if ( spx.size() < maxdim + 1 ) { spx.resize(maxdim + 1); }
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
	cell_ind _add_unsafe(const std::vector<size_t> &s) {
        size_t dim = s.size() - 1;

		// get index of simplex
        size_t ind = _ncells[dim]++;

        // determine faces
        if (dim > 0){
            // first we'll add faces
            size_t old_size = faces[dim-1].size();

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

                size_t ind = find_idx(__face);
                // return NO_IND if faces haven't been added
                if (ind == bats::NO_IND) {
                    // remove boundary placed so far
                    faces[dim-1].resize(old_size);
                    coeff[dim-1].resize(old_size);
					_ncells[dim]--; // we don't add the cell
                    return cell_ind(dim, bats::NO_IND);
                }
                faces[dim-1].emplace_back(ind);
                coeff[dim-1].emplace_back(c);
                c = -c;
            }
			// if we made it to this point, the simplex will be added
			// we should sort the faces and coeff
			bats::util::fill_sortperm(
				faces[dim-1].begin() + ((dim + 1) * ind),
				faces[dim-1].begin() + ((dim + 1) * (ind+1)),
				__perm
			);
			bats::util::apply_perm(
				faces[dim-1].begin() + ((dim + 1) * ind),
				__face,
				__perm
			);
			bats::util::apply_perm(
				coeff[dim-1].begin() + ((dim + 1) * ind),
				__tmpc,
				__perm
			);
        }

        // if we succeeded (faces were present), we add the simplex

        // add to map
        spx_to_idx.emplace(s, ind);
        for (auto v : s) {
            spx[dim].emplace_back(v);
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

	// adds a simplex to the complex without doing any checks
    // recurse on faces
	// put any added cells in ci
	cell_ind _add_recursive(std::vector<size_t> &s, std::vector<cell_ind> &ci) {
        size_t dim = s.size() - 1;

		// get index of simplex
        size_t ind = _ncells[dim]++;

        // determine faces
        if (dim > 0){
            // first we'll add faces

			std::vector<size_t> face; // newly allocated face

            int c = -1;
            // loop over faces in lexicographical order
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

                size_t ind = find_idx(face);
                // add face recursively
                if (ind == bats::NO_IND) {
					cell_ind fci = _add_recursive(face, ci);
					ind = fci.ind;
                }
                faces[dim-1].emplace_back(ind);
                coeff[dim-1].emplace_back(c);
                c = -c;
            }
			// if we made it to this point, the simplex will be added
			// we should sort the faces and coeff
			bats::util::fill_sortperm(
				faces[dim-1].begin() + ((dim + 1) * ind),
				faces[dim-1].begin() + ((dim + 1) * (ind+1)),
				__perm
			);
			bats::util::apply_perm(
				faces[dim-1].begin() + ((dim + 1) * ind),
				__face,
				__perm
			);
			bats::util::apply_perm(
				coeff[dim-1].begin() + ((dim + 1) * ind),
				__tmpc,
				__perm
			);
        }

        // if we succeeded (faces were present), we add the simplex

        // add to map
        spx_to_idx.emplace(s, ind);
        for (auto v : s) {
            spx[dim].emplace_back(v);
        }
		ci.emplace_back(cell_ind(dim, ind));

        return cell_ind(dim, ind);
    }

	// add simplex to complex with appropriate checks
    std::vector<cell_ind> add_safe_recursive(std::vector<size_t> &s) {
        size_t dim = s.size() - 1;
        // check that we have reserved memory for dimension
        reserve(dim);

        // ensure simplex is sorted
        std::sort(s.begin(), s.end());
        // check if simplex is already in complex
        if (spx_to_idx.count(s) > 0) {
            // simplex is already in complex
            return std::vector{cell_ind(dim, spx_to_idx[s])};
        }
        // we're clear to add
		std::vector<cell_ind> ci;
		_add_recursive(s, ci);
        return ci;
    }


public:


    SimplicialComplex() { reserve(0); }

	/**
	Initialization up to a certain dimension
	@param[in] maxdim - the maximum dimension simplex expected to be added
	*/
    SimplicialComplex(size_t maxdim) { reserve(maxdim); }

	/**
	Initialization on a certain vertex
	@param[in] n - the maximum expected vertex index
	@param[in] maxdim - the maximum dimension simplex expected to be added
	*/
	SimplicialComplex(size_t n, size_t maxdim) {
		reserve(maxdim);
		reserve(0, n);
	}

    // constructor that reserves space for given number of simplices in each dimension
    SimplicialComplex(const std::vector<size_t> &dim) {
        for (size_t d = 0; d < dim.size(); d++) {
            reserve(d, dim[d]);
        }
    }

    SimplicialComplex(const std::string &&fname) {
        reserve(0);
        std::ifstream file (fname, std::ios::in);
        if (file.is_open()) {
            std::string line;
            getline(file, line); // TODO: should be "SimplicialComplex"
            std::vector<size_t> spx;
            while (getline(file, line)) {
                bats::util::read_simplex(line, spx);
                this->add_safe(spx);
            }
            file.close();
        } else {
            std::cerr << "unable to read simplicial complex from " << fname << std::endl;
        }
    }

    // // read from file
    // SimplicialComplex(const std::string fname) {
    //
    // }

    // inline spx_map& trie() { return spx_to_idx; }

	/**
	Find the index of a simplex
	@param[in]	s	A vector containing the simplex
	@return The index associated with the simplex.
	Returns bats::NO_IND if the simplex is not in the complex.
	*/
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

	inline std::vector<cell_ind> add_recursive(
		std::vector<size_t> &s
	) { return add_safe_recursive(s); }
	inline std::vector<cell_ind> add_recursive(
		std::vector<size_t> &&s
	) { return add_safe_recursive(s); }
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

	/**
	fill s with simple i in dimension dim
	*/
	void get_simplex(size_t dim, size_t i, std::vector<size_t>& s) const {
		s.clear();
		for (auto it = simplex_begin(dim, i); it != simplex_end(dim, i); it++) {
			s.emplace_back(*it);
		}
		return;
	}

	/**
	return simplex i in dimension dim
	*/
	std::vector<size_t> get_simplex(size_t dim, size_t i) const {
		std::vector<size_t> s;
		s.reserve(dim+1);
		get_simplex(dim, i, s);
		return s;
	}

	inline auto get_cell(size_t dim, size_t i, std::vector<size_t>& s) const {return get_simplex(dim, i, s);}
	inline auto get_cell(size_t dim, size_t i) const {return get_simplex(dim, i);}

	// return simplices in dimension dim
	std::vector<std::vector<size_t>> get_simplices(const size_t dim) const {
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
	std::vector<std::vector<size_t>> get_simplices() const {
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


	// take union of simplicial complex with Y
	void union_add(const SimplicialComplex& Y) {
		std::vector<size_t> s;
		for (size_t dim = 0; dim < Y.maxdim() + 1; dim++){
			for (auto i : Y.spx[dim]) {
				s.emplace_back(i);
				if (s.size() == dim + 1) {
					add_safe(s);
					s.clear();
				}
			}
		}
	}

	friend SimplicialComplex simplicial_union(
		const SimplicialComplex& X,
		const SimplicialComplex& Y
	) {
		SimplicialComplex XY;
		for (size_t k = 0; k <= X.maxdim(); k++) {
			for (size_t nk = 0; nk < X.ncells(k); nk++) {
				auto s = X.get_simplex(k, nk);
				XY.add_safe(s);
			}
		}
		for (size_t k = 0; k <= Y.maxdim(); k++) {
			for (size_t nk = 0; nk < Y.ncells(k); nk++) {
				auto s = Y.get_simplex(k, nk);
				XY.add_safe(s);
			}
		}
		return XY;
	}

	// get indices in subcomplex A in dimension dim
	std::vector<size_t> get_indices(const SimplicialComplex& A, size_t dim) const {
		if (dim > A.maxdim()) { return std::vector<size_t>(); }
		std::vector<size_t> inds(A.ncells(dim));
		for (size_t i = 0; i < A.ncells(dim); i++) {
			auto s = A.get_simplex(dim, i);
			inds[i] = find_idx(s);
		}
		return inds;
	}

	std::vector<std::vector<size_t>> get_indices(const SimplicialComplex& A) const {
		std::vector<std::vector<size_t>> inds;
		for (size_t dim = 0; dim < maxdim() + 1; dim++) {
			inds.emplace_back(get_indices(A, dim));
		}
		return inds;
	}

	// return intersection of simplicial complexes X and Y
	friend SimplicialComplex intersection(
		const SimplicialComplex& X,
		const SimplicialComplex& Y
	) {
		SimplicialComplex XY(X.maxdim());
		for (size_t k = 0; k <= X.maxdim(); k++) {
			for (size_t nk = 0; nk < X.ncells(k); nk++) {
				auto s = X.get_simplex(k, nk);
				if (Y.find_idx(s) != bats::NO_IND) {
					XY._add_unsafe(s);
				}
			}
		}
		return XY;
	}

	// return inds, vals in boundary
	auto boundary(const size_t dim, const size_t k) const {
		size_t indbegin = (dim+1) * k;
		size_t indend = indbegin + dim + 1;
		std::vector<size_t> ind;
		std::vector<int> val;
		for (size_t i = indbegin; i < indend; i++) {
			ind.emplace_back(faces[dim-1][i]);
			val.emplace_back(coeff[dim-1][i]);
		}
		return std::tuple(ind, val);
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

    void save(std::string &fname) const {
        std::ofstream file (fname, std::ios::out);
        if (file.is_open()) {
            file << "SimplicialComplex\n";
            for (size_t dim = 0; dim < maxdim() + 1; dim++) {
                for (size_t i =0; i < ncells(dim); i++) {
                    bats::util::write_simplex(file, simplex_begin(dim, i), simplex_end(dim, i));
                }
            }
            file.close();
        } else {
            std::cerr << "unable to write simplicial complex to " << fname << std::endl;
        }
    }

	void print_summary() const {
		std::cout << "SimplicialComplex, maxdim = " << maxdim() << std::endl;
		for (size_t k = 0; k < maxdim() + 1; k++) {
			std::cout << "\tdim " << k << " : " << ncells(k) << " cells" << std::endl;
		}
		std::cout << ncells() << " total" << std::endl;
	}

	void print_cells() const {
		std::cout << "SimplicialComplex, maxdim = " << maxdim() << std::endl;
		for (auto& sx : get_simplices()) {
			for (auto& s : sx) {
				std::cout << s << ",";
			}
			std::cout << std::endl;
		}
	}

};




// helper functions for TriangulatedProduct
// prod_ind turns product of 0-simplices to 0-simplex
// ind_prod performs inverse operation
inline size_t  prod_ind(size_t i, size_t j, size_t n) {
	return i + n*j;
}
// std::tuple<size_t, size_t> ind_prod(size_t k, size_t n) {
// 	return std::tie(k % n, k / n);
// }

// loop over all possible simplices
template <typename CpxT, typename itT>
void product_paths(
	CpxT &XY,
	itT xit,
	const itT xend,
	itT yit,
	const itT yend,
	std::vector<size_t> &s,
	const size_t n
) {
	// put on vertex
	s.emplace_back(prod_ind(*xit, *yit, n));
	++xit; ++yit;
	if (xit == xend && yit == yend) {
		XY.add_recursive(s);
		s.pop_back();
		return;
	}
	if (xit != xend) {
		// increment in x direction
		product_paths(XY, xit, xend, --yit, yend, s, n);
		++yit;
	}
	if (yit != yend) {
		// increment in y direction
		product_paths(XY, --xit, xend, yit, yend, s, n);
	}
	// pop off vertex
	s.pop_back();
	return;
}

// loop over all possible simplices
template <typename CpxT, typename itT>
void product_paths(
	CpxT &XY,
	itT xit,
	const itT xend,
	itT yit,
	const itT yend,
	std::vector<size_t> &s,
	const size_t n,
	std::vector<cell_ind> &ci // holds all top-level simplices added
) {
	// put on vertex
	s.emplace_back(prod_ind(*xit, *yit, n));
	++xit; ++yit;
	if (xit == xend && yit == yend) {
		ci.emplace_back((XY.add_recursive(s)).back());
		s.pop_back();
		return;
	}
	if (xit != xend) {
		// increment in x direction
		product_paths(XY, xit, xend, --yit, yend, s, n, ci);
		++yit;
	}
	if (yit != yend) {
		// increment in y direction
		product_paths(XY, --xit, xend, yit, yend, s, n, ci);
	}
	// pop off vertex
	s.pop_back();
	return;
}


// Triangulated product of X and Y up to simplex dimension k
// use translator
// n is number of possible vertices on X
template <typename CpxT>
CpxT TriangulatedProduct(
	const CpxT& X,
	const CpxT& Y,
	const size_t maxdim,
	const size_t n
) {

	CpxT XY(n * Y.ncells(0), maxdim);

	std::vector<size_t> s; // placeholder simplex
	// loop over simplex dimension of product
	for (size_t dim = 0; dim <=maxdim; dim++) {
		// dimension of simplices we'll take a product of
		for (size_t dX = 0; dX <= dim; dX++) {
			size_t dY = dim - dX;
			if (dX > X.maxdim() || dY > Y.maxdim()) {continue;}

			// loop over simplices of dimension dX
			for (size_t iX = 0; iX < X.ncells(dX); iX++) {
				// loop over simplices of dimension dY
				for (size_t iY = 0; iY < Y.ncells(dY); iY++) {
					product_paths(
						XY,
						X.simplex_begin(dX, iX),
						X.simplex_end(dX, iX),
						Y.simplex_begin(dY, iY),
						Y.simplex_end(dY, iY),
						s,
						n
					);
				}
			}
		}
	}
	return XY;
}

template <typename CpxT>
inline CpxT TriangulatedProduct(
	const CpxT& X,
	const CpxT& Y,
	const size_t maxdim
) { return TriangulatedProduct(X, Y, maxdim, X.ncells(0)); }

template <typename CpxT>
inline CpxT TriangulatedProduct(
	const CpxT& X,
	const CpxT& Y
) { return TriangulatedProduct(X, Y, X.maxdim() + Y.maxdim(), X.ncells(0)); }

} // namespace bats
