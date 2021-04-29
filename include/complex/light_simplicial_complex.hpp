#pragma once
/*
A lightweight simplicial complex
based on k-dimensional simplicies on a vertex set of size n
every possible simplex is given a unique key
and a hash table is used to store the simplices in the complex

This avoids overhead of vectors and tries for the most part

the idea comes from Ripser
but this implementation is independent
and more flexible

the binomial numbering scheme is more compact
but there are inefficiencies, e.g. in binary search.

We'll just use a base-n numbering system
e.g. for simplex (x_0, ... x_k), we use sum x_i n^i
*/

#include <iostream>
#include <unordered_map>
#include <vector>
#include <linalg/csc_matrix.hpp>

namespace bats {

template <typename index_type=size_t, typename hash_table=std::unordered_map<index_type, size_t>>
class LightSimplicialComplex {
private:
    const index_type _n; // fixed size of vertex set
    const index_type _k; // fixed size of maximum simplex dimension
    index_type* _offset = nullptr;     // table of powers of n
    std::vector<std::vector<index_type>> index_to_key;
	std::vector<hash_table> key_to_index;


    // get offset for index i
    inline index_type& offset(size_t i) {
        return *(_offset + i);
    }
    inline const index_type& offset(size_t i) const {
        return *(_offset + i);
    }

    void _initialze() {
        _offset = new index_type[_k+1];

		offset(0) = 1;
		for (size_t i = 1; i < _k+1; i++) {
			offset(i) = _n * offset(i-1);
		}
    }

public:

	// get largest vertex from s of dimension dim
	// s is sum over j=0, j <= dim, x_i n^i where j is sorted
	// we can use integer division to do this
	inline index_type max_vertex(index_type s, size_t dim) const {
		return s / offset(dim); // integer division
	}

private:

	struct simplex_boundary_iterator {
		const LightSimplicialComplex* p; // parent complex
		size_t dim;
		int c;
		int i;
		index_type before; // sum of simplices before
		index_type after; // sum of vertices after



		simplex_boundary_iterator(
			index_type s,
			size_t dim,
			const LightSimplicialComplex* p
		) :  p(p), dim(dim), before(s), after(0) {
			// index for removal - remove last index first
			i = (int) dim;
			// determine coefficient (-1)^i = -1 if i is odd, 1 if even
			c = (i & 0x1) ? -1 : 1;
		}

		std::tuple<index_type, int> next() {
			// face and coefficient
			int coeff = c;
			// get max vertex and remainder
			index_type v = before / p->offset(i);
			before = before % p->offset(i);
			index_type face = before + after;
			after += v * p->offset(i-1);

			// prepare for next iteration
			c = -c;
			i--;
			return std::make_tuple(face, coeff);
		}

		// returns true if i >= 0 (iterator not finished)
		// returns fals if finished
		operator bool() const {
		    return !(i < 0);
		}
	};

public:
    LightSimplicialComplex() : _n(0), _k(0), _offset(nullptr) {}

    LightSimplicialComplex(
        const index_type n,
        const index_type k
    ) : _n(n), _k(k) {
        _initialze();
        index_to_key.resize(_k+1);
		key_to_index.resize(_k+1);
    }

    // copy constructor
    LightSimplicialComplex(
        const LightSimplicialComplex& other
    ) : _n(other._n), _k(other._k),
     index_to_key(other.index_to_key), key_to_index(other.key_to_index) {
        _initialze();
    }

    // move constructor
    // transfer ownership of binomial table
    LightSimplicialComplex(
        LightSimplicialComplex&& other
    ) : _n(other._n), _k(other._k), _offset(other._offset) {
        other._offset = nullptr;
        index_to_key = std::move(index_to_key);
        key_to_index = std::move(key_to_index);
    }

    ~LightSimplicialComplex() {
        if (_offset!= nullptr) {
            delete[] _offset;
        }
    }

    inline index_type maxdim() const { return _k; }
    inline size_t ncells(size_t dim) const {
        return index_to_key[dim].size();
    }
    size_t ncells() const {
      size_t ct = 0;
      for (size_t dim = 0; dim < maxdim() + 1; dim++) {
        ct += ncells(dim);
      }
      return ct;
    }

    // get simplex key
    // assume simplex is sorted
    index_type simplex_key(
        const std::vector<index_type>& s
    ) const {
        index_type key = 0;
        for (size_t j = 0; j < s.size(); j++) {
            key += s[j] * offset(j);
        }
        return key;
    }

	// get simplex corresponding to a particular key
	std::vector<index_type> key_to_simplex(
		const size_t dim,
		index_type key
	) const {
		std::vector<index_type> s(dim+1);
		for (size_t d = dim; d > 0; d--) {
			s[d] = max_vertex(key, d);
			key -= s[d] * offset(d);
		}
		s[0] = key;
		return s;
	}

	// get simplex at index i
	inline std::vector<index_type> get_simplex(const size_t dim, const size_t i) const {
		return key_to_simplex(dim, index_to_key[dim][i]);
	}

	std::vector<std::vector<index_type>> get_simplices(const size_t dim) const {
		std::vector<std::vector<index_type>> simplices(ncells(dim));
		for (size_t i = 0; i < ncells(dim); i++) {
			simplices[i] = get_simplex(dim, i);
		}
		return simplices;
	}

	// find index of simplex s
	// returns bats::NO_IND if s can not be found
	inline size_t find_idx(const size_t dim, const index_type key) const {
		auto loc = key_to_index[dim].find(key);
		return (loc == key_to_index[dim].end()) ? bats::NO_IND : loc->second;
	}

	// find index of simplex s
	// returns bats::NO_IND if s can not be found
	size_t find_idx(const std::vector<index_type> &s) const {
		auto key = simplex_key(s);
		size_t dim = s.size() - 1;
		return find_idx(dim, key);
	}

    // add simplex
    // assumptions:
    // s is sorted
    // boundary is already added
    // s has not already been added
    void add_unsafe(
        const std::vector<index_type>& s
    ) {
        auto k = simplex_key(s);
		size_t dim = s.size() - 1;

		// add to lookup
		key_to_index[dim].emplace(k, index_to_key[dim].size());

        // add to list of simplicies
		index_to_key[dim].emplace_back(k);
    }

	// TODO: make this safe
	inline void add(const std::vector<index_type>& s) {return add_unsafe(s);}

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

		// iterate over simplices in dimension dim
		for (auto s : index_to_key[dim]) {
			// iterate over boundary of simplex
			// for face of simplex
			auto bdry = simplex_boundary_iterator(s, dim, this);
			while (bdry) {
				auto [f, c] = bdry.next();
				rowind.emplace_back(key_to_index[dim-1].at(f));
				ival.emplace_back(c);
			}
		}

        return CSCMatrix<int, size_t>(
            m,
            n,
            colptr,
            rowind,
            ival
        );
    }


    void print_summary() const {
		std::cout << "LightSimplicialComplex, maxdim = " << maxdim() << std::endl;
		for (size_t k = 0; k < maxdim() + 1; k++) {
			std::cout << "\tdim " << k << " : " << ncells(k) << " cells" << std::endl;
		}
		std::cout << ncells() << " total" << std::endl;
	}

};



// // binomial simplicial complex
// template <typename index_type=size_t, typename hash_table=std::unordered_map<index_type, size_t>>
// class BinomialSimplicialComplex {
// private:
//     const index_type n; // fixed size of vertex set
//     const index_type k; // fixed size of maximum simplex dimension
//     index_type* _binom = nullptr;     // binomial coefficient table
//     std::vector<std::vector<index_type>> index_to_key;
// 	std::vector<hash_table> key_to_index;
//
//
//     // get binomial coefficient for _n choose _k
//     // assume _n < n and _k < k
//     inline index_type& binom(size_t _n, size_t _k) {
//         return *(_binom + (_k - 1)  * (n+1) + _n);
//     }
//     inline const index_type& binom(size_t _n, size_t _k) const {
//         return *(_binom + (_k - 1)  * (n+1) + _n);
//     }
//
//     void _initialze() {
//         _binom = new index_type[(n+1)*(k+1)];
//
//         // store so single value of k is in contiguous blocks
//         // set (n 1)
//         for (size_t _n = 0; _n < n+1; _n++) {
//             binom(_n, 1) = _n;
//         }
//         // set (0, k) = 0 for all k
//         for (size_t _k  = 1; _k < k+2; _k++) {
//             binom(0, _k) = 0;
//         }
//         // use identity (n k) = (n-1 k) + (n-1 k-1)
//         for (size_t _k = 2; _k < k+2; _k++) {
//             for (size_t _n = 1; _n < n+1; _n++) {
//                 binom(_n, _k) = binom(_n-1, _k) + binom(_n-1, _k-1);
//             }
//         }
//     }
//
// public:
//
// 	// get largest vertex from s of dimension dim
// 	// s is sum over j=0, j <= dim, (x_j j+1) where j is sorted
// 	// we look for x_dim by searching for
// 	// binom(x_dim, dim+1) <= s < binom(x_dim+1, dim+1)
// 	// use binary search
// 	index_type max_vertex(index_type s, size_t dim) const {
// 		index_type low(0), high(n);
// 		size_t d1 = dim + 1;
// 		while (s < binom(high, d1)) {
// 			index_type diff = high - low;
// 			index_type mid = high - (diff >> 1);
// 			if (s < binom(mid, d1)) {
// 				high = mid - 1;
// 			} else {
// 				low = mid;
// 			}
// 		}
// 		return high;
// 	}
//
// private:
//
// 	struct simplex_boundary_iterator {
// 		const LightSimplicialComplex* p; // parent complex
// 		size_t dim;
// 		int c;
// 		int i;
// 		index_type before; // sum of simplices before
// 		index_type after; // sum of vertices after
//
//
//
// 		simplex_boundary_iterator(
// 			index_type s,
// 			size_t dim,
// 			const LightSimplicialComplex* p
// 		) :  p(p), dim(dim), before(s), after(0) {
// 			// index for removal - remove last index first
// 			i = (int) dim;
// 			// determine coefficient (-1)^i = -1 if i is odd, 1 if even
// 			c = (i & 0x1) ? -1 : 1;
// 		}
//
// 		std::tuple<index_type, int> next() {
// 			// face and coefficient
// 			int coeff = c;
// 			index_type v = p->max_vertex(before, i);
// 			before -= p->binom(v, i+1);
// 			index_type face = before + after;
// 			after += p->binom(v, i);
//
// 			// prepare for next iteration
// 			c = -c;
// 			i--;
// 			return std::make_tuple(face, coeff);
// 		}
//
// 		// returns true if i >= 0 (iterator not finished)
// 		// returns fals if finished
// 		operator bool() const {
// 		    return !(i < 0);
// 		}
// 	};
//
// public:
//     LightSimplicialComplex() : n(0), k(0), _binom(nullptr) {}
//
//     LightSimplicialComplex(
//         const index_type _n,
//         const index_type _k
//     ) : n(_n), k(_k) {
//         _initialze();
//         index_to_key.resize(k+1);
// 		key_to_index.resize(k+1);
//     }
//
//     // copy constructor
//     LightSimplicialComplex(
//         const LightSimplicialComplex& other
//     ) : n(other.n), k(other.k),
//      index_to_key(other.index_to_key), key_to_index(other.key_to_index) {
//         _initialze();
//     }
//
//     // move constructor
//     // transfer ownership of binomial table
//     LightSimplicialComplex(
//         LightSimplicialComplex&& other
//     ) : n(other.n), k(other.k), _binom(other._binom) {
//         other._binom = nullptr;
//         index_to_key = std::move(index_to_key);
//         key_to_index = std::move(key_to_index);
//     }
//
//     ~LightSimplicialComplex() {
//         if (_binom != nullptr) {
//             delete[] _binom;
//         }
//     }
//
//     inline index_type maxdim() const { return k; }
//     inline size_t ncells(size_t dim) const {
//         return index_to_key[dim].size();
//     }
//     size_t ncells() const {
//       size_t ct = 0;
//       for (size_t dim = 0; dim < maxdim() + 1; dim++) {
//         ct += ncells(dim);
//       }
//       return ct;
//     }
//
//     // get simplex key
//     // assume simplex is sorted
//     index_type simplex_key(
//         const std::vector<index_type>& s
//     ) const {
//         index_type key = 0;
//         for (size_t j = 0; j < s.size(); j++) {
//             key += binom(s[j], j+1);
//         }
//         return key;
//     }
//
// 	// TODO
// 	std::vector<index_type> key_to_simplex(
// 		const size_t dim,
// 		index_type key
// 	) {
// 		std::vector<index_type> s(dim+1);
// 		for (size_t d = dim; d > 0; d--) {
// 			s[d] = max_vertex(key, d);
// 			std::cout << "(" << s[d] << "\t" << d << ") =  " <<  binom(s[d], d) << "\n";
// 			key -= binom(s[d], d);
// 		}
// 		s[0] = key;
// 		return s;
// 	}
//
//     // add simplex
//     // assumptions:
//     // s is sorted
//     // boundary is already added
//     // s has not already been added
//     void add_unsafe(
//         const std::vector<index_type>& s
//     ) {
//         auto k = simplex_key(s);
// 		size_t dim = s.size() - 1;
//
// 		// add to lookup
// 		key_to_index[dim][k] = index_to_key[dim].size();
//
//         // add to list of simplicies
// 		index_to_key[dim].emplace_back(k);
//     }
//
// 	// TODO: make this safe
// 	inline void add(const std::vector<index_type>& s) {return add_unsafe(s);}
//
// 	// get CSC integer matrix boundary in dimension dim
//     CSCMatrix<int, size_t> boundary_csc(const size_t dim) const {
//         size_t m = ncells(dim-1);
//         size_t n = ncells(dim);
//         // create colptr
//         std::vector<size_t> colptr(n+1);
//         auto vbeg = colptr.begin();
//         auto vend = colptr.end();
//         size_t val = 0;
//         while (vbeg != vend) {
//             *vbeg++ = val;
//             val += (dim + 1);
//         }
//
// 		// create rowind
//         std::vector<size_t> rowind;
//         rowind.reserve((dim+1)*n);
//
//         // create val
//         std::vector<int> ival;
//         ival.reserve((dim+1)*n);
//
// 		// iterate over simplices in dimension dim
// 		for (auto s : index_to_key[dim]) {
// 			// iterate over boundary of simplex
// 			// for face of simplex
// 			auto bdry = simplex_boundary_iterator(s, dim, this);
// 			while (bdry) {
// 				auto [f, c] = bdry.next();
// 				rowind.emplace_back(key_to_index[dim-1].at(f));
// 				ival.emplace_back(c);
// 			}
// 		}
//
//         return CSCMatrix<int, size_t>(
//             m,
//             n,
//             colptr,
//             rowind,
//             ival
//         );
//     }
//
//
//
//     // mostly for debugging
//     void print_binom() const {
//         for (size_t _n = 0; _n < n+1; _n++) {
//             for (size_t _k = 1; _k < k+2; _k++) {
//                 std::cout << "(" << _n << "\t" << _k << ") = " << binom(_n, _k) << std::endl;
//             }
//         }
//     }
//
//     void print_summary() const {
// 		std::cout << "LightSimplicialComplex, maxdim = " << maxdim() << std::endl;
// 		for (size_t k = 0; k < maxdim() + 1; k++) {
// 			std::cout << "\tdim " << k << " : " << ncells(k) << " cells" << std::endl;
// 		}
// 		std::cout << ncells() << " total" << std::endl;
// 	}
//
// };

}
