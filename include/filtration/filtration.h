#pragma once

#include <vector>
#include <algorithm>
#include <numeric>
#include <complex/abstract_complex.h>
#include <linalg/sparse_vector.h>
#include <linalg/col_matrix.h>
#include <util/sorted.h>
#include <morse/pairing.h>
#include <iostream>

// template over
//  TC - complex type
//  TF - filtration type
template <typename TF, class CpxT>
class Filtration
{
private:
    MorsePairing<CpxT> P;
    std::vector<std::vector<TF>> val; // filtration value for each cell

    inline void reserve(size_t maxdim) {
        while(val.size() < maxdim+1) {val.emplace_back(std::vector<TF>()); }
        return;
    }
    void reserve(size_t dim, size_t n) {
        reserve(dim);
        if (val[dim].size() < n) {
          val[dim].resize(n);
        }
    }

	// add to underlying complex
	template <class ...Ts>
	inline cell_ind _add_unsafe(const TF t, Ts (&...args)) {
		cell_ind ret = P._add_unsafe(args...);
		val[ret.dim][ret.ind] = t;
		return ret;
	}
	// add to underlying complex
	template <class ...Ts>
	inline cell_ind _add_unsafe_reserve(const TF t, Ts (&...args)) {
		cell_ind ret = P._add_unsafe_reserve(args...);
		reserve(ret.dim, ret.ind+1);
		val[ret.dim][ret.ind] = t;
		return ret;
	}



public:

    // TODO: complete pairs in MorsePairing
    // TODO: extract pairs in MorseParing to get barcode
	Filtration() {};

    // initialize on complex
    Filtration(CpxT &C) : P(C) {
        for (size_t dim = 0; dim < P.maxdim(); dim++){
			reserve(dim, P.ncells(dim));
		}
    };

	inline size_t maxdim() const { return val.size() - 1; }
	inline size_t ncells(const size_t dim) const { return val[dim].size(); }

	// add to underlying complex
	template <class ...Ts>
	inline cell_ind add(TF t, Ts (&...args)) {
		cell_ind ret = P.add(args...);
		reserve(ret.dim, ret.ind+1);
		val[ret.dim][ret.ind] = t;
		return ret;
	}

	// add to underlying complex and if possible pair with face
	template <class ...Ts>
	cell_ind add_pair(TF t, Ts (&...args)) {
		// first add to complex
		cell_ind ret = add(t, args...);
		// now iterate over boundary
		if (ret.dim > 1) {
			for (auto k = P.faces_begin(ret.dim, ret.ind); k < P.faces_end(ret.dim, ret.ind); k++) {
				if (val[ret.dim-1][*k] == t) {
					// only return if pair set successfully
					if (P.set_pair(ret.dim-1, *k, ret.ind)) {return ret;}
				}
			}
		} else if (ret.dim == 1) {
			// TODO: check if boundary is zero e.g. for cell complexes
			auto k = P.faces_begin(ret.dim, ret.ind);
			P.set_pair_edge(*k, *(k+1), ret.ind, val[0]);
		}
		return ret;
	}

	// add to underlying complex and if possible pair with face
	template <class ...Ts>
	cell_ind _add_pair_unsafe(TF t, Ts (&...args)) {
		// first add to complex
		cell_ind ret = _add_unsafe_reserve(t, args...);
		// now iterate over boundary
		if (ret.dim > 1) {
			for (auto k = P.faces_begin(ret.dim, ret.ind); k < P.faces_end(ret.dim, ret.ind); k++) {
				if (val[ret.dim-1][*k] == t) {
					// only return if pair set successfully
					if (P.set_pair(ret.dim-1, *k, ret.ind)) {return ret;}
				}
			}
		} else if (ret.dim == 1) {
			// TODO: check if boundary is zero e.g. for cell complexes
			auto k = P.faces_begin(ret.dim, ret.ind);
			P.set_pair_edge(*k, *(k+1), ret.ind, val[0]);
		}
		return ret;
	}

	MorsePairing<CpxT>& pairing() {
		return P;
	}

	inline std::vector<size_t> sortperm(const size_t dim) const { return bats::sortperm(val[dim]); }

	// friend void add_dimension_recursive_flag_unsafe(
	//     Filtration&,
	//     const std::vector<std::vector<size_t>>&,
	//     const size_t,
	//     const size_t,
	//     const std::vector<size_t>&,
	//     std::vector<size_t>&,
	//     const TF
	// );

	~Filtration() {
		//std::cout << "in Filtration destructor" << std::endl;
	}

};

// TODO:
// functions that extend filtrations on a pre-existing complex


// // empty initializer
// Filtration() : cpx(TC()) {}
//
// // initialize with maxdim
// Filtration(size_t maxdim) : cpx(TC(maxdim)) {
//     val = std::vector<std::vector<TF>>(maxdim+1);
// }
//
// // initialize with only complex
// Filtration(TC cpx) : cpx(cpx) {}
//
// Filtration(TC cpx,
//     std::vector<std::vector<TF>> &val) : cpx(cpx), val(val) {}




// bool add_unsafe(std::vector<size_t> &c, const TF t);
//
// bool add(std::vector<size_t> &c, const TF t);
// bool add(std::vector<size_t> &&c, const TF t);

// void add_dimension_recursive_flag_unsafe(std::vector<std::vector<size_t>> &nbrs,
//     size_t d, size_t maxd,
//     std::vector<size_t> &iter_idxs,
//     std::vector<size_t> &spx_idxs,
//     const TF t
// );

// // assume filtration value has been set for each cell.  Find sort permutaion in each dimension
// void sort() {
//   filt_perm.resize(cpx.maxdim() + 1);
//   inv_perm.resize(cpx.maxdim() + 1);
//   // loop over dimension
//   for (size_t dim = 0; dim < cpx.maxdim() + 1; dim++) {
//     // std::cout << "sorting dim = " << dim << std::endl;
//     size_t ncells_dim = cpx.ncells(dim);
//     // std::cout << "ncells dim = " << ncells_dim << std::endl;
//     filt_perm[dim] = std::vector<size_t>(ncells_dim);
//     // std::cout << "created vector" << std::endl;
//     std::iota(filt_perm[dim].begin(), filt_perm[dim].end(), 0);
//     // std::cout << "filled iota" << std::endl;
//     std::sort(filt_perm[dim].begin(), filt_perm[dim].end(),
//       [&](const int& a, const int& b) {
//         return (val[dim][a] < val[dim][b]);
//       }
//     );
//     // std::cout << "sorted" << std::endl;
//     // now filt_perm[dim] holds the sort permutation on cells for the filtration
//     inv_perm[dim] = std::vector<size_t>(ncells_dim);
//     for (size_t i = 0; i < ncells_dim; i++) {
//       inv_perm[dim][filt_perm[dim][i]] = i;
//     }
//     // std::cout << "inv sorted" << std::endl;
//     // now inv_perm[dim] holds the inverse permutation
//   }
//   return;
// }
//
// inline size_t ncells() { return cpx.ncells(); }
// inline size_t maxdim() { return cpx.maxdim(); }
//
// void print(size_t dim) {
//   for (auto i : filt_perm[dim]) {
//     std::cout << val[dim][i] << " | ";
//     cpx.print_cell(dim, i);
//   }
// }
//
// void print() {
//   for (size_t dim = 0; dim < cpx.maxdim() + 1; dim ++) {
//     std::cout << "dim " << dim << " : " << cpx.ncells(dim) << " cells" << std::endl;
//     print(dim);
//   }
// }
//
// // create boundary in dimension dim
// // template over column type
// template <class TVec>
// ColumnMatrix<TVec> boundary(size_t dim) {
//   // first get column matrix of boundary
//   ColumnMatrix<TVec> B = cpx.template boundary<TVec>(dim);
//   // // TODO: permute rows and columns
//   B.permute(filt_perm[dim-1], filt_perm[dim]);
//   return B;
// }

// // add cell c to complex at value t
// // Special method fo Simplicial Complexes
// // template <class TC, typename TF> class Filtration; // primary template
// // template <typename TF>
// template <>
// bool Filtration<SimplicialComplex, float>::add(std::vector<size_t> &&c, const float t) {
//   bool added = cpx.add(c);
//   if (added) {
//     size_t dim = c.size() - 1;
//     while (val.size() < dim + 1) {
//       val.emplace_back(std::vector<float>());
//     }
//     val[dim].emplace_back(t);
//   }
//   return added;
// }
//
// template <>
// bool Filtration<SimplicialComplex, float>::add_unsafe(std::vector<size_t> &c, const float t) {
//   bool added = cpx.add_unsafe(c);
//   if (added) {
//     size_t dim = c.size() - 1;
//     val[dim].emplace_back(t);
//   }
//   return added;
// }
//
// template<>
// void Filtration<SimplicialComplex, float>::add_dimension_recursive_flag_unsafe(
//   std::vector<std::vector<size_t>> &nbrs,
//   size_t d, size_t maxd,
//   std::vector<size_t> &iter_idxs,
//   std::vector<size_t> &spx_idxs,
//   const float t
// ) {
//   cpx.add_dimension_recursive_flag_unsafe(nbrs, d, maxd, iter_idxs, spx_idxs);
//   // add filtration time until val arrays have same number of elements as complex
//   for (size_t dim = d; dim < maxd + 1; dim++) {
//     //val[dim].reserve(cpx.ncells(dim));
//     while (val[dim].size() < cpx.ncells(dim)) {
//       val[dim].emplace_back(t);
//     }
//   }
// }
