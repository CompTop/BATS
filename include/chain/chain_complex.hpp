#pragma once
/*
class to store a chain complex
*/
#include <vector>
#include <util/permutation.hpp>

// template over matrix type
template <typename MT>
struct ChainComplex {
	std::vector<size_t> dim;
	std::vector<MT> boundary;

	ChainComplex() {}

	ChainComplex(const std::vector<size_t> &dim, const std::vector<MT> &boundary) : dim(dim), boundary(boundary) {}

	// produce a chain complex from a simplicial or cell complex
	template <typename CpxT>
	ChainComplex(const CpxT& X) {
		dim.resize(X.maxdim() + 1);
		boundary.resize(X.maxdim() + 1);
		for (size_t k = 0; k < X.maxdim() + 1; k++) {
			dim[k] = X.ncells(k);
			if (k == 0) {
				boundary[k] = MT(1, dim[k]);
			} else {
				boundary[k] = MT(X.boundary_csc(k));
			}
		}
	}

	// ChainComplex(const ChainComplex &C) : dim(C.dim), boundary(C.boundary) {}

	inline size_t maxdim() const { return dim.size() - 1; }
	//inline size_t dim(size_t k) const { return dim[k]; }
	//inline MT& boundary(size_t k) { return boundary[k]; }

	MT& operator[](size_t k) {
		if (k >= dim.size()) {
			// throw error
		}
		return boundary[k];
	}

	// permute basis in dimension k
	void permute_basis(size_t k, const std::vector<size_t> &perm) {
		if (k == 0) {
			// only worry about boundary[1]
			boundary[1].permute_rows(perm);
		} else if (k == maxdim()) {
			boundary[maxdim()].permute_cols(perm);
			// only worry about boundary[maxdim()]
		} else {
			// need to handle boundary[k] and boundary[k+1]
			boundary[k].permute_cols(perm);
			boundary[k+1].permute_rows(bats::util::inv_perm(perm));
		}
	}

	// permute basis in all dimensions
	void permute_basis(const std::vector<std::vector<size_t>> &perm) {
		for (size_t k = 0; k < perm.size(); k++) {
			permute_basis(k, perm[k]);
		}
	}

};

// defualt return
template <typename T, typename CpxT>
inline auto __ChainComplex(const CpxT &X, T) {
	using VT = SparseVector<T, size_t>;
	using MT = ColumnMatrix<VT>;

	return ChainComplex<MT>(X);
}
