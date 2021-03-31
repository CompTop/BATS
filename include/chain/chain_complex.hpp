#pragma once
/*
class to store a chain complex
*/
#include <vector>
#include <util/permutation.hpp>

namespace bats {


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

	// produce a relative chain complex C(X, A)
	template <typename CpxT>
	ChainComplex(const CpxT& X, const CpxT& A) {
		dim.resize(X.maxdim() + 1);
		boundary.resize(X.maxdim() + 1);
		auto inds = X.get_indices(A);
		std::vector<std::vector<size_t>> cinds;
		for (size_t k = 0; k < X.maxdim() + 1; k++) {
			std::sort(inds[k].begin(), inds[k].end());
			cinds.emplace_back(bats::util::sorted_complement(inds[k], X.ncells(k)));
			if (k == 0) {
				boundary[k] = MT(1, X.ncells(k) - A.ncells(k));
			} else {
				CSCMatrix<int, size_t> B = X.boundary_csc(k);
				boundary[k] = MT(B.submatrix(cinds[k-1], cinds[k]));
			}
			dim[k] = boundary[k].ncol();
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

	// tensor product of chain complexes A\otimes B
	// construct tensor product up to dimension dmax
	friend ChainComplex tensor_product(const ChainComplex &A, const ChainComplex &B, size_t dmax) {
		ChainComplex C;
		for (size_t n = 0; n <= dmax; n++) {
			for (size_t Adim = 0; Adim <= n; Adim++) {
				size_t Bdim = n - Adim;
				if (Adim > A.maxdim() || Bdim > B.maxdim()) {continue;}
				// add to tensor product
			}
		}
		return C;
	}

};

// defualt return
template <typename T, typename CpxT>
inline auto __ChainComplex(const CpxT &X, T) {
	using VT = SparseVector<T, size_t>;
	using MT = ColumnMatrix<VT>;

	return ChainComplex<MT>(X);
}
template <typename T, typename CpxT>
inline auto __ChainComplex(const CpxT &X, const CpxT &A, T) {
	using VT = SparseVector<T, size_t>;
	using MT = ColumnMatrix<VT>;

	return ChainComplex<MT>(X, A);
}

} // namespace bats
