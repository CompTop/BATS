#pragma once
/*
class to store a chain complex
*/
#include <vector>
#include <util/permutation.hpp>

namespace bats {


// template over matrix type
template <typename MT>
struct CochainComplex {
	std::vector<size_t> dim;
	std::vector<MT> coboundary;

	CochainComplex() {}

	CochainComplex(const std::vector<size_t> &dim, const std::vector<MT> &coboundary) : dim(dim), coboundary(coboundary) {}

	// produce a cochain complex from a simplicial or cell complex
	template <typename CpxT>
	CochainComplex(const CpxT& X) {
		dim.resize(X.maxdim() + 1);
		coboundary.resize(X.maxdim() + 1);
		for (size_t k = 0; k < X.maxdim() + 1; k++) {
			dim[k] = X.ncells(k);
			if (k == 0) {
				coboundary[k] = MT(1, dim[k]).transpose();
			} else {
				coboundary[k] = MT(X.boundary_csc(k)).transpose();
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
		return coboundary[k];
	}

	// permute basis in dimension k
	void permute_basis(size_t k, const std::vector<size_t> &perm) {
		if (k == 0) {
			// only worry about boundary[1]
			coboundary[1].permute_rows(perm);
		} else if (k == maxdim()) {
			coboundary[maxdim()].permute_cols(perm);
			// only worry about boundary[maxdim()]
		} else {
			// need to handle boundary[k] and boundary[k+1]
			coboundary[k].permute_cols(perm);
			coboundary[k+1].permute_rows(bats::util::inv_perm(perm));
		}
	}

	// permute basis in all dimensions
	void permute_basis(const std::vector<std::vector<size_t>> &perm) {
		for (size_t k = 0; k < perm.size(); k++) {
			permute_basis(k, perm[k]);
		}
	}

	void ipermute_basis(size_t k, const std::vector<size_t> &perm) {
		if (k == 0) {
			// only worry about boundary[1]
			coboundary[1].permute_rows(bats::util::inv_perm(perm));
		} else if (k == maxdim()) {
			coboundary[k].ipermute_cols(perm);
			// only worry about boundary[maxdim()]
		} else {
			// need to handle boundary[k] and boundary[k+1]
			coboundary[k].ipermute_cols(perm);
			coboundary[k+1].permute_rows(bats::util::inv_perm(perm));
		}
	}

};



// defualt return
template <typename T, typename CpxT>
inline auto __CochainComplex(const CpxT &X, T) {
	using VT = SparseVector<T, size_t>;
	using MT = ColumnMatrix<VT>;

	return CochainComplex<MT>(X);
}

} // namespace bats
