#pragma once
/*
compute homology-revealing bases for a chain complex
*/
#include <vector>
#include <set>
#include <chain/chain_complex.h>
#include "reduction.h"

template <typename MT>
struct ReducedChainComplex {

	using chain_type = typename MT::col_type;

	std::vector<size_t> dim;
	std::vector<MT> U; // basis matrices
	std::vector<MT> R; // reduced matrix
	std::vector<std::set<size_t>> I;
	std::vector<p2c_type> p2c;

	ReducedChainComplex() {}

	// compute reduced chain complex from chain complex
	ReducedChainComplex(const ChainComplex<MT> &C) {
		size_t dmax = C.maxdim() + 1;
		dim = C.dim;
		U.resize(dmax);
		R.resize(dmax);
		I.resize(dmax);
		p2c.resize(dmax);

		// TODO: can parallelize this
		for (size_t k = 0; k < dmax; k++) {
			U[k] = MT::identity(C.dim[k]);
			R[k] = C.boundary[k];
			p2c[k] = reduce_matrix(R[k], U[k]);
		}
		// TODO: can parallelize this
		for (size_t k = 0; k < dmax-1; k++) {
			I[k] = extract_basis_indices(R[k], p2c[k+1]);
		}
		I[dmax-1] = extract_basis_indices(R[dmax-1]);
	}


	// size of homology vector space in dimension k
	inline size_t hdim(size_t k) const { return I[k].size(); }
	inline size_t maxdim() const { return R.size() - 1; }

	// get preferred representative i in dimension k
	chain_type get_preferred_representative(
		size_t i,
		const size_t k
	) const {
		auto Iptr = I[k].begin();
		while (i > 0) {
			Iptr++;
			i--;
		}
		return U[k][*Iptr];
	}

	// modify y in-place to be preferred representative for homology class in dimension k
	void find_preferred_representative(
		typename MT::col_type &y,
		const size_t k
	) const {
		if (k == maxdim()) {
			// all cycles generate homology, so nothing to do
			return;
		}
		// else we need to find the preferred representative
		// const ColumnMatrix<TVec> &bdry,
		// const p2c_type &p2c
		size_t j = R[k+1].nrow();
		auto yit = y.upper_bound(j); // find iterator past last nonzero
		while (yit != y.nzbegin()) {
			--yit;
			j = yit->ind;
			// find if j is a pivot of bdry
			if (p2c[k+1].count(j) > 0) {
				auto i = p2c[k+1].at(j);

				// form column i of boundary in homology revealing basis
				auto bdri = u_solve(U[k], R[k+1][i]);
				auto ipiv = bdri.lastnz();
				auto a = (yit->val) / ipiv.val;
				y.axpy(-a, bdri);

				// get next nonzero
				yit = y.upper_bound(j-1);
			}
			// else j is a preferred representative
			// so we do nothing
		}
	}
};


// defualt return
template <typename T, typename CpxT>
inline auto __ReducedChainComplex(const CpxT &F, T) {
	using VT = SparseVector<T, size_t>;
	using MT = ColumnMatrix<VT>;

	return ReducedChainComplex(ChainComplex<MT>(F));
}
