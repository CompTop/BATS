#pragma once
/*
compute cohomology-revealing bases for a cochain complex
*/
#include <vector>
#include <set>
#include <chain/cochain_complex.hpp>
#include "reduction.hpp"
#include "basis.hpp"

namespace bats {


template <typename MT>
class ReducedCochainComplex {
public:

	using cochain_type = typename MT::col_type;

	std::vector<size_t> dim;
	std::vector<MT> U; // basis matrices
	std::vector<MT> R; // reduced matrix
	std::vector<std::vector<size_t>> I;
	std::vector<p2c_type> p2c;

	// size of homology vector space in dimension k
	inline size_t hdim(size_t k) const { return I[k].size(); }
	inline size_t maxdim() const { return R.size() - 1; }

	MT& operator[](size_t k) {
		if (k >= dim.size()) {
			// throw error
		}
		return R[k];
	}

private:
	// set cohomology indices after reduction is completed
	void set_indices() {
		// TODO: can parallelize this
		for (size_t k = 0; k < maxdim(); k++) {
			I[k] = extract_basis_indices(R[k], p2c[k+1]);
		}
		I[maxdim()] = extract_basis_indices(R[maxdim()]);
	}

public:

	ReducedCochainComplex() {}

	// compute reduced chain complex from chain complex
	ReducedCochainComplex(const CochainComplex<MT> &C) {
		size_t dmax = C.maxdim() + 1;
		dim = C.dim;
		U.resize(dmax);
		R.resize(dmax);
		I.resize(dmax);
		p2c.resize(dmax);

		// TODO: can parallelize this
		for (size_t k = 0; k < dmax; k++) {
			U[k] = (k == 0) ? MT::identity(1) : MT::identity(C.dim[k-1]);
			R[k] = C.coboundary[k];
			p2c[k] = reduce_matrix(R[k], U[k]);
		}

		set_indices();

	}

	// compute reduced boundary matrices only with flags
	template <typename algflag>
	ReducedCochainComplex(const CochainComplex<MT> &C, algflag) {
		size_t dmax = C.maxdim() + 1;
		dim = C.dim;
		R.resize(dmax);
		p2c.resize(dmax);
		I.resize(dmax);

		// TODO: can parallelize this
		for (size_t k = 0; k < dmax; k++) {
			R[k] = C.coboundary[k];
			p2c[k] = reduce_matrix(R[k], algflag());
		}

		set_indices();
	}

	// compute reduced boundary matrices and basis with flags
	template <typename algflag>
	ReducedCochainComplex(
		const CochainComplex<MT> &C,
		algflag,
		bats::compute_basis_flag
	) {
		size_t dmax = C.maxdim() + 1;
		dim = C.dim;
		R.resize(dmax);
		U.resize(dmax);
		p2c.resize(dmax);
		I.resize(dmax);

		// TODO: can parallelize this
		for (size_t k = 0; k < dmax; k++) {
			R[k] = C.coboundary[k];
			U[k] = (k == 0) ? MT::identity(1) : MT::identity(C.dim[k-1]);
			p2c[k] = reduce_matrix(R[k], U[k], algflag());
		}

		set_indices();
	}

	template <typename algflag>
	ReducedCochainComplex(
		const CochainComplex<MT> &C,
		algflag,
		bats::clearing_flag
	) {
		size_t dmax = C.maxdim() + 1;
		dim = C.dim;
		R.resize(dmax);
		p2c.resize(dmax);
		I.resize(dmax);

		// do bottom dimension normally
		R[0] = C.coboundary[0];
		p2c[0] = reduce_matrix(R[0], algflag());
		std::vector<size_t> clear_inds = get_clearing_inds(p2c[0]);
		for (ssize_t k = 0; k < dmax; k++) {
			R[k] = C.coboundary[k];
			p2c[k] = reduce_matrix_clearing(R[k], clear_inds, algflag());
			clear_inds = get_clearing_inds(p2c[k]);
		}

		set_indices();
	}

	template <typename algflag>
	ReducedCochainComplex(
		const CochainComplex<MT> &C,
		algflag,
		bats::compression_flag
	) {
		size_t dmax = C.maxdim() + 1;
		dim = C.dim;
		R.resize(dmax);
		p2c.resize(dmax);
		I.resize(dmax);

		// do top dimension normally
		R[dmax-1] = C.coboundary[dmax-1];
		p2c[dmax-1] = reduce_matrix(R[dmax-1], algflag());
		std::vector<bool> comp_inds = get_compression_inds(R[dmax-1]);

		for (ssize_t k = dmax-2; k >= 0; k--) {
			R[k] = C.coboundary[k];
			p2c[k] = reduce_matrix_compression(R[k], comp_inds, algflag());
			comp_inds = get_compression_inds(R[k]);
		}

		set_indices();
	}

	// TODO: FIX THIS
	// compute reduced boundary matrices and basis with flags
	// template <typename algflag>
	// ReducedCochainComplex(
	// 	const CochainComplex<MT> &C,
	// 	algflag,
	// 	bats::compression_flag,
	// 	bats::compute_basis_flag
	// ) {
	// 	size_t dmax = C.maxdim() + 1;
	// 	dim = C.dim;
	// 	R.resize(dmax);
	// 	U.resize(dmax);
	// 	p2c.resize(dmax);
	// 	I.resize(dmax);
	//
	// 	// do bottom dimension normally
	// 	R[0] = C.coboundary[0];
	// 	U[0] = MT::identity(C.dim[0]);
	// 	p2c[0] = reduce_matrix(R[0], U[0], algflag());
	// 	std::vector<bool> comp_inds = get_compression_inds(R[0]);
	// 	for (ssize_t k = 1; k < dmax; k++) {
	// 		R[k] = C.coboundary[k];
	// 		U[k] = MT::identity(C.dim[k]);
	// 		p2c[k] = reduce_matrix_compression(R[k], U[k], comp_inds, algflag());
	// 		comp_inds = get_compression_inds(R[k]);
	// 	}
	//
	// 	set_indices();
	// }

	// put vector/matrix in homology-revealing basis in dimension k
	template <typename TV>
	inline TV to_hom_basis(const TV &v, size_t k) const {
		return u_solve(U[k], v);
	}

	template <typename TV>
	inline TV from_hom_basis(const TV &v, size_t k) const {
		return U[k] * v;
	}



	// get preferred representative i in dimension k
	cochain_type get_preferred_representative(
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
	// assumes y is in homology-revealing basis
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
			// if (p2c[k+1].count(j) > 0) {
			if (p2c[k+1][j] != bats::NO_IND) {
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

	// find the preferred representative for a k-chain
	cochain_type chain_preferred_representative(
		const cochain_type &c, size_t k
	) const {
		auto x = to_hom_basis(c, k);
		find_preferred_representative(x, k);
		return from_hom_basis(x, k);
	}
};


// defualt return
template <typename T, typename CpxT, typename... Args>
inline auto __ReducedCochainComplex(const CpxT &F, T, Args... args) {
	using VT = SparseVector<T, size_t>;
	using MT = ColumnMatrix<VT>;

	return ReducedCochainComplex(CochainComplex<MT>(F), args...);
}

} // namespace bats
