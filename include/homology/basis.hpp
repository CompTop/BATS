#pragma once
/*
compute homology-revealing bases for a chain complex
*/
#include <vector>
#include <set>
#include <chain/chain_complex.hpp>
#include "reduction.hpp"

namespace bats {
// flag to indicate that basis should be computed
struct compute_basis_flag {};


/*
Maintains factorizations
B_k U_k = R_k
or B_k = R_k U_k^{-1}
*/
template <typename MT>
class ReducedChainComplex {
public:

	using chain_type = typename MT::col_type;

	//std::vector<size_t> dim;
	std::vector<MT> U; // basis matrices
	std::vector<MT> R; // reduced matrix
	std::vector<std::vector<size_t>> I;
	std::vector<p2c_type> p2c;

	// size of homology vector space in dimension k
	inline size_t hdim(size_t k) const { return I[k].size(); }
	inline size_t betti(size_t k) const { return I[k].size(); }
	inline size_t maxdim() const { return R.size() - 1; }
	inline size_t dim(size_t k) const { return R[k].ncol(); }

	MT& operator[](size_t k) {
		if (k >= R.size()) {
			// throw error
		}
		return R[k];
	}

private:
	// set homology indices after reduction is completed
	void set_indices() {
		// TODO: can parallelize this
		for (size_t k = 0; k < maxdim(); k++) {
			I[k] = extract_basis_indices(R[k], p2c[k+1]);
		}
		I[maxdim()] = extract_basis_indices(R[maxdim()]);
	}

public:

	ReducedChainComplex() {}

	// compute reduced chain complex from chain complex
	ReducedChainComplex(const ChainComplex<MT> &C) {
		size_t dmax = C.maxdim() + 1;
		//dim = C.dim;
		U.resize(dmax);
		R.resize(dmax);
		I.resize(dmax);
		p2c.resize(dmax);

		// TODO: can parallelize this
		for (size_t k = 0; k < dmax; k++) {
			U[k] = MT::identity(C.dim(k));
			R[k] = C.boundary[k];
			p2c[k] = reduce_matrix(R[k], U[k]);
		}

		set_indices();

	}

	// compute reduced boundary matrices only with flags
	template <typename algflag>
	ReducedChainComplex(const ChainComplex<MT> &C, algflag) {
		size_t dmax = C.maxdim() + 1;
		//dim = C.dim;
		R.resize(dmax);
		p2c.resize(dmax);
		I.resize(dmax);

		// TODO: can parallelize this
		for (size_t k = 0; k < dmax; k++) {
			R[k] = C.boundary[k];
			p2c[k] = reduce_matrix(R[k], algflag());
		}

		set_indices();
	}

	// compute reduced boundary matrices and basis with flags
	template <typename algflag>
	ReducedChainComplex(
		const ChainComplex<MT> &C,
		algflag,
		bats::compute_basis_flag
	) {
		size_t dmax = C.maxdim() + 1;
		//dim = C.dim;
		R.resize(dmax);
		U.resize(dmax);
		p2c.resize(dmax);
		I.resize(dmax);

		// TODO: can parallelize this
		for (size_t k = 0; k < dmax; k++) {
			R[k] = C.boundary[k];
			U[k] = MT::identity(C.dim(k));
			p2c[k] = reduce_matrix(R[k], U[k], algflag());
		}

		set_indices();
	}

	template <typename algflag>
	ReducedChainComplex(
		const ChainComplex<MT> &C,
		algflag,
		bats::clearing_flag
	) {
		size_t dmax = C.maxdim() + 1;
		// dim = C.dim;
		R.resize(dmax);
		p2c.resize(dmax);
		I.resize(dmax);

		// do top dimension normally
		R[dmax-1] = C.boundary[dmax-1];
		p2c[dmax-1] = reduce_matrix(R[dmax-1], algflag());
		std::vector<size_t> clear_inds = get_clearing_inds(p2c[dmax-1]);
		for (ssize_t k = dmax-2; k >= 0; k--) {
			R[k] = C.boundary[k];
			p2c[k] = reduce_matrix_clearing(R[k], clear_inds, algflag());
			clear_inds = get_clearing_inds(p2c[k]);
		}

		set_indices();
	}

	template <typename algflag>
	ReducedChainComplex(
		const ChainComplex<MT> &C,
		algflag,
		bats::compression_flag
	) {
		size_t dmax = C.maxdim() + 1;
		// dim = C.dim;
		R.resize(dmax);
		p2c.resize(dmax);
		I.resize(dmax);

		// do top dimension normally
		R[0] = C.boundary[0];
		p2c[0] = reduce_matrix(R[0], algflag());
		std::vector<bool> comp_inds = get_compression_inds(R[0]);
		// for (auto i : comp_inds) {
		// 	std::cout << i << std::endl;
		// }
		for (size_t k = 1; k < dmax; k++) {
			R[k] = C.boundary[k];
			p2c[k] = reduce_matrix_compression(R[k], comp_inds, algflag());
			comp_inds = get_compression_inds(R[k]);
		}

		set_indices();
	}

	// compute reduced boundary matrices and basis with flags
	template <typename algflag>
	ReducedChainComplex(
		const ChainComplex<MT> &C,
		algflag,
		bats::compression_flag,
		bats::compute_basis_flag
	) {
		size_t dmax = C.maxdim() + 1;
		//dim = C.dim;
		R.resize(dmax);
		U.resize(dmax);
		p2c.resize(dmax);
		I.resize(dmax);

		// do bottom dimension normally
		R[0] = C.boundary[0];
		U[0] = MT::identity(C.dim(0));
		p2c[0] = reduce_matrix(R[0], U[0], algflag());
		std::vector<bool> comp_inds = get_compression_inds(R[0]);
		for (size_t k = 1; k < dmax; k++) {
			R[k] = C.boundary[k];
			U[k] = MT::identity(C.dim(k));
			p2c[k] = reduce_matrix_compression(R[k], U[k], comp_inds, algflag());
			comp_inds = get_compression_inds(R[k]);
		}

		set_indices();
	}

	template <typename... Args>
	void update_reduction2(size_t k, Args ...args) {
		// step 1: make U[k] upper-triangular.  Similar to UQL factorization
		// but we apply updates to R[k] instead of factorizing
		// we can actually use the standard reduction algorithm on U
		// then put it in increasing pivot order
		p2c[k] = reduce_matrix_standard(U[k], R[k]); // we want standard reduction
		// p2c[k] can be used since it will just be updated later.

		// sort columns of U - apply same operations to R
		for (size_t j = 0; j < dim(k); j++) {
			// this loop will eventually terminate, since every pivot occurs
			while (p2c[k][j] != j) {
				// get pivot - we'll swap so this pivot is in correct location
				size_t p = p2c[k][j];
				U[k].swap_cols(p, j);
				R[k].swap_cols(p, j);
				std::swap(p2c[k][p], p2c[k][j]);
			}
		}

		// step 2: finish reduction of matrix R[k]
		p2c[k] = reduce_matrix(R[k], U[k], args...);
	}

	// update reduction B[k] U[k] = R[k]
	template <typename... Args>
	void update_reduction(size_t k, Args ...args) {
		// step 1: ensure U is upper triangular
		auto F = UQL(U[k]); // UQL factorization
		// step 2: move L and Q factors to R
		R[k] = R[k] * l_inv(F.L) * F.E.transpose();
		// step 3: update R
		U[k] = F.U;
		p2c[k] = reduce_matrix(R[k], U[k], args...);
	}

	// permute basis in dimension k
	// B_k U_k = R_k, so when we permute columns of B_k, we must permute rows of U_k
	// we also permute rows of R_{k+1}
	void permute_basis(size_t k, const std::vector<size_t> &perm) {
		auto iperm = bats::util::inv_perm(perm);
		if (k == 0) {
			// only worry about boundary[1]
			R[1].permute_rows(iperm);
		} else if (k == maxdim()) {
			// only need to worry about rows of U[k]
			U[k].permute_rows(iperm); // iperm?
			U[k].permute_cols(perm);
			R[k].permute_cols(perm);
			// U[k].permute_cols(perm);
		} else {
			// need to handle boundary[k] and boundary[k+1]
			U[k].permute_rows(iperm); // iperm?
			U[k].permute_cols(perm);
			R[k].permute_cols(perm);
			// U[k].permute_cols(perm);
			R[k+1].permute_rows(iperm);
		}
		// at end of this, homology classes are invalidated
	}

	// update basis in all dimensions
	template <typename... Args>
	void permute_basis(const std::vector<std::vector<size_t>> &perm, Args ...args) {
		for (size_t k = 0; k < perm.size(); k++) {
			permute_basis(k, perm[k]);
		}
		// next we update the factorizations
		for (size_t k = 0; k < perm.size(); k++) {
			update_reduction2(k, args...);
		}
	}

	// put vector/matrix in homology-revealing basis in dimension k
	template <typename TV>
	inline TV to_hom_basis(const TV &v, size_t k) const {
		return u_solve(U[k], v);
	}

	template <typename TV>
	inline TV from_hom_basis(const TV &v, size_t k) const {
		return U[k] * v;
	}



	// get preferred representative j in dimension k
	chain_type get_preferred_representative(
		const size_t j,
		const size_t k
	) const {
		return U[k][I[k][j]];
	}

	// modify y in-place to be preferred representative for homology class in dimension k
	// assumes y is in homology-revealing basis
	void find_preferred_representative(
		chain_type &y,
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
	chain_type chain_preferred_representative(
		const chain_type &c, size_t k
	) const {
		auto x = to_hom_basis(c, k);
		find_preferred_representative(x, k);
		return from_hom_basis(x, k);
	}

	void print_summary() const {
		std::cout << "ReducedChainComplex with " << maxdim() << " dimensions:" << std::endl;
		for (size_t k = 0; k < maxdim() + 1; k++) {
			std::cout << "\tdim " << k << ": " << dim(k)
			<< ", betti_" << k << ": " << hdim(k) << "\n";
		}
	}
};


// defualt return
template <typename T, typename CpxT, typename... Args>
inline auto __ReducedChainComplex(const CpxT &F, T, Args ...args) {
	using VT = SparseVector<T, size_t>;
	using MT = ColumnMatrix<VT>;

	return ReducedChainComplex(ChainComplex<MT>(F), args...);
}

template <typename T, typename CpxT, typename... Args>
inline auto Reduce(const CpxT &F, T, Args ...args) {
	return __ReducedChainComplex(F, T(), args...);
}

template <typename MT, typename... Args>
inline auto Reduce(const ChainComplex<MT>& C, Args ...args) {
	return ReducedChainComplex(C, args...);
}

} // namespace bats
