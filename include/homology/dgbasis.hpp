#pragma once

/*
compute homology-revealing bases for a chain complex
*/
#include <vector>
#include <set>
#include <dgvs/dgvs.hpp>
#include "reduction.hpp"
#include "parallel.hpp"
#include <chrono>

namespace bats {
// flag to indicate that basis should be computed
// struct compute_basis_flag {};

/*
Maintains factorizations
B_k U_k = R_k
or B_k = R_k U_k^{-1}
*/
template <typename MT>
struct ReducedDGVectorSpace {
public:

	using vect_type = typename MT::col_type;

	//std::vector<size_t> dim;
	int degree; // degree on the differential
	std::vector<MT> U; // basis matrices
	std::vector<MT> R; // reduced matrix
	std::vector<std::vector<size_t>> I;
	std::vector<p2c_type> p2c;


	// size of homology vector space in dimension k
	inline size_t hdim(size_t k) const {
		return (degree == -1) ? I[k].size() : I[k+1].size();
	}
	inline size_t betti(size_t k) const { return hdim(k); }
	inline size_t maxdim() const { return R.size() - 2; }
	inline size_t dim(size_t k) const {
		return (degree == -1) ? R[k].ncol() : R[k+1].ncol();
	}

	MT& operator[](size_t k) {
		if (k >= R.size()) {
			// throw error
		}
		return R[k];
	}

	/**
	initialize with chain complex C, but do not do reduction
	*/
	void initialize(const DGVectorSpace<MT>& C) {
		size_t dmax = C.maxdim() + 1;
		U.resize(dmax);
		R.resize(dmax);
		I.resize(dmax);
		p2c.resize(dmax);

		for (size_t k = 0; k < dmax; ++k) {
			U[k] = MT::identity(C.dim(k));
			R[k] = C[k];
			p2c[k].resize(C[k].nrow(), bats::NO_IND);
		}
	}

	// set homology indices after reduction is completed
	void set_indices() {
		if (degree == -1) {
			// homological type
			for (size_t k = 0; k < maxdim()+1; k++) {
				I[k] = extract_basis_indices(R[k], p2c[k+1]);
			}
			I[maxdim()+1] = extract_basis_indices(R[maxdim()]);
		} else if (degree == +1) {
			// cohomological type
			for (size_t k = 1; k < maxdim()+2; k++) {
				I[k] = extract_basis_indices(R[k], p2c[k-1]);
			}
		}

	}


	ReducedDGVectorSpace() {}

	// compute reduced chain complex from chain complex
	ReducedDGVectorSpace(const DGVectorSpace<MT> &C) : degree(C.degree) {
		size_t dmax = C.differential.size();
		//dim = C.dim;
		U.resize(dmax);
		R.resize(dmax);
		I.resize(dmax);
		p2c.resize(dmax);

		// TODO: can parallelize this
		for (size_t k = 0; k < dmax; k++) {
			size_t dimk = C.differential[k].ncol();
			U[k] = MT::identity(dimk);
			R[k] = C.differential[k];
			// partial_reduce_parallel(R[k], U[k], 1024);
			p2c[k] = reduce_matrix(R[k], U[k]);
		}

		set_indices();

	}


	/**
	put vector/matrix in homology-revealing basis in dimension k
	*/
	template <typename TV>
	inline TV to_hom_basis(const TV &v, size_t k) const {
		return u_solve(U[k], v);
	}

	/**
	put vector/matrix back in original basis in dimension k
	*/
	template <typename TV>
	inline TV from_hom_basis(const TV &v, size_t k) const {
		return U[k] * v;
	}



	// get preferred representative j in dimension k
	vect_type get_preferred_representative(
		const size_t j,
		const size_t k
	) const {
		return U[k][I[k][j]];
	}

	// modify y in-place to be preferred representative for homology class in dimension k
	// assumes y is in homology-revealing basis
	void find_preferred_representative(
		vect_type &y,
		size_t k
	) const {
		int dg_offset = (degree == -1) ? 0 : 1;
		k += dg_offset;
		if (k == R.size()-1) {
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
	vect_type chain_preferred_representative(
		const vect_type &c, size_t k
	) const {
		auto x = to_hom_basis(c, k);
		find_preferred_representative(x, k);
		return from_hom_basis(x, k);
	}

	void print_summary(bool print_nnz=false) const {
		std::cout << "ReducedDGVectorSpace, maxdim =  " << maxdim() << std::endl;
		for (size_t k = 0; k < maxdim() + 1; k++) {
			std::cout << "\tdim " << k << ": " << dim(k)
			<< ", betti_" << k << ": " << hdim(k);
			if (print_nnz) {
				std::cout << " nnz(R): " << R[k].nnz();
				if (U.size() > k) { // handle if basis not computed
					std::cout << " nnz(U): " << U[k].nnz();
				}
			}
			std::cout << "\n";
		}
	}
};


} // namespace bats
