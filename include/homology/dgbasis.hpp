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

	// compute reduced boundary matrices only with flags
	template <typename algflag>
	ReducedDGVectorSpace(const DGVectorSpace<MT> &C, algflag) : degree(C.degree) {
		size_t dmax = C.differential.size();
		//dim = C.dim;
		R.resize(dmax);
		p2c.resize(dmax);
		I.resize(dmax);

		// TODO: can parallelize this
		for (size_t k = 0; k < dmax; k++) {
			R[k] = C.differential[k];
			// partial_reduce_parallel(R[k], 1024);
			p2c[k] = reduce_matrix(R[k], algflag());
		}

		set_indices();
	}

	// compute reduced boundary matrices with basis and flags
	template <typename algflag>
	ReducedDGVectorSpace(
		const DGVectorSpace<MT> &C,
		 algflag,
		 bats::compute_basis_flag
	 ) : degree(C.degree) {
		size_t dmax = C.differential.size();
		//dim = C.dim;
		U.resize(dmax);
		R.resize(dmax);
		p2c.resize(dmax);
		I.resize(dmax);

		// TODO: can parallelize this
		for (size_t k = 0; k < dmax; k++) {
			size_t dimk = C.differential[k].ncol();
			U[k] = MT::identity(dimk);
			R[k] = C.differential[k];
			// partial_reduce_parallel(R[k], 1024);
			p2c[k] = reduce_matrix(R[k], U[k], algflag());
		}

		set_indices();
	}

	template <typename algflag>
	ReducedDGVectorSpace(
		const DGVectorSpace<MT> &C,
		algflag,
		bats::clearing_flag
	) : degree(C.degree) {
		size_t dmax = C.differential.size();
		// dim = C.dim;
		R.resize(dmax);
		p2c.resize(dmax);
		I.resize(dmax);

		if (degree == -1) {
			// do top dimension normally
			R[dmax-1] = C.differential[dmax-1];
			p2c[dmax-1] = reduce_matrix(R[dmax-1], algflag());
			std::vector<size_t> clear_inds = get_clearing_inds(p2c[dmax-1]);
			for (ssize_t k = dmax-2; k >= 0; --k) {
				R[k] = C.differential[k];
				p2c[k] = reduce_matrix_clearing(R[k], clear_inds, algflag());
				clear_inds = get_clearing_inds(p2c[k]);
			}
		} else { // degree == +1
			// do bottom dimension normally
			R[0] = C.differential[0];
			p2c[0] = reduce_matrix(R[0], algflag());
			std::vector<size_t> clear_inds = get_clearing_inds(p2c[0]);
			for (ssize_t k = 1; k < dmax; ++k) {
				R[k] = C.differential[k];
				p2c[k] = reduce_matrix_clearing(R[k], clear_inds, algflag());
				clear_inds = get_clearing_inds(p2c[k]);
			}
		}


		set_indices();
	}

	template <typename algflag>
	ReducedDGVectorSpace(
		const DGVectorSpace<MT> &C,
		algflag,
		bats::clearing_flag,
		bats::compute_basis_flag
	) : degree(C.degree) {
		size_t dmax = C.differential.size();
		// dim = C.dim;
		U.resize(dmax);
		R.resize(dmax);
		p2c.resize(dmax);
		I.resize(dmax);

		if (degree == -1) {
			// do top dimension normally
			U[dmax-1] = MT::identity(C.differential[dmax-1].ncol());
			R[dmax-1] = C.differential[dmax-1];
			p2c[dmax-1] = reduce_matrix(R[dmax-1], U[dmax-1], algflag());
			for (ssize_t k = dmax-2; k >= 0; --k) {
				size_t dimk = C.differential[k].ncol();
				U[k] = MT::identity(dimk);
				R[k] = C.differential[k];
				p2c[k] = reduce_matrix_clearing(R[k], U[k], R[k+1], p2c[k+1], algflag());
			}
		} else { // degree == +1
			// do bottom dimension normally
			U[0] = MT::identity(C.differential[0].ncol());
			R[0] = C.differential[0];
			p2c[0] = reduce_matrix(R[0], U[0], algflag());
			for (ssize_t k = 1; k < dmax; ++k) {
				size_t dimk = C.differential[k].ncol();
				U[k] = MT::identity(dimk);
				R[k] = C.differential[k];
				p2c[k] = reduce_matrix_clearing(R[k], U[k], R[k-1], p2c[k-1], algflag());
			}
		}


		set_indices();
	}

	template <typename algflag>
	ReducedDGVectorSpace(
		const DGVectorSpace<MT> &C,
		algflag,
		bats::compression_flag
	) : degree(C.degree) {
		size_t dmax = C.differential.size();
		// dim = C.dim;
		R.resize(dmax);
		p2c.resize(dmax);
		I.resize(dmax);

		if (degree == -1) {
			// do bottom dimension normally
			R[0] = C.differential[0];
			p2c[0] = reduce_matrix(R[0], algflag());
			std::vector<bool> comp_inds = get_compression_inds(R[0]);
			for (size_t k = 1; k < dmax; k++) {
				R[k] = C.differential[k];
				p2c[k] = reduce_matrix_compression(R[k], comp_inds, algflag());
				comp_inds = get_compression_inds(R[k]);
			}
		} else { // degree == +1
			// do top dimension normally
			R[dmax-1] = C.differential[dmax-1];
			p2c[dmax-1] = reduce_matrix(R[dmax-1], algflag());
			std::vector<bool> comp_inds = get_compression_inds(R[dmax-1]);
			for (ssize_t k = dmax-2; k >= 0; --k) {
				R[k] = C.differential[k];
				p2c[k] = reduce_matrix_compression(R[k], comp_inds, algflag());
				comp_inds = get_compression_inds(R[k]);
			}
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

	template <typename... Args>
	void update_reduction2(size_t k, Args ...args) {

		size_t k_ind = degree == -1 ? k : k + 1;
		// step 1: make U[k] upper-triangular.  Similar to UQL factorization
		// but we apply updates to R[k] instead of factorizing
		// we can actually use the standard reduction algorithm on U
		// then put it in increasing pivot order
		p2c[k_ind] = reduce_matrix_standard(U[k_ind], R[k_ind]); // we want standard reduction
		// p2c[k] can be used since it will just be updated later.

		// sort columns of U - apply same operations to R
		for (size_t j = 0; j < dim(k); j++) {
			// swap correct column if necessary
			if (p2c[k_ind][j] != j) {
				// get pivot - we'll swap so this pivot is in correct location
				size_t pj = U[k_ind][j].lastnz().ind; // pivot of column j
				size_t p = p2c[k_ind][j]; // location of desired pivot
				U[k_ind].swap_cols(p, j);
				R[k_ind].swap_cols(p, j);
				p2c[k_ind][j] = j; // we've put the pivot in the right place
				p2c[k_ind][pj] = p; // set pivot lookup to j
			}
		}

		// step 2: finish reduction of matrix R[k]
		p2c[k_ind] = reduce_matrix(R[k_ind], U[k_ind], args...);
	}

	// permute basis in dimension k
	// B_k U_k = R_k, so when we permute columns of B_k, we must permute rows of U_k
	// we also permute rows of R_{k+1}
	void permute_matrices(size_t k, const std::vector<size_t> &perm) {
		auto iperm = bats::util::inv_perm(perm);
		if (degree == -1) {
			if (k == 0) {
				// only worry about boundary[1]
				R[1].permute_rows(iperm);
			} else if (k == maxdim()) {
				// only need to worry about rows of U[k]
				U[k].permute_rows(iperm); // inverse perm here
			} else {
				// need to handle boundary[k] and boundary[k+1]
				U[k].permute_rows(iperm); // inverse perm here
				R[k+1].permute_rows(iperm); // regular perm here
			}
		} else { // degree == +1
			// permute rows of U^k = U[k+1]
			// permute rows of R^{k-1} = R[k]
			if (k == 0) {
				// only worry about differential[1]
				U[1].permute_rows(iperm);
			} else if (k == maxdim()) {
				// only need to worry about differential[k]
				R[k].permute_rows(iperm); // inverse perm here
			} else {
				// need to handle boundary[k] and boundary[k+1]
				U[k+1].permute_rows(iperm); // inverse perm here
				R[k].permute_rows(iperm); // regular perm here
			}
		}

		// at end of this, homology classes are invalidated
	}

	// update basis in all dimensions
	template <typename... Args>
	void permute_basis(const std::vector<std::vector<size_t>> &perm, Args ...args) {
		for (size_t k = 0; k < perm.size(); k++) {
			permute_matrices(k, perm[k]);
		}
		// next we update the factorizations
		for (size_t k = 0; k < perm.size(); k++) {
			update_reduction2(k, args...);
		}
		set_indices(); // update homology indices
	}

		// update basis in all dimensions by new filtration information
	template <typename Information_type, typename... Args>
	void update_basis_general(const Information_type &UI, Args ...args){
		using VectT = typename MT::col_type; // column vector type
		using ValT = typename VectT::val_type;

		// step 1. compute the permutation that
		// (a) move the deleted rows/columns to the end
		// (b) permute the intersection part of two filtrations
		size_t max_dim = UI.max_dim;

		// we want to extract permutation and deletion information first
		for (size_t k = 0; k < max_dim +1; ++k) {
			// step 1.1 find the permutation used to permute simplices
			// that is going to be deleted to the end
			std::vector<size_t> perm_deletion;
			if(!UI.deletion_indices[k].empty()){
				perm_deletion = perm_to_the_end(UI.deletion_indices[k], R[k].ncol());
			}else{
				perm_deletion = identity_perm(R[k].ncol());
			}

			// step 1.2 find the permutation used to permute
			// intersection of simplices and leave the ones in the end unmoved
			std::vector<size_t> perm_intersect = extension_perm(UI.permutations[k], R[k].ncol());

			// step 1.3 Combine the above 2 permutations together, i.e.,
			// compute the permutation that first do deletion permutation,
			// and then do inersection permutation.
			// We will use the perm_deletion to store the result.
			bats::util::apply_perm(perm_deletion, perm_intersect);

			// 1.4 permute U and R
			permute_matrices(k, perm_deletion);
		}

		// step 2: next we update the factorizations
		for (size_t k = 0; k < max_dim+1; ++k) { // for each dimension
			// step 2.1 make U reduced
			p2c[k] = reduce_matrix_standard(U[k], R[k]);

			// step 2.2 sort columns of U - apply same operations to R
			// to make U upper-triangular
			for (size_t j = 0; j < dim(k); j++) {
				// swap correct column if necessary
				if (p2c[k][j] != j) {
					// get pivot - we'll swap so this pivot is in correct location
					size_t pj = U[k][j].lastnz().ind; // pivot of column j
					size_t p = p2c[k][j]; // location of desired pivot
					U[k].swap_cols(p, j);
					R[k].swap_cols(p, j);
					p2c[k][j] = j; // we've put the pivot in the right place
					p2c[k][pj] = p; // set pivot look up to j
				}
			}

			// std::cout << "\nstep 2.3 delete columns in the end:" << std::endl;
			// step 2.3 delete columns in the end
			for (size_t i = 0; i < UI.deletion_indices[k].size(); i++){
				U[k].erase_column();
				U[k].erase_row_unsafe(); // after column deletion, safe to erase
				R[k].erase_column();
			}

			// step 2.4 delete rows in the end
			if(k!=0){
				// find the # of deletion
				size_t count = UI.deletion_indices[k-1].size();
				for (size_t i = 0; i < count; i++){
					R[k].erase_row_unsafe();
				}
			}


			// step 2.5 add rows to the specified location
			if(k!=0){
				// for(auto& ind: UI.addition_indices[k-1]){
				// 	R[k].insert_row(ind); // insert zero rows
				// }
				R[k].insert_rows(UI.addition_indices[k-1]);
			}


			// step 2.6 addition of columns
			if(k==0){
				// find the index of the rightmost column of U
				size_t final_ind_of_U = U[k].ncol() - 1;

				// find the # of addition
				size_t count = UI.addition_indices[k].size();

				for (size_t i = 0; i < count; i++){
					// change R
					R[k].append_column(); // insert zero columns
					// change U
					U[k].append_row();
					U[k].append_column(VectT(size_t(final_ind_of_U)));
					final_ind_of_U++;
				}
			}else{
				// boundary information
				const std::vector<std::vector<size_t>>& bd_info = UI.boundary_indices[k];

				// addition information
				const std::vector<size_t>& add_inds = UI.addition_indices[k];

				U[k].insert_rows(add_inds);

				std::vector<VectT> Rcol(add_inds.size());
				std::vector<VectT> Ucol(add_inds.size());
				// Then change U and R coming from the effect of
				// adding a column to the boundary matrix.
				// loop over each simplex that needs to be added
				for(size_t i = 0; i < add_inds.size(); i++){
					// the indices of its boundaries
					auto simplex_bd_ind = bd_info[i];
					// std::sort(simplex_bd_ind.begin(), simplex_bd_ind.end());
					// simplex index
					auto ind = add_inds[i];

					// create a vector of ones with length
					// equal to the boundary size  (Field F_2 for now)
					std::vector<ValT> vect_one(simplex_bd_ind.size(), 1);
					// creat the column vector
					Rcol[i] = VectT(simplex_bd_ind, vect_one);
					Ucol[i] = VectT(size_t(ind));
				}
				R[k].insert_columns(add_inds, Rcol);
				U[k].insert_columns(add_inds, Ucol);

				// // Then change U and R coming from the effect of
				// // adding a column to the boundary matrix.
				// // loop over each simplex that needs to be added
				// for(size_t i = 0; i < add_inds.size(); i++){
				// 	// the indices of its boundaries
				// 	auto simplex_bd_ind = bd_info[i];
				// 	// std::sort(simplex_bd_ind.begin(), simplex_bd_ind.end());
				// 	// simplex index
				// 	auto ind = add_inds[i];
				//
				// 	// create a vector of ones with length
				// 	// equal to the boundary size  (Field F_2 for now)
				// 	std::vector<ValT> vect_one(simplex_bd_ind.size(), 1);
				// 	// creat the column vector
				// 	auto vect = VectT(simplex_bd_ind, vect_one);
				// 	// insert the column at postion in B (position in new filtration)
				// 	R[k].insert_column(ind, vect);
				//
				// 	// correspondingly, U[k] will add a zero block
				// 	// U[k].insert_row(ind);
				// 	U[k].insert_column(ind, VectT({size_t(ind)}, {1}));
				// }
			}

			// step 2.7, make R[k] reduced
			p2c[k] = reduce_matrix(R[k], U[k], args...);
		}

		set_indices(); // update homology indices
	}


};


} // namespace bats
