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
			for (size_t k = 1; k < dmax; ++k) {
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
			for (size_t k = 1; k < dmax; ++k) {
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

	/**
	helper function for updates

	Reduces U, permutes to be upper triangular
	applies same column operations to R

	WARNING: after this, p2c[k] will contain pivots of U
	*/
	void _make_U_upper_triangular(size_t k) {
		// we can actually use the standard reduction algorithm on U
		// then put it in increasing pivot order
		p2c[k] = reduce_matrix_standard(U[k], R[k]); // we want standard reduction

		// sort columns of U - apply same operations to R
		for (size_t j = 0; j < U[k].ncol(); j++) {
			// swap correct column if necessary
			if (p2c[k][j] != j) {
				// get pivot - we'll swap so this pivot is in correct location
				size_t pj = U[k][j].lastnz().ind; // pivot of column j
				size_t p = p2c[k][j]; // location of desired pivot
				U[k].swap_cols(p, j);
				R[k].swap_cols(p, j);
				p2c[k][j] = j; // we've put the pivot in the right place
				p2c[k][pj] = p; // set pivot lookup to j
			}
		}
		return;
	}

	template <typename... Args>
	void update_reduction2(size_t k, Args ...args) {

		size_t k_ind = degree == -1 ? k : k + 1;
		// step 1: make U[k] upper-triangular.  Similar to UQL factorization
		// but we apply updates to R[k] instead of factorizing
		_make_U_upper_triangular(k_ind);

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

	// find the reverse index of simplex in persistent cohomology from homology
	// where n is the column size
	size_t find_reverse_index(const size_t& i, const size_t& n){
		return n - i - 1;
	}

	// find the reverse indices of several simplices in persistent cohomology from homology
	std::vector<size_t> find_reverse_index(const std::vector<size_t>& index_list, const size_t& n){
		std::vector<size_t> v;
		v.reserve(index_list.size());
		for(size_t ele: index_list){
			v.emplace_back(find_reverse_index(ele,n));
		}
		return v;
	}

	// find the reverse indices of boundary simplices
	std::vector<std::vector<size_t>> find_reverse_index(const std::vector<std::vector<size_t>> & index_list, const size_t& n){
		std::vector<std::vector<size_t>> v;
		v.reserve(index_list.size());
		for(std::vector<size_t> ele: index_list){
			v.emplace_back(find_reverse_index(ele,n));
		}
		return v;
	}

	/**
	delete cells in dimension k
	*/
	void _delete_cells(size_t k, const UpdateInfo2& UI) {
		if (UI.ndeletions[k] == 0) { return; } // quick exit

		if (degree == -1) {
			if (k < UI.perm.size()-1) {
				R[k+1].erase_final_rows_unsafe(UI.ndeletions[k]);
			}
			U[k].erase_final_columns(UI.ndeletions[k]);
			U[k].erase_final_rows_unsafe(UI.ndeletions[k]);
			R[k].erase_final_columns(UI.ndeletions[k]);
		} else { // degree == +1
			// delete initial columns of R[k+1]
			R[k+1].erase_initial_cols(UI.ndeletions[k]);
			// delete initial columns of U[k+1]
			U[k+1].erase_initial_cols(UI.ndeletions[k]);
			// delete initial rows of U[k+1]
			U[k+1].erase_initial_rows(UI.ndeletions[k]);

			// delete rows one dimension down
			R[k].erase_initial_rows(UI.ndeletions[k]);
		}
		return;
	}

	void _insert_cells(size_t k, const UpdateInfo2& UI) {
		using VectT = typename MT::col_type; // column vector type
		using ValT = typename VectT::val_type; // field type
		if (UI.insertion_indices[k].size() == 0) { return; } // quick exit
		std::cout << "\t\t_insert_cells " << k << std::endl;
		if (degree == -1) {
			if ( k < UI.perm.size() - 1 ) {
				// insert rows one dimension up
				R[k+1].insert_rows(UI.insertion_indices[k]);
			}

			// create vectors of columns to insert
			std::vector<VectT> Ucols;
			Ucols.reserve(UI.insertion_indices[k].size());
			std::vector<VectT> Rcols;
			Rcols.reserve(UI.insertion_indices[k].size());
			for (size_t i = 0; i < UI.insertion_indices[k].size(); ++i) {
				Ucols.emplace_back(VectT(size_t(UI.insertion_indices[k][i])));
				Rcols.emplace_back(UI.insertion_cols[k][i].cast_values(ValT()));
			}

			R[k].insert_columns(UI.insertion_indices[k], Rcols);
			U[k].insert_rows(UI.insertion_indices[k]); //zero rows
			U[k].insert_columns(UI.insertion_indices[k], Ucols);

		} else { // degree == +1

			// insert zero columns into R[k+1]
			auto t0 = std::chrono::steady_clock::now();
			std::vector<VectT> Rcols(UI.insertion_indices[k].size(), VectT());
			R[k+1].insert_columns(UI.insertion_indices[k], Rcols);
			auto t1 = std::chrono::steady_clock::now();
	        std::cout << "\t\t  insert Rcol "
	            << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
	            << "ms" << std::endl;

			// row insertion tends to be really slow for cohomology
			// for now, just do transpose
			t0 = std::chrono::steady_clock::now();
			U[k+1].insert_rows(UI.insertion_indices[k]);
			t1 = std::chrono::steady_clock::now();
	        std::cout << "\t\t  insert Urow "
	            << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
	            << "ms" << std::endl;
			// U[k+1] = U[k+1].transpose();
			// std::vector<VectT> Urows(UI.insertion_indices[k].size(), VectT());
			// U[k+1].insert_columns(UI.insertion_indices[k], Urows);
			// U[k+1] = U[k+1].transpose();

			// create vectors of columns to insert
			t0 = std::chrono::steady_clock::now();
			std::vector<VectT> Ucols;
			Ucols.reserve(UI.insertion_indices[k].size());
			for (size_t i = 0; i < UI.insertion_indices[k].size(); ++i) {
				Ucols.emplace_back(VectT(size_t(UI.insertion_indices[k][i])));
			}
			U[k+1].insert_columns(UI.insertion_indices[k], Ucols); // identity columns
			t1 = std::chrono::steady_clock::now();
			std::cout << "\t\t  insert Ucol "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
				<< "ms" << std::endl;

			// actually need to insert boundary transpose
			// if we want to insert row di, need to form ri = di V
			// => ri.T = V.T di.T
			t0 = std::chrono::steady_clock::now();
			auto Vt = U[k].transpose();
			t1 = std::chrono::steady_clock::now();
			std::cout << "\t\t  Vt = U[k].transpose(): "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
				<< "ms" << std::endl;
			// vector of rows to insert
			// we'll actually insert as columns of the transpose
			// std::vector<VectT> Rcols;
			t0 = std::chrono::steady_clock::now();
			Rcols.clear();
			Rcols.reserve(UI.insertion_indices[k].size());
			typename VectT::tmp_type tmp;
			for (size_t i = 0; i < UI.insertion_indices[k].size(); ++i) {
				auto Di = UI.insertion_cols[k][i].cast_values(ValT());
				Rcols.emplace_back(Vt.gemv(Di, tmp));
			}
			t1 = std::chrono::steady_clock::now();
			std::cout << "\t\t  form Rrows : "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
				<< "ms" << std::endl;

			// t0 = std::chrono::steady_clock::now();
			// R[k] = R[k].transpose(); // take transpose before insertion
			// t1 = std::chrono::steady_clock::now();
			// std::cout << "\t\t  transpose 1 : "
			// 	<< std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
			// 	<< "ms" << std::endl;
			// t0 = std::chrono::steady_clock::now();
			// R[k].insert_columns(UI.insertion_indices[k], Rcols);
			// t1 = std::chrono::steady_clock::now();
			// std::cout << "\t\t  insert Rrows : "
			// 	<< std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
			// 	<< "ms" << std::endl;
			// t0 = std::chrono::steady_clock::now();
			// R[k] = R[k].transpose(); // transpose again so we have inserted rows
			// t1 = std::chrono::steady_clock::now();
			// std::cout << "\t\t  transpose 2 : "
			// 	<< std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
			// 	<< "ms" << std::endl;
			t0 = std::chrono::steady_clock::now();
			R[k].insert_rows(UI.insertion_indices[k], Rcols);
			t1 = std::chrono::steady_clock::now();
			std::cout << "\t\t  insert Rrows : "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
				<< "ms" << std::endl;
		}
	}
	/**
	General update factorization
	*/
	void update_basis(
		const UpdateInfo2& UI
	) {
		// using VectT = typename MT::col_type; // column vector type
		// using ValT = typename VectT::val_type; // field type

		// step 1: perform permutations
		for (size_t k = 0; k < UI.perm.size(); ++k) {
			permute_matrices(k, UI.perm[k]);
		}

		for (size_t k = 0; k < UI.perm.size(); ++k) {
			size_t k_ind = degree == -1 ? k : k + 1;
			// step 2: correct U matrices
			_make_U_upper_triangular(k_ind);
			// step 3: process deletions
			_delete_cells(k, UI);
			// step 4: process insertions
			_insert_cells(k, UI);
		}

		// finalize reduction
		if (degree == -1) {
			for (ssize_t k = UI.perm.size()-1; k >= 0; --k) {
				if (k < (ssize_t) UI.perm.size() - 1) {
					p2c[k] = reduce_matrix_clearing(R[k], U[k], R[k+1], p2c[k+1], bats::standard_reduction_flag());
				} else {
					p2c[k] = reduce_matrix(R[k], U[k]);
				}
			}
		} else { // degree == +1
			p2c[0] = reduce_matrix(R[0], U[0]);
			for (size_t k = 1; k < UI.perm.size() + 1; ++k) {
				p2c[k] = reduce_matrix_clearing(R[k], U[k], R[k-1], p2c[k-1], bats::standard_reduction_flag());
			}
		}


		set_indices(); // update homology indices
	}

	/*
	General Update factorization DU = R by Updating Information
	Note:
	1) when degree = +1 ,i.e., cohomology, we need to
	 modify deletion and addition indices: i <-> n-i-1 , where n is column size
	2) add/remove rows of (co)boundary matrix do not need to modify U,
	while add/remove columns need. However, adding a row requires a mutlplication of it with U matrix.
	3) There are types of permutations used in BATs and I call them:
			for a vector perm and a vector v
			i) new_ind_perm: perm[i] is the new index of v[i]
	 	ii) math_perm: inverse of above
	, in the following implementation, permutations use math_perm
	(I know they are confusing, but changing it is hard...).
	**/
	template <typename Information_type, typename... Args>
	void update_basis_general(const Information_type &UI, Args ...args){
		using VectT = typename MT::col_type; // column vector type
		using ValT = typename VectT::val_type;

		// step 1. compute the permutation that
		// (a) move the deleted rows/columns to the end
		// (b) permute the intersection part of two filtrations
		size_t max_dim = UI.max_dim;
		// std::cout << "before permutation R[1] is" << std::endl;
		// R[1].print();
		// std::cout << "and U[1] is"  << std::endl;
		// U[1].print();

		// we want to extract permutation and deletion information first
		auto t0 = std::chrono::steady_clock::now();
		for (size_t k = 0; k < max_dim + 1; ++k) {
			size_t m = R[k].nrow();
			size_t n = R[k].ncol();
			// step 1.1 find the permutation used to permute simplices
			// that is going to be deleted to the end
			std::vector<size_t> perm_deletion;
			if (degree == +1){
				perm_deletion = identity_perm(m); // cohomology need to permute row
			}else{
				perm_deletion = identity_perm(n);
			}

			if(!UI.deletion_indices[k].empty()){
				if (degree == +1) {
					std::vector<size_t> deletion_inds = find_reverse_index(UI.deletion_indices[k], m);
					std::sort(deletion_inds.begin(), deletion_inds.end()); // sort it
					perm_deletion = perm_to_the_end(deletion_inds, m);
				}else{
					perm_deletion = perm_to_the_end(UI.deletion_indices[k], n);
				}
			}

			// step 1.2 find the permutation used to permute
			// intersection of simplices and leave the ones in the end unmoved
			std::vector<size_t> perm_intersect;
			if (degree == +1) {
				// auto perm_extended = extension_perm(UI.permutations[k], m);
				// reverse
				std::vector<size_t> perm_rever = find_reverse_index(UI.permutations[k], UI.permutations[k].size());
				std::reverse(perm_rever.begin(), perm_rever.end());
				perm_intersect = extension_perm(perm_rever, m);
			}else{
				perm_intersect = extension_perm(UI.permutations[k], n);
			}


			// step 1.3 Combine the above 2 permutations together, i.e.,
			// compute the permutation that first do deletion permutation,
			// and then do inersection permutation (it is valid since they are disjoint).
			// We will use the perm_deletion to store the result.
			bats::util::apply_perm(perm_deletion, perm_intersect);

			// 1.4 permute U and R
			permute_matrices(k, perm_deletion);
		}
		auto t1 = std::chrono::steady_clock::now();
        std::cout << "\tStep 1. Permuting Basis Matrices";
        std::cout << " takes "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
            << "ms" << std::endl;

		// std::cout << "after permutation, R[1] is" << std::endl;
		// R[1].print();
		// std::cout << "and U[1] is"  << std::endl;
		// U[1].print();
		// std::cout << "check R[1] U[1]^{-1} is " << std::endl;

		// step 2: next we update the factorizations
		for (size_t k = 0; k < max_dim+1; ++k) { // for each dimension
			size_t k_ind = degree == -1 ? k : k + 1;

			_make_U_upper_triangular(k_ind);

			t0 = std::chrono::steady_clock::now();
			// step 2.3 delete k-simplices in the end(column for homology and row for cohomology)
			if (degree == -1){	// cohomology remove block
				U[k].erase_final_columns(UI.deletion_indices[k].size());
				U[k].erase_final_rows_unsafe(UI.deletion_indices[k].size());
				R[k].erase_final_columns(UI.deletion_indices[k].size());
			}
			if(degree == +1){ // cohomology remove rows
				R[k].erase_rows(UI.deletion_indices[k].size());
			}
			t1 = std::chrono::steady_clock::now();
			std::cout << "\tStep 2.3 Removing "<< k << "-simplices";
			std::cout << " takes "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
				<< "ms" << std::endl;

			// if(k==1){
			// 	std::cout << "step2.3 after removing the last "<< UI.deletion_indices[k].size() << " rows of R[1]" << std::endl;
			// 	R[1].print();
			// }


			// step 2.4 delete (k-1)-simplices in the end(row for homology and column for cohomology)
			if(k!=0 and UI.deletion_indices[k-1].size() > 0){
				t0 = std::chrono::steady_clock::now();
				if (degree == -1){
					R[k].erase_final_rows_unsafe(UI.deletion_indices[k-1].size());

				}else{
					// std::cout << "before removing the last columns of R[k]" << std::endl;
					// R[k].print();
					// std::cout << "and U[k] is"  << std::endl;
					// U[k].print();

					U[k].erase_final_columns(UI.deletion_indices[k-1].size());
					U[k].erase_final_rows_unsafe(UI.deletion_indices[k-1].size());
					R[k].erase_final_columns(UI.deletion_indices[k-1].size());
					// if(k == 1){
						// std::cout << "step2.4 Delete last "<< UI.deletion_indices[k-1].size()<< " columns of R[1] " << std::endl;
						// std::cout << "after removing the last columns of R[k]" << std::endl;
						// R[k].print();
						// std::cout << "and U[k] is"  << std::endl;
						// U[k].print();
					// }
				}
				t1 = std::chrono::steady_clock::now();
				std::cout << "\tStep 2.4 Removing "<< k-1 << "-simplices";
				std::cout << " takes "
					<< std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
					<< "ms" << std::endl;

			}

			// step 2.5 add (k-1)-simplices (zero rows for homology and columns for cohomology)
			// to the specified locations in updating information
			if(k!=0 and  UI.addition_indices[k-1].size() > 0){
				t0 = std::chrono::steady_clock::now();
				std::vector<size_t> add_inds = UI.addition_indices[k-1];
				if (degree == +1){ // cohomology add zero columns in R and 1 block in U
					// std::cout << "\nstep2.5 need to insert new 1 block to U[k]" << std::endl;
					add_inds = find_reverse_index(add_inds, R[k].ncol() + add_inds.size()); // reverse it
					std::sort(add_inds.begin(), add_inds.end()); // sort indices before insert

					// std::cout << "addition indices are" << std::endl;
					// print_1D_vectors(add_inds);

					std::vector<vect_type> zero_cols(add_inds.size());
					R[k].insert_columns(add_inds, zero_cols);

					std::vector<VectT> Ucols(add_inds.size());
					for(size_t i = 0; i < add_inds.size(); i++){
						// U[k] needs a new one block
						U[k].insert_row(add_inds[i]); //insert zero row of U
						Ucols[i] = VectT(size_t(add_inds[i]));
					}
					U[k].insert_columns(add_inds, Ucols);
				}else{ // homology insert zero rows at specified locations
					R[k].insert_rows(add_inds);
				}
				// R[1].print(); U[1].print();
				t1 = std::chrono::steady_clock::now();
				std::cout << "\tStep 2.5 Add "<<k-1<<"-simplices takes";
				std::cout << " takes "
					<< std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
					<< "ms" << std::endl;
			}


			t0 = std::chrono::steady_clock::now();
			// step 2.6 add k-simplices (columns for homology, row for cohomology)
			if (!UI.addition_indices[k].empty()){
				if(k==0){
					// find the # of addition
					size_t count = UI.addition_indices[k].size();
					if (degree == -1){ // homology add columns
						// find the index of the rightmost column of U
						size_t final_ind_of_U = U[k].ncol() - 1;
						// append zero columns
						for (size_t i = 0; i < count; i++){
							R[k].append_column();
							U[k].append_row();
							// create an indicator (one) column vector at a given index
							U[k].append_column(VectT(size_t(final_ind_of_U)));
							final_ind_of_U++;
						}
					}else{ // cohomology append zero rows
						for (size_t i = 0; i < count; i++){
							R[k].append_row();
						}
					}
				}else{ // k >= 1
					std::vector<std::vector<size_t>> bd_info; // boundary information
					std::vector<size_t> add_inds; // addition indices
					if (degree == -1){ // homology add columns
						bd_info = UI.boundary_indices[k];
						add_inds = UI.addition_indices[k]; // ascending order in Y

						U[k].insert_rows(add_inds); //zero rows

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
					}
					else{ // cohomology add rows (!Notice: need to assert in ascending order of indices)
						bd_info = find_reverse_index(UI.boundary_indices[k], R[k].ncol());
						add_inds = find_reverse_index(UI.addition_indices[k], R[k].nrow() + UI.addition_indices[k].size());
						// we need to sort addition indices and corresponding permute their boundary indices
						std::vector<size_t> p = bats::util::sortperm(add_inds);
						bats::util::apply_perm_swap(add_inds, bats::util::inv_perm(p));
						bats::util::apply_perm_swap(bd_info, bats::util::inv_perm(p));

                        auto U_transpose = U[k].T();
                        std::vector<VectT> new_rows_in_R;
                        for(size_t i = 0; i < add_inds.size(); i++){
                            // the indices of its boundaries
                            auto simplex_bd_ind = bd_info[i];
                            // boundary indices are reversed but inserting rows require ascending order of indices
                            std::sort(simplex_bd_ind.begin(), simplex_bd_ind.end());

                            // create sparse vector of new row in D
                            std::vector<ValT> bd_values(simplex_bd_ind.size());
                            std::fill(bd_values.begin(),bd_values.end(), ValT(1)); //F2 for now
                            SparseVector<ValT, size_t> new_row_in_D(simplex_bd_ind, bd_values);
							// compute new rows in R
                            auto new_row_in_R = U_transpose.gemv(new_row_in_D);
                            new_rows_in_R.emplace_back(new_row_in_R);
                        }
                        R[k].insert_sparse_rows(add_inds, new_rows_in_R);

                        //  // if (k== 1){
                        //  //  std::cout << "after addition" << std::endl;
                        //  //  std::cout << "R[1]" << std::endl;
                        //  //  R[1].print();
                        //  // }
                        // }

						// if (k == 1){
						// 	std::cout << "\nWhen k = "<< k << std::endl;
						// 	std::cout << "Addition column indices" << std::endl;
						// 	print_1D_vectors(UI.addition_indices[k]);
						// 	std::cout << "Transformed to Need to add addition row indices" << std::endl;
						// 	print_1D_vectors(add_inds);
						// 	std::cout << "Need to add addition boudaries indices(column indices)" << std::endl;
						// 	print_2D_vectors(bd_info);
						// 	std::cout << "before addition" << std::endl;
						// 	std::cout << "R[k]" << std::endl;
						// 	R[1].print();
						// }

					}
				}
			}
			t1 = std::chrono::steady_clock::now();
			std::cout << "\tStep 2.6 Add "<<k<<"-simplices takes";
			std::cout << " takes "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
				<< "ms" << std::endl;
			// if(k == 1){
			// 	std::cout << "\nWhen k = "<< k << std::endl;
			// 	std::cout << "after addition" << std::endl;
			// 	std::cout << "R[k]" << std::endl;
			// 	R[k].print();
			// }

			// step 2.7, make R[k] reduced
			// TODO: reduce matrix by reduce_matrix_clearing()
			// might need to change the update dimension order!
			t0 = std::chrono::steady_clock::now();
			p2c[k] = reduce_matrix(R[k], U[k], args...);
			t1 = std::chrono::steady_clock::now();
			std::cout << "\tStep 2.7 Make final redution on dim "<<k<<" takes";
			std::cout << " takes "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
				<< "ms" << std::endl;
		}
		// std::cout << "\nafter final reduction R[1] = " << std::endl;
		// R[1].print();

		//cohomology the highest dimension need add/remove zero columns(simplices at the highest dimension)
		if(degree == +1){
			// TODO: modify U[max_dim+1] as well, but not necessary because R[max_dim+1] is a zero matrix
			for (size_t i = 0; i < UI.addition_indices[max_dim].size(); i++){
				R[max_dim+1].append_column();
			}
			R[max_dim+1].erase_final_columns(UI.deletion_indices[max_dim].size());
		}

		t0 = std::chrono::steady_clock::now();
		set_indices(); // update homology indices
		t1 = std::chrono::steady_clock::now();
		std::cout << "\tStep 3 Set Barcode takes";
			std::cout << " takes "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
				<< "ms" << std::endl;
	}

	/*
	General Update factorization DU = R by Updating Information
	Note:
	1) when degree = +1 ,i.e., cohomology, we need to
	 modify deletion and addition indices: i <-> n-i-1 , where n is column size
	2) add/remove rows of (co)boundary matrix do not need to modify U,
	while add/remove columns need. However, adding a row requires a mutlplication of it with U matrix.
	3) There are types of permutations used in BATs and I call them:
			for a vector perm and a vector v
			i) new_ind_perm: perm[i] is the new index of v[i]
	 	ii) math_perm: inverse of above
	, in the following implementation, permutations use math_perm
	(I know they are confusing, but changing it is hard...).
	**/
	template <typename Information_type, typename... Args>
	void update_basis_general_clearing(const Information_type &UI, Args ...args){
		using VectT = typename MT::col_type; // column vector type
		using ValT = typename VectT::val_type;

		// step 1. compute the permutation that
		// (a) move the deleted rows/columns to the end
		// (b) permute the intersection part of two filtrations
		size_t max_dim = UI.max_dim;
		// std::cout << "before permutation R[1] is" << std::endl;
		// R[1].print();
		// std::cout << "and U[1] is"  << std::endl;
		// U[1].print();

		// we want to extract permutation and deletion information first
		auto t0 = std::chrono::steady_clock::now();
		for (size_t k = 0; k < max_dim + 1; ++k) {
			size_t m = R[k].nrow();
			size_t n = R[k].ncol();
			// step 1.1 find the permutation used to permute simplices
			// that is going to be deleted to the end
			std::vector<size_t> perm_deletion;
			if (degree == +1){
				perm_deletion = identity_perm(m); // cohomology need to permute row
			}else{
				perm_deletion = identity_perm(n);
			}

			if(!UI.deletion_indices[k].empty()){
				if (degree == +1) {
					std::vector<size_t> deletion_inds = find_reverse_index(UI.deletion_indices[k], m);
					std::sort(deletion_inds.begin(), deletion_inds.end()); // sort it
					perm_deletion = perm_to_the_end(deletion_inds, m);
				}else{
					perm_deletion = perm_to_the_end(UI.deletion_indices[k], n);
				}
			}

			// step 1.2 find the permutation used to permute
			// intersection of simplices and leave the ones in the end unmoved
			std::vector<size_t> perm_intersect;
			if (degree == +1) {
				// auto perm_extended = extension_perm(UI.permutations[k], m);
				// reverse
				std::vector<size_t> perm_rever = find_reverse_index(UI.permutations[k], UI.permutations[k].size());
				std::reverse(perm_rever.begin(), perm_rever.end());
				perm_intersect = extension_perm(perm_rever, m);
			}else{
				perm_intersect = extension_perm(UI.permutations[k], n);
			}


			// step 1.3 Combine the above 2 permutations together, i.e.,
			// compute the permutation that first do deletion permutation,
			// and then do inersection permutation (it is valid since they are disjoint).
			// We will use the perm_deletion to store the result.
			bats::util::apply_perm(perm_deletion, perm_intersect);

			// 1.4 permute U and R
			permute_matrices(k, perm_deletion);
		}
		auto t1 = std::chrono::steady_clock::now();
        std::cout << "\tStep 1. Permuting Basis Matrices";
        std::cout << " takes "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
            << "ms" << std::endl;

		// std::cout << "after permutation, R[1] is" << std::endl;
		// R[1].print();
		// std::cout << "and U[1] is"  << std::endl;
		// U[1].print();
		// std::cout << "check R[1] U[1]^{-1} is " << std::endl;

		// step 2: next we update the factorizations
		for (size_t k = 0; k < max_dim+1; ++k) { // for each dimension
			// step 2.1 make U reduced
			t0 = std::chrono::steady_clock::now();
			p2c[k] = reduce_matrix_standard(U[k], R[k]);
			t1 = std::chrono::steady_clock::now();
			std::cout << "\tStep 2.1 Making U["<<k<<"] Reduced takes";
			std::cout << " takes "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
				<< "ms" << std::endl;
			// if (k == 1){
			// 	std::cout << "step2.1 after reduce U[1], R[1] is" << std::endl;
			// 	R[1].print();
			// 	std::cout << "and U[1] is"  << std::endl;
			// 	U[1].print();
			// }


			// step 2.2 sort columns of U - apply same operations to R
			// to make U upper-triangular
			t0 = std::chrono::steady_clock::now();
			if (!p2c[k].empty()){
				for (size_t j = 0; j < U[k].ncol(); j++) {
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
			}
			t1 = std::chrono::steady_clock::now();
			std::cout << "\tStep 2.2 Making U["<<k<<"] Upper triangle takes";
			std::cout << " takes "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
				<< "ms" << std::endl;

			// if (k == 1){
			// 	std::cout << "step2.2 after make U[1] upper-triangle, R[1] is" << std::endl;
			// 	R[1].print();
			// 	std::cout << "and U[1] is"  << std::endl;
			// 	U[1].print();
			// }

			t0 = std::chrono::steady_clock::now();
			// step 2.3 delete k-simplices in the end(column for homology and row for cohomology)
			if (degree == -1){	// cohomology remove block
				U[k].erase_final_columns(UI.deletion_indices[k].size());
				U[k].erase_final_rows_unsafe(UI.deletion_indices[k].size());
				R[k].erase_final_columns(UI.deletion_indices[k].size());
			}
			if(degree == +1){ // cohomology remove rows
				R[k].erase_rows(UI.deletion_indices[k].size());
			}
			t1 = std::chrono::steady_clock::now();
			std::cout << "\tStep 2.3 Removing "<< k << "-simplices";
			std::cout << " takes "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
				<< "ms" << std::endl;

			// if(k==1){
			// 	std::cout << "step2.3 after removing the last "<< UI.deletion_indices[k].size() << " rows of R[1]" << std::endl;
			// 	R[1].print();
			// }


			// step 2.4 delete (k-1)-simplices in the end(row for homology and column for cohomology)
			if(k!=0 and UI.deletion_indices[k-1].size() > 0){
				t0 = std::chrono::steady_clock::now();
				if (degree == -1){
					R[k].erase_final_rows_unsafe(UI.deletion_indices[k-1].size());

				}else{
					// std::cout << "before removing the last columns of R[k]" << std::endl;
					// R[k].print();
					// std::cout << "and U[k] is"  << std::endl;
					// U[k].print();

					U[k].erase_final_columns(UI.deletion_indices[k-1].size());
					U[k].erase_final_rows_unsafe(UI.deletion_indices[k-1].size());
					R[k].erase_final_columns(UI.deletion_indices[k-1].size());
					// if(k == 1){
						// std::cout << "step2.4 Delete last "<< UI.deletion_indices[k-1].size()<< " columns of R[1] " << std::endl;
						// std::cout << "after removing the last columns of R[k]" << std::endl;
						// R[k].print();
						// std::cout << "and U[k] is"  << std::endl;
						// U[k].print();
					// }
				}
				t1 = std::chrono::steady_clock::now();
				std::cout << "\tStep 2.4 Removing "<< k-1 << "-simplices";
				std::cout << " takes "
					<< std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
					<< "ms" << std::endl;

			}

			// step 2.5 add (k-1)-simplices (zero rows for homology and columns for cohomology)
			// to the specified locations in updating information
			if(k!=0 and  UI.addition_indices[k-1].size() > 0){
				t0 = std::chrono::steady_clock::now();
				std::vector<size_t> add_inds = UI.addition_indices[k-1];
				if (degree == +1){ // cohomology add zero columns in R and 1 block in U
					// std::cout << "\nstep2.5 need to insert new 1 block to U[k]" << std::endl;
					add_inds = find_reverse_index(add_inds, R[k].ncol() + add_inds.size()); // reverse it
					std::sort(add_inds.begin(), add_inds.end()); // sort indices before insert

					// std::cout << "addition indices are" << std::endl;
					// print_1D_vectors(add_inds);

					std::vector<vect_type> zero_cols(add_inds.size());
					R[k].insert_columns(add_inds, zero_cols);

					std::vector<VectT> Ucols(add_inds.size());
					for(size_t i = 0; i < add_inds.size(); i++){
						// U[k] needs a new one block
						U[k].insert_row(add_inds[i]); //insert zero row of U
						Ucols[i] = VectT(size_t(add_inds[i]));
					}
					U[k].insert_columns(add_inds, Ucols);
				}else{ // homology insert zero rows at specified locations
					R[k].insert_rows(add_inds);
				}
				// R[1].print(); U[1].print();
				t1 = std::chrono::steady_clock::now();
				std::cout << "\tStep 2.5 Add "<<k-1<<"-simplices takes";
				std::cout << " takes "
					<< std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
					<< "ms" << std::endl;
			}


			t0 = std::chrono::steady_clock::now();
			// step 2.6 add k-simplices (columns for homology, row for cohomology)
			if (!UI.addition_indices[k].empty()){
				if(k==0){
					// find the # of addition
					size_t count = UI.addition_indices[k].size();
					if (degree == -1){ // homology add columns
						// find the index of the rightmost column of U
						size_t final_ind_of_U = U[k].ncol() - 1;
						// append zero columns
						for (size_t i = 0; i < count; i++){
							R[k].append_column();
							U[k].append_row();
							// create an indicator (one) column vector at a given index
							U[k].append_column(VectT(size_t(final_ind_of_U)));
							final_ind_of_U++;
						}
					}else{ // cohomology append zero rows
						for (size_t i = 0; i < count; i++){
							R[k].append_row();
						}
					}
				}else{ // k >= 1
					std::vector<std::vector<size_t>> bd_info; // boundary information
					std::vector<size_t> add_inds; // addition indices
					if (degree == -1){ // homology add columns
						bd_info = UI.boundary_indices[k];
						add_inds = UI.addition_indices[k]; // ascending order in Y

						U[k].insert_rows(add_inds); //zero rows

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
					}
					else{ // cohomology add rows (!Notice: need to assert in ascending order of indices)
						bd_info = find_reverse_index(UI.boundary_indices[k], R[k].ncol());
						add_inds = find_reverse_index(UI.addition_indices[k], R[k].nrow() + UI.addition_indices[k].size());
						// we need to sort addition indices and corresponding permute their boundary indices
						std::vector<size_t> p = bats::util::sortperm(add_inds);
						bats::util::apply_perm_swap(add_inds, bats::util::inv_perm(p));
						bats::util::apply_perm_swap(bd_info, bats::util::inv_perm(p));

						auto U_transpose = U[k].T();
                        std::vector<VectT> new_rows_in_R;
                        for(size_t i = 0; i < add_inds.size(); i++){
                            // the indices of its boundaries
                            auto simplex_bd_ind = bd_info[i];
                            // boundary indices are reversed but inserting rows require ascending order of indices
                            std::sort(simplex_bd_ind.begin(), simplex_bd_ind.end());

                            // create sparse vector of new row in D
                            std::vector<ValT> bd_values(simplex_bd_ind.size());
                            std::fill(bd_values.begin(),bd_values.end(), ValT(1)); //F2 for now
                            SparseVector<ValT, size_t> new_row_in_D(simplex_bd_ind, bd_values);
							// compute new rows in R
                            auto new_row_in_R = U_transpose.gemv(new_row_in_D);
                            new_rows_in_R.emplace_back(new_row_in_R);
                        }
                        R[k].insert_sparse_rows(add_inds, new_rows_in_R);


						// for(size_t i = 0; i < add_inds.size(); i++){
						// 	// the indices of its boundaries
						// 	auto simplex_bd_ind = bd_info[i];
						// 	// boundary indices are reversed but inserting rows require ascending order of indices
						// 	std::sort(simplex_bd_ind.begin(), simplex_bd_ind.end());
						// 	// index of new added row
						// 	auto ind = add_inds[i];

						// 	// create sparse vector of new row in D
						// 	std::vector<ValT> bd_values(simplex_bd_ind.size());
						// 	std::fill(bd_values.begin(),bd_values.end(), ValT(1)); //F2 for now
						// 	SparseVector<ValT, size_t> new_row_in_D(simplex_bd_ind, bd_values);

						// 	std::vector<ValT> new_row_in_R((R[k].ncol()));
						// 	for (size_t j = 0; j < R[k].ncol(); j++){
						// 		// U[k]'s column vectors
						// 		auto U_cols = U[k].cols();
						// 		// dot product of sparse vectors
						// 		new_row_in_R[j] = new_row_in_D * U_cols[j];
						// 	}

						// 	R[k].insert_row(ind, new_row_in_R);
						// }
					}
				}
			}
			t1 = std::chrono::steady_clock::now();
			std::cout << "\tStep 2.6 Add "<<k<<"-simplices takes";
			std::cout << " takes "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
				<< "ms" << std::endl;
			// if(k == 1){
			// 	std::cout << "\nWhen k = "<< k << std::endl;
			// 	std::cout << "after addition" << std::endl;
			// 	std::cout << "R[k]" << std::endl;
			// 	R[k].print();
			// }

		}
		t0 = std::chrono::steady_clock::now();
		if (degree == -1) {
			// do top dimension normally
			p2c[max_dim-1] = reduce_matrix(R[max_dim-1], U[max_dim-1], bats::standard_reduction_flag());
			for (ssize_t k = max_dim-2; k >= 0; --k) {
				p2c[k] = reduce_matrix_clearing(R[k], U[k], R[k+1], p2c[k+1], bats::standard_reduction_flag());
			}
		} else { // degree == +1
			// do bottom dimension normally
			p2c[0] = reduce_matrix(R[0], U[0], bats::standard_reduction_flag());
			for (ssize_t k = 1; k < max_dim; ++k) {
				p2c[k] = reduce_matrix_clearing(R[k], U[k], R[k-1], p2c[k-1], bats::standard_reduction_flag());
			}
		}
		t1 = std::chrono::steady_clock::now();
		std::cout << "\tStep 3 Make final redution on all dims takes";
		std::cout << " takes "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
			<< "ms" << std::endl;

		// std::cout << "\nafter final reduction R[1] = " << std::endl;
		// R[1].print();

		//cohomology the highest dimension need add zero columns
		if(degree == +1){
			for (size_t i = 0; i < UI.addition_indices[max_dim].size(); i++){
				R[max_dim+1].append_column();
			}
		}

		t0 = std::chrono::steady_clock::now();
		set_indices(); // update homology indices
		t1 = std::chrono::steady_clock::now();
		std::cout << "\tStep 4 Set Barcode takes";
			std::cout << " takes "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
				<< "ms" << std::endl;
	}


};


} // namespace bats
