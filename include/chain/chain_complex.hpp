#pragma once
/*
class to store a chain complex
*/
#include <vector>
#include <util/permutation.hpp>
#include <filtration/update_information.hpp>

namespace bats {


// template over matrix type
template <typename MT>
struct ChainComplex {
	//std::vector<size_t> dim;
	std::vector<MT> boundary;

	ChainComplex() {}

	// chain complex on maxd
	ChainComplex(size_t maxd) : boundary(maxd+1) {}

	ChainComplex(const std::vector<MT> &boundary) : boundary(boundary) {}

	// produce a chain complex from a simplicial or cell complex
	template <typename CpxT>
	ChainComplex(const CpxT& X) {
		//dim.resize(X.maxdim() + 1);
		boundary.resize(X.maxdim() + 1);
		for (size_t k = 0; k < X.maxdim() + 1; k++) {
			//dim[k] = X.ncells(k);
			if (k == 0) {
				boundary[k] = MT(1, X.ncells(k));
			} else {
				boundary[k] = MT(X.boundary_csc(k));
			}
		}
	}

	// produce a relative chain complex C(X, A)
	template <typename CpxT>
	ChainComplex(const CpxT& X, const CpxT& A) {
		//dim.resize(X.maxdim() + 1);
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
			//dim[k] = boundary[k].ncol();
		}
	}

	inline size_t maxdim() const { return boundary.size() - 1; }
	inline size_t dim(size_t k) const { return boundary[k].ncol(); }
	size_t dim() const {
		size_t ret = 0;
		for (size_t k = 0; k < boundary.size(); ++k) { ret += dim(k); }
		return ret;
	}

	// check that ChainComplex is a valid complex
	// checks that composition of boundary maps is 0
	bool is_valid_complex() const {
		for (size_t k = 0; k < maxdim(); k++) {
			if ( !((boundary[k]*boundary[k+1]).is_zero()) ) {return false;}
		}
		return true;
	}

	// construct subcomplex defind by inds[k] in dimension k
	// assume this is a valid subcomplex
	ChainComplex subcomplex(std::vector<std::vector<size_t>>& inds) const {
		std::vector<MT> newboundary(boundary.size());
		for (size_t k = 0; k < maxdim() + 1; k++) {
			std::sort(inds[k].begin(), inds[k].end());
			if (k == 0) {
				newboundary[k] = MT(1, inds[k].size());
			} else {
				newboundary[k] = dim[k].submatrix(inds[k-1], inds[k]);
			}
		}
		return ChainComplex(newboundary);
	}

	// construct relative complex
	// by quotienting out inds in each dimension
	// assume inds defines a subcomplex
	ChainComplex relative_complex(std::vector<std::vector<size_t>>& inds) const {
		std::vector<MT> newboundary(boundary.size());
		std::vector<std::vector<size_t>> cinds;
		for (size_t k = 0; k < maxdim() + 1; k++) {
			std::sort(inds[k].begin(), inds[k].end());
			cinds.emplace_back(bats::util::sorted_complement(inds[k], dim(k)));
			if (k == 0) {
				newboundary[k] = MT(1, dim(k) - inds[k].size());
			} else {
				newboundary[k] = boundary[k].submatrix(cinds[k-1], cinds[k]);
			}
		}
		return ChainComplex(newboundary);
	}


	/**
	reference to k-dimensional boundary

	@param k	dimension
	*/
	MT& operator[](size_t k) {
		return boundary[k];
	}
	const MT& operator[](size_t k) const {
		return boundary[k];
	}

	// permute basis in dimension k
	void permute_basis(size_t k, const std::vector<size_t> &perm) {
		if (k == 0) {
			// only worry about boundary[1]
			boundary[1].permute_rows(perm);
		} else if (k == maxdim()) {
			boundary[k].permute_cols(perm);
			// only worry about boundary[maxdim()]
		} else {
			// need to handle boundary[k] and boundary[k+1]
			boundary[k].permute_cols(perm);
			boundary[k+1].permute_rows(perm);
			// boundary[k+1].ipermute_rows(bats::util::inv_perm(perm));
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
			boundary[1].permute_rows(bats::util::inv_perm(perm));
		} else if (k == maxdim()) {
			boundary[k].ipermute_cols(perm);
			// only worry about boundary[maxdim()]
		} else {
			// need to handle boundary[k] and boundary[k+1]
			boundary[k].ipermute_cols(perm);
			boundary[k+1].permute_rows(bats::util::inv_perm(perm));
		}
	}

	void ipermute_basis(const std::vector<std::vector<size_t>> &perm) {
		for (size_t k = 0; k < perm.size(); k++) {
			ipermute_basis(k, perm[k]);
		}
	}

template <typename Information_type>
	void update_basis_general(size_t k, const Information_type &UI) {
		using VectT = typename MT::col_type; // column vector type of boundary matrix
		using ValT = typename VectT::val_type;


		std::vector<size_t> PI_row = UI.permutations[k]; // row permutation information
		std::vector<size_t> PI_col = UI.permutations[k+1]; // column permutation information
		std::vector<size_t> del_ind_row = UI.deletion_indices[k]; // row deletion indices
		std::vector<size_t> del_ind_col = UI.deletion_indices[k+1]; // column addition indices
		std::vector<size_t> add_ind_row = UI.addition_indices[k]; // row addition indices
		std::vector<size_t> add_ind_col = UI.addition_indices[k+1]; // column addition indices
		
		// if k == 0, boundary[0] are all zeros
		if (k == 0) { 
			// step 1. delete columns
			if (!del_ind_row.empty()){
				// std::cout << "do deletion" << std::endl;
				for (size_t i = 0; i < del_ind_row.size(); i++){
					boundary[0].erase_column();
				}
			}
			// step 2. permute columns(no need here since are all zeros)

			// step 3. addition
			if (!add_ind_row.empty()){
				// std::cout << "do addition" << std::endl;
				for (size_t i = 0; i < add_ind_row.size(); i++){
					boundary[0].append_column();
				}
			}
		}
		// std::cout << "\nafter modification of boundary[0]" << std::endl;
		// boundary[0].print();

		// Start for boundary[k+1]:
		// Step 1. delete columns(we assume the deletion indices in A are in ascending order)
		if (!del_ind_col.empty()){
			// std::cout << "do deletion of column: ";
			// print_1D_vectors(del_ind_col);

			// find the permutation that will permute 
			// the column that is going to remove to the end
			std::vector<size_t> perm = perm_to_the_end(del_ind_col, boundary[k+1].ncol());
			boundary[k+1].permute_cols(perm); // permute 
			for(size_t i = 0; i < del_ind_col.size(); i++){
				// size_t current_index = ind - i;
				boundary[k+1].erase_column(); // erase the last column
				// i++;
			}
			// std::cout << "\n";
		}
		// std::cout << "\nafter deletion of columns" << std::endl;
		// boundary[k+1].print();

		// 2. delete rows(we assume the deletion indices in A are in ascending order)
		if (!del_ind_row.empty()){
			// std::cout << "do deletion of row: ";
			// print_1D_vectors(del_ind_row);

			std::vector<size_t> perm = perm_to_the_end(del_ind_row, boundary[k+1].nrow());
			boundary[k+1].permute_rows(perm); // permute 
			for(size_t i = 0; i < del_ind_row.size(); i++){
				boundary[k+1].erase_row(); // erase the last row
			}
		}
		// std::cout << "\nafter deletion of rows" << std::endl;
		// boundary[k+1].print();

		// 3. permute rows of boundary[k+1]
		if(!PI_row.empty() && PI_row != identity_perm(PI_row.size())){
			// std::cout << "permute rows with permutation" << std::endl;
			// print_1D_vectors(PI_row);
			boundary[k+1].ipermute_rows(bats::util::inv_perm(PI_row));
		}

		// std::cout << "\nafter permutation of rows" << std::endl;
		// boundary[k+1].print();

		// 4. permute columns of boundary[k+1]
		if(!PI_col.empty() && PI_col != identity_perm(PI_col.size())){
			// std::cout << "permute columns with permutation" << std::endl;
			// print_1D_vectors(PI_col);

			boundary[k+1].ipermute_cols(PI_col); // permutation defintion confict with bats
		}
		// std::cout << "\nafter permutation of columns" << std::endl;
		// boundary[k+1].print();

		// 5. add rows
		if(!add_ind_row.empty()){
			// std::cout << "do addition of row ";
			// print_1D_vectors(add_ind_row);
			for(auto ind: add_ind_row){
				boundary[k+1].insert_row(ind);
			}
		}
		// std::cout << "\nafter addition of row "<< std::endl;
		// boundary[k+1].print();

		// 6. add columns
		if(!add_ind_col.empty()){
			// std::cout << "do addition of column\n";

			auto bd_info = UI.boundary_indices[k+1];
			for(size_t i = 0; i < add_ind_col.size(); i++){
				// we will add a column to the boundary matrix
				auto simplex_bd_inds = bd_info[i]; // the indices of its boundaries
                std::sort(simplex_bd_inds.begin(), simplex_bd_inds.end()); // sort indices for now

				// std::cout << "at index "<< add_ind_col[i] << std::endl;
				// std::cout << "and its boundary index is ";
				// print_1D_vectors(simplex_bd_inds);

				// create a vector of ones
				std::vector<ValT> vect_one(simplex_bd_inds.size(), 1);
				auto vect = VectT(simplex_bd_inds, vect_one); // the added vector for the column
				boundary[k+1].insert_column(add_ind_col[i], vect);
			}
			// std::cout << "\n";
		}
		// std::cout << "after addition of column "<< std::endl;
		// boundary[k+1].print();
	}

	template <typename Information_type>
	void update_basis_general(const Information_type &UI_fast) {
		for (size_t k = 0; k < UI_fast.max_dim; k++) {
			update_basis_general(k, UI_fast);	
		}
	}

	// tensor product of chain complexes A\otimes B
	// construct tensor product up to dimension dmax
	friend ChainComplex tensor_product(const ChainComplex &A, const ChainComplex &B, size_t dmax) {
		ChainComplex C(dmax);
		using VT = typename MT::col_type;
		using T = typename VT::val_type;

		size_t ind_shift[dmax+1][dmax+1]; // keep track of index shifts in direct sum
		for (size_t n = 0; n <= dmax; n++) {
			if (n== 0) {
				// just do this explicitly
				ind_shift[0][0] = 0;
				C[0] = MT(1, A.dim(0) * B.dim(0));
				continue;
			}
			C[n] = MT(C[n-1].ncol(), 0);

			size_t shift = 0; // keep track of index shift
			for (size_t Adim = 0; Adim <= n; Adim++) {
				size_t Bdim = n - Adim;
				ind_shift[Adim][Bdim] = shift;
				if (Adim > A.maxdim() || Bdim > B.maxdim()) {continue;}
				shift += A.dim(Adim) * B.dim(Bdim); // how much we'll shift for next pair of dimensions
				// add to tensor product
				// determine shift from dimensions
				for (size_t iA = 0; iA < A.dim(Adim); iA++) {
					// loop over simplices of dimension dY
					for (size_t iB = 0; iB < B.dim(Bdim); iB++) {
						auto dxy = (A[Adim][iA].kron(VT(iB), B.dim(Bdim))).shift_inds(ind_shift[Adim-1][Bdim])\
						 + VT(iA, T(Adim & 0x1 ? -1 : 1)).kron(B[Bdim][iB], B[Bdim].nrow()).shift_inds(ind_shift[Adim][Bdim-1]);
						C[n].append_column(dxy);
					}
				}
			}
		}
		return C;
	}

	/**
	A preprocessing step for computing homology using the reduction algorithm.

	@warning using this function will invalidate any basis used by
	a homology algorithm since no basis vector will be obtained for
	the cleared columns
	*/
	void clear_compress_apparent_pairs() {

		std::vector<bool> comp_inds;
		std::vector<bool> clear_inds;
		for (size_t dim = 1; dim < maxdim() + 1; ++dim) {
			// prepare clearing and compression indices
			clear_inds.resize(boundary[dim].nrow());
			std::fill(clear_inds.begin(), clear_inds.end(), false);
			comp_inds.resize(boundary[dim].ncol());
			std::fill(comp_inds.begin(), comp_inds.end(), false);
			// look at pivots to determine clearing/compression inds
			for (size_t j = 0; j < boundary[dim].ncol(); ++j) {
				auto it = boundary[dim][j].nzend();
				// check that nnz > 0
				if (it != boundary[dim][j].nzbegin()) {
					--it; // point to last nonzero
					size_t i = it->ind; // pivot index
					if (!clear_inds[i]) {
						// this is the first time we have seen this pivot
						// mark j as a compression ind
						comp_inds[j] = true;
						// mark i as a clearing ind
						clear_inds[i] = true;
					}
				}
			}
			// we can now clear indices one dimension lower
			if (dim > 1) {
				boundary[dim-1].clear_cols(clear_inds);
			}
			// we can now compress indices one dimension higher
			if (dim < maxdim()) {
				boundary[dim+1].clear_rows(comp_inds);
			}
		}




	}

	inline friend ChainComplex tensor_product(
		const ChainComplex &A,
		const ChainComplex& B
	) {return tensor_product(A, B, A.maxdim() + B.maxdim()); }

};

// defualt return
template <typename T, typename CpxT>
inline auto __ChainComplex(const CpxT &X, T) {
	using VT = SparseVector<T, size_t>;
	using MT = ColumnMatrix<VT>;

	return ChainComplex<MT>(X);
}

template <typename T, typename CpxT>
inline auto Chain(const CpxT& X, T) {return __ChainComplex(X, T());}

template <typename T, typename CpxT>
inline auto __ChainComplex(const CpxT &X, const CpxT &A, T) {
	using VT = SparseVector<T, size_t>;
	using MT = ColumnMatrix<VT>;

	return ChainComplex<MT>(X, A);
}

template <typename T, typename CpxT>
inline auto Chain(const CpxT& X, const CpxT &A, T) {return __ChainComplex(X, A, T());}

} // namespace bats
