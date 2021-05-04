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

	MT& operator[](size_t k) {
		if (k >= boundary.size()) {
			// throw error
		}
		return boundary[k];
	}
	const MT& operator[](size_t k) const {
		if (k >= boundary.size()) {
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
