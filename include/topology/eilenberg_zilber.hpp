#pragma once
/*
Construction of the Eilenberg-Zilber chain map
C_*(X) \otimes C_*(Y) \to C_*(X \times Y)
Input is 2 simplicial complexes, X, Y
*/

#include <complex/simplicial_complex.hpp>
#include <chain/chain_complex.hpp>
#include <chain/chain_map.hpp>
#include <tuple>

namespace bats {

/*
Eilenberg-Zilberg map
returns:
	C_*(X) \otimes C_*(Y)
	F_* (Eilenberg-Zilber map)
	X \times Y
*/
template <typename CpxT, typename T>
auto EilenbergZilber(
    const CpxT& X,     // first simplicial complex
    const CpxT& Y,     // second simplicial complex
    const size_t maxdim,       		// maximum constructed dimension
	T								// field type
) {

	using VT = SparseVector<T, size_t>;
	using MT = ColumnMatrix<VT>;

	const size_t n = X.ncells(0); // number of 0-cells in X
	CpxT XY( n*n, maxdim); // product simplicial complex

	ChainComplex<MT> CX(X); // chain complex of X
	ChainComplex<MT> CY(Y); // chain complex of Y
	ChainComplex<MT> CXCY(maxdim); // tensor product chain complex

	ChainMap<MT> F(maxdim); // chain map


	std::vector<size_t> s; // placeholder simplex
	std::vector<cell_ind> ci; // list of added simplices

	size_t ind_shift[maxdim+1][maxdim+1]; // keep track of index shifts in direct sum

	// loop over simplex dimension of product
	for (size_t dim = 0; dim <=maxdim; dim++) {
		if (dim == 0) {
			// just do this explicitly
			ind_shift[0][0] = 0;
			CXCY[0] = MT(1, X.ncells(0) * Y.ncells(0));
			F[0] = MT::identity(X.ncells(0) * Y.ncells(0));
			continue;
		}
		CXCY[dim] = MT(CXCY[dim-1].ncol(), 0);
		// dimension of simplices we'll take a product of
		size_t shift = 0; // keep track of index shift
		for (size_t dX = 0; dX <= dim; dX++) {
			size_t dY = dim - dX;
			ind_shift[dX][dY] = shift;
			if (dX > X.maxdim() || dY > Y.maxdim()) {continue;}
			shift += X.ncells(dX) * Y.ncells(dY); // how much we'll shift for next pair of dimensions

			// loop over simplices of dimension dX
			for (size_t iX = 0; iX < X.ncells(dX); iX++) {
				auto sX = X.get_simplex(dX, iX);
				// loop over simplices of dimension dY
				for (size_t iY = 0; iY < Y.ncells(dY); iY++) {
					auto sY = Y.get_simplex(dY, iY);
					ci.clear();
					// create target simplices
					product_paths(
						XY,
						sX.begin(), //X.simplex_begin(dX, iX),
						sX.end(), //X.simplex_end(dX, iX),
						sY.begin(), //Y.simplex_begin(dY, iY),
						sY.end(), //Y.simplex_end(dY, iY),
						s,
						n,
						ci
					);
					// create column in tensor product
					// d(x \otimes y) = (dx \otimes y) + (-1)^deg(x) (x \otimes dy)
					auto dxy = (CX[dX][iX].kron(VT(iY), Y.ncells(dY))).shift_inds(ind_shift[dX-1][dY])\
					 + VT(iX, T(dX & 0x1 ? -1 : 1)).kron(CY[dY][iY], CY[dY].nrow()).shift_inds(ind_shift[dX][dY-1]);
					CXCY[dim].append_column(dxy);

					// create column in map
					// map this basis vector in tensor product
					// to linear combination of simplices in XY
					// TODO: dump inds in ci to column vector
					std::vector<size_t> ind;
					std::vector<T> val;

					for (auto& c: ci) {
						ind.emplace_back(c.ind);
						val.emplace_back(T(1));
					}
					F[dim].append_column(ind, val);
				}
			}
		}
	}
	// need to do this loop at end because some simplices added
	// recursively for higher dimensions
	for (size_t dim = 0; dim <=maxdim; dim++) {
		F[dim].set_nrow(XY.ncells(dim));
	}
    return std::tuple(CXCY, F, XY);
}

// compute shift for kronecker product of
// chain of dimension dx in CX with
// chain of dimension dy in CY
template <typename ChainCpx>
size_t kron_chain_shift(
	const size_t dx,
	const ChainCpx& CX,
	const size_t dy,
	const ChainCpx& CY
) {
	size_t dim = dx + dy;

	// compute index shift
	size_t shift = 0;
	for (size_t _dx = 0; _dx < dx; _dx++) {
		size_t _dy = dim - _dx;
		if (_dx > CX.maxdim() || _dy > CY.maxdim()) {continue;}
		shift += CX.dim(_dx) * CY.dim(_dy); // how much we'll shift for next pair of dimensions
	}
	return shift;
}

/*
construct the Kronecker product of two chains
*/
template <typename VT, typename ChainCpx>
auto kron_chains(
	const VT& cx,  // chain in CX
	const size_t dx, // dimension of chain
	const ChainCpx& CX, // full chain complex
	const VT& cy, // chain in CY
	const size_t dy, // dimension of chain
	const ChainCpx& CY // full chain complex
) {
	size_t shift = kron_chain_shift(dx, CX, dy, CY);
	// shift now holds appropriate index shift for dx, dy
	auto cxcy = cx.kron(cy, CY.dim(dy)).shift_inds(shift);
	return cxcy;
}

/*
construct Kronecker product of
homology representatives of dimension dx in RX with
homology representatives of dimension dy in RY
find representative in RXRY (reducted tensor product)
return this as a map on homology
*/
template <typename ReducedCpx>
auto kron_homology(
	const size_t dx,
	const ReducedCpx& RX,
	const size_t dy,
	const ReducedCpx& RY,
	const ReducedCpx& RXRY
) {
	using chain_type = typename ReducedCpx::chain_type;
	std::vector<chain_type> col;
	size_t dim = dx + dy;

	size_t shift = kron_chain_shift(dx, RX, dy, RY);
	for (size_t i = 0; i < RX.hdim(dx); i++ ) {
		auto cx = RX.get_preferred_representative(i, dx);
		for (size_t j = 0; j < RY.hdim(dy); j++) {
			auto cy = RY.get_preferred_representative(j, dy);
			// form kronecker product of chains
			auto cxcy = cx.kron(cy, RY.dim(dy)).shift_inds(shift);
			// find preferred representative in RXRY
			RXRY.find_preferred_representative(cxcy, dim);
			col.emplace_back(cxcy[RXRY.I[dim]]);
		}
	}

	return ColumnMatrix<chain_type>(RXRY.hdim(dim), RX.hdim(dx)*RY.hdim(dy), col);
}

// obtain kronecker product indices
template <typename CpxT>
std::vector<std::vector<size_t>> kron_index(
	const CpxT& X,
	std::vector<std::vector<size_t>>& Ainds,
	const CpxT& Y,
	std::vector<std::vector<size_t>>& Binds,
	size_t maxdim
) {
	std::vector<std::vector<size_t>> ABinds(maxdim+1);

	// loop over simplex dimension of product
	for (size_t dim = 0; dim <=maxdim; dim++) {
		if (dim == 0) {
			// just do this explicitly
			ABinds[0].reserve(Ainds.size() * Binds.size());
			for (auto& iA : Ainds[0]) {
				for (size_t iB = 0; iB < Y.ncells(0); iB++) {
					ABinds[0].emplace_back(iA * Y.ncells(0) + iB);
				}
			}
			for (auto iA : bats::util::sorted_complement(Ainds[0], X.ncells(0))) {
				for (auto iB : Binds[0]) {
					ABinds[0].emplace_back(iA * Y.ncells(0) + iB);
				}
			}
			std::sort(ABinds[0].begin(), ABinds[0].end());
			continue;
		}

		// dimension of simplices we'll take a product of
		size_t shift = 0; // keep track of index shift
		for (size_t dX = 0; dX <= dim; dX++) {
			size_t dY = dim - dX;
			if (dX > X.maxdim() || dY > Y.maxdim()) {continue;}

			// A x Y
			for (auto& iA : Ainds[dX]) {
				for (size_t iB = 0; iB < Y.ncells(dY); iB++) {
					ABinds[dim].emplace_back(iA * Y.ncells(dY) + iB + shift);
				}
			}

			// X\A x B (don't do cells already inserted)
			for (auto iA : bats::util::sorted_complement(Ainds[dX], X.ncells(dX))) {
				for (auto iB : Binds[dY]) {
					ABinds[dim].emplace_back(iA * Y.ncells(dY) + iB + shift);
				}
			}
			shift += X.ncells(dX) * Y.ncells(dY); // how much we'll shift for next pair of dimensions
		}
		std::sort(ABinds[dim].begin(), ABinds[dim].end());
	}
	return ABinds;
}

/*
Relative Eilenberg-Zilberg map
C_*(X, A) \otimes C_*(Y, B) \to C_*(X \times Y, X \times B \cup A \times Y)
returns:
	C_*(X, A) \otimes C_*(Y, B)
	F_* (Relative Eilenberg-Zilber map)
	C_*(X \times Y, X \times B \cup A \times Y)
*/
template <typename CpxT, typename T>
auto EilenbergZilber(
    const CpxT& X,     // first simplicial complex
	const CpxT& A,
    const CpxT& Y,     // second simplicial complex
	const CpxT& B,
    const size_t maxdim,       		// maximum constructed dimension
	T								// field type
) {

	using VT = SparseVector<T, size_t>;
	using MT = ColumnMatrix<VT>;

	auto [CXCY, F, XY] = EilenbergZilber(X, Y, maxdim, T());

	auto Ainds = X.get_indices(A);
	auto Binds = Y.get_indices(B);
	auto ABinds = kron_index(X, Ainds, Y, Binds, maxdim);
	auto RCXCY = CXCY.relative_complex(ABinds);

	CpxT R = TriangulatedProduct(X, B, maxdim, X.ncells(0));
	R.union_add(TriangulatedProduct(A, Y, maxdim, X.ncells(0)));

	auto Rinds = XY.get_indices(R);
	ChainComplex<MT> RCXY(XY, R); // relative chain complex

	auto RF = F.relative_map(Rinds, ABinds);

    return std::tuple(RCXCY, RF, RCXY);
}

} // namespace bats
