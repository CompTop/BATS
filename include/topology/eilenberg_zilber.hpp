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
template <typename T>
auto EilenbergZilber(
    const SimplicialComplex& X,     // first simplicial complex
    const SimplicialComplex& Y,     // second simplicial complex
    const size_t maxdim,       		// maximum constructed dimension
	T								// field type
) {

	using VT = SparseVector<T, size_t>;
	using MT = ColumnMatrix<VT>;

	const size_t n = X.ncells(0); // number of 0-cells in X
	SimplicialComplex XY; // product simplicial complex

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
				// loop over simplices of dimension dY
				for (size_t iY = 0; iY < Y.ncells(dY); iY++) {
					ci.clear();
					// create target simplices
					product_paths(
						XY,
						X.simplex_begin(dX, iX),
						(X.simplex_end(dX, iX) - 1),
						Y.simplex_begin(dY, iY),
						(Y.simplex_end(dY, iY) - 1),
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

} // namespace bats
