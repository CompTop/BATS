#pragma once
/*
class to store a chain complex
*/
#include <vector>

// template over matrix type
template <typename MT>
struct ChainComplex {
	std::vector<size_t> dim;
	std::vector<MT> boundary;


	ChainComplex(std::vector<size_t> &dim, std::vector<MT> &boundary) : dim(dim), boundary(boundary) {}

	// produce a chain complex from a simplicial or cell complex
	template <typename CpxT>
	ChainComplex(CpxT& X) {
		dim.resize(X.maxdim() + 1);
		boundary.resize(X.maxdim() + 1);
		for (int k = 0; k < X.maxdim() + 1; k++) {
			dim[k] = X.ncells(k);
			if (k == 0) {
				boundary[k] = MT(1, dim[k]);
			} else {
				boundary[k] = MT(X.boundary_csc(k));
			}
		}
	}

	inline size_t maxdim() const { return dim.size() - 1; }
	//inline size_t dim(size_t k) const { return dim[k]; }
	//inline MT& boundary(size_t k) { return boundary[k]; }

	MT& operator[](size_t k) {
		if (k >= dim.size()) {
			// throw error
		}
		return boundary[k];
	}

};
