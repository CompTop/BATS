#pragma once
/*
class to store a differential graded vector space
*/
#include <vector>
#include <util/permutation.hpp>
#include <filtration/update_information.hpp>

namespace bats {

/**
@brief a class for a differential graded vector space.

This encapsulates both chain and cochain complex constructions

degree is the degree of the differential:
-1 for chain complexes (default)
+1 for cochain complexes

differential holds the differential maps

We store maps starting on the edge (-1,0)
(-1) -- (0) -- (1) -- (2) -- ...
This map is 0 in the case of standard chain/cochain complexes
but can be non-zero for augmented chain/cochain complexes

TODO: need to handle +1 boundary in maxdim
*/
template <typename MT>
struct DGVectorSpace {
    int degree; // degree on the differential
    std::vector<MT> differential;

	/**
	Access maps in various dimensions

	if degree is +1 (cohomolgical type), then lowest map index is -1
	* -[-1]-> * -[0]-> * -[1]-> ...

	if degree is -1 (homological type) then lowest map index is 0
	* <-[0]- * <-[1]- * <-[2]- ...
	*/
	MT& operator[](ssize_t k) {
		return (degree == -1) ? differential[k] : differential[k+1];
	}
	const MT& operator[](ssize_t k) const {
		return (degree == -1) ? differential[k] : differential[k+1];
	}

	// default to degree -1
    DGVectorSpace() : degree(-1) {}

	/**
	Construct a DGVector space with maxd dimensions
	*/
	DGVectorSpace(size_t maxd, int deg=-1) : degree(deg), differential(maxd+1) {}

	/**
	Construct a DGVector space explicitly from differentials
	*/
	DGVectorSpace(const std::vector<MT> &diff, int deg=-1) : degree(deg), differential(diff) {}

	// produce a chain complex from a simplicial or cell complex
	template <typename CpxT>
	DGVectorSpace(const CpxT& X, const int deg=-1, const bool augmented=false) : degree(deg) {
		differential.resize(X.maxdim() + 2); // extra maps for augmentation and +1 beyond
		if (degree == -1) {
			// homological type
			// first, handle augmentation
			differential[0] = augmented ? MT(1, X.ncells(0), 1) : MT(0, X.ncells(0));
			for (size_t k = 1; k < X.maxdim() + 1; k++) {
				differential[k] = MT(X.boundary_csc(k));
			}
			// add dummy map one dimension up
			differential[X.maxdim() + 1] = MT(X.ncells(X.maxdim()), 0);
		} else if (degree == +1) {
			// cohomological type
			// first handle augmentation
			differential[0] = augmented ? MT(X.ncells(0), 1, 1) : MT(X.ncells(0), 0);
			for (size_t k = 1; k < X.maxdim() + 1; k++) {
				differential[k] = MT(X.boundary_csc(k)).T();
			}
			// add dummy map one dimension up
			differential[X.maxdim() + 1] = MT(0, X.ncells(X.maxdim()));
		} else {
			throw std::runtime_error("Degree of DGVectorSpace should be +1 or -1");
		}

	}



	inline ssize_t maxdim() const { return differential.size()-2; }
	inline size_t dim(ssize_t k) const {
		return (this->operator[](k)).ncol();
	}

	// permute basis in dimension k
	void permute_basis(size_t k, const std::vector<size_t> &perm) {
		// TODO: this assumes that degree = -1
		if (k == 0) {
			// only worry about boundary[1]
			differential[1].permute_rows(perm);
		} else if (k == maxdim()) {
			differential[k].permute_cols(perm);
			// only worry about boundary[maxdim()]
		} else {
			// need to handle boundary[k] and boundary[k+1]
			differential[k].permute_cols(perm);
			differential[k+1].permute_rows(perm);
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
		// TODO: this assumes that degree = -1
		if (degree == -1) {
			if (k == 0) {
				// only worry about boundary[1]
				differential[1].permute_rows(bats::util::inv_perm(perm));
			} else if (k == maxdim()) {
				differential[k].ipermute_cols(perm);
				// only worry about boundary[maxdim()]
			} else {
				// need to handle boundary[k] and boundary[k+1]
				differential[k].ipermute_cols(perm);
				differential[k+1].permute_rows(bats::util::inv_perm(perm));
			}
		} else { // degree == +1
			// differential[0] : * -> C^0
			// differential[k] : C^k -> C^{k+1}
			if (k == 0) {
				// only worry about boundary[1]
				differential[1].ipermute_cols(perm);
				// std::cout << "d_k = " << k+1 << ": cols " << differential[1].ncol() << ", " << perm.size() << std::endl;
			} else if (k == maxdim()) {
				differential[k].permute_rows(bats::util::inv_perm(perm));
				// std::cout << "d_k = " << k << ": rows " << differential[k].nrow() << ", " << perm.size() << std::endl;
				// only worry about boundary[maxdim()]
			} else {
				// need to handle boundary[k] and boundary[k+1]
				// std::cout << "d_k = " << k << ": rows " << differential[k].nrow() << ", " << perm.size() << std::endl << std::endl;
				// std::cout << "d_k = " << k+1 << ": cols " << differential[k+1].ncol() << ", " << perm.size() << std::endl;
				differential[k+1].ipermute_cols(perm);
				differential[k].permute_rows(bats::util::inv_perm(perm));
			}
		}

	}

	void ipermute_basis(const std::vector<std::vector<size_t>> &perm) {
		for (size_t k = 0; k < perm.size(); k++) {
			ipermute_basis(k, perm[k]);
		}
	}

	/**
	checks that d_{k-degree} d_k = 0
	*/
	bool is_differential(size_t k) {
		return (differential[k+degree] * differential[k]).is_zero();
	}

};


} // namespace bats
