#pragma once

#include <vector>

// class that wraps pre-pairing of complex
// no template type
class MorsePairing
{
private:
	// keeps track of whether a cell is paired
	std::vector<std::vector<bool>> ispaired;
	// cells that are paired up - dimensions 0 -- maxdim-1
	std::vector<std::vector<size_t>> up;
	// cells that are paired down - dimensions 1 -- maxdim
	std::vector<std::vector<size_t>> down;
	// maximum dimension of cells
	size_t _maxdim;

	// vertex side information for cluster
	std::vector<size_t> parent;

 	// find root node of tree containing 0-cell i
	size_t find_parent(size_t i) {
		if (i != parent[i]) {
			parent[i] = find_parent(parent[i]);
		}
		return parent[i];
	}

	// make sure that we can store pairs up to cells in dimension dim
	void reserve(size_t dim) {
		while (ispaired.size() < dim+1) {
			ispaired.emplace_back(std::vector<bool>());
		}
		while (up.size() < dim) {
			up.emplace_back(std::vector<size_t>());
		}
		while (down.size() < dim) {
			down.emplace_back(std::vector<size_t>());
		}
		_maxdim = std::max(_maxdim, dim);
	}

	void reserve(size_t dim, size_t n) {
		// first make sure there are enough dimensions available
		reserve(dim);

		// check ispaired vector
		if ( ispaired[dim].size() < n ) {
			ispaired[dim].resize(n, false);
		}

		// check parent vector
		if (dim == 0) {
			// reserve parent vector
			if (parent.size() < n) {
				size_t st = parent.size();
				parent.reserve(n);
				// set all new 0-cells to have themselves as parent
				for (size_t i = st; i < n; i++) {
					parent[i] = i;
				}
			}
		}

		// check up vector
		if (dim < maxdim()) {
			if (up[dim].size() < n) {
				up[dim].reserve(n);
			}
		}

		// check down vector
		if (dim > 0) {
			if (down[dim-1].size() < n) {
				down[dim-1].reserve(n);
			}
		}
	}

	// set pair (i,j) in dimension dim
	bool _set_pair_unsafe(size_t dim, size_t i, size_t j) {
		up[dim].emplace_back(i);
		down[dim].emplace_back(j);
		ispaired[dim][i] = true;
		ispaired[dim+1][j] = true;
		return true;
	}

	// set pair with reservation checks
	bool _set_pair_supersafe(size_t dim, size_t i, size_t j) {
		if (maxdim() < dim + 1) {
			reserve(dim + 1);
		}
		// check that ispaired has enough entries
		if ( size(dim) < i ) {
			reserve(dim, i);
		}
		if ( size(dim + 1) < j ) {
			reserve(dim + 1, j);
		}
		// check if already paired
		if (ispaired[dim][i] || ispaired[dim+1][j]) {
			// was already paired
			return false;
		}
		return _set_pair_unsafe(dim, i, j);
	}

	// set pair with checks
	bool _set_pair_safe(size_t dim, size_t i, size_t j) {
		// check if already paired
		if (ispaired[dim][i] || ispaired[dim+1][j]) {
			// was already paired
			return false;
		}
		return _set_pair_unsafe(dim, i, j);
	}

	// add an edge (i, j) with index ei
	bool _set_edge(size_t i, size_t j, size_t ei) {
		size_t pi = find_parent(i);
		size_t pj = find_parent(j);
		if (pi == pj) { return false; } // no pair set

		// merge components with i and j
		if (pi < pj) {
			// component with pi was born first
			parent[pj] = pi;
			// homology born at pj dies with ei
			return _set_pair_unsafe(0, pj, ei);
		} else {
			// component with pj was born first
			parent[pi] = pj;
			// homology born at pi dies with ei
			return _set_pair_unsafe(0, pi, ei);
		}
		return true;
	}


public:

	// morse pairing for complex up to maxdim cells
	MorsePairing(size_t maxdim) { reserve(maxdim); }

	// morse pairing on complex with given number of cells in each dim
	MorsePairing(std::vector<size_t> ncells) {
		for (size_t dim = 0; dim < ncells.size(); dim++) {
			reserve(dim, ncells[dim]);
		}
	}

	inline size_t maxdim() const { return _maxdim; }
	inline size_t size(size_t dim) const { return ispaired[dim].size(); }

	// return whether or not cell i in dimension dim is paired
	inline bool is_paired(size_t dim, size_t i) const { return ispaired[dim][i]; }

	// set pair between cell i in dimension dim, and cell j in dimension dim+1
	inline bool set_pair(size_t dim, size_t i, size_t j) { return _set_pair_safe(dim, i, j); }

	// set pair between vertices and edge
	inline bool set_pair_edge(size_t i, size_t j, size_t ei) {return _set_edge(i, j, ei); }

	inline const std::vector<size_t>& up_paired(size_t dim) const {return up[dim]; }
	inline const std::vector<size_t>& down_paired(size_t dim) const {return down[dim-1]; }
	std::vector<size_t> unpaired(size_t dim) const {
		std::vector<size_t> ind;
		for (size_t i = 0; i < ispaired[dim].size(); i++ ){
			if (!ispaired[dim][i]) {
				ind.emplace_back(i);
			}
		}
		return ind;
	}


};
