#pragma once

#include <vector>
#include <complex/abstract_complex.hpp>
#include <linalg/csc_matrix.hpp>

// class that wraps pre-pairing of complex
// template over complex type
template <class CpxT>
class MorsePairing
{
private:
	// pointer to complex
	CpxT* cpx;

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
	size_t find_parent(const size_t i) {
		if (i != parent[i]) {
			parent[i] = find_parent(parent[i]);
		}
		return parent[i];
	}

	// make sure that we can store pairs up to cells in dimension dim
	void reserve(const size_t dim) {
		while (ispaired.size() < dim+1) {
			ispaired.emplace_back(std::vector<bool>());
		}
		while (up.size() < dim+1) {
			up.emplace_back(std::vector<size_t>());
		}
		while (down.size() < dim) {
			down.emplace_back(std::vector<size_t>());
		}
		_maxdim = std::max(_maxdim, dim);
	}

	void reserve(const size_t dim, const size_t n) {
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
				parent.resize(n);
				// set all new 0-cells to have themselves as parent
				for (size_t i = st; i < n; i++) {
					parent[i] = i;
				}
			}
		}

		// don't reserve up/down - just let them grow dynamically
		// // check up vector
		// if (dim < maxdim()) {
		// 	if (up[dim].size() < n) {
		// 		up[dim].reserve(n);
		// 	}
		// }
		//
		// // check down vector
		// if (dim > 0) {
		// 	if (down[dim-1].size() < n) {
		// 		down[dim-1].reserve(n);
		// 	}
		// }
	}

	// set pair (i,j) in dimension dim
	bool _set_pair_unsafe(const size_t dim, const size_t i, const size_t j) {
		up[dim].emplace_back(i);
		down[dim].emplace_back(j);
		ispaired[dim][i] = true;
		ispaired[dim+1][j] = true;
		return true;
	}

	// set pair with reservation checks
	bool _set_pair_supersafe(const size_t dim, const size_t i, const size_t j) {
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
	inline bool _set_pair_safe(const size_t dim, const size_t i, const size_t j) {
		return (ispaired[dim][i] || ispaired[dim+1][j]) ? false : _set_pair_unsafe(dim, i, j);
	}

	// add an edge (i, j) with index ei
	bool _set_edge(const size_t i, const size_t j, const size_t ei) {
		size_t pi = find_parent(i);
		size_t pj = find_parent(j);
		if (pi == pj) { return false; } // no pair set

		// merge components with i and j
		// TODO: should be able to call a rank vector
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

	// add an edge (i, j) with index ei
	// use rank vector to compare
	template <typename TR>
	bool _set_edge(const size_t i, const size_t j, const size_t ei, const std::vector<TR> &rank) {
		size_t pi = find_parent(i);
		size_t pj = find_parent(j);
		if (pi == pj) { return false; } // no pair set

		// merge components with i and j
		if ( (rank[pi] < rank[pj]) || ((rank[pi] == rank[pj]) && (pi < pj)) ) {
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

	// add to underlying complex
	// assume reserved
	template <class ...Ts>
	inline cell_ind _add_unsafe(Ts (&...args)) { return cpx->_add_unsafe(args...); }

	template <class ...Ts>
	inline cell_ind _add_unsafe_reserve(Ts (&...args)) {
		cell_ind ret = cpx->_add_unsafe(args...);
		reserve(ret.dim, ret.ind+1);
		return ret;
	}



	// add to underlying complex and if possible pair with face
	template <class ...Ts>
	cell_ind _add_pair_unsafe(Ts (&...args)) {
		// first add to complex
		cell_ind ret = _add_unsafe(args...);
		// now iterate over boundary
		if (ret.dim > 1) {
			for (auto k = cpx->faces_begin(ret.dim, ret.ind); k < cpx->faces_end(ret.dim, ret.ind); k++) {
				if (!ispaired[ret.dim-1][*k]) {
					// set the pair then exit
					_set_pair_unsafe(ret.dim - 1, *k, ret.ind);
					return ret;
				}
			}
		} else if (ret.dim == 1) {
			auto k = cpx->faces_begin(ret.dim, ret.ind);
			_set_edge(*k, *(k+1), ret.ind);
		}
		return ret;
	}


public:

	MorsePairing() : _maxdim(0) {};

	MorsePairing(CpxT &C) : cpx(&C), _maxdim(0) {
		for (size_t dim = 0; dim < C.maxdim() + 1; dim++){
			reserve(dim, C._reserved[dim]);
			// reserve(dim, C.ncells(dim));
		}
	};

	// morse pairing for complex up to maxdim cells
	MorsePairing(size_t maxdim) : _maxdim(0) { reserve(maxdim); }

	// morse pairing on complex with given number of cells in each dim
	MorsePairing(std::vector<size_t> ncells) : _maxdim(0) {
		for (size_t dim = 0; dim < ncells.size(); dim++) {
			reserve(dim, ncells[dim]);
		}
	}

	// clear all pairings, but keep size and memory allocated
	void clear() {
		for (size_t i = 0; i < ispaired.size(); i++) {
			ispaired[i] = std::vector<bool>(ispaired[i].size(), false);
		}
		for (size_t i = 0; i < up.size(); i++) {
			up[i].clear();
		}
		for (size_t i = 0; i < down.size(); i++) {
			down[i].clear();
		}
	}

	inline size_t maxdim() const { return _maxdim; }
	inline size_t size(const size_t dim) const { return ispaired[dim].size(); }
	inline size_t ncells(const size_t dim) const { return ispaired[dim].size(); }

	// return whether or not cell i in dimension dim is paired
	inline bool is_paired(size_t dim, size_t i) const { return ispaired[dim][i]; }

	// set pair between cell i in dimension dim, and cell j in dimension dim+1
	inline bool set_pair(size_t dim, size_t i, size_t j) { return _set_pair_safe(dim, i, j); }

	// set pair between vertices and edge
	inline bool set_pair_edge(const size_t i, const size_t j, const size_t ei) {return _set_edge(i, j, ei); }
	template <typename TR>
	inline bool set_pair_edge(const size_t i, const size_t j, const size_t ei, const std::vector<TR> &rank) {
		return _set_edge(i, j, ei, rank);
	}

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

	// add to underlying complex
	template <class ...Ts>
	inline cell_ind add(Ts (&...args)) {
		cell_ind ret = cpx->add(args...);
		reserve(ret.dim, ret.ind+1);
		return ret;
	}

	// add to underlying complex and if possible pair with face
	template <class ...Ts>
	cell_ind add_pair(Ts (&...args)) {
		// first add to complex
		cell_ind ret = add(args...);
		// now iterate over boundary
		if (ret.dim > 1) {
			for (auto k = cpx->faces_begin(ret.dim, ret.ind); k < cpx->faces_end(ret.dim, ret.ind); k++) {
				if (!ispaired[ret.dim-1][*k]) {
					// set the pair then exit
					_set_pair_unsafe(ret.dim - 1, *k, ret.ind);
					return ret;
				}
			}
		} else if (ret.dim == 1) {
			auto k = cpx->faces_begin(ret.dim, ret.ind);
			_set_edge(*k, *(k+1), ret.ind);
		}
		return ret;
	}

	template <class ...Ts>
	inline auto faces_begin(Ts (&...args)) {
		return cpx->faces_begin(args...);
	}
	template <class ...Ts>
	inline auto faces_end(Ts (&...args)) {
		return cpx->faces_end(args...);
	}

	inline CSCMatrix<int, size_t> boundary_csc(const size_t dim) const { return cpx->boundary_csc(dim); }

	template <typename T1, typename T2>
	friend class Filtration;

};
