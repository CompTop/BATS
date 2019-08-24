#pragma once

#include <vector>
#include <iostream>
#include <algorithm>
#include "field.h"


// template over type of values
template <typename TI, typename TV>
class SparseVector
{
private:
	// store index-value pairs
	std::vector<std::pair<TI, TV>> indval;

	// sort in-place
	void sort() {
		// sort in-place
		std::sort(indval.begin(), indval.end());
	}

public:

	SparseVector() {}

	// copy constructor
	SparseVector(const SparseVector &x) : indval(x.indval) {}

	SparseVector(const std::vector<std::pair<TI, TV>> indval) : indval(indval) {}

	SparseVector(const std::vector<TI> &ind, const std::vector<TV> &val) {
		size_t nz = ind.size();
		for (size_t i = 0; i < nz; i++) {
			indval.push_back(std::make_pair(ind[i], val[i]));
		}
	}

	// constructor that takes in integers for val
	SparseVector(const std::vector<TI> &ind, const std::vector<int> &val) {
		size_t nz = ind.size();
		for (size_t i = 0; i < nz; i++) {
			indval.push_back(std::make_pair(ind[i], TV(val[i])));
		}
	}

	// cosntructor that returns indicator in given index
	SparseVector(const TI i) {
		indval.push_back(std::make_pair(i, TV(1)));
	}
	// get index and set index

	// get const iterator through nzs
	auto nzbegin() const {
		return indval.cbegin();
	}

	// get const iterator end of nzs
	auto nzend() const {
		return indval.cend();
	}

	// find nonzero index of last element with index < i
	auto find_last_nz(TI i) {
		auto it = std::lower_bound(
			indval.cbegin(),
			indval.cend(),
			std::make_pair(i, TV(0))
		);
		// if there is no element with index < i, return pointer to end
		return it == indval.cend() ? it : --it;
		//return it--;
	}


	// return last nonzero
	std::pair<TI, TV> last() const {
		return indval.back();
	}

	TI last_nzind() const {
		return indval.back().first;
	}

	// return ith nonzero value
	TV nzval(size_t i) const {
		return indval[i].second;
	}

	// return ith nonzero index
	TI nzind(size_t i) const {
		return indval[i].first;
	}

	// return ith nonzero pair
	std::pair<TI, TV> nzpair(size_t i) const {
		return indval[i];
	}

	// nnz
	inline size_t nnz() const {
		return indval.size();
	}


	// permute in-place
	void permute(const std::vector<size_t>  &perm) {
		for (size_t i = 0; i < nnz(); i++) {
			indval[i].first = perm[indval[i].first];
		}
		sort();
	}

	// return self + ax
	SparseVector caxpy(const TV &a, const SparseVector &x, size_t xoffset=0) const {

		auto xend = x.indval.cend() - xoffset;

		// where to put new vector
		std::vector<std::pair<TI, TV>> tmp;
		auto i1 = indval.cbegin();
		auto i2 = x.indval.cbegin();
		while (i1 < indval.cend() && i2 < xend) {
			if ((*i1).first == (*i2).first) {
				TV val = (a * ((*i2).second)) + (*i1).second;
				// std::cout << "a: " << a << " x: " << ((*i2).second)) << " y: " << (*i1).second << std::endl;
				if (!val.iszero()) {
					tmp.push_back(std::make_pair((*i1).first, val));
				}
				++i1;
				++i2;
			} else if ((*i1).first < (*i2).first) {
				tmp.push_back(*i1);
				++i1;
			} else {
				tmp.push_back(std::make_pair((*i2).first, a * (*i2).second));
				++i2;
			}
		}
		// run through rest of entries and dump in
		// at most one of the loops does anything
		while (i1 < indval.cend()) {
			tmp.push_back(*i1);
			++i1;
		}
		while (i2 < xend) {
			tmp.push_back(std::make_pair((*i2).first, a * (*i2).second));
			++i2;
		}

		return SparseVector(tmp);
	}

	// set
	// y <- ax + y
	void axpy(const TV &a, const SparseVector &x, size_t xoffset=0) {

		// check if there's anything to do
		auto xend = x.indval.cend() - xoffset;

		// if we won't update, return
		if (x.indval.cbegin() >= xend) { return; }

		const SparseVector res = caxpy(a, x, xoffset);
		indval = res.indval;
		return;
	}

	// zeros out pivot in x
	void eliminate_pivot(const SparseVector &x) {
		auto piv = last();
		auto pivx = x.last();
		TV alpha = - piv.second / pivx.second;
		// std::cout << "alpha = " << alpha << std::endl;
		axpy(alpha, x);
	}
	// scal - in place

	// add, subtract, multiply by scalar
	inline SparseVector operator+(const SparseVector &x) {return caxpy(TV(1), x); }
	inline SparseVector operator-(const SparseVector &x) {return caxpy(TV(-1), x); }

	void print() {
		for (size_t i = 0; i < nnz(); i++) {
			std::cout << indval[i].first << " : " << indval[i].second << std::endl;
		}
	}

	void print_row() {
		for (size_t i = 0; i < nnz(); i++) {
			std::cout << "(" << indval[i].first << "," << indval[i].second << ") ";
		}
		std::cout << std::endl;
	}
};

// specialized F2 implementation
template <typename TI, typename IntT>
class SparseVector<TI, ModP<IntT, 2>> {
private:
	using TV = ModP<IntT, 2>;

	// store nonzero indices
	std::vector<TI> ind;

	// sort in-place
	void sort() {
		// sort in-place
		std::sort(ind.begin(), ind.end());
	}

public:

	SparseVector() {}

	// copy constructor
	SparseVector(const SparseVector &x) : ind(x.ind) {}

	// construct from index vector
	SparseVector(const std::vector<TI> ind): ind(ind) {}

	SparseVector(const std::vector<TI> &inds, const std::vector<TV> &val) {
		size_t nz = inds.size();
		ind.reserve(nz);
		for (size_t i = 0; i < nz; i++) {
			if (val[i] != 0) {
				ind.push_back(inds[i]);
			}
		}
	}

	// constructor that takes in integers for val
	SparseVector(const std::vector<TI> &inds, const std::vector<int> &val) {
		size_t nz = inds.size();
		ind.reserve(nz);
		for (size_t i = 0; i < nz; i++) {
			if ((val[i] & 0x1) == 1) {
				ind.push_back(inds[i]);
			}
		}
	}

	// cosntructor that returns indicator in given index
	SparseVector(const TI i) {
		ind.push_back(i);
	}

	// nnz
	inline size_t nnz() const {
		return ind.size();
	}

	// get const iterator through nzs
	auto nzbegin() const {
		return ind.cbegin();
	}

	// get const iterator end of nzs
	auto nzend() const {
		return ind.cend();
	}

	// return ith nonzero value
	constexpr TV nzval(size_t i) const {
		return TV(1);
	}

	// return ith nonzero index
	TI nzind(size_t i) const {
		return ind[i];
	}

	// find nonzero index of last element with index < i
	auto find_last_nz(TI i) {
		auto it = std::lower_bound(
			ind.cbegin(),
			ind.cend(),
			i
		);
		// if there is no element with index < i, return pointer to end
		return it == ind.cend() ? it : --it;
		//return it--;
	}


	// permute in-place
	void permute(const std::vector<size_t>  &perm) {
		for (size_t i = 0; i < nnz(); i++) {
			ind[i] = perm[ind[i]];
		}
		sort();
	}

	// return result of self + ax
	SparseVector caxpy(const SparseVector &x, size_t xoffset=0) const {
		// first check if there is anything to do.
		auto xend = x.ind.cend() - xoffset;

		// where to put new vector
		std::vector<TI> tmp;
		auto i1 = ind.cbegin();
		auto i2 = x.ind.cbegin();
		while (i1 < ind.cend() && i2 < xend) {
			if ((*i1) == (*i2)) {
				++i1;
				++i2;
			} else if ((*i1) < (*i2)) {
				tmp.push_back(*i1);
				++i1;
			} else {
				tmp.push_back(*i2);
				++i2;
			}
		}
		// run through rest of entries and dump in
		// at most one of the loops does anything
		while (i1 < ind.cend()) {
			tmp.push_back(*i1);
			++i1;
		}
		while (i2 < xend) {
			tmp.push_back(*i2);
			++i2;
		}
		return SparseVector(tmp);
	}

	// set
	// y <- ax + y
	void axpy(const SparseVector &x, size_t xoffset=0) {

		// check if there's anything to do
		auto xend = x.ind.cend() - xoffset;
		// if we won't update, return
		if (x.ind.cbegin() >= xend) { return; }
		// if all nzs are in x, dump
		if (nnz() == 0) {
			std::copy(x.ind.cbegin(), xend, std::back_inserter(ind));
			//ind = x.ind;
			return;
		}

		const SparseVector res = caxpy(x, xoffset);
		ind = res.ind;
		return;
	}
	void axpy(const TV &a, const SparseVector &x, size_t xoffset=0) {
		if (a == 0) { return; }
		axpy(x, xoffset);
		return;
	}

	// zeros out pivot in x
	void eliminate_pivot(const SparseVector &x) {
		axpy(x);
	}

	// add, subtract, multiply by scalar
	inline SparseVector operator+(const SparseVector &x) {return caxpy(x); }
	inline SparseVector operator-(const SparseVector &x) {return caxpy(x); }


	void print() {
		for (size_t i = 0; i < nnz(); i++) {
			std::cout << ind[i] << std::endl;
		}
	}

	void print_row() {
		for (size_t i = 0; i < nnz(); i++) {
			std::cout << ind[i] << " ";
		}
		std::cout << std::endl;
	}

};

// // TODO: sparse F2 vector implementation
// template <typename TI>
// class SparseF2Vector
// {
// private:
// 	// non-zero indices
// 	std::vector<TI> ind;
// public:
//
//
// 	// get index and set index
// };
