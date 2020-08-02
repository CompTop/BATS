#pragma once

/*
An implementation of a sparse vector
maintains indices in sorted order
*/


#include <cstddef>
#include <vector>
#include <iostream>
#include <algorithm>
#include <string>
#include <sstream>
#include "field.h"
#include "abstract_vector.h"


// template over type of values
template <typename TV, typename TI=size_t>
class SparseVector
{
private:

	using key_type = nzpair<TI, TV>; // std::pair<TI, TV>

	std::vector<key_type> indval;

	// sort in-place
	void sort() {
		// sort in-place
		std::sort(indval.begin(), indval.end());
	}

public:
	using val_type = TV;

	SparseVector() {}

	SparseVector(const std::vector<key_type> &indval) : indval(indval) {}

	SparseVector(const std::vector<TI> &ind, const std::vector<TV> &val) {
		size_t nz = ind.size();
		indval.reserve(nz);
		for (size_t i = 0; i < nz; i++) {
			indval.emplace_back(key_type(ind[i], val[i]));
		}
	}

	// constructor from other sparse vector
	template<typename TI2>
	SparseVector(const SparseVector<int, TI2> &other) {
		indval.reserve(other.nnz());
		for (auto it = other.nzbegin(); it < other.nzend(); ++it) {
			indval.emplace_back(key_type((*it).ind, (*it).val));
		}
	}

	// constructor that loops over index and value iterators
	// can be iterators over different type
	template <typename IT1, typename IT2>
	SparseVector(IT1 indit, IT2 valit, size_t n) {
		indval.reserve(n);
		for (size_t i = 0; i < n; i++) {
			indval.emplace_back(key_type(TI(*indit++), TV(*valit++)));
		}
	}


	// cosntructor that returns indicator in given index
	SparseVector(const TI i) {
		indval.push_back(key_type(i, TV(1)));
	}

	// construct from line string
	SparseVector(std::string &line) {
		std::string token;
		std::istringstream iss(line);
		while (getline(iss, token, ',')) {
			// std::cout << token << ',' << std::endl;
			if (token.size() > 0) {
				indval.emplace_back(key_type(token));
			}
		}
	}

	// SparseVector& operator=(const SparseVector &other) {
	// 	indval = other.indval;
	// 	return *this;
	// }

	template<typename T>
	bool operator==(const T &other) const {
		auto it1 = nzbegin();
		auto it2 = other.nzbegin();
		while (it1 != nzend() && it2 != other.nzend()) {
			if (*it1++ != *it2++) {return false;}
		}
		if (it1 != nzend() || it2 != other.nzend()) { return false;}
		return true;
	}

	// extract indices by iterator
	template<typename T>
	SparseVector operator[](const T &indset) const {
		auto it = indset.cbegin();
		auto end = indset.cend();
		std::vector<key_type> indval2;
		indval2.reserve(indval.size());
		auto selfit = indval.cbegin();
		auto selfend = indval.cend();
		size_t ct = 0;
		while (it != end && selfit != selfend) {
			if (*it < selfit->ind) {
				// zero entry
				// indval2.emplace_back(key_type(ct, TV(0)));
				++it;
				++ct;
			} else if (*it == selfit->ind) {
				indval2.emplace_back(key_type(ct, selfit->val));
				++selfit;
				++it;
				++ct;
			} else {
				++selfit;
			}
		}
		// only enter this loop if index set is longer than nonzeros
		// do not want to put zero values in
		// while (it != end) {
		// 	indval2.emplace_back(key_type(ct, TV(0)));
		// 	++it;
		// 	++ct;
		// }
		return SparseVector(indval2);
	}

	// get const iterator through nzs
	inline auto nzbegin() const { return indval.cbegin(); }
	inline auto nzend() const { return indval.cend(); }
	// get nonconst iterator through nzs
	inline auto nzbegin() { return indval.begin(); }
	inline auto nzend(){ return indval.end(); }

	inline void clear() {indval.clear();}

	void clear_zeros() {
		std::vector<key_type> indval2;
		for (auto iv : indval) {
			if (iv.val != TV(0)) {
				indval2.emplace_back(iv);
			}
		}
		indval = indval2;
	}

	// nnz
	inline size_t nnz() const {return indval.size(); }

	std::vector<size_t> nzinds() const {
		std::vector<size_t> ind;
		ind.reserve(nnz());
		for (auto it = nzbegin(); it != nzend(); it++) {
			ind.emplace_back((*it).ind);
		}
		return ind;
	}

	// returns iterator pointing to first element that is not less than i
	inline auto lower_bound(const TI &i) {
		return std::lower_bound(nzbegin(), nzend(), key_type(i, TV(0)));
	}
	inline auto lower_bound(const TI &i) const {
		return std::lower_bound(nzbegin(), nzend(), key_type(i, TV(0)));
	}
	// returns iterator pointing to first element that is greater than i
	inline auto upper_bound(const TI &i) {
		return std::upper_bound(nzbegin(), nzend(), key_type(i, TV(0)));
	}
	inline auto upper_bound(const TI &i) const {
		return std::upper_bound(nzbegin(), nzend(), key_type(i, TV(0)));
	}

	// replace the item at the iterator location with k
	template <typename itT>
	auto replace(itT &it, const key_type &k) {
		*it = k;
		return it;
	}
	template <typename itT>
	inline auto replace(itT &it, const TI ind, const TV val) { return replace(it, key_type(ind, val));}
	template <typename itT>
	inline auto replace(itT &it, const TV val) { return replace(it, key_type((*it).ind, val));}

	TV getval(const size_t i) const {
		auto it = lower_bound(i);
		if (it == nzend() || (*it).ind != i) {
			return TV(0);
		}
		return (*it).val;
	}

	TV operator[](size_t i) const {
		return getval(i);
	}

	inline auto emplace_back(const key_type &k) {
		return indval.emplace_back(k);
	}
	inline auto emplace_back(TI ind, TV val) {
		return emplace_back(key_type(ind, val));
	}

	// find nonzero index of last element with index < i
	auto find_last_nz(TI i) {
		auto it = std::lower_bound(
			indval.cbegin(),
			indval.cend(),
			key_type(i, TV(0))
		);
		// if there is no element with index < i, return pointer to end
		return it == indval.cend() ? it : --it;
		//return it--;
	}


	// return last nonzero
	inline const key_type& lastnz() const {
		return indval.back();
	}

	inline const key_type& firstnz() const {
		return indval[0];
	}



	// TI last_nzind() const {
	// 	return indval.back().ind;
	// }
	//
	// // return ith nonzero value
	// TV nzval(size_t i) const {
	// 	return indval[i].val;
	// }
	//
	// // return ith nonzero index
	// TI nzind(size_t i) const {
	// 	return indval[i].ind;
	// }

	// // return ith nonzero pair
	// const key_type& nzpair(size_t i) const {
	// 	return indval[i];
	// }



	// permute in-place
	void permute(const std::vector<size_t>  &perm) {
		for (size_t i = 0; i < nnz(); i++) {
			indval[i].ind = perm[indval[i].ind];
		}
		sort();
	}

	// apply J-matrix transformation
	// index i -> (m-1) - i
	void J(const size_t m) {
		for (size_t i = 0; i < nnz(); i++) {
			indval[i].ind = (m - 1) - indval[i].ind;
		}
		sort();
	}

	// set
	// y <- ax + y
	template <class SVT>
	void axpy(
		const TV &a,
		const SVT &x
	) {

		if (a == TV(0)) {return;}

		// set i2 to find first ind >= firstind
		auto i2 = x.nzbegin();
		// where to put new vector
		if (i2 == x.nzend()) { return; } // nothing to do
		// something to do...
		auto i1 = indval.cbegin();
		std::vector<key_type> tmp;
		while (i1 != indval.cend() && i2 != x.nzend()) {
			if ( i1->ind == i2->ind) {
				TV val = (a * ((*i2).val)) + (*i1).val;
				// std::cout << "a: " << a << " x: " << ((*i2).val)) << " y: " << (*i1).val << std::endl;
				if (val != TV(0)) {
					tmp.push_back(key_type(i1->ind, val));
				}
				++i1;
				++i2;
			} else if ((i1->ind) < (i2->ind)) {
				tmp.push_back(*i1);
				++i1;
			} else {
				tmp.push_back(key_type(i2->ind, a * (i2->val)));
				++i2;
			}
		}
		// run through rest of entries and dump in
		// at most one of the loops does anything
		while (i1 != indval.cend()) {
			tmp.emplace_back(*i1);
			++i1;
		}
		while (i2 != x.nzend()) {
			tmp.emplace_back(key_type(i2->ind, a * (i2->val)));
			++i2;
		}

		// copy temp vector to indval
		// indval.resize(tmp.size());
		// std::copy(tmp.cbegin(), tmp.cend(), indval.begin());
		indval = tmp;

		clear_zeros();

		return;
	}

	// return self + ax[firstind:lastind]
	// template over sparse vector type
	template <class SVT>
	void axpy(
		const TV &a,
		const SVT &x,
		const TI &firstind,
		const TI &lastind
	) {

		if (a == TV(0)) {return;}

		// set i2 to find first ind >= firstind
		auto i2 = std::lower_bound(
			x.nzbegin(),
			x.nzend(),
			key_type(firstind, TV(0))
		);
		// where to put new vector
		if (i2 == x.nzend() || !(i2->ind < lastind)) { return; } // nothing to do
		// something to do...
		auto i1 = indval.cbegin();
		std::vector<key_type> tmp;
		while (i1 != indval.cend() && i2 != x.nzend()) {
			if (i1->ind == i2->ind) {
				TV val = (a * (i2->val)) + i1->val;
				if (val != TV(0)) {
					tmp.push_back(key_type(i1->ind, val));
				}
				++i1;
				++i2;
				if (!(i2->ind < lastind)) { break; }
			} else if (i1->ind < i2->ind) {
				tmp.push_back(*i1);
				++i1;
			} else {
				tmp.push_back(key_type(i2->ind, a * (i2->val)));
				++i2;
				if (!(i2->ind < lastind)) { break; }
			}
		}
		// run through rest of entries and dump in
		// at most one of the loops does anything
		while (i1 != indval.cend()) {
			tmp.push_back(*i1);
			++i1;
		}
		while (i2 != x.nzend() && (i2->ind < lastind)) {
			tmp.push_back(key_type(i2->ind, a * (i2->val)));
			++i2;
		}
		// copy temp vector to indval
		indval = tmp;
		return;
	}

	// return self + coeff*x[inds]
	// template over sparse vector type
	template <class SVT>
	void axpy(
		const SVT &x,
		const std::vector<TV> &coeff,
		const std::vector<TI> &inds
	) {

		// set i2 to find first ind >= firstind
		auto i2 = x.nzbegin(); // iterate over x
		// where to put new vector
		if (i2 == x.nzend()) { return; } // nothing to do
		// something to do...
		auto i1 = indval.cbegin(); // iterate over self
		size_t ii = 0; // index for inds

		std::vector<key_type> tmp;
		while (i1 != indval.cend() && i2 != x.nzend() && ii < inds.size()) {
			if (ii < i1->ind  && inds[ii] < i2->ind) {
				// need to increment ii until something happens
				++ii;
			} else if (i1->ind == ii && inds[ii] < i2->ind) {
				// no match from nonzero in x, but match from nonzero in self
				tmp.push_back(*i1);
				++i1;
				++ii;
			} else if (ii < i1->ind && inds[ii] == i2->ind) {
				// no match from nonzero in self, but match from x
				TV val = (coeff[ii] * (i2->val));
				if (val != TV(0)) {
					tmp.push_back(key_type(ii, val));
				}
				++ii;
				++i2;
			} else if (i1->ind == ii && i2->ind == inds[ii]) {
				// match from both x and self
				TV val = (coeff[ii] * (i2->val)) + i1->val;
				// std::cout << "a: " << a << " x: " << ((*i2).val)) << " y: " << (*i1).val << std::endl;
				if (val != TV(0)) {
					tmp.push_back(key_type(i1->ind, val));
				}
				++i1;
				++i2;
				++ii;
			} else {
				throw std::runtime_error("Unexpected condition in sparse axpy!");
			}
		}
		// run through rest of entries and dump in
		// at most one of the loops does anything
		while (i1 != indval.cend()) {
			tmp.push_back(*i1);
			++i1;
		}
		while (i2 != x.nzend() && ii < inds.size()) {
			if (inds[ii] < i2->ind) {
				++ii;
			} else {
				tmp.push_back(key_type(ii, coeff[ii] * (i2->val)));
				++i2;
				++ii;
			}
		}

		// copy temp vector to indval
		indval = tmp;

		return;
	}

	// scale each entry i by coeff[i]
	// assume coeff[i] is non-zero
	void scale_inplace(const std::vector<TV> &coeff) {
		auto it = nzbegin();
		while (it != nzend()) {
			it->val = (it->val * coeff[it->ind]);
			++it;
		}
	}



	// zeros out pivot in x
	// void eliminate_pivot(const SparseVector &x) {
	// 	auto piv = last();
	// 	auto pivx = x.last();
	// 	TV alpha = - piv.val / pivx.val;
	// 	// std::cout << "alpha = " << alpha << std::endl;
	// 	axpy(alpha, x);
	// }
	// scal - in place

	// add, subtract, multiply by scalar
	//inline SparseVector operator+(const SparseVector &x) {return caxpy(TV(1), x); }
	//inline SparseVector operator-(const SparseVector &x) {return caxpy(TV(-1), x); }

	void print() const {
		auto it = indval.cbegin();
		while (it != indval.cend()) {
			std::cout << *it << std::endl;
			++it;
		}
	}

	void print_row() const {
		auto it = indval.cbegin();
		while (it != indval.cend()) {
			std::cout << *it << ' ';
			++it;
		}
		std::cout << std::endl;
	}

	template <typename IO>
	void write(IO &io) const {
		auto it = indval.cbegin();
		while (it != indval.cend()) {
			io << (*it).ind << ':' << (*it).val << ',';
			++it;
		}
		io << '\n';
	}
};

// // specialized F2 implementation
// template <typename IntT, typename TI>
// class SparseVector<ModP<IntT, 2>, TI> {
// private:
// 	using TV = ModP<IntT, 2>;
//
// 	// store nonzero indices
// 	std::vector<TI> ind;
//
// 	// sort in-place
// 	void sort() {
// 		// sort in-place
// 		std::sort(ind.begin(), ind.end());
// 	}
//
// public:
//
// 	SparseVector() {}
//
// 	// copy constructor
// 	SparseVector(const SparseVector &x) : ind(x.ind) {}
//
// 	// construct from index vector
// 	SparseVector(const std::vector<TI> ind): ind(ind) {}
//
// 	SparseVector(const std::vector<TI> &inds, const std::vector<TV> &val) {
// 		size_t nz = inds.size();
// 		ind.reserve(nz);
// 		for (size_t i = 0; i < nz; i++) {
// 			if (!(val[i] == 0)) {
// 				ind.push_back(inds[i]);
// 			}
// 		}
// 	}
//
// 	// constructor that takes in integers for val
// 	SparseVector(const std::vector<TI> &inds, const std::vector<int> &val) {
// 		size_t nz = inds.size();
// 		ind.reserve(nz);
// 		for (size_t i = 0; i < nz; i++) {
// 			if ((val[i] & 0x1) == 1) {
// 				ind.push_back(inds[i]);
// 			}
// 		}
// 	}
//
// 	// constructor that loops over index and value iterators
// 	// can be iterators over different type
// 	template <typename IT1, typename IT2>
// 	SparseVector(IT1 indit, IT2 valit, size_t n) {
// 		ind.reserve(n);
// 		for (size_t i = 0; i < n; i++) {
// 			ind.emplace_back(TI(*indit++));
// 		}
// 	}
//
// 	// cosntructor that returns indicator in given index
// 	SparseVector(const TI i) {
// 		ind.push_back(i);
// 	}
//
// 	// nnz
// 	inline size_t nnz() const {
// 		return ind.size();
// 	}
//
// 	// get const iterator through nzs
// 	auto nzbegin() const {
// 		return ind.cbegin();
// 	}
//
// 	// get const iterator end of nzs
// 	auto nzend() const {
// 		return ind.cend();
// 	}
//
// 	// return ith nonzero value
// 	constexpr TV nzval(size_t i) const {
// 		return TV(1);
// 	}
//
// 	// return ith nonzero index
// 	TI nzind(size_t i) const {
// 		return ind[i];
// 	}
//
// 	// find nonzero index of last element with index < i
// 	auto find_last_nz(TI i) {
// 		auto it = std::lower_bound(
// 			ind.cbegin(),
// 			ind.cend(),
// 			i
// 		);
// 		// if there is no element with index < i, return pointer to end
// 		return it == ind.cend() ? it : --it;
// 		//return it--;
// 	}
//
//
// 	// permute in-place
// 	void permute(const std::vector<size_t>  &perm) {
// 		for (size_t i = 0; i < nnz(); i++) {
// 			ind[i] = perm[ind[i]];
// 		}
// 		sort();
// 	}
//
// 	// return result of self + ax
// 	SparseVector caxpy(const SparseVector &x, size_t xoffset=0) const {
// 		// first check if there is anything to do.
// 		auto xend = x.ind.cend() - xoffset;
//
// 		// where to put new vector
// 		std::vector<TI> tmp;
// 		tmp.reserve(nnz() + x.nnz());
// 		auto i1 = ind.cbegin();
// 		auto i2 = x.ind.cbegin();
// 		while (i1 < ind.cend() && i2 < xend) {
// 			if ((*i1) == (*i2)) {
// 				++i1;
// 				++i2;
// 			} else if ((*i1) < (*i2)) {
// 				tmp.push_back(*i1);
// 				++i1;
// 			} else {
// 				tmp.push_back(*i2);
// 				++i2;
// 			}
// 		}
// 		// run through rest of entries and dump in
// 		// at most one of the loops does anything
// 		while (i1 < ind.cend()) {
// 			tmp.push_back(*i1);
// 			++i1;
// 		}
// 		while (i2 < xend) {
// 			tmp.push_back(*i2);
// 			++i2;
// 		}
// 		return SparseVector(tmp);
// 	}
//
// 	// set
// 	// y <- ax + y
// 	void axpy(const SparseVector &x, size_t xoffset=0) {
//
// 		// check if there's anything to do
// 		auto xend = x.ind.cend() - xoffset;
// 		// if we won't update, return
// 		if (x.ind.cbegin() >= xend) { return; }
// 		// if all nzs are in x, dump
// 		if (nnz() == 0) {
// 			std::copy(x.ind.cbegin(), xend, std::back_inserter(ind));
// 			//ind = x.ind;
// 			return;
// 		}
//
// 		const SparseVector res = caxpy(x, xoffset);
// 		ind = res.ind;
// 		return;
// 	}
// 	void axpy(const TV &a, const SparseVector &x, size_t xoffset=0) {
// 		if (a == 0) { return; }
// 		axpy(x, xoffset);
// 		return;
// 	}
//
// 	// zeros out pivot in x
// 	void eliminate_pivot(const SparseVector &x) {
// 		axpy(x);
// 	}
//
// 	// add, subtract, multiply by scalar
// 	inline SparseVector operator+(const SparseVector &x) {return caxpy(x); }
// 	inline SparseVector operator-(const SparseVector &x) {return caxpy(x); }
//
//
// 	void print() {
// 		for (size_t i = 0; i < nnz(); i++) {
// 			std::cout << ind[i] << std::endl;
// 		}
// 	}
//
// 	void print_row() {
// 		for (size_t i = 0; i < nnz(); i++) {
// 			std::cout << ind[i] << " ";
// 		}
// 		std::cout << std::endl;
// 	}
//
// };

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
