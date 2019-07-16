#pragma once

#include <vector>
#include <iostream>
#include <algorithm>

// template over type of values
template <typename TI, typename TV>
class SparseVector
{
private:
  // store index-value pairs
  std::vector<std::pair<TI, TV>> indval;
  // // non-zero indices
  // std::vector<TI> ind;
  // // non-zero values
  // std::vector<TV> val;
public:

	SparseVector() {}

  SparseVector(std::vector<TI> ind, std::vector<TV> val) {
    size_t nz = ind.size();
    for (size_t i = 0; i < nz; i++) {
      indval.push_back(std::make_pair(ind[i], val[i]));
    }
  }

	// constructor that takes in integers for val
	SparseVector(std::vector<TI> ind, std::vector<int> val) {
		size_t nz = ind.size();
		for (size_t i = 0; i < nz; i++) {
      indval.push_back(std::make_pair(ind[i], TV(val[i])));
    }
	}

	// cosntructor that returns indicator in given index
	SparseVector(TI i) {
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

	// find last element with index < i
	auto find_last_nz(TI i) {
		auto it = std::lower_bound(
			indval.cbegin(),
			indval.cend(),
			std::make_pair(i, TV(0))
		);
		return it--;
	}


  // return last nonzero
  std::pair<TI, TV> last() const {
    return indval.back();
  }

	// return ith nonzero value
	TV nzval(size_t i) const {
		return indval[i].second;
	}

	// return ith nonzero index
	TV nzind(size_t i) const {
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

  // sort in-place
  void sort() {
    // sort in-place
    std::sort(indval.begin(), indval.end());
  }

  // permute in-place
  void permute(const std::vector<size_t>  &perm) {
    for (size_t i = 0; i < nnz(); i++) {
      indval[i].first = perm[indval[i].first];
    }
    sort();
  }

  // set
  // y <- ax + y
  void axpy(const TV &a, const SparseVector &x, size_t xoffset=0) {
    // first check if there is anything to do.
    if (x.nnz() == 0) { return; }
		auto xend = x.indval.cend() - xoffset;
    if (nnz() == 0) {
			std::copy(x.indval.begin(), xend, std::back_inserter(indval));
			//indval = x.indval;
			return;
		}

    // where to put new vector
    std::vector<std::pair<TI, TV>> tmp;
    auto i1 = indval.cbegin();
		auto i2 = x.indval.cbegin();
		do {
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
		} while (i1 < indval.cend() && i2 < xend);
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
		indval = tmp;
		return;


  }

  // axpy - in place
  // scal - in place

  // add, subtract, multiply by scalar

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

// TODO: sparse F2 vector implementation
template <typename TI>
class SparseF2Vector
{
private:
  // non-zero indices
  std::vector<TI> ind;
public:


  // get index and set index
};
