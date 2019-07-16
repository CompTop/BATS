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



  // return last nonzero
  std::pair<TI, TV> last() const {
    return indval.back();
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
  void axpy(const TV &a, const SparseVector &x) {
    // first check if there is anything to do.
    if (x.nnz() == 0) { return; }
    if (nnz() == 0) { indval = x.indval; return; }

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
		} while (i1 < indval.cend() && i2 < x.indval.cend());
		// run through rest of entries and dump in
		// at most one of the loops does anything
		while (i1 < indval.cend()) {
			tmp.push_back(*i1);
			++i1;
		}
		while (i2 < x.indval.cend()) {
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
