#pragma once

/*
An implementation similar to sparse_vector that uses
std::set instead of maintaining sort explicitly
*/

#include <cstddef>
#include <cassert>
#include <set>
#include <vector>
#include <iostream>
#include <algorithm>
#include "field.hpp"
#include "abstract_vector.hpp"

// template over index type and type of value
template <typename TV, typename TI=size_t>
class SetVector {
private:
	using key_type = nzpair<TI, TV>;

	std::set<key_type> indval;

public:

	SetVector() {}

	SetVector(const std::set<key_type> indval) : indval(indval) {}

	SetVector(const std::vector<TI> &ind, const std::vector<TV> &val) {
		assert (ind.size() == val.size());
		for (size_t i = 0; i < ind.size(); i++) {
			indval.emplace(key_type(ind[i], val[i]));
		}
	}

	// constructor that loops over index and value iterators
	// can be iterators over different type
	template <typename IT1, typename IT2>
	SetVector(IT1 indit, IT2 valit, size_t n) {
		for (size_t i = 0; i < n; i++) {
			indval.emplace_hint(indval.end(), key_type(TI(*indit++), TV(*valit++)));
		}
	}

	// constructor that returns indicator in given index
	SetVector(const TI i) {
		indval.emplace_hint(indval.end(), key_type(i, TV(1)));
	}

	// get const iterator through nzs
	inline auto nzbegin() const { return indval.cbegin(); }
	inline auto nzend() const { return indval.cend(); }
	// get non const iterator through nzs
	inline auto nzbegin() { return indval.cbegin(); }
	inline auto nzend() { return indval.cend(); }

    // nnz
    inline size_t nnz() const { return indval.size(); }

	// returns iterator pointing to first element that is not less than i
	inline auto lower_bound(const TI &i) {
		return indval.lower_bound(key_type(i, TV(0)));
	}
	inline auto lower_bound(const TI &i) const {
		return indval.lower_bound(key_type(i, TV(0)));
	}
	// returns iterator pointing to first element that greater than i
	inline auto upper_bound(const TI &i) {
		return indval.upper_bound(key_type(i, TV(0)));
	}
	inline auto upper_bound(const TI &i) const {
		return indval.upper_bound(key_type(i, TV(0)));
	}


	// set index
	auto set(const key_type &k) {
		auto ret = indval.emplace(k); // key value pair
		auto it = ret.first;
		if (!ret.second) {
			// delete this entry
			it = indval.erase(it);
			// put in the new key
			it = indval.emplace_hint(it, k);
		}
		return it;
	}
	inline auto set(const TI ind, const TV val) {return set(key_type(ind, val));}
	// set index with hint
	// replace the item at the iterator location with k
	template <typename itT>
	auto replace(itT &it, const key_type &k) {
		// delete this entry
		it = indval.erase(it);
		// put in the new key
		it = indval.emplace_hint(it, k);
		return it;
	}
	template <typename itT>
	inline auto replace(itT &it, const TI ind, const TV val) { return replace(it, key_type(ind, val));}
	template <typename itT>
	inline auto replace(itT &it, const TV val) { return replace(it, key_type((*it).ind, val));}

	inline const key_type& lastnz() const {
		return *(--nzend());
	}

	// get element with index ind.  If no such element, return zero
	TV get(const TI ind) const {
		auto it = indval.find(key_type(ind));
		if (it == indval.end()) {
			return TV(0);
		}
		return (*it).val;
	}

	void permute(const std::vector<size_t> &perm) {
		std::set<key_type> indval_new;
		for (auto it = nzbegin(); it != nzend(); ++it) {
			indval_new.emplace(key_type(perm[(*it).ind], (*it).val));
		}
		indval = indval_new;
	}

    // self + ax
	// template over sparse vector type
	template <class SVT>
    void axpy(
        const TV &a,
        const SVT &x
    ) {

        auto i1 = indval.begin();
        auto i2 = x.nzbegin();
        while (i1 != indval.end() && i2 != x.nzend()) {
            if ((*i1).ind == (*i2).ind) {
                TV val = (a * (*i2).val) + (*i1).val;
				// delete the current key
				i1 = indval.erase(i1);
                if (!(val == 0)) {
                    // put new value with this index
                    i1 = indval.emplace_hint(i1, nzpair((*i2).ind, val));
					++i1;
                }
                ++i2;
            } else if (*i1 < *i2) {
                // need to catch up i1
                ++i1;
            } else {
                // i2 has lower index
                i1 = indval.emplace_hint(i1, nzpair((*i2).ind, a * (*i2).val) );
                ++i2;
            }
        }
        // run through rest of entries in i2
        while (i2 != x.nzend()) {
            i1 = indval.emplace_hint(i1, nzpair((*i2).ind, a * (*i2).val) );
            ++i2;
        }
        return;
    }

	// self + ax[firstind:lastind]
	// used in triangular solves
	// template over sparse vector type
	template <class SVT>
	void axpy(
		const TV &a,
		const SVT &x,
		const TI &firstind,
		const TI &lastind
	) {

		// set i1 to find first ind >= firstind
		auto i1 = indval.lower_bound(key_type(firstind, TV(0)));

		// set i2 to find first ind >= firstind
		auto i2 = x.lower_bound(firstind);
		if (!((*i2).ind < lastind) || i2 == x.nzend()) { return; } // nothing to do
		while (i1 != indval.end() && i2 != x.nzend()) {
			if ((*i1).ind == (*i2).ind) {
				TV val = (a * (*i2).val) + (*i1).val;
				// delete the current key
				i1 = indval.erase(i1);
				if (!(val == 0)) {
					// put new value with this index
					i1 = indval.emplace_hint(i1, key_type((*i2).ind, val));
					++i1;
				}
				++i2;
				if (!((*i2).ind < lastind)) { break; }
			} else if (*i1 < *i2) {
				// need to catch up i1
				++i1;
			} else {
				// i2 has lower index
				i1 = indval.emplace_hint(i1, key_type((*i2).ind, a * (*i2).val) );
				++i2;
				if (!((*i2).ind < lastind)) { break; }
			}
		}
		// run through rest of entries in i2
		while (i2 != x.nzend() && !((*i2).ind < lastind)) {
			i1 = indval.emplace_hint(i1, nzpair((*i2).ind, a * (*i2).val) );
			++i2;
		}
		return;
	}

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
};
