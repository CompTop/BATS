#pragma once

/*
An implementation similar to sparse_vector that uses
std::set instead of maintaining sort explicitly
*/

#include <cstddef>
#include <cassert>
#include <set>
#include <iostream>
#include <algorithm>
#include "field.h"
#include "abstract_vector.h"

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

	// SetVector(const std::vector<TI> &ind, const std::vector<int> &val) {
	// 	assert (ind.size() == val.size());
	// 	for (size_t i = 0; i < ind.size(); i++) {
	// 		indval.emplace(nzpair(ind[i], TV(val[i])));
	// 	}
	// }

    // nnz
    inline size_t nnz() const { return indval.size(); }

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
	template <typename itT>
	inline auto set_hint(itT &it, const key_type &k) { return indval.emplace_hint(it, k);}
	template <typename itT>
	inline auto set_hint(itT &it, const TI ind, const TV val) { return indval.emplace_hint(it, key_type(ind, val));}

	// get element with index ind.  If no such element, return zero
	TV get(const TI ind) const {
		auto it = indval.find(key_type(ind));
		if (it == indval.end()) {
			return TV(0);
		}
		return (*it).val;
	}

    // self + ax
    void axpy(
        const TV &a,
        const SetVector &x
    ) {

        auto i1 = indval.begin();
        auto i2 = x.indval.cbegin();
        while (i1 != indval.end() && i2 != x.indval.cend()) {
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
        while (i2 != x.indval.cend()) {
            i1 = indval.emplace_hint(i1, nzpair((*i2).ind, a * (*i2).val) );
            ++i2;
        }
        return;
    }

	// self + ax[:lastind]
	// used in upper triangular solve
	void axpy(
		const TV &a,
		const SetVector &x,
		const TI &lastind
	) {

		auto i1 = indval.begin();
		auto i2 = x.indval.cbegin();
		while (i1 != indval.end() && i2 != x.indval.cend()) {
			if (!(*i2.ind < lastind)) { break; }
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
		while (i2 != x.indval.cend()) {
			if (!(*i2.ind < lastind)) { break; }
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
};
