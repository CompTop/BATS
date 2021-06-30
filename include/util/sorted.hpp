#pragma once
/*
Some operations on sorted containers
*/
#include <cstddef>
#include <vector>
#include <set>
#include <iostream>
#include <algorithm>
#include "common.hpp"
#include "permutation.hpp"

namespace bats {
namespace util {

/*
sets c = intersection(a, b)
over-writes c
TODO: template over container type as well
*/
template <typename T, typename C1, typename C2>
void intersect_sorted(const C1 &a, const C2 &b, std::vector<T> &c) {
    c.clear();
    auto ia = a.cbegin();
    auto ib = b.cbegin();
    while (ia != a.cend() && ib != b.cend()) {
        if (*ia < *ib) {
            ++ia;
        } else if (*ib < *ia) {
            ++ib;
        } else {

            c.emplace_back(*ia);
            ++ia;
            ++ib;
        }
    }
}

template <typename T, typename C1, typename C2>
void intersect_sorted(const C1 &a, const C2 &b, std::set<T> &c) {
    c.clear();
    auto ia = a.cbegin();
    auto ib = b.cbegin();
    auto ic = c.begin();
    while (ia != a.cend() && ib != b.cend()) {
        if (*ia < *ib) {
            ++ia;
        } else if (*ib < *ia) {
            ++ib;
        } else {

            ic = c.emplace_hint(ic, *ia);
            ++ia;
            ++ib;
        }
    }
}


// return true if intersection intersect_sorted would have something
// template over container types
template <typename C1, typename C2>
bool has_intersect_sorted(const C1 &a, const C2 &b) {
    auto ia = a.cbegin();
    auto ib = b.cbegin();
    while (ia != a.cend() && ib != b.cend()) {
        if (*ia < *ib) {
            ++ia;
        } else if (*ib < *ia) {
            ++ib;
        } else {
            // *ia == *ib
            return true;
        }
    }
    return false;
}

/*
sets c = intersection(a, b, (-inf, maxval))
over-writes c
TODO: template over container type as well
*/
template <typename T, typename C1, typename C2>
void intersect_sorted_lt(const C1 &a, const C2 &b, const T maxval, std::vector<T> &c) {
    c.clear();
    auto ia = a.cbegin();
    auto ib = b.cbegin();
    while (ia != a.cend() && ib != b.cend()) {
        if (*ia < *ib) {
            ++ia;
            if (ia == a.cend() || !(*ia < maxval)) {break;}
        } else if (*ib < *ia) {
            ++ib;
            if (ib == b.cend() || !(*ib < maxval)) {break;}
        } else {
            // *ia == *ib
            c.emplace_back(*ia);
            ++ia;
            ++ib;
          if (ia == a.cend() || ib == b.cend() || !(*ia < maxval) || !(*ib < maxval)) {break;}
        }
    }
}

// return true if intersection intersect_sorted_lt would have something
// template over container types
template <typename T, typename C1, typename C2>
bool has_intersect_sorted_lt(const C1 &a, const C2 &b, const T maxval) {
    auto ia = a.cbegin();
    auto ib = b.cbegin();
    while (ia != a.cend() && ib != b.cend()) {
        if (*ia < *ib) {
            ++ia;
            if (ia == a.cend() || !(*ia < maxval)) {break;}
        } else if (*ib < *ia) {
            ++ib;
            if (ib == b.cend() || !(*ib < maxval)) {break;}
        } else {
            // *ia == *ib
            return true;
        }
    }
    return false;
}


// copy x into y, then sort y
template <typename T>
void sort_into(const std::vector<T> &x, std::vector<T> &y) {
    // copy x into y
    y.resize(x.size());
    std::copy(x.cbegin(), x.cend(), y.begin());
    // sort y
    std::sort(y.begin(), y.end());
}




// inline void sort_ind_by_perm(
//     std::vector<size_t> &ind1,
//     std::vector<size_t> &ind2,
//     const std::vector<size_t> &perm
// )

// sort and reduce ind and val
template <typename TI, typename TV>
void sort_sum_reduce(
    std::vector<TI> &ind,
    std::vector<TV> &val,
    const size_t offset,
    std::vector<size_t> &perm,
    std::vector<TI> &tmpind,
    std::vector<TV> &tmpval
) {
    // sort ind and val simultaneously
    fill_sortperm(ind.cbegin() + offset, ind.cend(), perm);
    apply_perm(ind.data() + offset, tmpind, perm);
    apply_perm(val.data() + offset, tmpval, perm);

    size_t fill_ptr = offset;
    size_t look_ptr = offset + 1;
    while (ind.cbegin() + look_ptr < ind.cend()) {
        if (ind[look_ptr] == ind[fill_ptr]) {
            // keep adding to index at fill_ptr
            val[fill_ptr] += val[look_ptr];
            ++look_ptr;
        } else {
            // if accumulation was zero, start accumulating at same spot
            if (!(val[fill_ptr] == 0)) {
                fill_ptr++;
            }
            // start accumulating at new next fill_ptr
            ind[fill_ptr] = ind[look_ptr];
            val[fill_ptr] = val[look_ptr];
            ++look_ptr;
        }
    }

    // resize to fill_ptr
    // if last accumulation is zero, ignore
    if (!(val[fill_ptr] == 0)) {
        fill_ptr++;
    }
    ind.resize(fill_ptr);
    val.resize(fill_ptr);

}


// find nonzero index of last element with index < i
template<typename T, typename TI>
size_t find_sorted_lt(const TI &begin, const TI &end, const T &v) {
    if (begin == end) { return bats::NO_IND; }
    TI it = std::lower_bound(
        begin,
        end,
        v
    );
    it--;
    if (*it < v) {
        return it - begin;
    }
    return bats::NO_IND;
    //return it--;
}

// return complement of sorted vector in range [0,n)
// assume vector ind is sorted
template <typename T>
std::vector<T> sorted_complement(const std::vector<T> &ind, size_t n) {
    std::vector<T> cind;
    cind.reserve(n - ind.size());
    size_t i = 0;
    auto ip = ind.begin();
    while (ip != ind.end()) {
        if (*ip == i) {
            // we do not append i to the complement
            i++;
            ip++;
        } else {
            // we append i to the complement
            cind.emplace_back(i++);
        }
    }
    while (i < n) {
        // fill in the rest
        cind.emplace_back(i++);
    }
    return cind;
}

} // namespace util
} // namespace bats
