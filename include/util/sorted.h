#pragma once
/*
Some operations on sorted containers
*/
#include <cstddef>
#include <vector>
#include <iostream>
#include <algorithm>
#include "common.h"
/*
sets c = intersection(a, b)
over-writes c
TODO: template over container type as well
*/
template <typename T>
void intersect_sorted(const std::vector<T> &a, const std::vector<T> &b, std::vector<T> &c) {
    c.clear();
    auto ia = a.cbegin();
    auto ib = b.cbegin();
    while (ia < a.cend() && ib < b.cend()) {
        if (*ia < *ib) {
            // std::cout << "a small" << std::endl;
            ++ia;
        } else if (*ib < *ia) {
            // std::cout << "b small" << std::endl;
            ++ib;
        } else {
            // std::cout << "agree" << std::endl;
            // std::cout << *ia << std::endl;
            // *ia == *ib
            c.emplace_back(*ia);
            ++ia;
            ++ib;
        }
    }
    // for (auto i : c) {
    //   std::cout << i << std::endl;
    // }
}

/*
sets c = intersection(a, b, (-inf, maxval))
over-writes c
TODO: template over container type as well
*/
template <typename T>
void intersect_sorted_lt(const std::vector<T> &a, const std::vector<T> &b, T maxval, std::vector<T> &c) {
    c.clear();
    auto ia = a.cbegin();
    auto ib = b.cbegin();
    while (ia < a.cend() && ib < b.cend()) {
        if (*ia < *ib) {
            ++ia;
            if (!(*ia < maxval)) {break;}
        } else if (*ib < *ia) {
            ++ib;
            if (!(*ib < maxval)) {break;}
        } else {
            // *ia == *ib
            c.emplace_back(*ia);
            ++ia;
            ++ib;
          if (!(*ia < maxval) || !(*ib < maxval)) {break;}
        }
    }
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


// fill a vector that will return a sort permutation on data
template <typename T>
void fill_sortperm(
    const T& begin,
    const T& end,
    std::vector<size_t> &perm
) {
    perm.resize(end - begin);
    std::iota(perm.begin(), perm.end(), 0);
    std::sort(
        perm.begin(),
        perm.end(),
        [&](const size_t& a, const size_t& b) {
            return *(begin+a) < *(begin+b);
        }
    );
    return;
}

// fill a vector that will return a sort permutation on data
template <typename T>
inline void fill_sortperm(
    const std::vector<T> &data,
    std::vector<size_t> &perm
) {
    return fill_sortperm(data.cbegin(), data.cend(), perm);
}

namespace bats{
template <typename T>
std::vector<size_t> sortperm(const std::vector<T>& data) {
    std::vector<size_t> perm(data.size());
    fill_sortperm(data, perm);
    return perm;
}
}


// reorders data by permutation
template <typename T>
void apply_perm(
    T* begin,
    std::vector<T> &tmp,
    const std::vector<size_t> &perm
) {
    // use allocated array to do permutation
    tmp.resize(perm.size());
    // put permuted array in tmp
    for (size_t i = 0; i < perm.size(); i++) {
        tmp[i] = *(begin + perm[i]);
    }
    // put tmp in range
    for (size_t i = 0; i < perm.size(); i++) {
        *(begin+i) = tmp[i];
    }
    return;
}


template <typename T>
void apply_perm(
    T* begin,
    const std::vector<size_t> &perm
) {
    std::vector<T> tmp;
    apply_perm(begin, tmp, perm);
}

template <typename T>
inline void apply_perm(
    std::vector<T> &data,
    const std::vector<size_t> &perm
) {
    assert (perm.size() == data.size());
    return apply_perm(data.data(), perm);
}


// partial perm
// given a vector of inds, and a total size n
// return a vector ret of length n so that
// ret[i] = index in inds if i in inds
// ret[i] = NO_IND otherwise
std::vector<size_t> partial_perm(const std::vector<size_t> &ind, const size_t n) {
    std::vector<size_t> ret(n, bats::NO_IND);
    for (size_t k = 0; k < ind.size(); k++) {
        ret[ind[k]] = k;
    }
    return ret;
}
// TODO: could also implement this with a map to save space

// sort indices so they appear in permutation order
inline void sort_ind_by_perm(std::vector<size_t> &ind, const std::vector<size_t> &perm) {
    std::sort(
        ind.begin(),
        ind.end(),
        [&](const size_t& a, const size_t& b) {
            return perm[a] < perm[b];
        }
    );
}

// fill indperm so that ind[ind_perm] is in perm order
void fill_partial_sortperm(
    const std::vector<size_t> &ind,
    const std::vector<size_t> &perm,
    std::vector<size_t> &indperm
) {
    indperm.resize(ind.size());
    std::iota(indperm.begin(), indperm.end(), 0);
    std::sort(
        indperm.begin(),
        indperm.end(),
        [&](const size_t& a, const size_t& b) {
            return perm[ind[a]] < perm[ind[b]];
        }
    );
    return;
}

// sort so that ind1 is ordered by perm
// apply same sort to ind2
void sort_ind_pair_by_perm(
    std::vector<size_t> &ind1,
    std::vector<size_t> &ind2,
    const std::vector<size_t> &perm
) {
    std::vector<size_t> indperm, tmp;
    fill_partial_sortperm(ind1, perm, indperm);
    apply_perm(ind1.data(), tmp, indperm);
    apply_perm(ind2.data(), tmp, indperm);
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
