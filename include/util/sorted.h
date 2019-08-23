#pragma once
/*
Some operations on sorted containers
*/
#include <cstddef>
#include <vector>
#include <iostream>
#include <algorithm>

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
