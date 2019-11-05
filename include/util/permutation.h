#pragma once
/*
utilities for permutations
*/
#include <vector>
#include <random>
#include <cstddef>
#include <iostream>
#include <algorithm>

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

template <typename T>
std::vector<size_t> stable_sortperm(const std::vector<T>& data) {
    std::vector<size_t> perm(data.size());
    std::iota(perm.begin(), perm.end(), 0);
    std::stable_sort(
        perm.begin(),
        perm.end(),
        [&](const size_t& a, const size_t& b) {
            return data[a] < data[b];
        }
    );
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


// random permutation on n indices
std::vector<size_t> rand_perm(const size_t n) {

    std::default_random_engine generator;
    std::uniform_real_distribution<float> distribution(0.0,1.0);

    std::vector<float> x(n);
    for (auto it = x.begin(); it < x.end(); ++it) {
        *it = distribution(generator);
    }

    return bats::sortperm(x);
}

// compute inverse permutation
inline std::vector<size_t> inv_perm(const std::vector<size_t> &p) {
    return bats::sortperm(p);
}
