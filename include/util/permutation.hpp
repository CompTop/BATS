#pragma once
/*
utilities for permutations
*/
#include <vector>
#include <random>
#include <cstddef>
#include <iostream>
#include <algorithm>
#include <cmath>

namespace bats {
namespace util {

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


template <typename T>
std::vector<size_t> sortperm(const std::vector<T>& data) {
    std::vector<size_t> perm(data.size());
    fill_sortperm(data, perm);
    return perm;
}

// fill a vector that will return a sort permutation on data
template <typename TI>
std::vector<size_t> sortperm(
    const TI& begin,
    const TI& end
) {
    std::vector<size_t> perm(end - begin);
    std::iota(perm.begin(), perm.end(), 0);
    std::sort(
        perm.begin(),
        perm.end(),
        [&](const size_t& a, const size_t& b) {
            return *(begin+a) < *(begin+b);
        }
    );
    return perm;
}

/**
filll a sortperm using custom comparator

@param first    random acces iterator at beginning of range to sort
@param last     random access iterator just past the last element of range to sort
@param comp     comparison function, comp(a, b) should return whether a should come before b
*/
template <typename RAI, class Compare>
std::vector<size_t> sortperm(
    RAI first,
    RAI last,
    Compare comp
) {
    std::vector<size_t> perm(std::distance(first, last));
    std::iota(perm.begin(), perm.end(), 0);
    std::sort(
        perm.begin(),
        perm.end(),
        [&](const size_t& a, const size_t& b) {
            return comp(*(first+a), *(first+b));
        }
    );
    return perm;
}

/**
filll a sortperm using custom comparator

@param first    random acces iterator at beginning of range to sort
@param last     random access iterator just past the last element of range to sort
@param comp     comparison function, comp(a, b) should return whether a should come before b
*/
template <class Compare>
std::vector<size_t> sortperm(
    const size_t first,
    const size_t last,
    Compare comp
) {
    std::vector<size_t> perm(last - first);
    std::iota(perm.begin(), perm.end(), 0);
    std::sort(
        perm.begin(),
        perm.end(),
        [&](const size_t& a, const size_t& b) {
            return comp((first+a), (first+b));
        }
    );
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

// fill a vector that will return indices of first k elements of range
template <typename TI>
std::vector<size_t> firstk(
    const TI& begin,
    const TI& end,
    const size_t k
) {
    std::vector<size_t> perm(end - begin);
    std::iota(perm.begin(), perm.end(), 0);
    std::partial_sort(
        perm.begin(),
        perm.begin() + k,
        perm.end(),
        [&](const size_t& a, const size_t& b) {
            return *(begin+a) < *(begin+b);
        }
    );
    perm.resize(k); // resize to first k elements
    return perm;
}

// get top k indices
// does not sort
template <typename T>
std::vector<size_t> top_k(
    const std::vector<T> &data,
    const size_t k
) {

    std::vector<size_t> perm(data.size());
    std::iota(perm.begin(), perm.end(), 0);

    auto begin = data.begin();

    std::nth_element(
        perm.begin(),
        perm.begin() + k,
        perm.end(),
        [&](const size_t& a, const size_t& b) {
            return *(begin+a) < *(begin+b);
        }
    );

    perm.resize(k); // resize to first k elements
    return perm;
}

// get top p indices
// does not sort
template <typename T>
std::vector<size_t> top_p(
    const std::vector<T> &data,
    const double p
) {
    size_t k = (size_t) std::round(data.size() * p);

    return top_k(data, k);
}



// reorders data by permutation
template <typename T, typename T2>
void apply_perm(
    T begin,
    std::vector<T2> &tmp,
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

/*
uses std::swap instead of copying to temporary vector
data[i] should map to data[perm[i]]
i.e. [2, 0, 1] applied to [0.0, 1.0, 2.0] should be [1.0, 2.0, 0.0]
*/
template <typename T>
inline void apply_perm_swap(
    std::vector<T> &data,
    const std::vector<size_t> &perm
) {
    assert (perm.size() == data.size());
    std::vector<T> tmp(data.size());
    for (size_t i = 0; i < perm.size(); i++) {
        std::swap(tmp[perm[i]], data[i]);
    }
    // put tmp in range
    for (size_t i = 0; i < perm.size(); i++) {
        std::swap(tmp[i], data[i]);
    }
}

/*
uses std::swap instead of copying to temporary vector
data[perm[i]] should map to data[i]
i.e. [2, 0, 1] applied to [0.0, 1.0, 2.0] should be [2.0, 0.0, 1.0]
*/
template <typename T>
inline void apply_iperm_swap(
    std::vector<T> &data,
    const std::vector<size_t> &perm
) {
    assert (perm.size() == data.size());
    std::vector<T> tmp(data.size());
    for (size_t i = 0; i < perm.size(); i++) {
        std::swap(tmp[i], data[perm[i]]);
    }
    // put tmp in range
    for (size_t i = 0; i < perm.size(); i++) {
        std::swap(tmp[i], data[i]);
    }
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

    return sortperm(x);
}

// compute inverse permutation
inline std::vector<size_t> inv_perm(const std::vector<size_t> &p) {
    // return sortperm(p);
    std::vector<size_t> invp(p.size());
    for (size_t j = 0; j < p.size(); j++) {
        invp[p[j]] = j;
    }
    return invp;
}

/*
count number of inversions
i<j, but p[i] > p[j]
short brute-force for-loop.  For small permutations this should be ok
*/
template <typename T>
size_t perm_inversions(const std::vector<T> &p) {
    size_t n = p.size();
    size_t ct = 0;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i+1; j < n; j++) {
            ct += (p[i] > p[j]);
        }
    }
    return ct;
}

// return sign of sort permutation of an array p
template <typename T>
inline int perm_sign(const std::vector<T> &p) {
    return (perm_inversions(p) & 0x1) == 1 ? -1 : 1;
}

} // namespace util

// The following functions when called will not need to use bats::util::function_name

// Return an identity permutation, k is its length
std::vector<size_t> identity_perm(const size_t& k){
    std::vector<size_t> l(k);
    std::iota(l.begin(), l.end(), 0);
    return l;
}

/*
If we sort vector v, this function will return the indices of sorted values 
*/
template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

    // initialize original index locations
    std::vector<size_t> idx(v.size());
    // fill the vector with increasing numbers starting from 0
    std::iota(idx.begin(), idx.end(), 0); 

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values 
    std::stable_sort(idx.begin(), idx.end(),
        [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

    return idx;
}

/*
Make indices list to permutation, e.g., [5 0 2] to [2 0 1]

Given a list of correpsonding old indices of permutations, find the permutation
e.g., given [5 0 2], which is the corresponding old indices of an old vector v_old 
for a new vector v: 
v[0] = v_old[5]
v[1] = v_old[0]
v[2] = v_old[2]
need to find the permutation [2 0 1]
*/
template <typename T>
std::vector<size_t> find_perm_from_vector(const std::vector<T> &v) {

    // initialize original index locations
    std::vector<size_t> idx(v.size());
    // fill the vector with increasing numbers starting from 0
    std::iota(idx.begin(), idx.end(), 0); 

    // sort indexes based on comparing values in v
    std::stable_sort(idx.begin(), idx.end(),
        [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

    // get the desired permutation
    std::vector<size_t> perm(v.size());
    for(size_t i = 0; i < v.size(); i++){
        perm[idx[i]] = i;
    }
    return perm;
}


// get the permutation for permuting i-th element to the end 
std::vector<size_t> perm_to_the_end(const size_t& index, const size_t& length){
    std::vector<size_t> v;
    v.reserve(length);

    for(size_t i = 0; i < length; i++){
        if(i!=index)v.emplace_back(i);
    }
    v.emplace_back(index);
    return v;
}

// get the permutation for permuting elements, with their indices in a list, to the end 
// the index_list should be sorted 
// eg. passing ([1,3],5) will return [0,2,4,1,3]
std::vector<size_t> perm_to_the_end(const std::vector<size_t>& index_list, const size_t& length){
    std::vector<size_t> v;
    v.reserve(length);
    auto it = index_list.begin();
    for(size_t i = 0; i < length; i++){
        if(i != *it){
            v.emplace_back(i);
        }else{
            it++;
        }
    }

    for(auto& element: index_list){
        v.emplace_back(element);
    }
    return v;
}

// extend a permutation to a desired length 
// with the elements appended unmoved, e.g.
// eg., passing ([2,0,1], 6) return [2,0,1,3,4,5]
std::vector<size_t> extension_perm(const std::vector<size_t>& perm, const size_t& length){
    if(!perm.empty()){
        std::vector<size_t> v;
        v.reserve(length);
        size_t i = 0;
        for(auto& element: perm){
            v.emplace_back(element);
            i++;
        }
        while (i<length)
        {
            v.emplace_back(i);
            i++;
        }
        return v;
    }else{
        return identity_perm(length);
    }
    
}

} // namespace bats
