#pragma once
#include <set>
#include <random>

namespace bats {
namespace util {

// return union of two sets
// template over containter types
template <typename CT1, typename CT2>
std::set<size_t> set_union(const CT1 &s1, const CT2 &s2) {
    std::set<size_t> s;
    for (auto x : s1) {
        s.insert(x);
    }
    for (auto x : s2) {
        s.insert(x);
    }
    return s;
}

// intersection of two sets
// template over container type
template <typename C1, typename C2>
std::set<size_t> set_intersection(const C1 &a, const C2 &b) {
	std::set<size_t> s;

	auto is = s.end();
    auto ia = a.cbegin();
    auto ib = b.cbegin();
    while (ia != a.cend() && ib != b.cend()) {
        if (*ia < *ib) {
            ++ia;
        } else if (*ib < *ia) {
            ++ib;
        } else {
            s.emplace_hint(is, *ia);
			is = s.end();
            ++ia;
            ++ib;
        }
	}
	return s;
}

// random subset of n numbers of size ns
std::set<size_t> random_subset(
	const size_t n,
	const size_t ns
) {
	std::default_random_engine generator;
	std::uniform_int_distribution<size_t> distribution(0,n-1);

	std::set<size_t> inds;

	while (inds.size() < ns) {
		inds.insert( distribution(generator) );
	}

	return inds;
}

// convert vector to set
template <typename T>
std::set<T> to_set(
	const std::vector<T> &v
) {
	std::set<T> s;

	for (auto x : v) {
		s.insert(x);
	}

	return s;
}

} // namespace util
} // namespace bats
