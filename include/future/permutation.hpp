#pragma once

#include <vector>
#include <list>
#include <numeric>
#include <iostream>

namespace bats {
namespace future {

// forward declarations
struct ElementaryPermutation;
class CompositePermutation;


/*
	Permutation types
*/

struct ElementaryPermutation {
/*
	elementary permutation to swap two elements of a set
*/
	size_t i;
	size_t j;

	ElementaryPermutation() {}
	ElementaryPermutation(size_t i, size_t j) : i(i), j(j) {}

	// apply elementary permutation to data
	// assume random access iterator type
	template <typename T>
	inline T& operator()(T& a) const { std::swap(a[i], a[j]); return a;}
	// handle lvalues
	template <typename T>
	inline T& operator()(T&& a) const { std::swap(a[i], a[j]); return a; }

	// specialization to CompositePermutation
	inline CompositePermutation& operator()(CompositePermutation &a) const;

	friend std::ostream& operator<<(std::ostream& os, const ElementaryPermutation& p);
};
std::ostream& operator<<(std::ostream& os, const ElementaryPermutation& p) {
	os << '(' << p.i << ',' << p.j << ')';
	return os;
}

class CompositePermutation {
/*
	Concatenation of elementary permutations
*/
	// typedef std::vector<ElementaryPermutation> container_type;
	typedef std::list<ElementaryPermutation> container_type;
private:
	container_type perm;

public:

	CompositePermutation() {};

	inline void append(const ElementaryPermutation& p) { perm.emplace_back(p); }
	inline void swap(size_t i, size_t j) { perm.emplace_back(ElementaryPermutation(i,j)); }

	// apply elementary permutation to data
	// assume random access iterator type
	template <typename T>
	T& operator()(T& a) const {
		for (auto p : perm) {
			a = p(a);
		}
		return a;
	}
	// handle lvalues
	template <typename T>
	T& operator()(T&& a) const {
		for (auto p : perm) {
			a = p(a);
		}
		return a;
	}

	inline CompositePermutation& operator()(const ElementaryPermutation &p) { perm.emplace(perm.begin(), p); return *this; }

	friend std::ostream& operator<<(std::ostream& os, const CompositePermutation& p);
};
std::ostream& operator<<(std::ostream& os, const CompositePermutation& p) {
	if (p.perm.size() == 0) {
		os << "()";
	} else {
		for (auto pi : p.perm) {
			os << pi;
		}
	}
	return os;
}

// specialization to CompositePermutation
inline CompositePermutation& ElementaryPermutation::operator()(CompositePermutation &a) const {
	a.append(*this);
	return a;
}


class Permutation {
private:
	std::vector<size_t> p;

	inline void fill_identity() {
		std::iota(p.begin(), p.end(), 0);
	}

public:

	Permutation() {}
	Permutation(size_t n) : p(n) {
		fill_identity();
	}
};

}
}
