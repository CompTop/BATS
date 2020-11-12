#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include <vector>
#include <set>

#include <linalg/field.hpp>
#include <linalg/sparse_vector.hpp>
#include <linalg/set_vector.hpp>

#define F2 ModP<int, 2>
#define F3 ModP<int, 3>
#define Q Rational<int>
#define F5 ModP<int, 5>

TEST_CASE_TEMPLATE("Sparse Vector Ops", F, F2, F3, F5, Q, int) {

	SUBCASE("Construction") {
		std::vector<size_t> ind;
	    std::vector<int> val_int;

		ind = {1,3,5};
		val_int = {-1,1,-1};

		SetVector<F, size_t> x(ind.cbegin(), val_int.cbegin(), ind.size());
		SparseVector<F, size_t> y(ind.cbegin(), val_int.cbegin(), ind.size());

		REQUIRE(x == y);
	}

	SUBCASE("axpy") {
		std::vector<size_t> ind;
		std::vector<F> val;
		ind = {1,3};
	    val = {-1,1};
	    SparseVector<F, size_t> a(ind, val);
	    ind = {0,3};
	    val = {-1,1};
	    SparseVector<F, size_t> b(ind, val);

	    a.axpy(-1, b);

		ind = {0, 1};
		val = {1, -1};
		SparseVector<F, size_t> t(ind, val); // target
		REQUIRE(a == t);
	}

}
