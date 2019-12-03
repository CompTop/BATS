#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include <vector>
#include <set>

#include <linalg/field.h>
#include <linalg/sparse_vector.h>
#include <linalg/set_vector.h>
#include <linalg/col_matrix.h>

#define F2 ModP<int, 2>
#define F3 ModP<int, 3>
#define Q Rational<int>
#define F5 ModP<int, 5>

TEST_CASE_TEMPLATE("Col Matrix Solve", F, F2, F3, F5, Q, int) {

	using VT = SparseVector<F, size_t>;
	using MatT = ColumnMatrix<VT>;

	MatT I = MatT::identity(5);

	std::vector<size_t> ind = {0,2,3};
    std::vector<F> val = {-1, 1, -1};

	VT y(ind, val);

	SUBCASE("U solve") {
		auto x1 = u_solve(I, y);

		REQUIRE(x1 == y);
	}

	SUBCASE("L solve") {
		auto x2 = l_solve(I, y);

		REQUIRE(x2 == y);
	}


}
