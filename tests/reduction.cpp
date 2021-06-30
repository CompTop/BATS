#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include <vector>
#include <set>
#include <random>

#include <bats.hpp>
using namespace bats;

#define F2 ModP<int, 2>
#define F3 ModP<int, 3>
#define Q Rational<int>
#define F5 ModP<int, 5>

// number of seeds for random checks
#define N_SEEDS 4

TEST_CASE_TEMPLATE("Identity Reduction", F, F2, F3, F5, Q) {
	using VT = SparseVector<F, size_t>;
	using MatT = ColumnMatrix<VT>;

	MatT B = MatT::identity(5);
	MatT U = MatT::identity(5);

	bats::reduce_matrix(B, U);
	CHECK(B == MatT::identity(5));
	CHECK(U == MatT::identity(5));
}

TEST_CASE_TEMPLATE("Ones Reduction", F, F2, F3, F5, Q) {
	using VT = SparseVector<F, size_t>;
	using MatT = ColumnMatrix<VT>;

	std::vector<VT> cols{VT({0,1},{1,1}), VT({0,1},{1,1})};
	MatT B(2, 2, cols);
	MatT R(B);
	CHECK(R == B);
	MatT U = MatT::identity(2);

	auto p2c = bats::reduce_matrix(R, U);
	// basic sanity checks
	CHECK(U.is_upper());
	CHECK(R.is_reduced());
	CHECK(R * u_inv(U) == B);

	// check known facts about this reduction
	CHECK(p2c[1] == 0);
	cols = {VT({0,1}, {1,1}), VT()};
	CHECK(R == MatT(2,2,cols));
	cols = {VT({0},{1}), VT({0,1},{-1,1})};
	CHECK(U == MatT(2,2,cols));

}

TEST_CASE_TEMPLATE("Random Seeds", T, F2, F3, F5, Q) {
	using VT = SparseVector<T, size_t>;
	using MatT = ColumnMatrix<VT>;

	SUBCASE("Square") {
		for ( unsigned seed = 0; seed < N_SEEDS; seed++) {
			std::default_random_engine generator(seed);
			auto A = MatT::random(10, 10, 0.4, 1, generator);
			MatT R(A);
			auto U = MatT::identity(A.ncol());

			bats::reduce_matrix(R, U);

			// sanity checks
			CHECK(U.is_upper());
			CHECK(R.is_reduced());
			CHECK(R * u_inv(U) == A);
		}
	}

	SUBCASE("Short") {
		for ( unsigned seed = 0; seed < N_SEEDS; seed++) {
			std::default_random_engine generator(seed);
			auto A = MatT::random(10, 20, 0.4, 1, generator);

			MatT R(A);
			auto U = MatT::identity(A.ncol());

			bats::reduce_matrix(R, U);

			// sanity checks
			CHECK(U.is_upper());
			CHECK(R.is_reduced());
			CHECK(R * u_inv(U) == A);
		}
	}

	SUBCASE("Tall") {
		for ( unsigned seed = 0; seed < N_SEEDS; seed++) {
			std::default_random_engine generator(seed);
			auto A = MatT::random(20, 10, 0.4, 1, generator);

			MatT R(A);
			auto U = MatT::identity(A.ncol());

			bats::reduce_matrix(R, U);

			// sanity checks
			CHECK(U.is_upper());
			CHECK(R.is_reduced());
			CHECK(R * u_inv(U) == A);
		}
	}

}
