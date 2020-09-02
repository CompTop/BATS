#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include <vector>
#include <set>
#include <random>

#include <bats.hpp>

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

TEST_CASE_TEMPLATE("Shapes", F, F2, F3, F5, Q) {
	using VT = SparseVector<F, size_t>;
	using MatT = ColumnMatrix<VT>;

	// TODO: generate some random matrices, set seeds
	MatT I = MatT::identity(5);

	CHECK(I.is_upper());
	CHECK(I.is_lower());
	CHECK(I.is_pivot_matrix());
	CHECK(I.is_EL());
	CHECK(I.is_ELhat());
	CHECK(I.is_EU());
	CHECK(I.is_EUhat());
}

#define CHECK_LEUP(F, A) \
CHECK(F.LEUP_prod() == A); \
CHECK(F.L.is_lower()); \
CHECK(F.U.is_upper()); \
CHECK(F.E.is_EL()); \
CHECK(F.P.is_pivot_matrix());

TEST_CASE_TEMPLATE("LEUP Factorization", T, F2, F3, F5, Q) {
	using VT = SparseVector<T, size_t>;
	using MatT = ColumnMatrix<VT>;

	SUBCASE("Square") {
		for ( unsigned seed = 0; seed < 3; seed++) {
			std::default_random_engine generator(seed);
			auto A = MatT::random(10, 10, 0.2, 1, generator);

			auto F = LEUP(A);
			CHECK_LEUP(F, A)
		}
	}

	SUBCASE("Short") {
		for ( unsigned seed = 0; seed < 3; seed++) {
			std::default_random_engine generator(seed);
			auto A = MatT::random(10, 20, 0.2, 1, generator);

			auto F = LEUP(A);
			CHECK_LEUP(F, A)
		}
	}

	SUBCASE("Tall") {
		for ( unsigned seed = 0; seed < 3; seed++) {
			std::default_random_engine generator(seed);
			auto A = MatT::random(20, 10, 0.2, 1, generator);

			auto F = LEUP(A);
			CHECK_LEUP(F, A)
		}
	}

}

#define CHECK_PUEL(F, A) \
CHECK(F.PUEL_prod() == A); \
CHECK(F.L.is_lower()); \
CHECK(F.U.is_upper()); \
CHECK(F.E.is_ELhat()); \
CHECK(F.P.is_pivot_matrix());

TEST_CASE_TEMPLATE("PUEL Factorization", T, F2, F3, F5, Q) {
	using VT = SparseVector<T, size_t>;
	using MatT = ColumnMatrix<VT>;

	SUBCASE("Square") {
		for ( unsigned seed = 0; seed < 3; seed++) {
			std::default_random_engine generator(seed);
			auto A = MatT::random(10, 10, 0.2, 1, generator);

			auto F = PUEL(A);
			CHECK_PUEL(F, A)
		}
	}

	SUBCASE("Short") {
		for ( unsigned seed = 0; seed < 3; seed++) {
			std::default_random_engine generator(seed);
			auto A = MatT::random(10, 20, 0.2, 1, generator);

			auto F = PUEL(A);
			CHECK_PUEL(F, A)
		}
	}

	SUBCASE("Tall") {
		for ( unsigned seed = 0; seed < 3; seed++) {
			std::default_random_engine generator(seed);
			auto A = MatT::random(20, 10, 0.2, 1, generator);

			auto F = PUEL(A);
			CHECK_PUEL(F, A)
		}
	}

}


#define CHECK_PLEU(F, A) \
CHECK(F.PLEU_prod() == A); \
CHECK(F.L.is_lower()); \
CHECK(F.U.is_upper()); \
CHECK(F.E.is_EU()); \
CHECK(F.P.is_pivot_matrix());

TEST_CASE_TEMPLATE("PLEU Factorization", T, F2, F3, F5, Q) {
	using VT = SparseVector<T, size_t>;
	using MatT = ColumnMatrix<VT>;

	SUBCASE("Square") {
		for ( unsigned seed = 0; seed < 3; seed++) {
			std::default_random_engine generator(seed);
			auto A = MatT::random(10, 10, 0.2, 1, generator);

			auto F = PLEU(A);
			CHECK_PLEU(F, A)
		}
	}

	SUBCASE("Short") {
		for ( unsigned seed = 0; seed < 3; seed++) {
			std::default_random_engine generator(seed);
			auto A = MatT::random(10, 20, 0.2, 1, generator);

			auto F = PLEU(A);
			CHECK_PLEU(F, A)
		}
	}

	SUBCASE("Tall") {
		for ( unsigned seed = 0; seed < 3; seed++) {
			std::default_random_engine generator(seed);
			auto A = MatT::random(20, 10, 0.2, 1, generator);

			auto F = PLEU(A);
			CHECK_PLEU(F, A)
		}
	}

}

#define CHECK_UELP(F, A) \
CHECK(F.UELP_prod() == A); \
CHECK(F.L.is_lower()); \
CHECK(F.U.is_upper()); \
CHECK(F.E.is_EUhat()); \
CHECK(F.P.is_pivot_matrix());

TEST_CASE_TEMPLATE("UELP Factorization", T, F2, F3, F5, Q) {
	using VT = SparseVector<T, size_t>;
	using MatT = ColumnMatrix<VT>;

	SUBCASE("Square") {
		for ( unsigned seed = 0; seed < 3; seed++) {
			std::default_random_engine generator(seed);
			auto A = MatT::random(10, 10, 0.2, 1, generator);

			auto F = UELP(A);
			CHECK_UELP(F, A)
		}
	}

	SUBCASE("Short") {
		for ( unsigned seed = 0; seed < 3; seed++) {
			std::default_random_engine generator(seed);
			auto A = MatT::random(10, 20, 0.2, 1, generator);

			auto F = UELP(A);
			CHECK_UELP(F, A)
		}
	}

	SUBCASE("Tall") {
		for ( unsigned seed = 0; seed < 3; seed++) {
			std::default_random_engine generator(seed);
			auto A = MatT::random(20, 10, 0.2, 1, generator);

			auto F = UELP(A);
			CHECK_UELP(F, A)
		}
	}

}

#define CHECK_LQU(F, A) \
CHECK(F.LQU_prod() == A); \
CHECK(F.L.is_lower()); \
CHECK(F.U.is_upper()); \
CHECK(F.E.is_pivot_matrix()); \


TEST_CASE_TEMPLATE("LQU Factorization", T, F2, F3, F5, Q) {
	using VT = SparseVector<T, size_t>;
	using MatT = ColumnMatrix<VT>;

	SUBCASE("Square") {
		for ( unsigned seed = 0; seed < 3; seed++) {
			std::default_random_engine generator(seed);
			auto A = MatT::random(10, 10, 0.2, 1, generator);

			auto F = LQU(A);
			CHECK_LQU(F, A)
		}
	}

	SUBCASE("Short") {
		for ( unsigned seed = 0; seed < 3; seed++) {
			std::default_random_engine generator(seed);
			auto A = MatT::random(10, 20, 0.2, 1, generator);

			auto F = LQU(A);
			CHECK_LQU(F, A)
		}
	}

	SUBCASE("Tall") {
		for ( unsigned seed = 0; seed < 3; seed++) {
			std::default_random_engine generator(seed);
			auto A = MatT::random(20, 10, 0.2, 1, generator);

			auto F = LQU(A);
			CHECK_LQU(F, A)
		}
	}

}
