#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <bats.hpp>
using namespace bats;

#define F2 ModP<int, 2>
#define F3 ModP<int, 3>
#define Q Rational<int>

TEST_CASE_TEMPLATE("Field Operations", F, F2, F3, Q, ModP<int, 5>) {
	F aa;
	F a0 = 0;
	F a1 = 1;
	F a2 = 2;
	aa = a0;

	SUBCASE("Check equalities") {
		REQUIRE(aa == a0);
		REQUIRE(a0 == 0);
		REQUIRE(a1 == 1);
		REQUIRE(a2 == 2);
	}

	SUBCASE("Check inequalities") {
		REQUIRE(a1 != a0);
		REQUIRE(a0 != 1);
		REQUIRE(a1 != 0);
	}

	SUBCASE("Unary Operations") {
		REQUIRE(-a0 == a0);
		REQUIRE(a1.inv() == a1);
	}

	SUBCASE("Binary Operations") {
		REQUIRE(a0 + a0 == 0);
		REQUIRE(a1 + a0 == 1);
		REQUIRE(a1 - a1 == 0);
		REQUIRE(a2 - a1 == a1);
		REQUIRE(a1 * a0 == 0);
		REQUIRE(a1 * a1 == 1);
		REQUIRE(a1 / a1 == a1);
		REQUIRE(a2 / a1 == a2);
		REQUIRE(a2 * a2 == 4);
		REQUIRE(a1 * a1.inv() == a1);
	}

	SUBCASE("Update operations") {
		aa = a0;
		aa += a1;
		REQUIRE(aa == 1);
		aa += a1;
		REQUIRE(aa == 2);
		aa -= a1;
		REQUIRE(aa == 1);
		aa *= a2;
		REQUIRE(aa == 2);
		aa /= a1;
		REQUIRE(aa == 2);
		aa -= a1;
		REQUIRE(aa == 1);
	}
}

TEST_CASE("F2 Operations") {
	F2 a0 = 0;
	F2 a1 = 1;
	F2 aa;

	SUBCASE("Basic Operations") {
		REQUIRE(a1 + a1 == 0);
		REQUIRE(a0 + a1 == a1);
		REQUIRE(a0 + a0 == a0);
		REQUIRE(a1 - a1 == a0);
		REQUIRE(a1 * a1 == a1);
		REQUIRE(a1 * a0 == a0);
		REQUIRE(a1 / a1 == a1);
		REQUIRE(a1.inv() == a1);
	}

	SUBCASE("Update Operations") {
		aa = a1;
		aa += a0;
		REQUIRE(aa == 1);
		aa += a1;
		REQUIRE(aa == 0);
		aa -= a1;
		REQUIRE(aa == 1);
		aa *= a1;
		REQUIRE(aa == 1);
		aa /= a1;
		REQUIRE(aa == 1);
	}

}

TEST_CASE("F3 Operations") {
	F3 a0 = 0;
	F3 a1 = 1;
	F3 a2 = 2;
	F3 aa;

	SUBCASE("Basic Operations") {
		REQUIRE(a1 + a1 == a2);
		REQUIRE(a0 + a1 == a1);
		REQUIRE(a0 + a0 == a0);
		REQUIRE(a1 - a1 == a0);
		REQUIRE(a1 * a2 == a2);
		REQUIRE(a1 * a0 == a0);
		REQUIRE(a2 / a1 == a2);
		REQUIRE(a1.inv() == a1);
		REQUIRE(a2.inv() == a2);
		REQUIRE(-a1 == a2);
		REQUIRE(!(a1 == 0));
	}

	SUBCASE("Update Operations") {
		aa = a1;
		aa += a0;
		REQUIRE(aa == 1);
		aa += a1;
		REQUIRE(aa == 2);
		aa -= a1;
		REQUIRE(aa == 1);
		aa *= a1;
		REQUIRE(aa == 1);
		aa /= a1;
		REQUIRE(aa == 1);
	}

}

TEST_CASE("Rational Operations") {
	Q a0 = 0;
	Q a1 = 1;
	Q a2 = 2;
	Q aa;

	SUBCASE("Basic Operations") {
		REQUIRE(a1 + a1 == a2);
		REQUIRE(a0 + a1 == a1);
		REQUIRE(a0 + a0 == a0);
		REQUIRE(a1 - a1 == a0);
		REQUIRE(a1 * a2 == a2);
		REQUIRE(a1 * a0 == a0);
		REQUIRE(a2 / a1 == a2);
		REQUIRE(a1.inv() == a1);
		REQUIRE(a2.inv() == Q(1,2));
		REQUIRE(-a1 == Q(-1));
	}

	SUBCASE("Update Operations") {
		aa = a1;
		aa += a0;
		REQUIRE(aa == 1);
		aa += a1;
		REQUIRE(aa == 2);
		aa -= a1;
		REQUIRE(aa == 1);
		aa *= a1;
		REQUIRE(aa == 1);
		aa /= a1;
		REQUIRE(aa == 1);
	}

}
