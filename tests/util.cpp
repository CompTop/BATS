
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include <vector>

#include <util/sorted.hpp>

TEST_CASE("sorting") {
	std::vector<int> a = {1,2,3};
  	std::vector<int> b = {1, 3, 5};
  	std::vector<int> c, targ;
  	bats::util::intersect_sorted(a, b, c);
	targ = {1,3};
	REQUIRE(c == targ);


  	bats::util::intersect_sorted_lt(a, b, 3, c);
	targ = {1};
	REQUIRE(c == targ);

}
