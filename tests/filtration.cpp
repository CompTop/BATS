#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include <vector>
#include <set>
#include <random>

#include <bats.hpp>
using namespace bats;

TEST_CASE("Filtered Simplicial Complex") {

    // first build a filtration
    bats::Filtration<double, bats::SimplicialComplex> F;

    std::vector<size_t> s;
    for (size_t i = 0; i < 3; i++) {
        s = {i}; F.add(0., s);
    }
    REQUIRE(F.maxdim() == 0);
    REQUIRE(F.ncells(0) == 3);

    s = {0, 1}; F.add(1., s);
    s = {2, 0}; F.add(1., s);
    s = {1, 2}; F.add(1., s);
    REQUIRE(F.maxdim() == 1);
    REQUIRE(F.ncells(1) == 3);

    // check that we use stable sort perms
    std::vector<size_t> perm = {0,1,2};
    REQUIRE(F.sortperm(0) == perm);
    REQUIRE(F.sortperm(1) == perm);





}
