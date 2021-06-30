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
