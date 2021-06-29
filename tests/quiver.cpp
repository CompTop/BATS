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

#define N_SEEDS 4

// dummy struct for node type
struct VectorSpace {
	size_t dim;
};

TEST_CASE_TEMPLATE("Square", T, F2, F3, F5) {
	using VT = SparseVector<T, size_t>;
	using MatT = ColumnMatrix<VT>;

	SUBCASE("-->") {

		size_t n = 64;
		size_t d = 10;
		for ( unsigned seed = 0; seed < N_SEEDS; seed++) {
			std::default_random_engine generator(seed);

			// create quiver
			bats::Diagram<VectorSpace, MatT> A(n, n-1);
			for (size_t i = 0; i < n-1; i++) {
				// persistence type quiver
				// edge i : node i -> node i+1
				A.set_edge(i, i, i+1, MatT::random(d, d, 0.5, 1, generator));
			}

			// compute barcode
			auto ps1 = bats::barcode_sparse_rightleft(A, 0);
			auto ps2 = bats::barcode_sparse_leftright(A, 0);
			auto ps3 = bats::barcode_sparse_divide_conquer(A, 0);

			CHECK(bats::barcode_equality(ps1, ps2));
			CHECK(bats::barcode_equality(ps1, ps3));
			CHECK(bats::barcode_equality(ps2, ps3));

		}
	}

	SUBCASE("<--") {

		size_t n = 64;
		size_t d = 10;
		for ( unsigned seed = 0; seed < N_SEEDS; seed++) {
			std::default_random_engine generator(seed);

			// create quiver
			bats::Diagram<VectorSpace, MatT> A(n, n-1);
			for (size_t i = 0; i < n-1; i++) {
				// persistence type quiver
				// edge i : node i+1 -> node i
				A.set_edge(i, i+1, i, MatT::random(d, d, 0.5, 1, generator));
			}

			// compute barcode
			auto ps1 = bats::barcode_sparse_rightleft(A, 0);
			auto ps2 = bats::barcode_sparse_leftright(A, 0);
			auto ps3 = bats::barcode_sparse_divide_conquer(A, 0);

			CHECK(bats::barcode_equality(ps1, ps2));
			CHECK(bats::barcode_equality(ps1, ps3));
			CHECK(bats::barcode_equality(ps2, ps3));

		}
	}

	SUBCASE("zigzag") {

		size_t n = 64;
		size_t d = 10;
		for ( unsigned seed = 0; seed < N_SEEDS; seed++) {
			std::default_random_engine generator(seed);

			// create quiver
			bats::Diagram<VectorSpace, MatT> A(n, n-1);
			for (size_t i = 0; i < n-1; i++) {
				// zigzag type quiver
				if ((i & 0x1) == 1) {
					// <--
					A.set_edge(i, i+1, i, MatT::random(d, d, 0.5, 1, generator));
				} else {
					// -->
					A.set_edge(i, i, i+1, MatT::random(d, d, 0.5, 1, generator));
				}


			}

			// compute barcode
			auto ps1 = bats::barcode_sparse_rightleft(A, 0);
			auto ps2 = bats::barcode_sparse_leftright(A, 0);
			auto ps3 = bats::barcode_sparse_divide_conquer(A, 0);

			CHECK(bats::barcode_equality(ps1, ps2));
			CHECK(bats::barcode_equality(ps1, ps3));
			CHECK(bats::barcode_equality(ps2, ps3));

		}
	}


}
