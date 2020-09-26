#include <vector>
#include <set>
#include <string>
#include <random>

#include <bats.hpp>
#include <util/set.h>

#include <omp.h>

#define FT ModP<int, 2>
// #define FT Rational<int>
#define VT SparseVector<FT>
#define MT ColumnMatrix<VT>

// dummy struct for node type
struct VectorSpace {};

// seed for deterministic number generation
#define SEED 0

// persistence-type quiver
Diagram<VectorSpace, MT> random_quiver(
	size_t n,
	size_t d,
	double p=0.5
) {
	std::default_random_engine generator(SEED);

	Diagram<VectorSpace, MT> Q(n, n-1);
	for (size_t i = 0; i < n-1; i++) {
		// persistence type quiver
		// edge i : node i -> node i+1
		// Q.set_edge(i, i+1, i, MT::random(d, d, 0.5, 1, generator));
		Q.set_edge(i, i+1, i, MT::random(d, d, p, 1, generator));
	}

	return Q;
}

template <typename T>
void time_barcode_dq(T& Q) {
	std::cout << "divide and conquer (";
	auto nthread = omp_get_max_threads();
	std::cout << nthread << " threads)" << std::endl;
	auto start = omp_get_wtime();
	auto ps = barcode_sparse_divide_conquer(Q, 0);
	auto end = omp_get_wtime();
	std::cout << "\ttime elapsed: " << end - start << "sec. " << std::endl;
}

template <typename T>
void time_barcode_seq(T& Q) {
	std::cout << "sequential" << std::endl;
	auto start = omp_get_wtime();
	auto ps = barcode_sparse_rightleft(Q, 0);
	auto end = omp_get_wtime();
	std::cout << "\ttime elapsed: " << end - start << "sec. " << std::endl;
}

int main() {

	for (size_t dim: {10, 100, 200})
	{
		size_t n = 64;
		std::cout << "\nn= " << n << " dim = " << dim << std::endl;
		auto Q = random_quiver(n, dim);

		time_barcode_seq(Q);

		omp_set_num_threads(1);
		time_barcode_dq(Q);
		omp_set_num_threads(2);
		time_barcode_dq(Q);
		omp_set_num_threads(4);
		time_barcode_dq(Q);
		omp_set_num_threads(8);
		time_barcode_dq(Q);
		omp_set_num_threads(16);
		time_barcode_dq(Q);
		omp_set_num_threads(24);
		time_barcode_dq(Q);
		omp_set_num_threads(32);
		time_barcode_dq(Q);
	}

	for (size_t n: {8, 16, 32, 64, 128})
	{
		size_t dim = 100;
		std::cout << "\nn= " << n << " dim = " << dim << std::endl;
		auto Q = random_quiver(n, dim);

		time_barcode_seq(Q);

		omp_set_num_threads(1);
		time_barcode_dq(Q);
		omp_set_num_threads(2);
		time_barcode_dq(Q);
		omp_set_num_threads(4);
		time_barcode_dq(Q);
		omp_set_num_threads(8);
		time_barcode_dq(Q);
		omp_set_num_threads(16);
		time_barcode_dq(Q);
		omp_set_num_threads(24);
		time_barcode_dq(Q);
		omp_set_num_threads(32);
		time_barcode_dq(Q);
	}


	return 0;
}
