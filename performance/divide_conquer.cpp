#include <vector>
#include <set>
#include <string>
#include <random>
#include <fstream> // write to file

#include <bats.hpp>
#include <util/set.h>

#include <omp.h>

#define FT2 ModP<int, 2>
#define VT2 SparseVector<FT2>
#define MT2 ColumnMatrix<VT2>

#define FT3 ModP<int, 3>
#define VT3 SparseVector<FT3>
#define MT3 ColumnMatrix<VT3>

// dummy struct for node type
struct VectorSpace {};

// seed for deterministic number generation
#define SEED 0

// persistence-type quiver
template <typename MatT>
Diagram<VectorSpace, MatT> random_quiver(
	size_t n,
	size_t d,
	MatT,
	double p=0.5
) {
	std::default_random_engine generator(SEED);

	Diagram<VectorSpace, MatT> Q(n, n-1);
	for (size_t i = 0; i < n-1; i++) {
		// persistence type quiver
		// edge i : node i -> node i+1
		// Q.set_edge(i, i+1, i, MT::random(d, d, 0.5, 1, generator));
		Q.set_edge(i, i+1, i, MatT::random(d, d, p, 1, generator));
	}

	return Q;
}

template <typename T, typename Stream>
void time_barcode_dq(
	T& Q,
	Stream &f
) {
	std::cout << "divide and conquer (";
	auto nthread = omp_get_max_threads();
	std::cout << nthread << " threads)" << std::endl;
	auto start = omp_get_wtime();
	auto ps = barcode_sparse_divide_conquer(Q, 0);
	auto end = omp_get_wtime();
	std::cout << "\ttime elapsed: " << end - start << "sec. " << std::endl;
	f << end - start << std::endl;
}

template <typename T, typename Stream>
void time_barcode_seq(
	T& Q,
	Stream &f
) {
	std::cout << "sequential" << std::endl;
	auto start = omp_get_wtime();
	auto ps = barcode_sparse_rightleft(Q, 0);
	auto end = omp_get_wtime();
	std::cout << "\ttime elapsed: " << end - start << "sec. " << std::endl;
	f << end - start << std::endl;
}

int main() {

	std::ofstream outfile1("dq1_f2.csv");
	// write header
	outfile1 << "n, dim, threads, time\n";


	for (size_t dim: {10, 100, 200})
	{
		size_t n = 64;
		std::cout << "\nn= " << n << " dim = " << dim << std::endl;
		auto Q = random_quiver(n, dim, MT2(), 0.5);

		// 0 threads = sequential
		outfile1 << n << ", " << dim << ", " << 0 << ", ";
		time_barcode_seq(Q, outfile1);

		for (size_t nthreads: {1, 2, 4, 8, 16, 24, 32}) {
			outfile1 << n << ", " << dim << ", " << nthreads << ", ";
			omp_set_num_threads(nthreads);
			time_barcode_dq(Q, outfile1);
		}
	}

	std::ofstream outfile2("dq2_f2.csv");
	// write header
	outfile2 << "n, dim, threads, time\n";

	for (size_t n: {8, 16, 32, 64, 128})
	{
		size_t dim = 100;
		std::cout << "\nn= " << n << " dim = " << dim << std::endl;
		auto Q = random_quiver(n, dim, MT2(), 0.5);

		// 0 threads = sequential
		outfile2 << n << ", " << dim << ", " << 0 << ", ";
		time_barcode_seq(Q, outfile2);

		for (size_t nthreads: {1, 2, 4, 8, 16, 24, 32}) {
			outfile2 << n << ", " << dim << ", " << nthreads << ", ";
			omp_set_num_threads(nthreads);
			time_barcode_dq(Q, outfile2);
		}
	}

	std::ofstream outfile3("dq1_f3.csv");
	// write header
	outfile3 << "n, dim, threads, time\n";


	for (size_t dim: {10, 100, 200})
	{
		size_t n = 64;
		std::cout << "\nn= " << n << " dim = " << dim << std::endl;
		auto Q = random_quiver(n, dim, MT3(), 0.5);

		// 0 threads = sequential
		outfile3 << n << ", " << dim << ", " << 0 << ", ";
		time_barcode_seq(Q, outfile3);

		for (size_t nthreads: {1, 2, 4, 8, 16, 24, 32}) {
			outfile3 << n << ", " << dim << ", " << nthreads << ", ";
			omp_set_num_threads(nthreads);
			time_barcode_dq(Q, outfile3);
		}
	}

	std::ofstream outfile4("dq2_f3.csv");
	// write header
	outfile4 << "n, dim, threads, time\n";

	for (size_t n: {8, 16, 32, 64, 128})
	{
		size_t dim = 100;
		std::cout << "\nn= " << n << " dim = " << dim << std::endl;
		auto Q = random_quiver(n, dim, MT3(), 0.5);

		// 0 threads = sequential
		outfile4 << n << ", " << dim << ", " << 0 << ", ";
		time_barcode_seq(Q, outfile4);

		for (size_t nthreads: {1, 2, 4, 8, 16, 24, 32}) {
			outfile4 << n << ", " << dim << ", " << nthreads << ", ";
			omp_set_num_threads(nthreads);
			time_barcode_dq(Q, outfile4);
		}
	}


	return 0;
}
