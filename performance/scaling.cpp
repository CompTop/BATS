#include <vector>
#include <set>
#include <string>
#include <random>
#include <fstream> // write to file
#include <sstream>
#include <iomanip>
#include <ctime>

#include <bats.hpp>
// #include <util/set.h>

#include <omp.h>
using namespace bats;

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

	// generate file with timestamp
	auto t = std::time(nullptr);
	auto tm = *std::localtime(&t);
	std::ostringstream oss1;
	oss1 << "scale_n_" << std::put_time(&tm, "%Y_%m_%d_%H") << ".csv";
	std::cout << "writing to " << oss1.str() << std::endl;
	std::ofstream outfile1(oss1.str());
	// write header
	outfile1 << "n, dim, threads, time\n";


	for (size_t dim: {128, 256})
	{
		for (size_t n: {4, 8, 16, 32, 48, 64, 96, 128, 144, 192, 240, 256, 288, 512, 768, 1024}) {
			std::cout << "\nn= " << n << " dim = " << dim << std::endl;
			auto Q = random_quiver(n, dim, MT2(), 0.5);

			// 0 threads = sequential
			outfile1 << n << ", " << dim << ", " << 0 << ", ";
			time_barcode_seq(Q, outfile1);

			auto nthread = omp_get_max_threads();
			outfile1 << n << ", " << dim << ", " << nthread << ", ";
			time_barcode_dq(Q, outfile1);
		}

	}

	std::ostringstream oss2;
	oss2 << "scale_d_" << std::put_time(&tm, "%Y_%m_%d_%H") << ".csv";
	std::cout << "writing to " << oss2.str() << std::endl;
	std::ofstream outfile2(oss2.str());
	// write header
	outfile2 << "n, dim, threads, time\n";

	for (size_t n: {256, 512})
	{
		for (size_t dim: {2, 4, 8, 16, 32, 64, 128, 256}) {
			std::cout << "\nn= " << n << " dim = " << dim << std::endl;
			auto Q = random_quiver(n, dim, MT2(), 0.5);

			// 0 threads = sequential
			outfile2 << n << ", " << dim << ", " << 0 << ", ";
			time_barcode_seq(Q, outfile2);

			auto nthread = omp_get_max_threads();
			outfile2 << n << ", " << dim << ", " << nthread << ", ";
			time_barcode_dq(Q, outfile2);
		}
	}


	return EXIT_SUCCESS;
}
