#include <vector>
#include <set>
#include <string>
#include <random>
#include <fstream> // write to file
#include <sstream>
#include <iomanip>
#include <ctime>

#include <bats.hpp>
#include <omp.h>
using namespace bats;

#define FT ModP<int, 2>
// #define FT Rational<int>
#define VT SparseVector<FT>
#define MT ColumnMatrix<VT>

template <typename Stream>
void run_problem(
	size_t nsets, // number of subsets
	size_t ns,    // number of points
	size_t nthread, //number of threads
	size_t hdim,	// homology dimension
	Stream &f
) {
	// sample circle
	size_t d = 2; // 2d points
	size_t n = 1000; // 100 points
	auto dist = Euclidean(); // Euclidean distance
	double rmax = 0.4;

	omp_set_num_threads(nthread);

	// generate data
	auto X = sample_sphere<double>(d, n);

	auto nthread2 = omp_get_max_threads();
	std::cout << "\nsubsets: " << nsets << "\npoints: " << ns <<  "\nthreads: " << nthread2 << std::endl;
	f << nsets << ", " << ns << ", " << nthread2 << ", ";

	auto start = omp_get_wtime();
		// landmark sets
		std::vector<std::set<size_t>> subset;
		for (size_t i =0; i < nsets; i++) {
			subset.emplace_back(bats::util::random_subset(n, ns));
		}
		// Diagram of Sets and inclusions
		auto SetDgm = linear_subset_union_diagram(subset);

		// Diagram of Spaces and maps
		auto TopDgm = Rips(SetDgm, X, dist, rmax, 2);
	auto end = omp_get_wtime();
	std::cout << "\tproblem setup: " << end - start << "sec. " << std::endl;
	f << end - start << ", ";

	start = omp_get_wtime();
		// diagram in Chain
		auto ChainDgm = Chain(TopDgm, FT());

		// diagram in Homology
		auto HkDgm = Hom(ChainDgm, hdim);
	end = omp_get_wtime();
	std::cout << "\tchain and hom: " << end - start << "sec. " << std::endl;
	f << end - start << ", ";

	start = omp_get_wtime();
		auto ps1 = barcode(HkDgm, hdim, flags::rightward());
	end = omp_get_wtime();
	std::cout << "\tbarcode sequential: " << end - start << "sec. " << std::endl;
	f << end - start << ", ";

	start = omp_get_wtime();
		auto ps2 = barcode(HkDgm, hdim, flags::divide_conquer());
	end = omp_get_wtime();
	std::cout << "\tbarcode divide conquer: " << end - start << "sec. " << std::endl;
	f << end - start << std::endl;

}

int main() {

	size_t nsets = 128;
	size_t ns = 100; // number of points in each subset
	size_t NREP = 10; // number of repetitions per experiment

	// generate file with timestamp
	auto t = std::time(nullptr);
	auto tm = *std::localtime(&t);
	std::ostringstream oss;
	oss << "zigzag_circle_" << std::put_time(&tm, "%Y_%m_%d_%H") << ".csv";
	std::cout << "writing to " << oss.str() << std::endl;
	std::ofstream outfile(oss.str());
	// std::ofstream outfile("zigzag_circle.csv");
	// write header
	outfile << "subsets,points,threads,setup,chain/hom,sequential,divide/conquer\n";


	for (size_t nthread: {1, 2, 4, 8, 16, 24, 32, 48, 64}) {
		for (size_t rep = 0; rep < NREP; ++rep) {
			run_problem(nsets, ns, nthread, 1, outfile);
		}
	}

	return 0;
}
