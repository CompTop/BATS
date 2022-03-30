#include <bats.hpp>
#include <util/io.hpp>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include <chrono>
#include <random>

using FT = ModP<int, 2>;
#define VT SparseVector<FT>
#define MT ColumnMatrix<VT>


using CpxT = bats::SimplicialComplex;
using namespace bats;

int main(int argc, char* argv[]) {
    // std::vector<std::vector<std::vector<size_t>>> Splx_3;
	size_t d = 2; // dimension of Euclidean Space
	size_t n = 150; // number of points to sample

	// maximum simplex dimension
    size_t maxdim = bats::util::io::parse_argv(argc, argv, "-maxdim", 3);
    double rmax = bats::util::io::parse_argv(argc, argv, "-rmax", 2);

    auto X = bats::sample_sphere<double>(d, n);
    // auto X_data = X.data;
    // X_data.swap_rows(0,1);
    // auto Y = bats::DataSet(X_data);
    auto Y = bats::sample_sphere<double>(d, n);

	auto dist = bats::Euclidean(); // metric

    auto FX = bats::RipsFiltration<CpxT>(X, dist, 0.5, maxdim);
	auto FY = bats::RipsFiltration<CpxT>(X, dist, 0.7, maxdim);

	auto t0 = std::chrono::steady_clock::now();
	auto U = bats::UpdateInfo2(FX, FY);
	auto t1 = std::chrono::steady_clock::now();
    std::cout << "\tconstruct update info takes "
        << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
        << "ms" << std::endl;

	for (int degree : {-1, +1}) {
		// int degree = -1;
		std::cout << "degree(" << degree << ")" << std::endl;

		t0 = std::chrono::steady_clock::now();
		auto CX = FilteredDGVectorSpace<double, MT>(FX, degree);
        auto RX = ReducedFilteredDGVectorSpace(CX, bats::standard_reduction_flag(),
                bats::clearing_flag(), bats::compute_basis_flag());
		t1 = std::chrono::steady_clock::now();
		std::cout << "\tinitial reduction takes "
	        << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
	        << "ms" << std::endl;

		t0 = std::chrono::steady_clock::now();
		RX.update_basis(U);
		t1 = std::chrono::steady_clock::now();
	    std::cout << "\tcompute update takes "
	        << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
	        << "ms" << std::endl;

		t0 = std::chrono::steady_clock::now();
		auto CY = FilteredDGVectorSpace<double, MT>(FY, degree);
        auto RY = ReducedFilteredDGVectorSpace(CY, bats::standard_reduction_flag(),
                bats::clearing_flag(), bats::compute_basis_flag());
		t1 = std::chrono::steady_clock::now();
		std::cout << "\treduction from scratch takes "
	        << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
	        << "ms" << std::endl;
	}






    return EXIT_SUCCESS;
}
