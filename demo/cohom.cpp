#include <iostream>
#include <vector>
#include <util/io.hpp>
#include <chrono>

#include <bats.hpp>

using FT=ModP<int, 3>;
using VT=SparseVector<FT>;
using MT=ColumnMatrix<VT>;

using CpxT = bats::LightSimplicialComplex<size_t, std::unordered_map<size_t, size_t>>;
// using CpxT = bats::SimplicialComplex;

int main(int argc, char* argv[]) {

	if constexpr (false) {
	    std::vector<size_t> spx;
		bats::SimplicialComplex X;
		spx = {0}; X.add(spx);
		spx = {1}; X.add(spx);
		spx = {2}; X.add(spx);
		spx = {0,1}; X.add(spx);
		spx = {0,2}; X.add(spx);
		spx = {1,2}; X.add(spx);

		auto f = bats::SimplicialMap(X, X);
		auto F = bats::ChainMap<MT>(f);

		auto CX = bats::ChainComplex<MT>(X);
		auto RX = bats::ReducedChainComplex(CX);

		for (size_t k = 0; k < 2; ++k) {
			std::cout << "Homology induced map in dim " << k << std::endl;
			auto FHk = induced_map(F, RX, RX, k);
			FHk.print();
		}

		auto CCX = bats::CochainComplex<MT>(X);
		auto RCCX = bats::ReducedCochainComplex(CCX);

		for (size_t k = 0; k < 2; ++k) {
			std::cout << "Cohomology induced map in dim " << k << std::endl;
			auto FHk = induced_map(F, RCCX, RCCX, k);
			FHk.print();
		}


	}
	if constexpr (true) {
		size_t d = 2; // dimension of Euclidean Space
	    size_t n = 100;

	    // maximum simplex dimension
	    size_t maxdim = bats::util::io::parse_argv(argc, argv, "-maxdim", 2);
	    double rmax = bats::util::io::parse_argv(argc, argv, "-rmax", 0.25);

	    bats::DataSet<double> x = bats::sample_sphere<double>(d, n);


	    auto X = bats::RipsComplex<CpxT>(x, bats::LInfDist(), rmax, maxdim);
		// for (size_t k = 0; k < X.maxdim() + 1; ++k) {
		// 	std::cout << "dim " << k << ": " X.ncells(k)
		// }
		X.print_summary();

		{
			auto start = std::chrono::steady_clock::now();
			auto CX = bats::ChainComplex<MT>(X);
			// for (auto k = 0; k < CX.maxdim()+1; k++){
			// 	CX[k].print();
			// }

			auto RCX = bats::ReducedChainComplex(CX,
				bats::extra_reduction_flag(),
				bats::clearing_flag());
			auto end = std::chrono::steady_clock::now();
			for (size_t k = 0; k < RCX.maxdim() + 1; k++){
				std::cout << RCX.hdim(k) << std::endl;
			}

	        std::cout << "homology reduction: "
	            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
	            << "ms" << std::endl;
		}

		{
			auto start = std::chrono::steady_clock::now();
		    auto CCX = bats::CochainComplex<MT>(X);
			std::cout << "formed cochain complex" << std::endl;
			// for (auto k = 0; k < CX.maxdim()+1; k++){
			// 	CCX[k].print();
			// }

			auto RCCX = bats::ReducedCochainComplex(CCX,
				bats::extra_reduction_flag(),
				bats::clearing_flag());
			std::cout << "formed reduced cochain complex" << std::endl;
			auto end = std::chrono::steady_clock::now();
			for (size_t k = 0; k < RCCX.maxdim() + 1; k++){
				std::cout << RCCX.hdim(k) << std::endl;
			}

	        std::cout << "cohomology reduction: "
	            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
	            << "ms" << std::endl;
		}

		{
			auto start = std::chrono::steady_clock::now();
			auto CX = bats::ChainComplex<MT>(X);
			auto RCX = bats::ReducedChainComplex(CX,
				bats::extra_reduction_flag(),
				bats::clearing_flag());
			auto end = std::chrono::steady_clock::now();

	        std::cout << "homology reduction: "
	            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
	            << "ms" << std::endl;
		}
		{
			auto start = std::chrono::steady_clock::now();
			auto CCX = bats::CochainComplex<MT>(X);
			auto RCCX = bats::ReducedCochainComplex(CCX,
				bats::extra_reduction_flag(),
				bats::clearing_flag());
			auto end = std::chrono::steady_clock::now();

			std::cout << "cohomology reduction: "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
				<< "ms" << std::endl;
		}
	}


    return EXIT_SUCCESS;

}
