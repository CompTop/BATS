#include <iostream>
#include <vector>
#include <util/io.hpp>
#include <chrono>

#include <bats.hpp>

#define FT ModP<int, 3>
#define VT SparseVector<FT>
#define MT ColumnMatrix<VT>

int main(int argc, char* argv[]) {

	if (false) {
	    std::vector<size_t> spx;
		bats::SimplicialComplex X;
		spx = {0}; X.add(spx);
		spx = {1}; X.add(spx);
		spx = {2}; X.add(spx);
		spx = {0,1}; X.add(spx);
		spx = {0,2}; X.add(spx);
		spx = {1,2}; X.add(spx);
	}
	size_t d = 3; // dimension of Euclidean Space
    size_t n = 1000;

    // maximum simplex dimension
    size_t maxdim = bats::util::io::parse_argv(argc, argv, "-maxdim", 3);
    double rmax = bats::util::io::parse_argv(argc, argv, "-rmax", 0.25);

    bats::DataSet<double> x = bats::sample_sphere<double>(d, n);


    auto X = bats::RipsComplex(x, bats::LInfDist(), rmax, maxdim);

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


    return 0;

}
