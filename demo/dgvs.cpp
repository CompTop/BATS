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

int main() {

	if constexpr (true) {
	    std::vector<size_t> spx;

		bats::SimplicialComplex X;
		spx = {0}; X.add(spx);
		spx = {1}; X.add(spx);
		spx = {2}; X.add(spx);
		spx = {3}; X.add(spx);
		spx = {0,1}; X.add(spx);
		spx = {1,2}; X.add(spx);
		spx = {2,3}; X.add(spx);
		spx = {0,3}; X.add(spx);

		bats::SimplicialComplex Y;
		spx = {0}; Y.add(spx);
		spx = {1}; Y.add(spx);
		spx = {2}; Y.add(spx);
		spx = {0,1}; Y.add(spx);
		spx = {0,2}; Y.add(spx);
		spx = {1,2}; Y.add(spx);

		// map from X -> Y collapses one edge
		std::vector<size_t> f0 = {0,1,2,2};
		// map from X -> Y collapses entire complex
		// std::vector<size_t> f0 = {0,0,0,0};
		auto f = bats::SimplicialMap(X, Y, f0);
		// auto F = bats::ChainMap<MT>(f);

		std::cout << "\nHomology-type\n";
		auto DGX = bats::DGVectorSpace<MT>(X, -1, false);
		std::cout << "\nReduced X\n";
		bats::ReducedDGVectorSpace RDGX(DGX);
		RDGX.print_summary();

		auto DGY = bats::DGVectorSpace<MT>(Y, -1, false);
		std::cout << "\nReduced Y\n";
		bats::ReducedDGVectorSpace RDGY(DGY);
		RDGY.print_summary();

		std::cout <<"\nMap\n";
		bats::DGLinearMap<MT> DF(f, -1);
		for (ssize_t k = 0; k < DF.maxdim()+1; ++k) {
			std::cout << "k = " << k << std::endl;
			DF[k].print();
		}
		for (ssize_t k = 0; k < DF.maxdim()+1; ++k) {
			auto DFtil = bats::induced_map(DF, RDGX, RDGY, k);
			std::cout << "induced map: " << k  << std::endl;
			DFtil.print();
		}

		std::cout << "\nCohomology-type\n";
		auto DGX2 = bats::DGVectorSpace<MT>(X, +1, false);
		std::cout << "\nReduced X\n";
		bats::ReducedDGVectorSpace RDGX2(DGX2);
		RDGX2.print_summary();

		auto DGY2 = bats::DGVectorSpace<MT>(Y, +1, false);
		std::cout << "\nReduced Y\n";
		bats::ReducedDGVectorSpace RDGY2(DGY2);
		RDGY2.print_summary();

		std::cout <<"\nMap\n";
		bats::DGLinearMap<MT> DF2(f, +1);
		for (ssize_t k = 0; k < DF2.maxdim()+1; ++k) {
			std::cout << "k = " << k << std::endl;
			DF2[k].print();
		}
		for (ssize_t k = 0; k < DF2.maxdim()+1; ++k) {
			auto DFtil = bats::induced_map(DF2, RDGY, RDGX, k);
			std::cout << "induced map: " << k  << std::endl;
			DFtil.print();
		}

		// std::cout <<"\nMap\n";
		// bats::DGLinearMap<MT> DF(f, +1);
		// for (ssize_t k = 0; k < DF.maxdim()+1; ++k) {
		// 	std::cout << "k = " << k << std::endl;
		// 	DF[k].print();
		// }
		// for (ssize_t k = 0; k < DF.maxdim()+1; ++k) {
		// 	auto DFtil = bats::induced_map(DF, RDG, RDG, k);
		// 	std::cout << "induced map: " << k  << std::endl;
		// 	DFtil.print();
		// }
		//
		//
		// std::cout <<"\nDual Map\n";
		// bats::DGLinearMap<MT> DF2(f, -1);
		// for (ssize_t k = 0; k < DF2.maxdim()+1; ++k) {
		// 	std::cout << "k = " << k << std::endl;
		// 	DF2[k].print();
		// }
		// for (ssize_t k = 0; k < DF.maxdim()+1; ++k) {
		// 	auto DFtil = bats::induced_map(DF2, RDG2, RDG2, k);
		// 	std::cout << "induced map: " << k  << std::endl;
		// 	DFtil.print();
		// }



	}



    return EXIT_SUCCESS;

}
