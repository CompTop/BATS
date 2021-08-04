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
		spx = {0,1}; X.add(spx);
		spx = {0,2}; X.add(spx);
		spx = {1,2}; X.add(spx);

		auto f = bats::SimplicialMap(X, X);
		auto F = bats::ChainMap<MT>(f);

		std::cout << "\nHomology-type\n";
		auto DG = bats::DGVectorSpace<MT>(X, -1, false);
		for (ssize_t k = 0; k < DG.maxdim()+2; k++) {
			std::cout << "k = " << k << std::endl;
			DG[k].print();
		}
		std::cout << "\nReduced\n";
		bats::ReducedDGVectorSpace RDG(DG);
		RDG.print_summary();

		std::cout << "\nCohomology-type\n";
		auto DG2 = bats::DGVectorSpace<MT>(X, +1, false);
		for (ssize_t k = -1; k < DG2.maxdim()+1; k++) {
			std::cout << "k = " << k << std::endl;
			DG2[k].print();
		}
		std::cout << "\nReduced\n";
		bats::ReducedDGVectorSpace RDG2(DG2);
		RDG2.print_summary();


		std::cout <<"\nMap\n";
		bats::DGLinearMap<MT> DF(f, +1);
		for (ssize_t k = 0; k < DF.maxdim()+1; ++k) {
			std::cout << "k = " << k << std::endl;
			DF[k].print();
		}
		for (ssize_t k = 0; k < DF.maxdim()+1; ++k) {
			auto DFtil = bats::induced_map(DF, RDG, RDG, k);
			std::cout << "induced map: " << k  << std::endl;
			DFtil.print();
		}


		std::cout <<"\nDual Map\n";
		bats::DGLinearMap<MT> DF2(f, -1);
		for (ssize_t k = 0; k < DF2.maxdim()+1; ++k) {
			std::cout << "k = " << k << std::endl;
			DF2[k].print();
		}
		for (ssize_t k = 0; k < DF.maxdim()+1; ++k) {
			auto DFtil = bats::induced_map(DF2, RDG2, RDG2, k);
			std::cout << "induced map: " << k  << std::endl;
			DFtil.print();
		}



	}



    return EXIT_SUCCESS;

}
