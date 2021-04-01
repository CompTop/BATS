#include <iostream>
#include <vector>
#include <util/io.hpp>
#include <chrono>

#include <bats.hpp>

#define FT ModP<int, 3>
#define VT SparseVector<FT>
#define MT ColumnMatrix<VT>

int main() {

	std::cout << "standard version:" << std::endl;
	{
		std::vector<size_t> spx;
		bats::SimplicialComplex X;
		spx = {0}; X.add(spx);
		spx = {1}; X.add(spx);
		spx = {2}; X.add(spx);
		spx = {0,1}; X.add(spx);
		spx = {0,2}; X.add(spx);
		spx = {1,2}; X.add(spx);

		X.print_summary();

		auto [CXCX, F, XX] = bats::EilenbergZilber(X, X, 2, FT());
		XX.print_summary();
		std::cout << "tensor valid complex = " << CXCX.is_valid_complex() << std::endl;

		std::cout << "\ntensor product dimensions: " << std::endl;
		for (size_t k = 0; k <= CXCX.maxdim(); k++) {
			std::cout << "dim " << k << ": " << CXCX.dim(k) << std::endl;
		}

		std::cout << "\ntensor product homology: " << std::endl;
		auto RCXCX = bats::ReducedChainComplex(CXCX);
		for (size_t k = 0; k <= RCXCX.maxdim(); k++) {
			std::cout << "betti " << k << ": " << RCXCX.betti(k) << std::endl;
		}

		std::cout << "\ntriangulated product: " << std::endl;
		auto CXX = bats::__ChainComplex(XX, FT());
		auto RCXX = bats::ReducedChainComplex(CXX);
		for (size_t k = 0; k <= RCXX.maxdim(); k++) {
			std::cout << "betti " << k << ": " << RCXX.betti(k) << std::endl;
		}

		std::cout << "\nEilenberg-Zilber induced maps" << std::endl;
		for (size_t k = 0; k <= F.maxdim(); k++) {
			//F[k].print();
			std::cout << "dimension " << k << std::endl;
			auto Fk = bats::induced_map(F, RCXCX, RCXX, k);
			Fk.print();
		}
	}
	std::cout << "\n\nrelative version:" << std::endl;
	{
		std::vector<size_t> spx;
		bats::SimplicialComplex X, A;
		spx = {0}; X.add(spx); A.add(spx);
		spx = {1}; X.add(spx);
		spx = {2}; X.add(spx); A.add(spx);
		spx = {0,1}; X.add(spx);
		spx = {1,2}; X.add(spx);


		std::cout << "H(X, A): " << std::endl;
		auto CXA = bats::__ChainComplex(X, A, FT());
		auto RXA = bats::ReducedChainComplex(CXA);
		for (size_t k = 0; k <= RXA.maxdim(); k++) {
			std::cout << "betti " << k << ": " << RXA.betti(k) << std::endl;
		}



		auto [RCXCX, RF, RCXX] = bats::EilenbergZilber(X, A, X, A, 2, FT());

		std::cout << "\ntensor valid complex = " << RCXCX.is_valid_complex() << std::endl;


		std::cout << "\ntensor product dimensions: " << std::endl;
		for (size_t k = 0; k <= RCXCX.maxdim(); k++) {
			std::cout << "dim " << k << ": " << RCXCX.dim(k) << std::endl;
		}

		std::cout << "\ntensor product homology: " << std::endl;
		auto RRCXCX = bats::ReducedChainComplex(RCXCX);
		for (size_t k = 0; k <= RRCXCX.maxdim(); k++) {
			std::cout << "betti " << k << ": " << RRCXCX.betti(k) << std::endl;
		}

		std::cout << "\ntriangulated product: " << std::endl;
		auto RRCXX = bats::ReducedChainComplex(RCXX);
		for (size_t k = 0; k <= RRCXX.maxdim(); k++) {
			std::cout << "betti " << k << ": " << RRCXX.betti(k) << std::endl;
		}

		std::cout << "\nEilenberg-Zilber induced maps" << std::endl;
		for (size_t k = 0; k <= RF.maxdim(); k++) {
			//F[k].print();
			std::cout << "dimension " << k << std::endl;
			auto Fk = bats::induced_map(RF, RRCXCX, RRCXX, k);
			Fk.print();
		}

	}




    return 0;

}
