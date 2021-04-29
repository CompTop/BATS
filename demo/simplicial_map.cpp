#include <iostream>
#include <vector>
#include <bats.hpp>

int main() {

	std::cout << "\nCreating simplicial complex" << std::endl;

	// SimplicialComplex X(1);
	bats::LightSimplicialComplex<size_t> X(2, 1);
	X.add({0});
	X.add({1});
	X.add({0,1});

	// for a simplicial map, we simply say where the vertices end up
	std::vector<size_t> f0 = {1, 0};

	auto f = bats::SimplicialMap(X, X, f0);

	for (size_t k = 0; k < f.maxdim() + 1; k++) {
		std::cout << "\nmap in dimension " << k << std::endl;
		f[k].print();
	}
	return 0;
}
