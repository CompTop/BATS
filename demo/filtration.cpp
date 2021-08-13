#include <bats.hpp>
#include <iostream>

using namespace bats;

#define FT ModP<int, 3>
#define VT SparseVector<FT>
#define MT ColumnMatrix<VT>

int main() {

	std::vector<size_t> spx;
	Filtration<double, SimplicialComplex> F;
	spx = {0}; F.add(0.0, spx);
	spx = {1}; F.add(0.0, spx);
	spx = {2}; F.add(0.0, spx);
	spx = {0,1}; F.add(1.0, spx);
	spx = {0,2}; F.add(1.0, spx);
	spx = {1,2}; F.add(1.0, spx);

	auto FC = FilteredChainComplex<double, MT>(F);

	auto RFC = ReducedFilteredChainComplex(FC);

	// persistence pairs for H1
	auto ps = RFC.persistence_pairs(1);

	for (auto p : ps) {
		std::cout << p.str() << std::endl;
	}

	auto X = F.sublevelset(0.5);
	X.print_summary();


	return 0;

}
