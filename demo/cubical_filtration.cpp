#include <iostream>
#include <vector>

#include <bats.hpp>

#define FT ModP<int, 3>
#define VT SparseVector<FT>
#define MT ColumnMatrix<VT>

int main() {

	std::vector<size_t> cbx;
	bats::Filtration<double, bats::CubicalComplex> F(2); // 2-dimensional cubical cpx
	// spx = {0}; F.add(0.0, spx);
	// spx = {1}; F.add(0.0, spx);
	cbx = {0,0,0,1}; F.add_recursive(1.0, cbx);

	auto FC = bats::FilteredChainComplex<double, MT>(F);

	auto RFC = bats::ReducedFilteredChainComplex(FC);

	// persistence pairs for H1
	auto ps = RFC.persistence_pairs(0);

	for (auto p : ps) {
		std::cout << p.str() << std::endl;
	}

	return 0;
}
