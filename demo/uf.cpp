#include <bats.hpp>
#include <iostream>
#include <vector>

using namespace bats;

#define FT ModP<int, 3>
#define VT SparseVector<FT>
#define MT ColumnMatrix<VT>

int main() {

	std::vector<size_t> ind = {0,1,2,3,4,5,6};
	std::vector<double> val = {1.0, 2.0, 3.0, 4,5,6};

	// persistence pairs for H1
	auto ps = rips_union_find_pairs2(ind, val);

	for (auto p : ps) {
		std::cout << p.str() << std::endl;
	}

	return 0;

}
