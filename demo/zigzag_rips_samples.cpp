#include <vector>
#include <set>
#include <string>

#include <bats.hpp>
#include <util/set.h>

#define FT ModP<int, 2>
// #define FT Rational<int>
#define VT SparseVector<FT>
#define MT ColumnMatrix<VT>

int main() {

	size_t nsets = 16;
	size_t ns = 100; // number of points in each subset

	// sample circle
	size_t d = 2; // 2d points
	size_t n = 1000; // 100 points
	auto dist = Euclidean(); // Euclidean distance
	double rmax = 0.5;

	// generate data
	auto X = sample_sphere<double>(d, n);

	// landmark sets
	std::vector<std::set<size_t>> subset;
	for (size_t i =0; i < nsets; i++) {
		subset.emplace_back(random_subset(n, ns));
	}

	// Diagram of Sets and inclusions
	auto SetDgm = linear_subset_union_diagram(subset);

	// Diagram of Spaces and maps
	auto TopDgm = Rips(SetDgm, X, dist, rmax, d);

	// std::string dname = "rss.dgm";
	// TopDgm.save(dname);

	// diagram in Chain
	auto ChainDgm = Chain<MT>(TopDgm);

	// diagram in Homology
	auto HkDgm = Hom(ChainDgm, 1);

	// for (auto M : HkDgm.edata) {
	//     M.print();
	// }

	auto ps = barcode_sparse_divide_conquer(HkDgm, 1);
	for (auto p : ps) {
		std::cout << p.str() << std::endl;
	}

	return 0;
}
