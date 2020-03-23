#include <bats.h>
#include <util/io.h>
#include <string>

#define FT ModP<int, 3>
#define VT SparseVector<FT>
#define MT ColumnMatrix<VT>

int main (int argc, char* argv[]) {

	size_t d = 2; // dimension of Euclidean Space
	size_t n = 100;

	// maximum simplex dimension
    size_t maxdim = parse_argv(argc, argv, "-maxdim", 3);
    double rmax = parse_argv(argc, argv, "-rmax", 0.5);

	auto X = sample_cube<double>(d, n);

	auto dist = Euclidean();

	// generate a cover
	auto L = greedy_landmarks(X, 10, dist);

	// auto F = RipsFiltration(X, dist, rmax, maxdim);
	auto F = WitnessFiltration(L, X, dist, rmax, maxdim);

	auto FC = FilteredChainComplex<double, MT>(F);

	auto RFC = ReducedFilteredChainComplex(FC);

	// persistence pairs for H1
	auto ps = RFC.persistence_pairs(1);

	for (auto p : ps) {
	    std::cout << p.str() << " " << p.death - p.birth << std::endl;
	}

	return 0;
}
