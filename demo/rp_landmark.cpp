#include <bats.h>
#include <util/io.h>
#include <string>

#define FT ModP<int, 2>
#define VT SparseVector<FT>
#define MT ColumnMatrix<VT>

int main (int argc, char* argv[]) {

	size_t d = 3; // dimension of Euclidean Space = RP^{d-1}
	size_t n = 200;

	// maximum simplex dimension
    size_t maxdim = parse_argv(argc, argv, "-maxdim", 3); // should be at least d
    double rmax = parse_argv(argc, argv, "-rmax", 0.2);

	//auto X = sample_cube<double>(d, n);
	auto X = sample_sphere<double>(d, n);

	force_repel_rp(X, 0.05);
	for (size_t i = 0; i < 10; i++) {
		force_repel_rp(X, 0.05);
	}


	auto dist = RPAngleDist(); //AngleDist();

	// generate a cover
	// auto L = greedy_landmarks(X, 10, dist);
	// auto cover = landmark_cover(X, L, dist, 3);

	auto F = RipsFiltration(X, dist, rmax, maxdim);
	//auto F = RipsFiltration(X, cover, dist, rmax, maxdim);

	//auto F = FlagFiltration(edges, ts, 3, 2, 0.);

	auto FC = FilteredChainComplex<double, MT>(F);

	auto RFC = ReducedFilteredChainComplex(FC);

	// persistence pairs for H1
	auto ps = RFC.persistence_pairs(2);

	for (auto p : ps) {
		if (p.death > p.birth)
			std::cout << p.str() << " " << p.death - p.birth << std::endl;
	}

	return 0;
}
