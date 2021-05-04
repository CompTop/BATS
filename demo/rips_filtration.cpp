#include <bats.hpp>
#include <util/io.hpp>
#include <string>

#define FT ModP<int, 2>
#define VT SparseVector<FT>
#define MT ColumnMatrix<VT>

using CpxT = bats::LightSimplicialComplex<size_t, std::unordered_map<size_t, size_t>>;
// using CpxT = bats::SimplicialComplex;

int main (int argc, char* argv[]) {

	size_t d = 2; // dimension of Euclidean Space
	size_t n = 250; // number of points to sample

	// maximum simplex dimension
    size_t maxdim = bats::util::io::parse_argv(argc, argv, "-maxdim", 3);
    double rmax = bats::util::io::parse_argv(argc, argv, "-rmax", 0.4);

	//auto X = sample_cube<double>(d, n);
	auto X = bats::sample_sphere<double>(d, n);

	// auto dist = RPAngleDist(); //AngleDist();
	auto dist = bats::Euclidean();

	{
		X.data.print();
		auto R = bats::RipsComplex<CpxT>(X, dist, rmax, maxdim);
		// auto C = bats::__ChainComplex(R, FT());
		bats::ChainComplex<MT> C(R);
		auto RC = bats::ReducedChainComplex(C);
		std::cout << "non-filtered homology: " << RC.hdim(1) << std::endl;
	}

	// generate a cover
	// auto L = greedy_landmarks(X, 10, dist);
	// auto cover = landmark_cover(X, L, dist, 3);

	auto F = bats::RipsFiltration<CpxT>(X, dist, rmax, maxdim);
	for (size_t i = 0; i <= F.maxdim(); i++) {
		std::cout << F.ncells(i) << " in dim " << i << std::endl;
	}
	//auto F = RipsFiltration(X, cover, dist, rmax, maxdim);

	//auto F = FlagFiltration(edges, ts, 3, 2, 0.);

	{
		auto FC = bats::__FilteredChainComplex(F, FT());

		auto RFC = bats::ReducedFilteredChainComplex(
			FC
			// bats::extra_reduction_flag(),
			// bats::compression_flag()
		);
		// auto RFC = __ReducedFilteredChainComplex(F, FT());

		std::cout << "hdim(1) = " << RFC.RC.hdim(1) << std::endl;

		// persistence pairs for H1
		auto ps = RFC.persistence_pairs(1);

		for (auto p : ps) {
			if (p.death > p.birth)
				std::cout << p.str() << " " << p.death - p.birth << std::endl;
		}
	}
	{
		std::cout << "\n\nauto reduction\n";
		std::cout << F.complex().ncells(0) << ", " << F.ncells(0) << std::endl;
		auto RFC = bats::__ReducedFilteredChainComplex(F, FT(), bats::extra_reduction_flag());

		std::cout << "hdim(1) = " << RFC.RC.hdim(1) << std::endl;

		// persistence pairs for H1
		auto ps = RFC.persistence_pairs(1);

		for (auto p : ps) {
			if (p.death > p.birth)
				std::cout << p.str() << " " << p.death - p.birth << std::endl;
		}
	}

	return 0;
}
