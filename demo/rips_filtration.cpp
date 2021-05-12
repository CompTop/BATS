#include <bats.hpp>
#include <util/io.hpp>
#include <string>
#include <chrono>

#define FT ModP<int, 2>
#define VT SparseVector<FT>
#define MT ColumnMatrix<VT>

using CpxT = bats::LightSimplicialComplex<size_t, std::unordered_map<size_t, size_t>>;
// using CpxT = bats::SimplicialComplex;

int main (int argc, char* argv[]) {

	// maximum simplex dimension
	size_t d = bats::util::io::parse_argv(argc, argv, "-dim", 2); // dimension of Euclidean Space
	size_t n = bats::util::io::parse_argv(argc, argv, "-npoints", 250);
    size_t maxdim = bats::util::io::parse_argv(argc, argv, "-maxdim", 2);
    double rmax = bats::util::io::parse_argv(argc, argv, "-rmax", 3.0);

	std::cout << "d = " << d << std::endl
			<< "n = " << n << std::endl
			<< "maxdim = " << maxdim << std::endl
			<< "rmax = " << rmax << "\n\n";

	//auto X = sample_cube<double>(d, n);
	auto X = bats::sample_sphere<double>(d, n);

	// auto dist = RPAngleDist(); //AngleDist();
	auto dist = bats::Euclidean();

	{

		auto start = std::chrono::steady_clock::now();
		auto R = bats::RipsComplex<CpxT>(X, dist, rmax, maxdim);
		auto end = std::chrono::steady_clock::now();
	    std::cout << "Construction of Complex: "
	        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
	        << "ms" << std::endl;
		// auto C = bats::__ChainComplex(R, FT());
		start = std::chrono::steady_clock::now();
		bats::ChainComplex<MT> C(R);
		auto RC = bats::ReducedChainComplex(C);
		end = std::chrono::steady_clock::now();
	    std::cout << "Chain and reduction: "
	        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
	        << "ms" << std::endl;
		RC.print_summary();
		// std::cout << "non-filtered homology: " << RC.hdim(1) << "\n\n";
	}

	// generate a cover
	// auto L = greedy_landmarks(X, 10, dist);
	// auto cover = landmark_cover(X, L, dist, 3);
	{
		auto start = std::chrono::steady_clock::now();
		auto F = bats::RipsFiltration<CpxT>(X, dist, rmax, maxdim);
		auto end = std::chrono::steady_clock::now();
	    std::cout << "Construction of Filtration: "
	        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
	        << "ms" << std::endl;
		for (size_t i = 0; i <= F.maxdim(); i++) {
			std::cout << F.ncells(i) << " in dim " << i << std::endl;
		}

		start = std::chrono::steady_clock::now();
		auto FC = bats::Chain(F, FT());

		auto RFC = bats::ReducedFilteredChainComplex(
			FC,
			bats::extra_reduction_flag(),
			// bats::standard_reduction_flag(),
			// bats::compression_flag(),
			bats::compute_basis_flag()
		);
		end = std::chrono::steady_clock::now();
	    std::cout << "Chain and reduction: "
	        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
	        << "ms" << std::endl;
		// auto RFC = __ReducedFilteredChainComplex(F, FT());

		RFC.print_summary();
		// std::cout << "hdim(1) = " << RFC.RC.hdim(1) << std::endl;

		// // persistence pairs for H1
		// auto ps = RFC.persistence_pairs(1);

		// for (auto p : ps) {
		// 	if (p.death > p.birth)
		// 		std::cout << p.str() << " " << p.death - p.birth << std::endl;
		// }
	}
	// if (false) {
	// 	std::cout << "\n\nauto reduction\n";
	// 	std::cout << F.complex().ncells(0) << ", " << F.ncells(0) << std::endl;
	// 	auto RFC = bats::Reduce(F, FT(), bats::extra_reduction_flag());
	//
	// 	std::cout << "hdim(1) = " << RFC.RC.hdim(1) << std::endl;
	//
	// 	// persistence pairs for H1
	// 	auto ps = RFC.persistence_pairs(1);
	//
	// 	// for (auto p : ps) {
	// 	// 	if (p.death > p.birth)
	// 	// 		std::cout << p.str() << " " << p.death - p.birth << std::endl;
	// 	// }
	// }

	return 0;
}
