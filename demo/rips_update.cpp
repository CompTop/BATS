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
	size_t n = bats::util::io::parse_argv(argc, argv, "-npoints", 100);
    size_t maxdim = bats::util::io::parse_argv(argc, argv, "-maxdim", 2);
    double rmax = bats::util::io::parse_argv(argc, argv, "-rmax", 3.0);

	std::cout << "d = " << d << std::endl
			<< "n = " << n << std::endl
			<< "maxdim = " << maxdim << std::endl
			<< "rmax = " << rmax << "\n\n";

	auto X = bats::sample_cube<double>(d, n, 0);
	// auto X = bats::sample_sphere<double>(d, n);

	// auto dist = RPAngleDist(); //AngleDist();
	auto dist = bats::Euclidean();


	{

		auto D = dist(X); // create distance matrix
		auto r_enc = bats::enclosing_radius(D);
		std::cout << "enclosing radius = " << r_enc << std::endl;

		auto start = std::chrono::steady_clock::now();
		auto F = bats::RipsFiltration<CpxT>(D, r_enc, maxdim);
		auto end = std::chrono::steady_clock::now();
	    std::cout << "Construction of Filtration: "
	        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
	        << "ms" << std::endl;
		for (size_t i = 0; i <= F.maxdim(); i++) {
			std::cout << F.ncells(i) << " in dim " << i << std::endl;
		}

		std::cout << "\nstandard reduction\n";
		start = std::chrono::steady_clock::now();
		auto FC = bats::Chain(F, FT());
		// FC.C.clear_compress_apparent_pairs();
		auto RFC = bats::ReducedFilteredChainComplex(
			FC,
			bats::standard_reduction_flag(),
			bats::compute_basis_flag()
		);
		end = std::chrono::steady_clock::now();
	    std::cout << "reduction: "
	        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
	        << "ms" << std::endl;
		// auto RFC = __ReducedFilteredChainComplex(F, FT());

		RFC.print_summary(true);


		// Now, let's perturb the data set
		bats::add_normal_noise(X, 0, 0.0, 0.05);
		D = dist(X); // create distance matrix
		r_enc = bats::enclosing_radius(D);
		std::cout << "\n\nenclosing radius = " << r_enc << std::endl;

		start = std::chrono::steady_clock::now();
		auto F2 = bats::RipsFiltration<CpxT>(D, r_enc, maxdim);
		end = std::chrono::steady_clock::now();
	    std::cout << "Construction of Filtration: "
	        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
	        << "ms" << std::endl;
		for (size_t i = 0; i <= F2.maxdim(); i++) {
			std::cout << F2.ncells(i) << " in dim " << i << std::endl;
		}

		// Compute update info
		start = std::chrono::steady_clock::now();
		auto uinfo = bats::Update_info(F, F2);
		end = std::chrono::steady_clock::now();
	    std::cout << "\nConstruction of Update info: "
	        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
	        << "ms" << std::endl;

		// Compute update
		start = std::chrono::steady_clock::now();
		RFC.update_filtration_general(uinfo, bats::standard_reduction_flag());
		end = std::chrono::steady_clock::now();
		std::cout << "\nApplication of update: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
			<< "ms" << std::endl;
	}



	return 0;
}
