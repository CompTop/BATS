#include <iostream>
#include <vector>
#include <chrono>

#include <bats.hpp>
#include <util/io.hpp>
#include <string>

#define FT ModP<int, 3>
#define VT SparseVector<FT>
#define MT ColumnMatrix<VT>



using namespace bats;

// using CpxT = SimplicialComplex;
using CpxT = LightSimplicialComplex<size_t, std::unordered_map<size_t, size_t>>;

int main(int argc, char* argv[]) {

    size_t d = 2; // dimension of Euclidean Space
    size_t n = 100;

    // maximum simplex dimension
    size_t maxdim = bats::util::io::parse_argv(argc, argv, "-maxdim", 3);
    double rmax = bats::util::io::parse_argv(argc, argv, "-rmax", 0.2);
    std::string fname = bats::util::io::parse_argv(argc, argv, "-file", std::string(""));

    DataSet<double> x;
    if (fname.empty()) {
        x = sample_cube<double>(d, n, 0); // seed with 0
    } else {
        x = DataSet(read_point_cloud(fname));
    }
    

    auto start1 = std::chrono::steady_clock::now();
    auto X = RipsComplex<CpxT>(x, LInfDist(), rmax, maxdim);
    auto end1 = std::chrono::steady_clock::now();
    std::cout << "Time to construct Rips Complex: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1).count()
        << "ms" << std::endl;
    X.print_summary();

    auto F = RipsFiltration<CpxT>(x, LInfDist(), rmax, maxdim);
    // F.print_summary();

    // Create new data y
    DataSet<double> y;
    y = sample_cube<double>(d, n, 0); 
    // auto Y = RipsComplex<CpxT>(y, LInfDist(), rmax, maxdim);
    auto F2 = RipsFiltration<CpxT>(x, LInfDist(), rmax, maxdim);
    // standard reduction of Y
    auto start = std::chrono::steady_clock::now();
    auto FC2 = bats::Chain(F2, FT());
    auto RFC2 = bats::ReducedFilteredChainComplex(
        FC2,
        bats::standard_reduction_flag(),
        bats::compute_basis_flag()
    );
    auto end = std::chrono::steady_clock::now();
    std::cout << "reduction: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
        << "ms" << std::endl;

    // // compare different options
	// {
    //     // standard reduction
	// 	start = std::chrono::steady_clock::now();
	// 	auto FC = bats::Chain(F, FT());
	// 	auto RFC = bats::ReducedFilteredChainComplex(
	// 		FC,
	// 		bats::standard_reduction_flag(),
	// 		bats::compute_basis_flag()
	// 	);
	// 	end = std::chrono::steady_clock::now();
	//     std::cout << "reduction: "
	//         << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
	//         << "ms" << std::endl;
	// 	RFC.print_summary(true);

	// 	// Compute update info
	// 	start = std::chrono::steady_clock::now();
	// 	auto uinfo = bats::Update_info(F, F2);
	// 	end = std::chrono::steady_clock::now();
	//     std::cout << "\nConstruction of Update info: "
	//         << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
	//         << "ms" << std::endl;

	// 	// Compute update
	// 	start = std::chrono::steady_clock::now();
	// 	RFC.update_filtration_general(uinfo, bats::standard_reduction_flag());
	// 	end = std::chrono::steady_clock::now();
	// 	std::cout << "\nUpdate persistence: "
	// 		<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
	// 		<< "ms" << std::endl;

    //     if(test_reduce_result(RFC , RFC2)){
    //         std::cout << "By compare two RFCC, update on RFCC success!!" << std::endl;
    //     }
	// }

    // // homology
    // { 
    //     std::cout << "DGVS(-1)" << std::endl;
    //     auto CX = FilteredDGVectorSpace<double, MT>(F, -1);
    //     auto RX = ReducedFilteredDGVectorSpace(CX);
    //     for (int k = 0; k < 3; ++k) {
    //         std::cout << "\t hdim " << k << ": " << RX.hdim(k) << std::endl;
    //     }

    //     // Compute update info
	// 	start = std::chrono::steady_clock::now();
	// 	auto uinfo = bats::Update_info(F, F2);
	// 	end = std::chrono::steady_clock::now();
	//     std::cout << "\nConstruction of Update info: "
	//         << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
	//         << "ms" << std::endl;

	// 	// Compute update
	// 	start = std::chrono::steady_clock::now();
	// 	RX.update_filtration_general(uinfo, bats::standard_reduction_flag());
	// 	end = std::chrono::steady_clock::now();
	// 	std::cout << "Update persistence: "
	// 		<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
	// 		<< "ms" << std::endl;

    //     // check result
    //     auto CX2 = FilteredDGVectorSpace<double, MT>(F2, -1);
    //     auto RX2 = ReducedFilteredDGVectorSpace(CX2);
    //     if(test_reduce_result(RX , RX2)){
    //         std::cout << "By compare two RFCC, update on RFCC success!!" << std::endl;
    //     }
    // }

    // cohomology
    {
        std::cout << "DGVS(+1)" << std::endl;
        auto CX = FilteredDGVectorSpace<double, MT>(F, +1);
        auto RX = ReducedFilteredDGVectorSpace(CX);
        for (int k = 0; k < 3; ++k) {
            std::cout << "\t hdim " << k << ": " << RX.hdim(k) << std::endl;
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
		RX.update_filtration_general(uinfo, bats::standard_reduction_flag());
		end = std::chrono::steady_clock::now();
		std::cout << "Update persistence: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
			<< "ms" << std::endl;

        // check result
        auto CX2 = FilteredDGVectorSpace<double, MT>(F2, -1);
        auto RX2 = ReducedFilteredDGVectorSpace(CX2);
        if(test_reduce_result(RX , RX2)){
            std::cout << "By compare two RFCC, update on RFCC success!!" << std::endl;
        }
    }



    return 0;
}
