#include <bats.hpp>
#include <util/io.hpp>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include <chrono>
#include <random>

using FT = ModP<int, 2>;
#define VT SparseVector<FT>
#define MT ColumnMatrix<VT>

using CpxT = bats::LightSimplicialComplex<size_t, std::unordered_map<size_t, size_t>>;
// using CpxT = bats::SimplicialComplex;
using namespace bats;

int main(int argc, char* argv[]) {
    // std::vector<std::vector<std::vector<size_t>>> Splx_3;
	size_t d = 2; // dimension of Euclidean Space
	size_t n = 100; // number of points to sample

	// maximum simplex dimension
    size_t maxdim = bats::util::io::parse_argv(argc, argv, "-maxdim", 2);
    double rmax_X = bats::util::io::parse_argv(argc, argv, "-rmax_X", 4.3);
    double rmax_Y = bats::util::io::parse_argv(argc, argv, "-rmax_Y", 0.0);
    auto X = bats::sample_sphere<double>(d, n);
    auto X_data = X.data;
    auto Y = bats::DataSet(X_data);
    add_normal_noise(Y, 0, 0.1);

	auto dist = bats::Euclidean(); // metric

    auto F_X = bats::RipsFiltration<CpxT>(X, dist, rmax_X, maxdim);
    
    // print_filtration_info(F_X);
    auto FCC = bats::Chain(F_X, FT());  //Build FilteredChainComplex
    auto RFCC = bats::Reduce(FCC); // Build ReducedFilteredChainComplex
    // RFCC.print_summary();

    /*
    Option 1, rebuild 
    */
    std::cout << "Rebuild everything"<<std::endl;
    auto start = std::chrono::steady_clock::now();
    auto F_Y = bats::RipsFiltration<CpxT>(Y, dist, rmax_Y, maxdim);
    // print_filtration_info(F_Y);

    auto t0 = std::chrono::steady_clock::now();
    auto FCC_Y = bats::Chain(F_Y, FT()); //FilteredChainComplex
    auto t1 = std::chrono::steady_clock::now();
    std::cout << "\tbuild FCC takes "
        << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
        << "ms" << std::endl;

    t0 = std::chrono::steady_clock::now();
    auto RFCC_Y = bats::Reduce(FCC_Y); // Build ReducedFilteredChainComplex
    t1 = std::chrono::steady_clock::now();
    std::cout << "\tbuild RFCC takes "
        << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
        << "ms" << std::endl;

    auto end = std::chrono::steady_clock::now();
    std::cout << "The whole process takes "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
        << "ms" << std::endl;


    // /*
    // Option 3, update ReducedFilteredChainComplex 
    // */ 
    // {
    //     FCC = bats::Chain(F_X, FT()); // necessary since last option will modify FCC ??
    //     RFCC = bats::Reduce(FCC);

    //     std::cout << "\nUpdate Reduced Filtered Chain Complex "<< std::endl;
    //     start = std::chrono::steady_clock::now();
        
    //     t0 = std::chrono::steady_clock::now();

    //     auto UI = bats::Update_info(F_X, F_Y); // get unfiltered information
    //     // get filtered info, this step uncessary if cells in filtration has been sorted by their filtration values
    //     // UI.filtered_info(FCC.perm); 
    //     t1 = std::chrono::steady_clock::now();
    //     std::cout << "\tbuild Updating Information success and";
    //     std::cout << " takes "
    //         << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
    //         << "ms" << std::endl;
        

    //     t0 = std::chrono::steady_clock::now();
    //     RFCC.update_filtration_general(UI);
    //     t1 = std::chrono::steady_clock::now();
    //     std::cout << "\tupdate RFCC success and";
    //     std::cout << " takes "
    //         << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
    //         << "ms" << std::endl;

    //     end = std::chrono::steady_clock::now();
    //     std::cout << "\tthe whole process takes "
    //         << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
    //         << "ms" << std::endl;
        
    //     // check test results
    //     if(test_reduce_result(RFCC_Y , RFCC)){
    //         std::cout << "By compare two RFCC, this method success!!"  << std::endl;
    //     }
        
    // }

    // /*
    // Option 4, update DGVS -1 (homology)
    // */ 
    // {
    //     std::cout << "\nDGVS(-1)" << std::endl;
    //     t0 = std::chrono::steady_clock::now();
    //     auto CX = FilteredDGVectorSpace<double, MT>(F_X, -1);
    //     auto RX = ReducedFilteredDGVectorSpace(CX, bats::standard_reduction_flag(), 
    //             bats::clearing_flag(), bats::compute_basis_flag());
    //     t1 = std::chrono::steady_clock::now();
    //     std::cout << "\tCompute RFCC of X takes";
    //     std::cout << " takes "
    //         << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
    //         << "ms" << std::endl;

    //     start = std::chrono::steady_clock::now();
    //     F_Y = bats::RipsFiltration<CpxT>(Y, dist, rmax - 0.5, maxdim);
        
    //     t0 = std::chrono::steady_clock::now();
    //     auto UI = bats::Update_info(F_X, F_Y);
    //     // get filtered info, this step uncessary if cells in filtration has been sorted by their filtration values
    //     // UI.filtered_info(FCC.perm); 
    //     t1 = std::chrono::steady_clock::now();
    //     std::cout << "\tbuild Updating Information success and";
    //     std::cout << " takes "
    //         << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
    //         << "ms" << std::endl;
        

    //     t0 = std::chrono::steady_clock::now();
    //     RX.update_filtration_general(UI, bats::standard_reduction_flag());
    //     // RFCC.update_filtration_general(UI);
    //     t1 = std::chrono::steady_clock::now();
    //     std::cout << "\tupdate RFCC success and";
    //     std::cout << " takes "
    //         << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
    //         << "ms" << std::endl;

    //     end = std::chrono::steady_clock::now();
    //     std::cout << "\tthe whole process takes "
    //         << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
    //         << "ms" << std::endl;
        
    //     // check result
    //     auto CX2 = FilteredDGVectorSpace<double, MT>(F_Y, -1);
    //     auto RX2 = ReducedFilteredDGVectorSpace(CX2);
    //     if(test_reduce_result(RX , RX2)){
    //         std::cout << "By comparing two RFCC, update on RFCC success!!" << std::endl;
    //     }
        
    // }

    // /*
    // Option 4, update DGVS -1 (homology)
    // */ 
    // {
    //     std::cout << "\nDGVS(-1) + clearing update" << std::endl;
    //     t0 = std::chrono::steady_clock::now();
    //     auto CX = FilteredDGVectorSpace<double, MT>(F_X, -1);
    //     auto RX = ReducedFilteredDGVectorSpace(CX, bats::standard_reduction_flag(), 
    //             bats::clearing_flag(), bats::compute_basis_flag());
    //     t1 = std::chrono::steady_clock::now();
    //     std::cout << "\tCompute RFCC of X takes";
    //     std::cout << " takes "
    //         << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
    //         << "ms" << std::endl;

    //     start = std::chrono::steady_clock::now();
    //     F_Y = bats::RipsFiltration<CpxT>(Y, dist, rmax - 0.5, maxdim);
        
    //     t0 = std::chrono::steady_clock::now();
    //     auto UI = bats::Update_info(F_X, F_Y);
    //     // get filtered info, this step uncessary if cells in filtration has been sorted by their filtration values
    //     // UI.filtered_info(FCC.perm); 
    //     t1 = std::chrono::steady_clock::now();
    //     std::cout << "\tbuild Updating Information success and";
    //     std::cout << " takes "
    //         << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
    //         << "ms" << std::endl;
        

    //     t0 = std::chrono::steady_clock::now();
    //     RX.update_filtration_general_clearing(UI, bats::standard_reduction_flag());
    //     // RFCC.update_filtration_general(UI);
    //     t1 = std::chrono::steady_clock::now();
    //     std::cout << "\tupdate RFCC success and";
    //     std::cout << " takes "
    //         << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
    //         << "ms" << std::endl;

    //     end = std::chrono::steady_clock::now();
    //     std::cout << "\tthe whole process takes "
    //         << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
    //         << "ms" << std::endl;
        
    //     // check result
    //     auto CX2 = FilteredDGVectorSpace<double, MT>(F_Y, -1);
    //     auto RX2 = ReducedFilteredDGVectorSpace(CX2);
    //     if(test_reduce_result(RX , RX2)){
    //         std::cout << "By comparing two RFCC, update on RFCC success!!" << std::endl;
    //     }
        
    // }

    /*
    Option 5, update DGVS +1 (cohomology)
    */ 
    {
        std::cout << "\nDGVS(+1)" << std::endl;
        t0 = std::chrono::steady_clock::now();
        auto CX = FilteredDGVectorSpace<double, MT>(F_X, +1);
        auto RX = ReducedFilteredDGVectorSpace(CX, bats::standard_reduction_flag(), 
                bats::clearing_flag(), bats::compute_basis_flag());
        t1 = std::chrono::steady_clock::now();
        std::cout << "\tCompute RFCC of X takes";
        std::cout << " takes "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
            << "ms" << std::endl;

        start = std::chrono::steady_clock::now();
        t0 = std::chrono::steady_clock::now();
        F_Y = bats::RipsFiltration<CpxT>(Y, dist, rmax_Y, maxdim);
        t1 = std::chrono::steady_clock::now();
        std::cout << "\tConstruct filtration Y takes";
        std::cout << " takes "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
            << "ms" << std::endl;

        t0 = std::chrono::steady_clock::now();
        auto UI = bats::Update_info(F_X, F_Y);
        // get filtered info, this step uncessary if cells in filtration has been sorted by their filtration values
        // UI.filtered_info(FCC.perm); 
        t1 = std::chrono::steady_clock::now();
        std::cout << "\tbuild Updating Information success and";
        std::cout << " takes "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
            << "ms" << std::endl;
        

        t0 = std::chrono::steady_clock::now();
        RX.update_filtration_general(UI);
        // RFCC.update_filtration_general(UI);
        t1 = std::chrono::steady_clock::now();
        std::cout << "\tupdate RFCC success and";
        std::cout << " takes "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
            << "ms" << std::endl;

        end = std::chrono::steady_clock::now();
        std::cout << "\tthe whole process takes "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
            << "ms" << std::endl;
        
        // check result
        auto CX2 = FilteredDGVectorSpace<double, MT>(F_Y, +1);
        auto RX2 = ReducedFilteredDGVectorSpace(CX2);
        if(test_reduce_result(RX , RX2)){
            std::cout << "By comparing two RFCC, update on RFCC success!!" << std::endl;
        }
        
    }

    /*
    Option 5, update DGVS +1 (cohomology) + clearing
    */ 
    {
        std::cout << "\nDGVS(+1) + clearning update" << std::endl;
        t0 = std::chrono::steady_clock::now();
        auto CX = FilteredDGVectorSpace<double, MT>(F_X, +1);
        auto RX = ReducedFilteredDGVectorSpace(CX, bats::standard_reduction_flag(), 
                bats::clearing_flag(), bats::compute_basis_flag());
        t1 = std::chrono::steady_clock::now();
        std::cout << "\tCompute RFCC of X takes";
        std::cout << " takes "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
            << "ms" << std::endl;

        start = std::chrono::steady_clock::now();
        t0 = std::chrono::steady_clock::now();
        F_Y = bats::RipsFiltration<CpxT>(Y, dist, rmax_Y, maxdim);
        t1 = std::chrono::steady_clock::now();
        std::cout << "\tConstruct filtration Y takes";
        std::cout << " takes "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
            << "ms" << std::endl;

        t0 = std::chrono::steady_clock::now();
        auto UI = bats::Update_info(F_X, F_Y);
        // get filtered info, this step uncessary if cells in filtration has been sorted by their filtration values
        // UI.filtered_info(FCC.perm); 
        t1 = std::chrono::steady_clock::now();
        std::cout << "\tbuild Updating Information success and";
        std::cout << " takes "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
            << "ms" << std::endl;
        

        t0 = std::chrono::steady_clock::now();
        RX.update_filtration_general_clearing(UI);
        // RFCC.update_filtration_general(UI);
        t1 = std::chrono::steady_clock::now();
        std::cout << "\tupdate RFCC success and";
        std::cout << " takes "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
            << "ms" << std::endl;

        end = std::chrono::steady_clock::now();
        std::cout << "\tthe whole process takes "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
            << "ms" << std::endl;
        
        // check result
        auto CX2 = FilteredDGVectorSpace<double, MT>(F_Y, +1);
        auto RX2 = ReducedFilteredDGVectorSpace(CX2);
        if(test_reduce_result(RX , RX2)){
            std::cout << "By comparing two RFCC, update on RFCC success!!" << std::endl;
        }
        
    }

    return EXIT_SUCCESS;
}