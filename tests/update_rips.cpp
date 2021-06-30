#include <bats.hpp>
#include <util/io.hpp>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include <chrono>
#include <random>

using FT = ModP<int, 2>;

using CpxT = bats::LightSimplicialComplex<size_t, std::unordered_map<size_t, size_t>>;
// using CpxT = bats::SimplicialComplex;


int main(int argc, char* argv[]) {
    // std::vector<std::vector<std::vector<size_t>>> Splx_3;
	size_t d = 2; // dimension of Euclidean Space
	size_t n = 50; // number of points to sample

	// maximum simplex dimension
    size_t maxdim = bats::util::io::parse_argv(argc, argv, "-maxdim", 3);
    double rmax = bats::util::io::parse_argv(argc, argv, "-rmax", 2);

    auto X = bats::sample_sphere<double>(d, n);
    auto X_data = X.data;
    X_data.swap_rows(0,1);
    auto Y = bats::DataSet(X_data);


	auto dist = bats::Euclidean(); // metric

    auto F_X = bats::RipsFiltration<CpxT>(X, dist, rmax, maxdim);
    
    // print_filtration_info(F_X);
    auto FCC = bats::Chain(F_X, FT());  //Build FilteredChainComplex
    auto RFCC = bats::Reduce(FCC); // Build ReducedFilteredChainComplex
    // RFCC.print_summary();

    /*
    Option 1, rebuild 
    */
    std::cout << "Rebuild everything"<<std::endl;
    auto start = std::chrono::steady_clock::now();
    auto F_Y = bats::RipsFiltration<CpxT>(Y, dist, rmax, maxdim);
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



    /*
    Option 2, update FilteredChainComplex with filtered information
    */ 
    {
        std::cout << "\nUpdate FilteredChainComplex"<< std::endl;
        FCC = bats::Chain(F_X, FT()); // necessary since last option will modify FCC ??
        
        start = std::chrono::steady_clock::now();
        F_Y = bats::RipsFiltration<CpxT>(Y, dist, rmax, maxdim);

        t0 = std::chrono::steady_clock::now();
        auto UI = bats::Update_info(F_X, F_Y); // unfiltered information
        // get filtered information, 
        // If we know the complex has been sorted by filtration value in advance
        // in the filtration construction step, then we do not need this step! 
        // UI.filtered_info(FCC.perm); // pass in old filtration permutation

        t1 = std::chrono::steady_clock::now();
        std::cout << "\tbuild Updating Information success and";
        std::cout << " takes "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
            << "ms" << std::endl;  

        t0 = std::chrono::steady_clock::now();
        FCC.update_filtration_general(UI);
        t1 = std::chrono::steady_clock::now();
        std::cout << "\tupdate FCC success and";
        std::cout << " takes "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
            << "ms" << std::endl;

        
        t0 = std::chrono::steady_clock::now();
        auto RFCC2 = bats::Reduce(FCC); //ReducedFilteredChainComplex
        t1 = std::chrono::steady_clock::now();
        std::cout << "\tbuild RFCC success and";
        std::cout << " takes "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
            << "ms" << std::endl;

        end = std::chrono::steady_clock::now();
        std::cout << "\tthe whole process takes "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
            << "ms" << std::endl;

        // check test results
        if(test_reduce_result(RFCC_Y , RFCC2)){
            std::cout << "By compare two RFCC, this method success!!" << std::endl;
        }
        

    }

    /*
    Option 3, update ReducedFilteredChainComplex 
    */ 
    {
        FCC = bats::Chain(F_X, FT()); // necessary since last option will modify FCC ??
        RFCC = bats::Reduce(FCC);

        std::cout << "\nUpdate Reduced Filtered Chain Complex "<< std::endl;
        start = std::chrono::steady_clock::now();
        F_Y = bats::RipsFiltration<CpxT>(Y, dist, rmax, maxdim);
        
        t0 = std::chrono::steady_clock::now();

        auto UI = bats::Update_info(F_X, F_Y); // get unfiltered information
        // UI.filtered_info(FCC.perm); // get filtered info, this step uncessary if cells in filtration has been sorted by their filtration values

        t1 = std::chrono::steady_clock::now();
        std::cout << "\tbuild Updating Information success and";
        std::cout << " takes "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
            << "ms" << std::endl;
        

        t0 = std::chrono::steady_clock::now();
        RFCC.update_filtration_general(UI);
        t1 = std::chrono::steady_clock::now();
        std::cout << "\tupdate RFCC success and";
        std::cout << " takes "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
            << "ms" << std::endl;

        end = std::chrono::steady_clock::now();
        std::cout << "\tthe whole process takes "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
            << "ms" << std::endl;
        
        // check test results
        if(test_reduce_result(RFCC_Y , RFCC)){
            std::cout << "By compare two RFCC, this method success!!"  << std::endl;
        }
        
    }

    return EXIT_SUCCESS;
}