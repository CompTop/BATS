#include <bats.hpp>
#include <util/io.hpp>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include <chrono>
#include <random>

// using FT = ModP<int, 2>;
#define FT ModP<int, 2>
#define VT SparseVector<FT>
#define MT ColumnMatrix<VT>
using namespace bats;

using CpxT = bats::LightSimplicialComplex<size_t, std::unordered_map<size_t, size_t>>;
// using CpxT = bats::SimplicialComplex;

int main() {

    /*
    First we build a simplicial complex which can be used to extend filtrations
    */
    CpxT X(5,2); // 3 is the size of vertex set, 2 is the size of maximum simplex dimension
    CpxT Y(5,2);

    std::vector<size_t> s; //store temporary simplices 

    s = {0}; X.add(s);
    s = {1}; X.add(s);
    s = {3}; X.add(s);
    s = {2}; X.add(s);
    s = {4}; X.add(s);

    s = {0,3}; X.add(s);
    s = {1,3}; X.add(s);
    s = {3,4}; X.add(s);
    s = {0,1}; X.add(s);
    s = {1,2}; X.add(s);
    s = {0,2}; X.add(s);
    s = {0,1,2}; X.add(s);
    s = {0,1,3}; X.add(s);
    //summary of simpilcial complex X
    std::cout << "Summary of simpilcial complex X" << std::endl;
    X.print_summary(); 

    // next we construct a simplicial complex with H1 = H2 = H3 = 1
    s = {0}; Y.add(s);
    s = {2}; Y.add(s);
    s = {3}; Y.add(s);
    s = {4}; Y.add(s);
    s = {1}; Y.add(s);

    s = {1,2}; Y.add(s); // intersect 4 - 1 in X, 0 - 1
    s = {0,4}; Y.add(s);
    s = {0,1}; Y.add(s); // intersect 3 - 3 in X, 2 - 3
    s = {1,4}; Y.add(s);
    s = {0,2}; Y.add(s); // intersect 5 - 2 in X, 4 - 5
    s = {0,3}; Y.add(s); // intersect 0 - 4 in X, 5 - 0
    s = {1,3}; Y.add(s); // intersect 1 - 0 in X, 6 - 2
    s = {2,3}; Y.add(s);

    s = {0,1,2}; Y.add(s);
    s = {0,1,3}; Y.add(s);
    s = {1,2,3}; Y.add(s);
    s = {0,2,3}; Y.add(s);
    std::cout << "Summary of simpilcial complex Y" << std::endl;
    Y.print_summary(); 

    // filtration value on each simplex
    std::vector<std::vector<double>> vals = {{0.1,0.3,0.1,0.2,0.5},{5,2,6,4,2,3},{7,9}}; 
    /*
    Now let's build a filtration 
    */
    auto F = bats::Filtration(X, vals);
    std::cout << "Filtration information of F_X" << std::endl;
    print_filtration_info(F);
    auto FCC = bats::Chain(F, FT());  //Build FilteredChainComplex
    auto RFCC = bats::Reduce(FCC); // Build ReducedFilteredChainComplex
    std::cout << "Summary of RFCC of X:" << std::endl;
    RFCC.print_summary();

    // new filtration value
    std::vector<std::vector<double>> vals_Y = {{0.1,0.0,0.3,0.1,0.2},
                                            {2,6,3,5,4,1,2,3},
                                            {8,10,9,7}}; 

    /*
    Option 1, rebuild
    */
    auto F_Y = bats::Filtration(Y, vals_Y);
    std::cout << "Filtration information of Y" << std::endl;
    print_filtration_info(F_Y);

    std::cout << "\nRebuild RFCC"<< std::endl;
    auto start = std::chrono::steady_clock::now();
    auto FCC_Y = bats::Chain(F_Y, FT()); //FilteredChainComplex
    auto RFCC_Y = bats::Reduce(FCC_Y); // Build ReducedFilteredChainComplex
    std::cout << "Summary of RFCC of Y:" << std::endl;
    RFCC_Y.print_summary();
	auto end = std::chrono::steady_clock::now();
    std::cout << "Time: "
    << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
    << "ms" << std::endl;

    // /*
    // Option 2, update FilteredChainComplex with filtered information
    // */ 
    // {
    //     std::cout << "\nUpdate FilteredChainComplex"<< std::endl;
    //     FCC = bats::Chain(F, FT()); 
    //     F_Y = bats::Filtration(Y, vals_Y); // build new filtration
    //     auto UI = bats::Update_info(F, F_Y); // get unfiltered information
    //     // if cells infiltration was not sorted by filtration values when construction,
    //     // then we need to get filtered information,
    //     UI.filtered_info(FCC.perm);
    //     FCC.update_filtration_general(UI);
    //     auto RFCC_2 = bats::Reduce(FCC); //ReducedFilteredChainComplex
    //     RFCC_2.print_summary();
    //     // check test results
    //     if(test_reduce_result(RFCC_Y , RFCC_2)){
    //         std::cout << "This method success!!" << std::endl;
    //     }

    //     // for(auto & mat: RFCC_2.RC.R){
    //     //     std::cout << "reduced matrix now is:" << std::endl;
    //     //     mat.print();
    //     // }
    // }

    // /*
    // Option 3, update ReducedFilteredChainComplex with filtered information
    // */ 
    // {
    //     std::cout << "\nUpdate Reduced Filtered Chain Complex"<< std::endl;
    //     FCC = bats::Chain(F, FT()); // necessary since last option might modify
    //     RFCC = bats::Reduce(FCC);

    //     start = std::chrono::steady_clock::now();
    //     F_Y = bats::Filtration(Y, vals_Y);
        
    //     auto UI = bats::Update_info(F, F_Y); // get unfiltered information
    //     // necessary to get filtered information,
    //     // if cells infiltration was not sorted by filtration values when construction
    //     UI.filtered_info(FCC.perm);  

    //     RFCC.update_filtration_general(UI);
    //     RFCC.print_summary();
    //     end = std::chrono::steady_clock::now();
        
    //     std::cout << "Time: "
    //         << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
    //         << "ms" << std::endl;
    //     // check test results
    //     if(test_reduce_result(RFCC_Y , RFCC)){
    //         std::cout << "This method success!!" << std::endl;
    //     }
    //     // for(auto & mat: RFCC.RC.R){
    //     //     std::cout << "reduced matrix now is:" << std::endl;
    //     //     mat.print();
    //     // }
        
    // }



    // homology
    { 
        std::cout << "\nDGVS(-1)" << std::endl;
        auto CX = FilteredDGVectorSpace<double, MT>(F, -1);
        auto RX = ReducedFilteredDGVectorSpace(CX);
        for (int k = 0; k < 3; ++k) {
            std::cout << "\t hdim " << k << ": " << RX.hdim(k) << std::endl;
        }

        // Compute update info
		start = std::chrono::steady_clock::now();
        F_Y = bats::Filtration(Y, vals_Y);
		auto uinfo = bats::Update_info(F, F_Y);
		end = std::chrono::steady_clock::now();
	    std::cout << "\nConstruction of Update info: "
	        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
	        << "ms" << std::endl;

		// Compute update
		start = std::chrono::steady_clock::now();
        uinfo.filtered_info(RX.perm); 
		RX.update_filtration_general(uinfo, bats::standard_reduction_flag());
		end = std::chrono::steady_clock::now();
		std::cout << "Update persistence: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
			<< "ms" << std::endl;

        // check result
        auto CX2 = FilteredDGVectorSpace<double, MT>(F_Y, -1);
        auto RX2 = ReducedFilteredDGVectorSpace(CX2);
        if(test_reduce_result(RX , RX2)){
            std::cout << "By comparing two RFCC, update on RFCC success!!" << std::endl;
        }
    }

    // // cohomology
    // {
    //     std::cout << "\nDGVS(+1)" << std::endl;
    //     std::cout << "RFCC of X" << std::endl;
    //     auto CX = FilteredDGVectorSpace<double, MT>(F, +1);
    //     auto RX = ReducedFilteredDGVectorSpace(CX);
    //     RX.RC.print_summary();
        
    //     std::cout << "RFCC of Y" << std::endl;
    //     auto CX2 = FilteredDGVectorSpace<double, MT>(F_Y, +1);
    //     auto RX2 = ReducedFilteredDGVectorSpace(CX2);
    //     RX2.RC.print_summary();

    //     // Compute update info
	// 	start = std::chrono::steady_clock::now();
	// 	auto uinfo = bats::Update_info(F, F_Y);
	// 	end = std::chrono::steady_clock::now();
	//     std::cout << "\nConstruction of Update info: "
	//         << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
	//         << "ms" << std::endl;

	// 	// Compute update
	// 	start = std::chrono::steady_clock::now();
    //     uinfo.filtered_info(RX.perm, +1); 
	// 	RX.update_filtration_general(uinfo, bats::standard_reduction_flag());
	// 	end = std::chrono::steady_clock::now();
	// 	std::cout << "Update persistence: "
	// 		<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
	// 		<< "ms" << std::endl;

    //     std::cout << "\nNew summary of RX" << std::endl;
    //     RX.RC.print_summary();

    //     // check result
        
    //     if(test_reduce_result(RX , RX2)){
    //         std::cout << "By comparing two RFCC, update on RFCC success!!" << std::endl;
    //     }
    // }


    return EXIT_SUCCESS;
}