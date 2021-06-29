#include <bats.hpp>
#include <util/io.hpp>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include <chrono>
#include <random>
#include <yuan/update_information.hpp>

using FT = ModP<int, 2>;

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
    std::cout << "\nsummary of simpilcial complex X" << std::endl;
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

    // filtration value on each simplex
    std::vector<std::vector<double>> vals = {{0.1,0.3,0.1,0.2,0.5},{5,2,6,4,2,3},{7,9}}; 
    /*
    Now let's build a filtration 
    */
    auto F = bats::Filtration(X, vals);
    std::cout << "\nFiltration information of F_X" << std::endl;
    print_filtration_info(F);
    auto FCC = bats::Chain(F, FT());  //Build FilteredChainComplex
    // std::cout << "\nOriginal FCC_X is " << std::endl;
    // for(auto & mat: FCC.C.boundary){
    //     mat.print();
    // }
    auto RFCC = bats::Reduce(FCC); // Build ReducedFilteredChainComplex
    // RFCC.print_summary();

    // new filtration value
    std::vector<std::vector<double>> vals_Y = {{0.1,0.0,0.3,0.1,0.2},
                                            {2,6,3,5,4,1,2,3},
                                            {8,10,9,7}}; 

    /*
    rebuild
    */
    std::cout << "\nRebuild everything";
    auto t0 = std::chrono::steady_clock::now();
    auto F_Y = bats::Filtration(Y, vals_Y);

    std::cout << "\nFiltration information of F_Y" << std::endl;
    print_filtration_info(F_Y);

    auto FCC_Y = bats::Chain(F_Y, FT()); //FilteredChainComplex
    auto RFCC_Y = bats::Reduce(FCC_Y); // Build ReducedFilteredChainComplex
    auto t1 = std::chrono::steady_clock::now();
            std::cout << " takes "
                << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
                << "ms" << std::endl;
    std::cout << "\nprint Summary of ReducedFilteredChainComplex:" << std::endl;
    RFCC_Y.print_summary();


    std::cout << "\nUpdate FilteredChainComplex";
    auto UI = Update_info(F, F_Y); // get unfiltered information

    // if cells infiltration was not sorted by filtration values when construction,
    // then we need to get filtered information.
    UI.filtered_info(FCC.perm);

    std::cout << "\nPrint the summary of filtration information" << std::endl;
    UI.print_detail();


    std::cout << "\nThe goal of FCC_Y is " << std::endl;
    for(auto & mat: FCC_Y.C.boundary){
        mat.print();
    }

    FCC.update_filtration_general(UI);

    std::cout << "\nAfter updating FCC_X is " << std::endl;
    for(auto & mat: FCC.C.boundary){
        mat.print();
    }

    // check if the FCC are the same
    for(size_t i = 0; i < FCC.C.boundary.size(); i++){
        if(FCC.C.boundary[i] == FCC_Y.C.boundary[i]){
            std::cout << "dimension "<< i << " match"<< std::endl;
        }
    }   

    return EXIT_SUCCESS;
}
