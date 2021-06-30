#include <bats.hpp>
#include <util/io.hpp>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include <chrono>
#include <random>
#include <fstream>

using FT = ModP<int, 2>;

using CpxT = bats::LightSimplicialComplex<size_t, std::unordered_map<size_t, size_t>>;
// using CpxT = bats::SimplicialComplex;

template<typename T>
auto time_elapse(T t1, T t0){
    return std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
}


int main(int argc, char* argv[]) {
    // Create an output filestream object
    std::ofstream myFile("results.csv");
    
    // Send data to the stream
    
    myFile << "build_Filtration,build_updating_information,build/update_FCC,build/update_RFCC,whole\n";
    
    std::vector<size_t> num_pts = {5, 10, 15, 20};
    for(auto& n: num_pts){
        myFile << "N = "<< n <<"\n";
        // std::vector<std::vector<std::vector<size_t>>> Splx_3;
        size_t d = 2; // dimension of Euclidean Space
        // size_t n = 50; // number of points to sample

        // maximum simplex dimension
        size_t maxdim = bats::util::io::parse_argv(argc, argv, "-maxdim", 3);
        double rmax = bats::util::io::parse_argv(argc, argv, "-rmax", 2);

        // create dataset
        auto X = bats::sample_sphere<double>(d, n);
        auto X_data = X.data;
        X_data.swap_rows(0,1);
        auto Y = bats::DataSet(X_data);

        // Compute Persistent Homology on old filtration
        auto dist = bats::Euclidean(); // metric
        auto F_X = bats::RipsFiltration<CpxT>(X, dist, rmax, maxdim);
        auto FCC = bats::Chain(F_X, FT());  //Build FilteredChainComplex
        auto RFCC = bats::Reduce(FCC); // Build ReducedFilteredChainComplex
        // RFCC.print_summary();


        /*
        Option 1, rebuild 
        */
        auto start = std::chrono::steady_clock::now();

        auto t0 = std::chrono::steady_clock::now();
        auto F_Y = bats::RipsFiltration<CpxT>(Y, dist, rmax, maxdim);
        auto t1 = std::chrono::steady_clock::now();
        myFile << time_elapse(t1,t0)<<",,";

        t0 = std::chrono::steady_clock::now();
        auto FCC_Y = bats::Chain(F_Y, FT()); //FilteredChainComplex
        t1 = std::chrono::steady_clock::now();
        myFile << time_elapse(t1,t0)<<",";
        

        t0 = std::chrono::steady_clock::now();
        auto RFCC_Y = bats::Reduce(FCC_Y); // Build ReducedFilteredChainComplex
        t1 = std::chrono::steady_clock::now();
        myFile << time_elapse(t1,t0)<<",";

        auto end = std::chrono::steady_clock::now();
        myFile << time_elapse(end,start)<<"\n";


        /*
        Option 2, update FilteredChainComplex with filtered information
        */ 
        {
            FCC = bats::Chain(F_X, FT()); // necessary since last option will modify FCC ??
            
            start = std::chrono::steady_clock::now();

            t0 = std::chrono::steady_clock::now();
            F_Y = bats::RipsFiltration<CpxT>(Y, dist, rmax, maxdim);
            t1 = std::chrono::steady_clock::now();
            myFile << time_elapse(t1,t0) <<",";

            t0 = std::chrono::steady_clock::now();
            auto UI = bats::Update_info(F_X, F_Y);
            t1 = std::chrono::steady_clock::now();
            myFile << time_elapse(t1,t0) <<",";

            t0 = std::chrono::steady_clock::now();
            FCC.update_filtration_general(UI);
            t1 = std::chrono::steady_clock::now();
            myFile << time_elapse(t1,t0) <<",";

            
            t0 = std::chrono::steady_clock::now();
            auto RFCC2 = bats::Reduce(FCC); //ReducedFilteredChainComplex
            t1 = std::chrono::steady_clock::now();
            myFile << time_elapse(t1,t0) <<",";

            end = std::chrono::steady_clock::now();
            myFile << time_elapse(end,start) <<"\n";
        }

        /*
        Option 3, update ReducedFilteredChainComplex 
        */ 
        {
            FCC = bats::Chain(F_X, FT()); // necessary since last option will modify FCC
            RFCC = bats::Reduce(FCC);

            start = std::chrono::steady_clock::now();

            t0 = std::chrono::steady_clock::now();
            F_Y = bats::RipsFiltration<CpxT>(Y, dist, rmax, maxdim);
            t1 = std::chrono::steady_clock::now();
            myFile << time_elapse(t1,t0) <<",";

            t0 = std::chrono::steady_clock::now();
            auto UI = bats::Update_info(F_X, F_Y); 
            t1 = std::chrono::steady_clock::now();
            myFile << time_elapse(t1,t0) <<",";
            
            t0 = std::chrono::steady_clock::now();
            RFCC.update_filtration_general(UI);
            t1 = std::chrono::steady_clock::now();
            myFile << time_elapse(t1,t0) <<",";

            end = std::chrono::steady_clock::now();
            myFile << time_elapse(end,start) <<"\n";
        }
    }

    // Close the file
    myFile.close();

    return EXIT_SUCCESS;
}