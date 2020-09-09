#include <iostream>
#include <vector>
#include <chrono>

#include <bats.hpp>
#include <util/io.h>
#include <string>


#define FT ModP<int, 3>
#define VT SparseVector<FT>
#define MT ColumnMatrix<VT>

int main(int argc, char* argv[]) {

    size_t d = 2; // dimension of Euclidean Space
    size_t n = 350;

    // maximum simplex dimension
    size_t maxdim = parse_argv(argc, argv, "-maxdim", 3);
    double rmax = parse_argv(argc, argv, "-rmax", 0.2);
    std::string fname = parse_argv(argc, argv, "-file", std::string(""));

    DataSet<double> x;
    if (fname.empty()) {
        x = sample_cube<double>(d, n);
    } else {
        x = DataSet(read_point_cloud(fname));
    }

    auto X = RipsComplex(x, LInfDist(), rmax, maxdim);

    auto CX = ChainComplex<MT>(X);


    {
        auto start = std::chrono::steady_clock::now();
        auto RX = ReducedChainComplex(
            CX,
            bats::standard_reduction_flag()
        );
        auto end = std::chrono::steady_clock::now();
        std::cout << "standard reduction: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
            << "ms" << std::endl;
    }
    {
        auto start = std::chrono::steady_clock::now();
        auto RX = ReducedChainComplex(
            CX,
            bats::standard_reduction_flag(),
            bats::compression_flag()
        );
        auto end = std::chrono::steady_clock::now();
        std::cout << "reduction with compression: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
            << "ms" << std::endl;
    }
    {
        auto start = std::chrono::steady_clock::now();
        auto RX = ReducedChainComplex(
            CX,
            bats::standard_reduction_flag(),
            bats::clearing_flag()
        );
        auto end = std::chrono::steady_clock::now();
        std::cout << "reduction with clearing: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
            << "ms" << std::endl;
    }
    {
        auto start = std::chrono::steady_clock::now();
        auto RX = ReducedChainComplex(
            CX,
            bats::extra_reduction_flag()
        );
        auto end = std::chrono::steady_clock::now();
        std::cout << "extra reduction: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
            << "ms" << std::endl;
    }
    {
        auto start = std::chrono::steady_clock::now();
        auto RX = ReducedChainComplex(
            CX,
            bats::extra_reduction_flag(),
            bats::compression_flag()
        );
        auto end = std::chrono::steady_clock::now();
        std::cout << "extra reduction with compression: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
            << "ms" << std::endl;
    }
    {
        auto start = std::chrono::steady_clock::now();
        auto RX = ReducedChainComplex(
            CX,
            bats::extra_reduction_flag(),
            bats::clearing_flag()
        );
        auto end = std::chrono::steady_clock::now();
        std::cout << "extra reduction with clearing: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
            << "ms" << std::endl;
    }

    return 0;
}
