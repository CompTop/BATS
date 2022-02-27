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
    size_t n = 350;

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

    {
        std::cout << "DGVS(-1)" << std::endl;
        auto CX = FilteredDGVectorSpace<double, MT>(F, -1);
        auto RX = ReducedFilteredDGVectorSpace(CX);
        for (int k = 0; k < 3; ++k) {
            std::cout << "\t hdim " << k << ": " << RX.hdim(k) << std::endl;
        }
    }

    {
        std::cout << "DGVS(+1)" << std::endl;
        auto CX = FilteredDGVectorSpace<double, MT>(F, +1);
        auto RX = ReducedFilteredDGVectorSpace(CX);
        for (int k = 0; k < 3; ++k) {
            std::cout << "\t hdim " << k << ": " << RX.hdim(k) << std::endl;
        }
    }


    // {
    //     std::cout << "DGVS(+1)" << std::endl;
    //     auto CX = DGVectorSpace<MT>(X, +1);
    //     auto RX = ReducedDGVectorSpace(CX);
    //     for (int k = 0; k < 2; ++k) {
    //         std::cout << "\t hdim " << k << ": " << RX.hdim(k) << std::endl;
    //     }
    // }

    auto CX = ChainComplex<MT>(X);

    {
        auto RX = __ReducedChainComplex(X, FT(), bats::extra_reduction_flag());
    }

    for (size_t k =0; k < 2; k++) {
        std::cout << "\ntesting flags, no basis" << std::endl;
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
    }

    for (size_t k =0; k < 2; k++) {
        std::cout << "\ntesting flags, with basis" << std::endl;
    {
        auto start = std::chrono::steady_clock::now();
        auto RX = ReducedChainComplex(
            CX,
            bats::standard_reduction_flag(),
            bats::compute_basis_flag()
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
            bats::compression_flag(),
            bats::compute_basis_flag()
        );
        auto end = std::chrono::steady_clock::now();
        std::cout << "standard reduction with compression: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
            << "ms" << std::endl;
    }
    {
        auto start = std::chrono::steady_clock::now();
        auto RX = ReducedChainComplex(
            CX,
            bats::extra_reduction_flag(),
            bats::compute_basis_flag()
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
            bats::compression_flag(),
            bats::compute_basis_flag()
        );
        auto end = std::chrono::steady_clock::now();
        std::cout << "extra reduction with compression: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
            << "ms" << std::endl;
    }

    }

    return 0;
}
