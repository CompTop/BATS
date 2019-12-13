#include <iostream>
#include <vector>
#include <chrono>

#include <bats.h>
#include <util/io.h>
#include <string>


#define FT ModP<int, 3>
#define VT SparseVector<FT>
#define MT ColumnMatrix<VT>

int main(int argc, char* argv[]) {

    size_t d = 2; // dimension of Euclidean Space
    size_t n = 100;

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

    auto RX = ReducedChainComplex(CX);

    for (size_t k = 0; k < RX.maxdim()+1; k++) {
        std::cout << "betti " << k << " = " << RX.hdim(k) << std::endl;
    }

    return 0;
}
