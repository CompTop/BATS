#include <iostream>
#include <vector>
#include <chrono>

#include <bats.h>


#define FT ModP<int, 3>
#define VT SparseVector<FT>
#define MT ColumnMatrix<VT>

int main() {

    size_t d = 3; // dimension of Euclidean Space
    // sphere will be sampled on S^{d-1}
    size_t n = 1000; // number of samples
    size_t m = 20; // number of landmarks

    // maximum simplex dimension
    size_t maxdim = d;

    auto dist = Euclidean();

    auto X = sample_sphere<double>(d, n);
    //X.data.print();
    auto L = greedy_landmarks(X, m, dist);
    //L.data.print();

    auto W = WitnessComplex(X, L, dist, maxdim);
    W.print_summary();
    //W.boundary_csc(1).print();

    auto CX = ChainComplex<MT>(W);

    auto RX = ReducedChainComplex(CX);

    for (size_t k = 0; k < RX.maxdim()+1; k++) {
        std::cout << "betti " << k << " = " << RX.hdim(k) << std::endl;
    }

    return 0;
}
