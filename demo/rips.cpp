
#include <topology/data_gen.h>
#include <topology/rips.h>

#include <linalg/field.h>
#include <linalg/sparse_vector.h>
#include <linalg/col_matrix.h>

#include <chain/chain_complex.h>
#include <homology/basis.h>

#include <iostream>
#include <vector>
#include <chrono>


#define FT ModP<int, 3>
#define VT SparseVector<FT>
#define MT ColumnMatrix<VT>

int main() {

    size_t d = 2; // dimension of Euclidean Space
    size_t n = 100;

    // maximum simplex dimension
    size_t maxdim = d+1;

    auto x = sample_cube<double>(d, n);

    auto X = RipsComplex(x, LInfDist(), 0.2, maxdim);

    auto CX = ChainComplex<MT>(X);

    auto RX = ReducedChainComplex(CX);

    for (size_t k = 0; k < RX.maxdim()+1; k++) {
        std::cout << "betti " << k << " = " << RX.hdim(k) << std::endl;
    }




    return 0;
}
