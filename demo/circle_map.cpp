#include <iostream>
#include <linalg/field.h>
#include <linalg/sparse_vector.h>
#include <linalg/col_matrix.h>

#include <complex/cell_complex.h>
#include <homology/reduction.h>
#include <complex/cell_map.h>

#include <chain/chain_complex.h>
#include <homology/basis.h>

#include <chain/chain_map.h>
#include <homology/induced_map.h>

#define FT ModP<int, 3>
#define VT SparseVector<FT, size_t>
#define Vint SparseVector<int, size_t>
#define MT ColumnMatrix<VT>
#define MTint ColumnMatrix<SparseVector<int, size_t>>

// Create a circle with 1 zero cell and 1 edge
CellComplex Circle() {
    CellComplex X(1);
    X.add_vertices(1);

    std::vector<size_t> b = {};
    std::vector<int> c = {};

    X.add(b, c, 1);
    return X;
}

// map that twists circle k times
CellularMap CircleTwist(int k) {
    CellularMap M(1);
    M[0] = MTint::identity(1);
    std::vector<size_t> ind = {0};
    std::vector<int> val = {k};
    std::vector<Vint> col = {Vint(ind, val)};
    M[1] = MTint(1,1, col);
    return M;
}

void run_example(int k) {
	std::cout << "\n\nexample: twist x " << k << std::endl;
	auto X = Circle();
	auto f = CircleTwist(k);

	auto CX = ChainComplex<MT>(X);
	auto RX = ReducedChainComplex(CX);
	auto F = ChainMap<MT>(f);

	auto FH0 = induced_map(F, RX, RX, 0);
	std::cout << "induced map in dim 0" << std::endl;
	FH0.print();

	auto FH1 = induced_map(F, RX, RX, 1);
	std::cout << "induced map in dimension 1" << std::endl;
	FH1.print();
}

int main() {

	run_example(1);
	run_example(2);
	run_example(3);
	run_example(4);
	run_example(-1);

	return 0;
}
