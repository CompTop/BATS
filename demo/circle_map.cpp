#include <iostream>
#include <bats.h>

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

SimplicialComplex CycleGraph_Simplicial(size_t n) {
    SimplicialComplex X(1); // 1 indicates maximum dimension of cell

    // add 0-cells
    for (size_t i = 0; i < n; i++) {
        //x = {i};
        X.add({i});
    }

    // add 1-cells
    for (size_t i = 0; i < n; i++) {
        X.add({i, (i+1) % n});
    }

    return X;
}

// Create a cycle graph on n vertices
CellComplex CycleGraph_Cellular(size_t n) {
    CellComplex X(1);
    X.add_vertices(n);

    std::vector<size_t> b;
    std::vector<int> c;

    for (size_t i = 0; i < n; i++) {
        b = {i, (i+1) % n};
        c = {-1, 1};
        X.add(b, c, 1);
    }
    return X;
}

// map that twists circle k times on simple cell complex
CellularMap CircleTwist(int k) {
    CellularMap M(1);
    M[0] = MTint::identity(1);
    std::vector<size_t> ind = {0};
    std::vector<int> val = {k};
    std::vector<Vint> col = {Vint(ind, val)};
    M[1] = MTint(1,1, col);
    return M;
}

// sends circle on n nodes to circle on m nodes twisting k times
CellularMap CircleTwist_Cellular(int k, size_t n, size_t m) {
    std::vector<size_t> ind;
    std::vector<int> val;
    size_t k2 = (k * m) / n; // stride length
    if( k2 * n !=  k * m ) {
        throw "incompatible sizes!";
    }

    CellularMap F(1);

    // map 0-cells
    std::vector<Vint> col0;
    for (size_t j = 0; j < n; j++) {
        ind = { (j * k2) % m };
        val = {1};
        col0.emplace_back(Vint(ind, val));
    }
    F[0] = MTint(m, n, col0);

    // map 1-cells
    std::vector<Vint> col1;
    for (size_t j = 0; j < n; j++) {
        ind.clear();
        val.clear();
        // each edge maps to k edges
        for (size_t i = 0; i < k2; i++) {
            ind.emplace_back(((j * k2) + i) % m);
            val.emplace_back(1);
        }
        col1.emplace_back(Vint(ind, val));
    }
    F[1] = MTint(m, n, col1);

    return F;
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

void run_example_cellular(int k) {
	std::cout << "\n\nexample: twist x " << k << std::endl;
	auto X = CycleGraph_Cellular(24);
	auto f = CircleTwist_Cellular(k, 24, 24);

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

	std::cout << "SIMPLE CELL STRUCTURE" << std::endl;
	run_example(1);
	run_example(2);
	run_example(3);
	run_example(4);
	run_example(-1);

	std::cout << "\n\nCYCLE GRAPH STRUCTURE" << std::endl;
	run_example_cellular(1);
	run_example_cellular(2);
	run_example_cellular(3);
	run_example_cellular(4);

	return 0;
}
