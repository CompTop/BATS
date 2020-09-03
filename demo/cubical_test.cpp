#include <iostream>
#include <vector>

#include <bats.hpp>

auto full_cube(const size_t N) {
	CubicalComplex X(3);

	for (size_t i = 0; i < N-1; i++) {
		for (size_t j = 0; j < N-1; j++) {
			for (size_t k = 0; k < N-1; k++) {
				X.add_recursive({i, i+1, j, j+1, k, k+1});
			}
		}
	}
	return X;
}

auto skeleton_filtration_diagram(const CubicalComplex &X) {
	Diagram<CubicalComplex, CellularMap> D(X.maxdim() + 1, X.maxdim());
	for (size_t d = 0; d < X.maxdim() + 1; d++) {
		D.set_node(d, X.skeleton(d));
	}
	for (size_t d = 0; d < X.maxdim(); d++) {
		D.set_edge(d, d, d+1, CubicalMap(D.node_data(d), D.node_data(d+1)));
	}

	return D;
}

int main() {

	std::cout << "\nCreating cubical complex" << std::endl;

	auto X = full_cube(2);

	X.print_summary();

	auto D = skeleton_filtration_diagram(X);

	// diagram in Chain
	std::cout << "Chain functor" << std::endl;
	auto ChainDgm = __Chain(D, ModP<int, 2>());

	// diagram in Homology
	std::cout << "Homology functor" << std::endl;
	for (size_t d = 0; d < 3; d++) {
		std::cout << "dimension " << d << std::endl;
		auto HkDgm = Hom(ChainDgm, d);

		auto ps = barcode_sparse_divide_conquer(HkDgm, d);
		for (auto p : ps) {
			std::cout << p.str() << std::endl;
		}
	}






	return 0;
}
