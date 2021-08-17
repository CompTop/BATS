#include <bats.hpp>
using namespace bats;

int main() {

	size_t n = 1000;
	auto X = sample_sphere<double>(2, n);

	auto dist = Euclidean(); //AngleDist();

	std::cout << "constructing zigzag" << std::endl;
	// auto [SetDgm, ripsvals] = DiscreteMorozovZigzagSets(X, dist, 10.0);
	auto [SetDgm, ripsvals] = OscillatingRipsZigzagSets(X, dist, 6.0, 4.0);

	if constexpr (true) {
		auto RipsDgm = Rips(SetDgm, X, dist, ripsvals, 2);

		// diagram in Chain
		std::cout << "Chain functor" << std::endl;
		auto ChainDgm = ChainFunctor(RipsDgm, ModP<int, 2>());

		// diagram in Homology
		std::cout << "Homology functor" << std::endl;
		auto HkDgm = Hom(ChainDgm, 1);
		std::cout << HkDgm.nnode() << "," << HkDgm.nedge() << std::endl;

		HkDgm.edge_data(1000).print();

		std::cout << "barcode" << std::endl;
		auto ps = barcode_sparse_rightleft(HkDgm, 1);
		std::cout << ps.size() << std::endl;

		std::cout << "success" << std::endl;
	}
	if constexpr (true) {
		// diagram in Homology
		std::cout << "RipsHom functor" << std::endl;
		auto HkDgm = RipsHom(SetDgm, X, dist, ripsvals, 1, ModP<int, 2>());
		std::cout << HkDgm.nnode() << "," << HkDgm.nedge() << std::endl;

		HkDgm.edge_data(1000).print();

		std::cout << "barcode" << std::endl;
		auto ps = barcode_sparse_rightleft(HkDgm, 1);
		std::cout << ps.size() << std::endl;

		std::cout << "success" << std::endl;
	}


	return 0;
}
