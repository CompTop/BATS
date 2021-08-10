#include <bats.hpp>
using namespace bats;

int main() {

	size_t n = 1000;
	auto X = sample_sphere<double>(2, n);

	auto dist = Euclidean(); //AngleDist();

	std::cout << "constructing zigzag" << std::endl;
	auto [RipsDgm, ripsvals] = DiscreteMorozovZigzag(X, dist, 10.0, 2);

	// diagram in Chain
	std::cout << "Chain functor" << std::endl;
	auto ChainDgm = __Chain(RipsDgm, ModP<int, 2>());

	// diagram in Homology
	std::cout << "Homology functor" << std::endl;
	auto HkDgm = Hom(ChainDgm, 1);
	std::cout << HkDgm.nnode() << "," << HkDgm.nedge() << std::endl;

	std::cout << "success" << std::endl;

	return 0;
}
