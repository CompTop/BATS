#include <iostream>
#include <vector>

#include <bats.hpp>

int main() {

	std::cout << "\nCreating cubical complex" << std::endl;

	CubicalComplex X(3); // 2 is max dimension
	std::cout << "\ndeclared complex" << std::endl;
	X.add_recursive({1,2,0,1,5,6});
	std::cout << "\nadded toplex" << std::endl;

	X.print_summary();

	auto B1 = X.boundary_csc(1);
	B1.print();

	auto B2 = X.boundary_csc(2);
	B2.print();

	(B1 * B2).print();

	auto B3 = X.boundary_csc(3);
	B3.print();

	(B2 * B3).print();

	auto cubes = X.get_cubes();
	std::cout << cubes.size() << std::endl;
	for (auto c: cubes) {
		std::cout << c.size() << ", ";
	}
	std::cout << '\n';

	CubicalComplex Y = X.skeleton(0); // 2 is max dimension
	std::cout << "Y cells " << Y.ncells() << std::endl;
	Y.add_recursive({1,2,0,1,5,5});

	auto M = CubicalMap(Y, X);
	M[0].print();

	return 0;
}
