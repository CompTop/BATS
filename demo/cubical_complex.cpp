#include <iostream>
#include <vector>

#include <bats.h>

int main() {

	std::cout << "\nCreating cubical complex" << std::endl;

	CubicalComplex X(3); // 2 is max dimension
	std::cout << "\ndeclared complex" << std::endl;
	X.add_toplex({0,1,0,1,0,1});
	std::cout << "\nadded toplex" << std::endl;

	X.print_summary();

	return 0;
}
