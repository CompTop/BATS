#include <vector>
#include <set>

#include <bats.h>

#define FT ModP<int, 3>
#define VT SparseVector<FT>
#define MT ColumnMatrix<VT>

int main() {

	// generate cylinder dataset

	// cylinder example
	// size_t d = 3;
	// double rmax = 1.0;
	// auto x = gen_cylinder(40, 30);
	// add_normal_noise(x, 0.0, 0.05);

	// random cube example
	size_t d = 3; // dimension
	double rmax = 0.15;
	auto x = sample_cube<double>(d, 1000);

	/// project onto first coordinate
	auto p = coordinate_projection(x, d, 0);

	// create open cover using projection
	auto cover = uniform_interval_cover(p, 20);

	// Diagram of Sets and inclusions
	auto SetDgm = linear_cover_union_diagram(cover);
	// auto SetDgm = linear_cover_intersection_diagram(cover);

	// Diagram of Spaces and maps
	auto TopDgm = Rips(SetDgm, x, d, rmax, 2);

	// diagram in Chain
	auto ChainDgm = Chain<MT>(TopDgm);

	// diagram in Homology
	auto HkDgm = Hom(ChainDgm, 1);

	for (auto M : HkDgm.edata) {
	    M.print();
	}

	// dump into Atype rep
	auto [data, mat, etype] = A_type_rep(HkDgm);

	for (size_t k = 0; k < mat.size(); k++) {
		std::cout << "*" << std::endl;
		std::cout << (etype[k] ? "v " : "^ ");
		mat[k].print();
	}
	std::cout << "*" << std::endl;

	return 0;
}
