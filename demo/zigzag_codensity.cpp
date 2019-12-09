#include <vector>
#include <set>

#include <bats.h>
#include <util/set.h>

#define FT ModP<int, 2>
// #define FT Rational<int>
#define VT SparseVector<FT>
#define MT ColumnMatrix<VT>

int main() {


	// sample circle
	size_t d = 2; // 2d points
	size_t n = 1000; // 100 points
	auto dist = Euclidean(); // Euclidean distance

	double rmax = 0.8;

	// generate data
	auto X = sample_sphere<double>(d, n);

	// add some noise
	add_normal_noise(X, 0.0, 0.05);

	std::vector<size_t> ks = {5, 10, 15, 20};
	std::vector<double> ps = {0.2, 0.2, 0.2, 0.2};


	std::vector<std::set<size_t>> subset;
	for (size_t i =0; i < ks.size(); i++) {
		size_t k = ks[i];
		double p = ps[i];
		auto dk = kdist(X, dist, k);
		subset.emplace_back(
			to_set(bats::top_p(dk, p))
		); // get top p percent of k-nn estimator
	}

	// Diagram of Sets and inclusions
	auto SetDgm = linear_subset_union_diagram(subset);

	// Diagram of Spaces and maps
	auto TopDgm = Rips(SetDgm, X, dist, rmax, d);

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

	(void) data; // to remove warnings

	// create quiver, factorize and check consistency
	auto taq = Type_A<FT>(mat,etype);
	taq.create_copy_of_mats();

	taq.forward_sweep();
	taq.backward_sweep();

	// print E mats
	for(size_t i=0;i<taq.n;i++){
		std::cout<<"Arrow Dir :"<<taq.arrow_dir[i]<<"\n";
		if(taq.arrow_dir[i]==0){
			taq.ELmats[i].print();
		}else{
			taq.ELHmats[i].print();
		}
	}

	//print barcodes
	std::cout<<"Barcodes:\n\n";
	taq.print_barcodes();

	std::cout<<"Compressed Barcodes:\n\n";
	taq.print_barcodes_compressed();

	bool cons =  taq.is_consistent();

	std::cout<<"Quiver factorization is"<<( cons?"":" NOT" )<<" consistent !\n\n";

	return 0;
	}
