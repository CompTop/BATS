#include <vector>
#include <set>
#include <string>

#include <bats.h>
#include <util/set.h>

#define FT ModP<int, 2>
// #define FT Rational<int>
#define VT SparseVector<FT>
#define MT ColumnMatrix<VT>

int main() {

	size_t nsets = 8;
	size_t ns = 100; // number of points in each subset

	// sample circle
	size_t d = 2; // 2d points
	size_t n = 1000; // 100 points
	auto dist = Euclidean(); // Euclidean distance
	double rmax = 0.5;

	// generate data
	auto X = sample_sphere<double>(d, n);

	// landmark sets
	std::vector<std::set<size_t>> subset;
	for (size_t i =0; i < nsets; i++) {
		subset.emplace_back(random_subset(n, ns));
	}

	// Diagram of Sets and inclusions
	auto SetDgm = linear_subset_union_diagram(subset);

	// Diagram of Spaces and maps
	auto TopDgm = Rips(SetDgm, X, dist, rmax, d);

	std::string dname = "rss.dgm";
	TopDgm.save(dname);

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
