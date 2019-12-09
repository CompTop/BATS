#include <vector>
#include <set>

#include <bats.h>

#define FT ModP<int, 2>
// #define FT Rational<int>
#define VT SparseVector<FT>
#define MT ColumnMatrix<VT>

// sample nlms landmark sets on n points for DataSet X
template<typename T>
std::vector<DataSet<T>> sample_random_landmarks(
	const DataSet<T> &X,
	const size_t nlms,
	const size_t n
) {
	std::vector<DataSet<T>> LM;
	for (size_t i = 0; i < nlms; i++) {
		LM.emplace_back(random_landmarks(X, n));
	}

	return LM;
}

// get witness neighborhoods for each landmark set
template <typename T, typename M>
std::vector<bats::Cover> witness_neighborhoods(
	const DataSet<T> &X,
	const std::vector<DataSet<T>> &LM,
	const M &dist,
    const size_t nu,
    const T rmax
) {

	std::vector<bats::Cover> cover(LM.size());

	// get neighborhoods
	#pragma omp parallel for
	for (size_t i = 0; i < LM.size(); i++) {
		cover[i] = witness_neighborhoods(X, LM[i], dist, nu, rmax);
	}

	return cover;
}

// create bivariate cover for cover[i], cover[i+1]
Diagram<bats::Cover, std::vector<size_t>> bivariate_cover_diagram(
	std::vector<bats::Cover> &cover
) {
	size_t n = 2 * cover.size() - 1; // number of nodes
	size_t m = n - 1; // number of edges
	Diagram<bats::Cover, std::vector<size_t>> D(n, m);

	// set covers
	for (size_t i = 0; i < cover.size(); i++) {
		D.set_node(2*i, cover[i]);
	}

	for (size_t i =0; i < cover.size() - 1; i++) {
		auto [coverij, fi, fj] = bivariate_cover(cover[i], cover[i+1]);
		D.set_node(2*i + 1, coverij); // set bivariate cover node
		D.set_edge(2*i, 2*i + 1, 2*i, fi); // set map to i
		D.set_edge(2*i+1, 2*i + 1, 2*i+2, fj); // set map to i+1
	}

	return D;
}


int main() {

	size_t nlms = 8;
	size_t nl = 30;

	// sample circle
	size_t d = 2; // 2d points
	size_t n = 2000; // 1000 points
	auto dist = Euclidean(); // Euclidean distance

	// generate data
	auto X = sample_sphere<double>(d, n);

	// get landmarks and covers
	auto LM = sample_random_landmarks(X, nlms, nl);
	auto cover = witness_neighborhoods(X, LM, dist, 2, 0.0);

	// Diagram of Covers
	auto CovDgm = bivariate_cover_diagram(cover);

	// Diagram of Nerves
	auto TopDgm = Nerve(CovDgm, d);

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
