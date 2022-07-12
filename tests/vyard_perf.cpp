#include <iostream>
#include <memory>
#include "vineyard.hpp"
#include <bats.hpp>
#include <util/io.hpp>
#include <chrono>
// C++ implementation of the approach
#include <bits/stdc++.h>


#define F ModP<int, 5>

// number of seeds for random checks
#define N_SEEDS 4


// Function to return the next random number
int getNum(vector<int>& v)
{
 
    // Size of the vector
    int n = v.size();
 
    // Generate a random number
    std::srand(time(NULL));
 
    // Make sure the number is within
    // the index range
    int index = rand() % n;
 
    // Get random number from the vector
    int num = v[index];
 
    // Remove the number from the vector
    std::swap(v[index], v[n - 1]);
    v.pop_back();
 
    // Return the removed number
    return num;
}
 
// Function to generate n non-repeating random numbers
std::vector<size_t> generateRandom(int n)
{
    std::vector<int> v(n);
 
    // Fill the vector with the values
    // 1, 2, 3, ..., n
    for (int i = 0; i < n; i++)
        v[i] = i;
    
    // While vector has elements
    // get a random number from the vector
    std::vector<size_t> perm;
    perm.reserve(n);
    while (v.size()) {
        perm.emplace_back(getNum(v));
    }
    return perm;
}

#define p2c_type std::vector<size_t>

template <class TVec>
p2c_type reduce_matrix_standard(ColumnMatrix<TVec> &M, ColumnMatrix<TVec> &U) {

	if (M.ncol() != U.ncol()) {throw std::runtime_error("Number of columns are not the same!");}

	// p2c_type pivot_to_col;
	p2c_type pivot_to_col(M.nrow(), bats::NO_IND);
	// create a temporary vector for use in axpys
	typename TVec::tmp_type tmp;

	// loop over columns
	// size_t ct = 0;
	for (size_t j = 0; j < M.ncol(); j++) {
		while(M[j].nnz() > 0) {
			// ++ct;
			// std::cout << j << " : ";
			// M[j].print_row();
			// piv is index-value nzpair
			auto piv = M[j].lastnz();
			// if (pivot_to_col.count(piv.ind) > 0) {
			if (pivot_to_col[piv.ind] != bats::NO_IND) {
				// eliminate pivot
				size_t k = pivot_to_col[piv.ind];
				auto a = piv.val / M[k].lastnz().val;
				M[j].axpy(-a, M[k], tmp);
				U[j].axpy(-a, U[k], tmp); // update change of basis
			} else {
				// new pivot
				pivot_to_col[piv.ind] = j;
				break;
			}
		}
	}
	// std::cout << "# iterations = " << ct << std::endl;
	return pivot_to_col;
}

// Vineyard Version
template <class TMat>
p2c_type reduce_matrix_standard(TMat &M, TMat &U) {

	if (M.ncol() != U.ncol()) {throw std::runtime_error("Number of columns are not the same!");}

	// p2c_type pivot_to_col;
	p2c_type pivot_to_col(M.nrow(), bats::NO_IND);
	// create a temporary vector for use in axpys
	// typename TVec::tmp_type tmp;
    typename TMat::col_type tmp;

	// loop over columns
	// size_t ct = 0;
	for (size_t j = 0; j < M.ncol(); j++) {
		while(M[j].nnz() > 0) {
			// ++ct;
			// std::cout << j << " : ";
			// M[j].print_row();
			// piv is index-value nzpair
			auto piv = M[j].lastnz();
			// if (pivot_to_col.count(piv.ind) > 0) {
			if (pivot_to_col[piv.ind] != bats::NO_IND) {
				// eliminate pivot
				size_t k = pivot_to_col[piv.ind];
				auto a = piv.val / M[k].lastnz().val;
				M[j].axpy(-a, M[k], tmp);
				U[j].axpy(-a, U[k], tmp); // update change of basis
			} else {
				// new pivot
				pivot_to_col[piv.ind] = j;
				break;
			}
		}
	}
	// std::cout << "# iterations = " << ct << std::endl;
	return pivot_to_col;
}

int main(int argc, char* argv[])
{   
    using MatT2 = VineyardMatrix<F>;
    using VT = SparseVector<F, size_t>;
	using MatT = ColumnMatrix<VT>;

    size_t n_rows = bats::util::io::parse_argv(argc, argv, "-m", 150);
    size_t n_cols = bats::util::io::parse_argv(argc, argv, "-n", 150);
   
    // size_t n_rows = 50;
    // size_t n_cols = 50;
    auto A = MatT::random(n_rows, n_cols, 0.4, 4);
    auto B = MatT2(n_rows, n_cols, A.cols());
    // A.print();
    // B.print();

    // Test permutation performance
    {
        auto perm = generateRandom(n_rows);
        // bats::print_1D_vectors(perm);
        auto start = std::chrono::steady_clock::now();
        A.permute_rows(perm);
        auto end = std::chrono::steady_clock::now();
        std::cout << "Permute rows of standard sparse matrix: "
                << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
                << "ms" << std::endl;
        
        start = std::chrono::steady_clock::now();
        B.permute_rows(perm);
        end = std::chrono::steady_clock::now();
        std::cout << "Permute rows of Vineyard matrix: "
                << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
                << "ms" << std::endl;
    }

    // Test Reduction Performance
    {
        auto start = std::chrono::steady_clock::now();
        auto U_A = MatT::identity(n_cols);
        auto p2cA = reduce_matrix_standard(A, U_A);
        auto end = std::chrono::steady_clock::now();
        std::cout << "Reduce standard sparse matrix: "
                << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
                << "ms" << std::endl;
        // std::cout << "A is" << std::endl;
        // A.print();
    }
    {
        auto start = std::chrono::steady_clock::now();
        auto U_B = MatT2::identity(n_cols);
        auto p2cB = reduce_matrix_standard(B, U_B);
        auto end = std::chrono::steady_clock::now();
        std::cout << "Reduce Vineyard matrix: "
                << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
                << "ms" << std::endl;
        // std::cout << "B is" << std::endl;
        // B.print();
    }
    return 0;
}