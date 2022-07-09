#include <iostream>
#include <memory>
#include "vineyard.hpp"
#include <bats.hpp>
#include <chrono>
// C++ implementation of the approach
#include <bits/stdc++.h>


#define F ModP<int, 5>

// perform reduction algorithm on a column matrix in-place
// apply change of basis to U
// invariant M * U = R
#define p2c_type std::vector<size_t>

template <class TMat>
p2c_type reduce_matrix_standard(TMat &M, TMat &U) {

	if (M.ncol() != U.ncol()) {throw std::runtime_error("Number of columns are not the same!");}

	// p2c_type pivot_to_col;
	p2c_type pivot_to_col(M.nrow(), bats::NO_IND);
	// create a temporary vector for use in axpys
	// typename TVec::tmp_type tmp;

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
				// M[j].axpy(-a, M[k], tmp);
				// U[j].axpy(-a, U[k], tmp); // update change of basis
                M[j].axpy(-a, M[k]);
				U[j].axpy(-a, U[k]); // update change of basis
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


int main()
{   
    using MatT = VineyardMatrix<F>;

    size_t n_rows = 4;
    size_t n_cols = 4;
    auto B = MatT::random(n_rows, n_cols, 0.4, 4);
    MatT R(B);
    auto U = MatT::identity(n_cols);
    auto p2c = reduce_matrix_standard(B, U);
    

    // A.print();
    // B.print();
    // bats::reduce_matrix(B, U);

    return 0;
}