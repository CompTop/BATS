#pragma once

#include <map>
#include <unordered_map>
#include <linalg/sparse_vector.h>
#include <linalg/col_matrix.h>

#define p2c_type std::unordered_map<size_t, size_t>

// perform reduction algorithm on a column matrix in-place
template <class TVec>
p2c_type reduce_matrix(ColumnMatrix<TVec> &M) {

	p2c_type pivot_to_col;

	// loop over columns
	for (size_t j = 0; j < M.ncol(); j++) {
		while(M[j].nnz() > 0) {
			// std::cout << j << " : ";
			// M[j].print_row();
			// piv is index-value nzpair
			auto piv = M[j].lastnz();
			if (pivot_to_col.count(piv.ind) > 0) {
				// eliminate pivot
				size_t k = pivot_to_col[piv.ind];
				auto a = piv.val / M[k].lastnz().val;
				M[j].axpy(-a, M[k]);
			} else {
				// new pivot
				pivot_to_col[piv.ind] = j;
				break;
			}
		}
	}
	return pivot_to_col;
}

// // perform reduction algorithm on a column matrix in-place
// // Maintains invariant M * U^{-1}
// // i.e if M <- M * A then U^{-1} <- A^{-1} U^{-1} = (U * A)^{-1}
// // thus Updates are M <- M * A, U <- U * A
// template <class TVec>
// std::map<size_t, size_t> reduce_matrix(ColumnMatrix<TVec> &M, ColumnMatrix<TVec> &U) {
//
// 	p2c_type pivot_to_col;
//
//   // loop over columns
// 	for (size_t j = 0; j < M.width(); j++) {
// 		while(M[j].nnz() > 0) {
// 			// std::cout << j << " : ";
// 			// M[j].print_row();
// 			// piv is index-value pair
// 			auto piv = M[j].last();
// 			if (pivot_to_col.count(piv.first) > 0) {
// 				size_t k = pivot_to_col[piv.first];
// 				// get coefficient
// 				auto pivk = M[k].last();
// 				auto alpha = - piv.second / pivk.second;
// 				// std::cout << "alpha = " << alpha << std::endl;
// 				M[j].axpy(alpha, M[k]);
// 				U[j].axpy(alpha, U[k]);
// 			} else {
// 				pivot_to_col[piv.first] = j;
// 				break;
// 			}
// 		}
// 	}
// 	return pivot_to_col;
// }



  /*
Function to reduce boundary matrix.
INPUT:
	B - boundary matrix
	pivot_to_col - map from pivot to column
OUTPUT:
	none - inputs are modified in-place.
*/
// void homology_reduction_alg(std::vector<SparseF2Vec<int>> &B, std::map<int, int> &pivot_to_col) {
// 	// loop over columns of boundary matrix
// 	for (size_t j = 0; j < B.size(); j++) {
// 		// if nnz = 0, the reduction is complete
// 		while (B[j].nnz() > 0) {
// 			int piv = B[j].last();
// 			if (pivot_to_col.count(piv) > 0) {
// 				int k = pivot_to_col[piv];
// 				// there is a column with that pivot
// 				B[j].add(B[k]);
// 			} else {
// 				// there is no column with that pivot
// 				pivot_to_col[piv] = j;
// 				break;
// 			}
// 		} // end column reduction
// 	} // end for loop
// 	return;
// }
