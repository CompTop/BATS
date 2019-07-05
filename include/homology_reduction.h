#pragma once

// perform reduction algorithm on a column matrix in-place
template <class TVec>
void reduce_matrix(ColumnMatrix<TVec> M) {

	std::map<size_t, size_t> pivot_to_col();
	
  // get pivots
	for (size_t j = 0; j < M.width(); j++) {
		while(M[j].nnz() > 0) {
			// piv is index-value pair
			auto piv = M[j].last();
			if (pivot_to_col.count(piv.first) > 0) {
				size_t k = pivot_to_col[piv.first];
				// get coefficient
				auto pivk = M[k].last;
				auto alpha = piv.second / pivk.second;
				M[j].axpy(alpha, B[k])
			} else {
				pivot_to_col[piv.first] = j;
				break;
			}
		}
	}
	return;
}

  /*
Function to reduce boundary matrix.
INPUT:
	B - boundary matrix
	pivot_to_col - map from pivot to column
OUTPUT:
	none - inputs are modified in-place.
*/
void homology_reduction_alg(std::vector<SparseF2Vec<int>> &B, std::map<int, int> &pivot_to_col) {
	// loop over columns of boundary matrix
	for (size_t j = 0; j < B.size(); j++) {
		// if nnz = 0, the reduction is complete
		while (B[j].nnz() > 0) {
			int piv = B[j].last();
			if (pivot_to_col.count(piv) > 0) {
				int k = pivot_to_col[piv];
				// there is a column with that pivot
				B[j].add(B[k]);
			} else {
				// there is no column with that pivot
				pivot_to_col[piv] = j;
				break;
			}
		} // end column reduction
	} // end for loop
	return;
}
