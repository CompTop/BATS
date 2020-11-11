#pragma once
/*
Create a cell map from simplicial map data

i.e. we take a map on 0-cells
and extend to a map on all simplices
*/
#include <iostream>
#include <vector>
#include <stdexcept>
#include <util/simplex.h>
#include "cell_map.h"
#include "simplicial_complex.h"

#include <linalg/sparse_vector.h>
#include <linalg/col_matrix.h>

// construct simplicial map from X to Y
// send vertex set of X to vertex set of Y using f0
CellularMap SimplicialMap(
	const SimplicialComplex &X,
	const SimplicialComplex &Y,
	const std::vector<size_t> &f0
) {
	std::vector<size_t> s; // simplex
	size_t maxd = X.maxdim();
	CellularMap f(maxd);
	for (size_t k = 0; k < maxd+1; k++) {
		std::vector<SparseVector<int, size_t>> col;
		for (size_t i = 0; i < X.ncells(k); i++) {
			// fill s with image of simplex i in dimension k
			s.clear();
			for (auto it = X.simplex_begin(k, i); it != X.simplex_end(k, i); ++it) {
				s.emplace_back(f0[*it]);
			}
			// for (auto si : s) {
			// 	std::cout << si << ',';
			// }
			// std::cout << std::endl;
			int sgn = simplex_sign(s);
			if (sgn != 0) {
				size_t j = Y.find_idx(s);
				if (j == bats::NO_IND) {
					throw std::out_of_range("No existing target simplex!");
				}
				col.emplace_back(SparseVector<int, size_t>({j}, {sgn}));
			} else {
				// simplex was degenerate
				col.emplace_back(SparseVector<int, size_t>());
			}
		}
		f[k] = ColumnMatrix<SparseVector<int, size_t>>(Y.ncells(k), X.ncells(k), col);
	}

	return f;
}


// construct simplicial map from X to Y
// no initial data is specified - assume inclusion map
CellularMap SimplicialMap(
	const SimplicialComplex &X,
	const SimplicialComplex &Y
) {
	std::vector<size_t> s; // simplex
	size_t maxd = X.maxdim();
	CellularMap f(maxd);
	for (size_t k = 0; k < maxd+1; k++) {
		std::vector<SparseVector<int, size_t>> col;
		for (size_t i = 0; i < X.ncells(k); i++) {
			// fill s with image of simplex i in dimension k
			s.clear();
			for (auto it = X.simplex_begin(k, i); it != X.simplex_end(k, i); ++it) {
				s.emplace_back(*it);
			}

			size_t j = Y.find_idx(s);
			if (j == bats::NO_IND) {
				throw std::out_of_range("No existing target simplex!");
			}
			col.emplace_back(SparseVector<int, size_t>({j}, {1}));
		}
		f[k] = ColumnMatrix<SparseVector<int, size_t>>(Y.ncells(k), X.ncells(k), col);
	}

	return f;
}
