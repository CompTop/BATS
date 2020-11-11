#pragma once

/*
function to produce a cubical map
*/

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
// Assumes that X \subseteq Y
CellularMap CubicalMap(
	const CubicalComplex &X,
	const CubicalComplex &Y
) {
	std::vector<size_t> s; // cell
	size_t maxd = X.maxdim();
	CellularMap f(maxd);
	for (size_t k = 0; k < maxd+1; k++) {
		std::vector<SparseVector<int, size_t>> col;
		for (size_t i = 0; i < X.ncells(k); i++) {
			// fill s with image of simplex i in dimension k
			s.clear();
			for (auto it = X.cell_begin(k, i); it != X.cell_end(k, i); ++it) {
				s.emplace_back(*it);
			}
			// for (auto si : s) {
			// 	std::cout << si << ',';
			// }
			// std::cout << std::endl;
			int sgn = 1;
			size_t j = Y.find_idx(s);
			if (j == bats::NO_IND) {
				throw std::out_of_range("No existing target cell!");
			}
			col.emplace_back(SparseVector<int, size_t>({j}, {sgn}));
		}
		f[k] = ColumnMatrix<SparseVector<int, size_t>>(Y.ncells(k), X.ncells(k), col);
	}

	return f;
}
