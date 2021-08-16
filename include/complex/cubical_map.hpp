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
#include <util/simplex.hpp>
#include "cell_map.hpp"
#include "simplicial_complex.hpp"

#include <linalg/sparse_vector.hpp>
#include <linalg/col_matrix.hpp>

namespace bats {

// construct simplicial map from X to Y
// Assumes that X \subseteq Y
inline CellularMap CubicalMap(
	const CubicalComplex &X,
	const CubicalComplex &Y
) {
	return CellularMap(X, Y);
}

} // namespace bats
