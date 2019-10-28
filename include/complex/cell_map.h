#pragma once
/*
Cellular map f: X->Y
data is how cells in X map to cells in Y
essentially a chain map
*/

#include <vector>
#include <linalg/sparse_vector.h>
#include <linalg/col_matrix.h>

class CellularMap
{
private:
	using vec_type = SparseVector<int, size_t>;
	using map_type = ColumnMatrix<vec_type>;

	// store image of k-cells as linear combination of k-cells
	std::vector<map_type> cell_map;

	inline void _resize(size_t dim) { cell_map.resize(dim); }
	inline void _safe_resize(size_t dim) { if (dim >= cell_map.size()) { _resize(dim); } }

public:

	CellularMap() {}

	// allocate for dim dimensions
	CellularMap(size_t dim) {
		cell_map.resize(dim+1);
	}

	inline size_t maxdim() const { return cell_map.size() - 1; }

	map_type& operator[](size_t k) {
		_safe_resize(k);
		return cell_map[k];
	}

};
