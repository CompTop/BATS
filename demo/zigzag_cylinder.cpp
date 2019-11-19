#include <vector>
#include <set>

#include <linalg/field.h>
#include <linalg/sparse_vector.h>
#include <linalg/col_matrix.h>

#include <complex/cell_complex.h>
#include <complex/simplicial_complex.h>
#include <homology/reduction.h>
#include <complex/cell_map.h>
#include <complex/simplicial_map.h>

#include <chain/chain_complex.h>
#include <homology/basis.h>

#include <chain/chain_map.h>
#include <homology/induced_map.h>

#include <topology/metric.h>
#include <topology/cover.h>
#include <topology/data_gen.h>
#include <topology/rips.h>
#include <topology/inclusion.h>

#include <multigraph/diagram.h>

#define FT ModP<int, 3>
#define VT SparseVector<FT>
#define MT ColumnMatrix<VT>

int main() {

	// generate cylinder dataset
	auto x = gen_cylinder(30, 20);

	/// project onto first coordinate
	auto p = coordinate_projection(x, 3, 0);

	// create open cover using projection
	auto cover = uniform_interval_cover(p, 10);

	// Diagram of Sets and inclusions
	Diagram<std::set<size_t>, std::vector<size_t>> SetDgm;
	for (size_t i = 0; i < cover.size(); i++) {
	    auto i1 = SetDgm.add_node(cover[i]);
	    if (i > 0) {
	        // add backward arrow
	        SetDgm.add_edge(i1, i1-1, vertex_inclusion_map(SetDgm.node[i1], SetDgm.node[i1-1]));
	    }
	    if (i == (cover.size() - 1)) { break; }
	    auto i2 = SetDgm.add_node(set_union(cover[i], cover[i+1]));
	    // map from i1 to i2
	    SetDgm.add_edge(i1, i2, vertex_inclusion_map(SetDgm.node[i1], SetDgm.node[i2]));
	}

	// Diagram of Rips complexes and inclusion maps
	size_t n = SetDgm.nnode();
	size_t m = SetDgm.nedge();

	Diagram<SimplicialComplex, CellularMap> TopDgm(n, m);

	#pragma omp parallel for
	for (size_t i = 0; i < n; i++) {
	    auto xi = get_subset(x, 3, SetDgm.node[i]);
	    TopDgm.set_node(i, RipsComplex(xi, 3, 1.0, 2));
	}

	#pragma omp parallel for
	for (size_t i = 0; i < m; i++) {
	    auto s = SetDgm.elist[i].src;
	    auto t = SetDgm.elist[i].targ;
		    TopDgm.set_edge(i, s, t,
		                    SimplicialMap(TopDgm.node[s] , TopDgm.node[t], SetDgm.edata[i])
		                   );
	}

	// Diagram of vector spaces and induced maps
	Diagram<ReducedChainComplex<MT>, MT> HkDgm(n, m);
	size_t k = 1; // homology dimension

	#pragma omp parallel for
	for (size_t i = 0; i < n; i++) {
	    auto CCX = ChainComplex<MT>(TopDgm.node[i]);
	    HkDgm.set_node(i, ReducedChainComplex(CCX));
	}

	#pragma omp parallel for
	for (size_t i = 0; i < m; i++) {
	    auto s = TopDgm.elist[i].src;
	    auto t = TopDgm.elist[i].targ;
	    auto F = ChainMap<MT>(TopDgm.edata[i]);
	    HkDgm.set_edge(i, s, t,
	                   induced_map(F, HkDgm.node[s], HkDgm.node[t], k)
	                   );
	}

	for (auto M : HkDgm.edata) {
	    M.print();
	}

	return 0;
}
