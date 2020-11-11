#pragma once
/*
Utility to dump a multigraph of vector spaces into something we can operate on
*/
#include <vector>

#include <multigraph/diagram.hpp>
#include <linalg/naive_dense.hpp>
#include <linalg/col_matrix.hpp>
#include "type_a.hpp"
#include <stdexcept>

namespace bats {

// return tuple of
//   1. vect of matrix data
//   2. vect of matrix wrappers
//   3. vect of etype (bool) = 0 = <-, 1 = ->
// template over node type and column type
template <typename NT, typename CT>
auto A_type_rep(Diagram<NT, ColumnMatrix<CT>> &D) {
    using FT = typename CT::val_type;
    using AD = A<Dense<FT,ColMaj>>;

    size_t m = D.nedge();
    // TODO: check m == n-1

    std::vector<FT*> data(m);
    std::vector<AD> mat(m);
    std::vector<bool> etype(m);

    // assume nodes and edges are sorted for a traversal
    // TODO: check this
    size_t i = 0;
    for (size_t k = 0; k < m; k++) {
        size_t mk = D.edata[k].nrow();
        size_t nk = D.edata[k].ncol();
        data[k] = D.edata[k].dump_dense();
        mat[k] = AD(mk, nk, data[k]);
        if (i == D.elist[k].src) {
            etype[k] = true;
        } else if (i == D.elist[k].targ) {
            etype[k] = false;
        } else {
            // throw error
        }
        i++;
    }

    return std::make_tuple(data, mat, etype);
}


template <typename MT>
std::vector<PersistencePair<size_t>> barcode(
    Diagram<ReducedChainComplex<MT>, MT> &D,
	size_t hdim
) {
    using FT = typename MT::val_type;

    auto [data, mat, etype] = A_type_rep(D);

    auto taq = Type_A<FT>(mat,etype);

    taq.create_copy_of_mats();

	taq.forward_sweep();
	taq.backward_sweep();

	if (!taq.is_consistent()) {
		throw std::runtime_error("quiver fact is not consistent!");
	}

    return taq.barcode_pairs(hdim);
}

} // namespace bats
