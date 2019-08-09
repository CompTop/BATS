#pragma once
// compute homology over a filtration

#include <linalg/field.h>
#include <linalg/sparse_vector.h>
#include <linalg/col_matrix.h>

#include <filtration/filtration.h>

#include <vector>

#define FF ModP<int, 3>
#define VecT SparseVector<size_t, FF>
#define MatT ColumnMatrix<VecT>


void extract_barcode(
    MatT& Rk,
    std::map<size_t, size_t>& p2ck1,
    std::vector<size_t>& inf_bars,
    std::vector<size_t>& births,
    std::vector<size_t>& deaths
) {
    // clear vectors
    inf_bars.clear();
    births.clear();
    deaths.clear();

    for (size_t j = 0; j < Rk.width(); j++) {
        if (Rk[j].nnz() == 0) {
            // homology born
            if (p2ck1.count(j) > 0) {
                births.emplace_back(j);
                deaths.emplace_back(p2ck1[j]); // killed by k+1 simplex
            } else {
                inf_bars.emplace_back(j);
            }
        }
    }

    return;
}



// template over type of object
// only needs to have method to get boundary
// compute barcode in dimension k
template <typename T>
void compute_barcode(
    T& F,
    size_t k,
    std::vector<size_t>& inf_bars,
    std::vector<size_t>& births,
    std::vector<size_t>& deaths
) {
    auto Bk = F.template boundary<VecT>(k);
    auto Bk1 = F.template boundary<VecT>(k+1);

    auto p2c = reduce_matrix(Bk);
    auto p2c1 = reduce_matrix(Bk1);

    extract_barcode(Bk, p2c1, inf_bars, births, deaths);
    return;

}
