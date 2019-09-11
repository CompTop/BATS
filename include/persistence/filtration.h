#pragma once
// compute homology over a filtration

#include <util/sorted.h>


#include <linalg/field.h>
#include <linalg/csc_matrix.h>
#include <linalg/sparse_vector.h>
#include <linalg/col_matrix.h>
#include <linalg/schur.h>

#include <morse/pairing.h>
#include <filtration/filtration.h>
#include "homology_reduction.h"

#include <vector>
#include <iostream>

// #define FF ModP<int, 3>
// #define VecT SparseVector<size_t, FF>
// #define MatT ColumnMatrix<VecT>

template <typename T, class CpxT, typename FT>
void standard_reduce(const Filtration<T, CpxT> &F, FT) {

    using CT = SparseVector<FT>;
    using MatT = ColumnMatrix<CT>;

    size_t maxd = F.maxdim(); // maxdim of filtration
    auto M = F.pairing();
    // for holding p2cs
    std::vector<p2c_type> p2c(maxd);
    std::vector<size_t> permd, permd1;
    permd1 = F.sortperm(0);
    for (size_t d = 1; d < maxd+1; d++) {
        auto B = M.boundary_csc(d);
        permd = F.sortperm(d);
        MatT BC(B);
        // permute rows and columns
        BC.permute(permd1, permd);
        p2c[d-1] = reduce_matrix(BC);
        permd1 = permd;
    }

}

// compute block matrices in dimension k
template <typename T, class CpxT>
void block_reduce(Filtration<T, CpxT> &F, size_t k) {
    using MatT = CSCMatrix<int, size_t>;

    // std::cout << "extracting pairing" << std::endl;
    auto M = F.pairing();
    std::vector<size_t> rind1, rind2, cind1, cind2, p1, p2;

    // std::cout << "getting permuations" << std::endl;
    rind1 = M.up_paired(k-1); // row indices that are paired with a column index
    rind2 = M.unpaired(k-1);  // row indices that are unpaired
    cind1 = M.down_paired(k); // column indices that are paired with a row index
    cind2 = M.unpaired(k);    // column indices that are unparied
    p1 = F.sortperm(k-1);
    p2 = F.sortperm(k);

    // std::cout << "sorting" << std::endl;
    sort_ind_pair_by_perm(rind1, cind1, p1); // this ensures that the rind1, cind1 block will be upper triangular
    sort_ind_by_perm(rind2, p1); // this ensures that rind2 will appear in filtration order
    sort_ind_by_perm(cind2, p2); // this ensures that cind2 will appear in filtration order

    // std::cout << "getting partial perms" << std::endl;
    std::vector<size_t> prind1 = partial_perm(rind1, F.ncells(k-1)); // partial permutation for rind1
    std::vector<size_t> prind2 = partial_perm(rind2, F.ncells(k-1)); // partial permutation for rind2

    MatT Bk = M.boundary_csc(k);

    MatT A, B, C, D;
    // std::cout << "block select 1" << std::endl;
    block_select(Bk, cind1,
        {&prind1, &prind2},
        {cind1.size(), cind2.size()},
        {&A, &C});
    // std::cout << "block select 2" << std::endl;
    block_select(Bk, cind2,
        {&prind1, &prind2},
        {cind1.size(), cind2.size()},
        {&B, &D});

    //MatT S = schur(A, B, C, D);

    return;

}


// void extract_barcode(
//     MatT& Rk,
//     std::map<size_t, size_t>& p2ck1,
//     std::vector<size_t>& inf_bars,
//     std::vector<size_t>& births,
//     std::vector<size_t>& deaths
// ) {
//     // clear vectors
//     inf_bars.clear();
//     births.clear();
//     deaths.clear();
//
//     for (size_t j = 0; j < Rk.width(); j++) {
//         if (Rk[j].nnz() == 0) {
//             // homology born
//             if (p2ck1.count(j) > 0) {
//                 births.emplace_back(j);
//                 deaths.emplace_back(p2ck1[j]); // killed by k+1 simplex
//             } else {
//                 inf_bars.emplace_back(j);
//             }
//         }
//     }
//
//     return;
// }



// template over type of object
// only needs to have method to get boundary
// compute barcode in dimension k
// template <typename T>
// void compute_barcode(
//     T& F,
//     size_t k,
//     std::vector<size_t>& inf_bars,
//     std::vector<size_t>& births,
//     std::vector<size_t>& deaths
// ) {
//     auto Bk = F.template boundary<VecT>(k);
//     auto Bk1 = F.template boundary<VecT>(k+1);
//
//     auto p2c = reduce_matrix(Bk);
//     auto p2c1 = reduce_matrix(Bk1);
//
//     extract_barcode(Bk, p2c1, inf_bars, births, deaths);
//     return;
//
// }
