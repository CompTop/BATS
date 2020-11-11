#pragma once
// compute homology over a filtration

#include <util/sorted.hpp>


#include <linalg/field.hpp>
#include <linalg/csc_matrix.hpp>
#include <linalg/sparse_vector.hpp>
#include <linalg/col_matrix.hpp>
#include <linalg/schur.hpp>

#include <morse/pairing.hpp>
#include <filtration/filtration.hpp>
#include "homology_reduction.hpp"

#include <vector>
#include <iostream>
#include <limits>

// #define FF ModP<int, 3>
// #define VecT SparseVector<size_t, FF>
// #define MatT ColumnMatrix<VecT>

template <typename T, class CpxT, typename FT>
std::vector<std::vector<T>> standard_reduce(const Filtration<T, CpxT> &F, FT) {

    using CT = SparseVector<FT>;
    using MatT = ColumnMatrix<CT>;

    size_t maxd = F.maxdim(); // maxdim of filtration
    auto M = F.pairing();
    // Copy the pairing
    // this will have pairs added as reduction occurs
    MorsePairing M2(M);
    M2.clear();
    // for holding p2cs
    p2c_type p2c;
    std::vector<size_t> perm_row, perm_col;
    //std::vector<size_t> iperm_row, iperm_col;
    perm_row = F.sortperm(0);
    // perm_col = F.sortperm(maxd);
    //iperm_row = bats::sortperm(perm_row);
    for (size_t d = 1; d < maxd+1; d++) {
    // for (size_t d = maxd; d > 0; d--) {
        //std::cout << "processing d = " << d << std::endl;
        auto B = M.boundary_csc(d);
        perm_col = F.sortperm(d);
        // perm_row = F.sortperm(d-1);
        //iperm_col = bats::sortperm(perm_col);
        MatT BC(B);
        // permute rows and columns
        BC.permute(perm_row, perm_col);
        p2c = reduce_matrix(BC);

        // put pairs into pairing in original order
        for (auto& [k, v] : p2c) {
            if (!M2.set_pair(d-1, perm_row[k], perm_col[v])) {
                // this shouldn't happen. If it does, there's a problem :(
                auto vald1 = F.get_val(d-1);
                auto vald = F.get_val(d);
                std::cout << "!!!failed to set pair " << perm_row[k] << ',' << perm_col[v];
                std::cout << "  in dimension " << d-1 << std::endl;
                std::cout << "  filtration times " << vald1[perm_row[k]] << ',' << vald[perm_col[v]] << std::endl;
            }
        }

        // in next dimension, column perm will be row perm.
        perm_row = perm_col;
        // perm_col = perm_row;
        //iperm_row = iperm_col;
    }

    // finally, we put all the completed pairs into a barcode
    std::vector<std::vector<T>> bars(maxd);
    for (size_t d = 1; d < maxd + 1; d++) {
        // fill bars by accessing filtration in original order
        auto up = M2.up_paired(d-1);
        auto down = M2.down_paired(d);
        auto unp = M2.unpaired(d-1);
        auto vald1 = F.get_val(d-1);
        auto vald = F.get_val(d);
        bars[d-1].reserve(2*(up.size() + unp.size()));
        for (size_t k = 0; k < up.size(); k++) {
            bars[d-1].emplace_back(vald1[up[k]]);
            bars[d-1].emplace_back(vald[down[k]]);
        }
        for (size_t k = 0; k < unp.size(); k++) {
            bars[d-1].emplace_back(vald1[unp[k]]);
            bars[d-1].emplace_back(std::numeric_limits<T>::infinity());
        }
    }

    return bars;
}

// compute block matrices in dimension k
template <typename T, class CpxT, typename FT>
void complete_pairs(
    Filtration<T, CpxT> &F,
    MorsePairing<CpxT> &M,
    size_t d,
    FT
) {
    using CSCM = CSCMatrix<int, size_t>; // csc matrix
    using CT = SparseVector<FT>; // column type
    using CM = ColumnMatrix<CT>; // column matrix

    std::vector<size_t> rind1, rind2, cind1, cind2, p1, p2;

    // std::cout << "getting permuations" << std::endl;
    rind1 = M.up_paired(d-1); // row indices that are paired with a column index
    rind2 = M.unpaired(d-1);  // row indices that are unpaired
    cind1 = M.down_paired(d); // column indices that are paired with a row index
    cind2 = M.unpaired(d);    // column indices that are unparied
    p1 = F.sortperm(d-1);
    p2 = F.sortperm(d);

    // std::cout << "sorting" << std::endl;
    sort_ind_pair_by_perm(rind1, cind1, p1); // this ensures that the rind1, cind1 block will be upper triangular
    sort_ind_by_perm(rind2, p1); // this ensures that rind2 will appear in filtration order
    sort_ind_by_perm(cind2, p2); // this ensures that cind2 will appear in filtration order

    // std::cout << "getting partial perms" << std::endl;
    std::vector<size_t> prind1 = partial_perm(rind1, F.ncells(d-1)); // partial permutation for rind1
    std::vector<size_t> prind2 = partial_perm(rind2, F.ncells(d-1)); // partial permutation for rind2

    CSCM Bk = M.boundary_csc(d);

    CSCM A, B, C, D;
    // std::cout << "block select 1" << std::endl;
    block_select(Bk, cind1,
        {&prind1, &prind2},
        {rind1.size(), rind2.size()},
        {&A, &C});
    // std::cout << "block select 2" << std::endl;
    block_select(Bk, cind2,
        {&prind1, &prind2},
        {rind1.size(), rind2.size()},
        {&B, &D});


    CM A1(A), B1(B), C1(C), D1(D);

    // get schur complement
    CM S = schur(A1, B1, C1, D1);
    //MatT S = schur(A, B, C, D);


    // run standard reduction alg on S
    auto p2c = reduce_matrix(S); // segfault here

    // put pairs into pairing in original order
    for (auto& [k, v] : p2c) {
        //M.set_pair(d-1, rind2[k], cind2[v]);
        if (!M.set_pair(d-1, rind2[k], cind2[v])) {
            // this shouldn't happen.  If it does, there's a problem :(
            std::cout << "!!!failed to set pair " << rind2[k] << ',' << cind2[v];
            std::cout << "  in dimension " << d-1 << std::endl;
        }
    }

    return;

}




// use pre-pairs in M to complete reduction over field FT
template <typename T, class CpxT, typename FT>
std::vector<std::vector<T>> block_reduce(Filtration<T, CpxT> &F, FT) {

    auto M = F.pairing();
    // Copy the pairing
    // this will have pairs added as reduction occurs
    MorsePairing M2(M);
    //std::vector<std::vector<T>> bars(maxd);

    size_t maxd = F.maxdim(); // maxdim of filtration
    // we assume that 0-1 pairs have been found already using union find.
    for (size_t d = maxd; d > 1; d--) {
        // complete the pairs between dims d and d-1
        complete_pairs(F, M2, d, FT());
    }

    // finally, we put all the completed pairs into a barcode
    std::vector<std::vector<T>> bars(maxd);
    for (size_t d = 1; d < maxd + 1; d++) {
        // fill bars by accessing filtration in original order
        auto up = M2.up_paired(d-1);
        auto down = M2.down_paired(d);
        auto unp = M2.unpaired(d-1);
        auto vald1 = F.get_val(d-1);
        auto vald = F.get_val(d);
        bars[d-1].reserve(2*(up.size() + unp.size()));
        for (size_t k = 0; k < up.size(); k++) {
            bars[d-1].emplace_back(vald1[up[k]]);
            bars[d-1].emplace_back(vald[down[k]]);
        }
        for (size_t k = 0; k < unp.size(); k++) {
            bars[d-1].emplace_back(vald1[unp[k]]);
            bars[d-1].emplace_back(std::numeric_limits<T>::infinity());
        }
    }

    return bars;
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
