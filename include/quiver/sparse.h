#pragma once

#include <vector>
#include <map>

#include <multigraph/diagram.h>
#include <linalg/col_matrix.h>
#include <linalg/sparse_fact.h>
#include <persistence/barcode.h>

template <typename Edge>
inline bool is_left_arrow(const Edge &e) {
    return e.targ < e.src;
}

template <typename NT, typename TM>
auto barcode_form_leftward(const Diagram<NT, TM> &dgm) {

    size_t m = dgm.nedge();

    using TC = typename TM::col_type;
    std::vector<SparseFact<TC>> facts(m);

    // copy matrices on edges
    std::vector<TM> mats = dgm.edata;

    // forward sweep
    // left to right
    for (size_t j = 0; j < m; j++) {
        if (is_left_arrow(dgm.elist[j])) {
            // Left arrow = LEUP

            // first apply change of basis from left
            if (j > 0) {
                if (is_left_arrow(dgm.elist[j-1])) {
                    // A = U * P * A
                    mats[j] = facts[j-1].U * facts[j-1].P * mats[j];
                } else {
                    // A = U^{-1} * P^{-1} * A
                    mats[j] = u_inv(facts[j-1].U) * facts[j-1].P.T() * mats[j];

                }

            }

             // Left arrow = LEUP
             facts[j] = LEUP(mats[j]);


        } else {
            // Right arrow = PLEU

            // first apply change of basis from left
            if (j > 0) {
                if (is_left_arrow(dgm.elist[j-1])) {
                         // A = A * P^{-1} * U^{-1}
                         mats[j] = mats[j] * facts[j-1].P.T() * u_inv(facts[j-1].U);
                    } else {
                         // A = A * P * U
                         mats[j] = mats[j] * facts[j-1].P * facts[j-1].U;
                    }
            }

            facts[j] = PUEL(mats[j]);

        }
    }
    // we don't do leftward sweep if we only want barcode form

    // dump E matrices to return
    for (size_t j = 0; j < m; j++) {
        mats[j] = facts[j].E;
    }


    return mats;
}

struct bar {
    size_t start;
    size_t start_ind;
    size_t end;
    size_t end_ind;

    bar() {}

    bar(size_t s, size_t si, size_t t, size_t ti) : start(s), start_ind(si), end(t), end_ind(ti) {}
};

template <typename TN, typename TM>
std::vector<bar> barcode_from_barcode_form(
    const std::vector<TM> &mat,
    const Diagram<TN, TM> &dgm
) {

    size_t m = dgm.nedge(); // number of edges

    std::vector<bar> bars;
    std::map<size_t, size_t> piv_to_ind, piv_to_ind2;

    // first edge
    if (is_left_arrow(dgm.elist[0])) {
        // <-
        // loop over rows
        for (size_t i = 0; i < mat[0].nrow(); i++) {
            bars.emplace_back(bar(0, i, 0, i));
            piv_to_ind[i] = i;
        }

    } else {
        // ->
        // loop over columns
        for (size_t j = 0; j < mat[0].ncol(); j++) {
            bars.emplace_back(bar(0, j, 0, j));
            piv_to_ind[j] = j;
        }
    }

    for (size_t k = 0; k < m; k++) {
        piv_to_ind2.clear(); // updated pivots

        if (is_left_arrow(dgm.elist[k])) {
            // <-
            // extend bar by row
            for (size_t j = 0; j < mat[k].ncol(); j++) {
                auto piv = mat[k][j].nzbegin();
                if (piv == mat[k][j].nzend()) {
                    // no extension to complete. start new bar
                    piv_to_ind2[j] = bars.size();
                    bars.emplace_back(bar(k+1, j, k+1, j));
                } else {
                    // extend bar that is in pivot row
                    size_t pi = piv_to_ind[piv->ind];
                    bars[pi].end = k+1;
                    bars[pi].end_ind = j;
                    piv_to_ind2[j] = pi;
                }
            }
        } else {
            // ->
            // extend bar by column
            // easier to work on transpose
            auto matT = mat[k].T();
            for (size_t j = 0; j < matT.ncol(); j++) {
                auto piv = matT[j].nzbegin();
                if (piv == matT[j].nzend()) {
                    // no extension to complete. start new bar
                    // i.e. row was not in range of matrix
                    piv_to_ind2[j] = bars.size();
                    bars.emplace_back(bar(k+1, j, k+1, j));
                } else {
                    // extend bar that is in pivot row
                    size_t pi = piv_to_ind[piv->ind];
                    bars[pi].end = k+1;
                    bars[pi].end_ind = j;
                    piv_to_ind2[j] = pi;
                }
            }

        }
        // swap roles of pivot lookups
        std::swap(piv_to_ind, piv_to_ind2);
    }


    return bars;
}

std::vector<PersistencePair<size_t>> bars_to_pairs(
    const std::vector<bar> &bars,
    size_t hdim
) {
    size_t n = bars.size();
    std::vector<PersistencePair<size_t>> pairs(n);
    for (size_t i = 0; i < n; i++) {
        pairs[i] = PersistencePair(hdim, bars[i].start_ind, bars[i].end_ind, bars[i].start, bars[i].end);
    }
    return pairs;
}

template <typename NT, typename TM>
std::vector<PersistencePair<size_t>> barcode_sparse(const Diagram<NT, TM> &dgm, size_t hdim) {
    auto mats = barcode_form_leftward(dgm);
    auto bars = barcode_from_barcode_form(mats, dgm);
    return bars_to_pairs(bars, hdim);
}
