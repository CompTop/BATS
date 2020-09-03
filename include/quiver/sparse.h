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

// pass facts[i].U, P to mats[i+1]
template <typename NT, typename TC, typename TM>
void pass_UP_right(
    const Diagram<NT, TM> &dgm,
    std::vector<SparseFact<TC>> &facts,
    std::vector<TM> &mats,
    ssize_t i
) {
    if (is_left_arrow(dgm.elist[i+1])) {
        if (is_left_arrow(dgm.elist[i])) {
            // A = U * P * A
            mats[i+1] = facts[i].U * facts[i].P * mats[i+1];
        } else {
            // A = U^{-1} * P^{-1} * A
            mats[i+1] = u_inv(facts[i].U) * facts[i].P.T() * mats[i+1];
        }
    } else {
        if (is_left_arrow(dgm.elist[i])) {
             // A = A * P^{-1} * U^{-1}
             mats[i+1] = mats[i+1] * facts[i].P.T() * u_inv(facts[i].U);
        } else {
             // A = A * P * U
             mats[i+1] = mats[i+1] * facts[i].P * facts[i].U;
        }
    }
    // TODO: clear facts[i].U and facts[i].P
}

// commute facts.L[i] through facts.E[i]
template <typename NT, typename TC, typename TM>
void commute_L_left(
    const Diagram<NT, TM> &dgm,
    std::vector<SparseFact<TC>> &facts,
    ssize_t i
) {
    if (is_left_arrow(dgm.elist[i])) {
        // EL - L commutation
        facts[i].L = EL_L_commute(facts[i].E, facts[i].L);
    } else {
        // L - ELhat commutation
        facts[i].L = L_EL_commute(facts[i].L, facts[i].E);
    }
}

// pass facts[i].L, to facts[i-1].L, with EL commute
template <typename NT, typename TC, typename TM>
void pass_L_left(
    const Diagram<NT, TM> &dgm,
    std::vector<SparseFact<TC>> &facts,
    ssize_t i
) {
    if (is_left_arrow(dgm.elist[i-1])) {
        if (is_left_arrow(dgm.elist[i])) {
            // L E_L L
            facts[i-1].L = facts[i-1].L * EL_L_commute(facts[i-1].E, facts[i].L);
        } else {
            // L E_L L^{-1}
            facts[i-1].L = facts[i-1].L * EL_L_commute(facts[i-1].E, l_inv(facts[i].L));
        }
    } else {
        if (is_left_arrow(dgm.elist[i])) {
            facts[i-1].L = L_EL_commute(l_inv(facts[i].L), facts[i-1].E) * facts[i-1].L;
        } else {
            facts[i-1].L = L_EL_commute(facts[i].L, facts[i-1].E) * facts[i-1].L;
        }
    }
}

// sweep 1 from left to right
// start at index j0
// end at index j1
template <typename NT, typename TC, typename TM>
void type_A_leftright_sweep1(
    const Diagram<NT, TM> &dgm,
    std::vector<SparseFact<TC>> &facts,
    std::vector<TM> &mats,
    ssize_t j0,
    ssize_t j1
) {
    // do first arrow
    facts[j0] = is_left_arrow(dgm.elist[j0]) ? LEUP(mats[j0]) : PUEL(mats[j0]);

    // left to right - apply change of basis, then do factorization
    for (ssize_t j = j0 + 1; j <= j1; j++) {
        // first apply change of basis from left
        pass_UP_right(dgm, facts, mats, j-1);
        // now perform factorization
        facts[j] = is_left_arrow(dgm.elist[j]) ? LEUP(mats[j]) : PUEL(mats[j]);
    }
}

// sweep 2 - commute L matrices back through quiver
// start at index j1
// end at index j0+1
template <typename NT, typename TC, typename TM>
void type_A_leftright_sweep2(
    const Diagram<NT, TM> &dgm,
    std::vector<SparseFact<TC>> &facts,
    ssize_t j0,
    ssize_t j1
) {
    // right to left
    for (ssize_t j = j1; j > j0; j--) {
        pass_L_left(dgm, facts, j);
    }
}

template <typename NT, typename TM>
auto barcode_form_leftright(const Diagram<NT, TM> &dgm) {

    size_t m = dgm.nedge();

    using TC = typename TM::col_type;
    std::vector<SparseFact<TC>> facts(m);

    // copy matrices on edges
    std::vector<TM> mats = dgm.edata;

    type_A_leftright_sweep1(dgm, facts, mats, 0, m-1);
    // we don't do second sweep if we only want barcode form
    // type_A_leftright_sweep2(dgm, facts, 0, m-1);
    // dump E matrices to return
    for (size_t j = 0; j < m; j++) {
        mats[j] = facts[j].E;
    }


    return mats;
}

// pass facts[i].U, P to mats[i+1]
template <typename NT, typename TC, typename TM>
void pass_PL_left(
    const Diagram<NT, TM> &dgm,
    std::vector<SparseFact<TC>> &facts,
    std::vector<TM> &mats,
    ssize_t i
) {
    if (is_left_arrow(dgm.elist[i-1])) {
        if (is_left_arrow(dgm.elist[i])) {
            // A = A * P * L
            mats[i-1] = mats[i-1] * facts[i].P * facts[i].L;
        } else {
            // A = A * P^{-1} * L^{-1}
            mats[i-1] = mats[i-1] * facts[i].P.T() * l_inv(facts[i].L);
        }
    } else {
        if (is_left_arrow(dgm.elist[i])) {
            // A = L^{-1} * P^{-1} * A
            mats[i-1] = l_inv(facts[i].L) * facts[i].P.T() * mats[i-1];
        } else {
            // A = L * P * A
            mats[i-1] = facts[i].L * facts[i].P * mats[i-1];
        }
    }
}

// commute facts.L[i] through facts.E[i]
template <typename NT, typename TC, typename TM>
void commute_U_right(
    const Diagram<NT, TM> &dgm,
    std::vector<SparseFact<TC>> &facts,
    ssize_t i
) {
    if (is_left_arrow(dgm.elist[i])) {
        // EU - U commutation // came from PLEU
        facts[i].U = EU_U_commute(facts[i].E, facts[i].U);
    } else {
        // L - ELhat commutation // came from UELP
        facts[i].U = U_EU_commute(facts[i].U, facts[i].E);
    }
}

// pass facts[i].U, to facts[i+1].U, with EL commute
template <typename NT, typename TC, typename TM>
void pass_U_right(
    const Diagram<NT, TM> &dgm,
    std::vector<SparseFact<TC>> &facts,
    ssize_t i
) {
    if (is_left_arrow(dgm.elist[i])) {
        if (is_left_arrow(dgm.elist[i+1])) {
            facts[i+1].U = U_EU_commute(facts[i].U, facts[i+1].E) * facts[i+1].U;
        } else {
            facts[i+1].U = facts[i+1].U * EU_U_commute(facts[i+1].E, u_inv(facts[i].U));
        }
    } else {
        if (is_left_arrow(dgm.elist[i+1])) {
            facts[i+1].U = U_EU_commute(u_inv(facts[i].U), facts[i+1].E) * facts[i+1].U;
        } else {
            facts[i+1].U = facts[i+1].U * EU_U_commute(facts[i+1].E, facts[i].U);
        }
    }
}

// sweep 1 from right to left
// start at index j0
// end at index j1
template <typename NT, typename TC, typename TM>
void type_A_rightleft_sweep2(
    const Diagram<NT, TM> &dgm,
    std::vector<SparseFact<TC>> &facts,
    ssize_t j0,
    ssize_t j1
) {
    // left to right
    for (ssize_t j = j0; j < j1 ; j++) {
        pass_U_right(dgm, facts, j);
    }
}

template <typename NT, typename TC, typename TM>
void type_A_rightleft_sweep1(
    const Diagram<NT, TM> &dgm,
    std::vector<SparseFact<TC>> &facts,
    std::vector<TM> &mats,
    ssize_t j0,
    ssize_t j1
) {
    // first factorization
    facts[j1] = is_left_arrow(dgm.elist[j1]) ? PLEU(mats[j1]) : UELP(mats[j1]);
    // right to left
    for (ssize_t j = j1-1; j >= j0; j--) {
        // change basis from the right
        pass_PL_left(dgm, facts, mats, j+1);
        // compute factorization
        facts[j] = is_left_arrow(dgm.elist[j]) ? PLEU(mats[j]) : UELP(mats[j]);
    }
}

template <typename NT, typename TM>
auto barcode_form_rightleft(const Diagram<NT, TM> &dgm) {

    size_t m = dgm.nedge();

    using TC = typename TM::col_type;
    std::vector<SparseFact<TC>> facts(m);

    // copy matrices on edges
    std::vector<TM> mats = dgm.edata;

    type_A_rightleft_sweep1(dgm, facts, mats, 0, m-1);
    // we don't do second sweep if we only want barcode form
    type_A_rightleft_sweep2(dgm, facts, 0, m-1);
    // dump E matrices to return
    for (size_t j = 0; j < m; j++) {
        mats[j] = facts[j].E;
    }

    return mats;
}

/*
Divide and conquer declarations
*/
template <typename NT, typename TC, typename TM>
void type_A_dq_EL(
    const Diagram<NT, TM> &dgm,
    std::vector<SparseFact<TC>> &facts,
    std::vector<TM> &mats,
    ssize_t j0,
    ssize_t j1
);

template <typename NT, typename TC, typename TM>
void type_A_dq_EU(
    const Diagram<NT, TM> &dgm,
    std::vector<SparseFact<TC>> &facts,
    std::vector<TM> &mats,
    ssize_t j0,
    ssize_t j1
);

template <typename NT, typename TC, typename TM>
void type_A_dq_common(
    const Diagram<NT, TM> &dgm,
    std::vector<SparseFact<TC>> &facts,
    std::vector<TM> &mats,
    ssize_t j0,
    ssize_t j1
) {
    // phase 1 - divide into two sub-quivers
    size_t j2 = j0 + j1 / 2; //edge at j2 will be LQU factorization
    size_t j0b = j2 - 1;
    size_t j1a = j2 + 1;


    // left side will be EL-type, right-side will be EU-type
    type_A_dq_EL(dgm, facts, mats, j0, j0b);
    type_A_dq_EU(dgm, facts, mats, j1a, j1);

    // now do the LQU factorization in the middle
    // first pass the UP and PL terms from left and right


    // take LQU factorization

    // commute L and U factors out

}


template <typename NT, typename TC, typename TM>
void type_A_dq_EU(
    const Diagram<NT, TM> &dgm,
    std::vector<SparseFact<TC>> &facts,
    std::vector<TM> &mats,
    ssize_t j0,
    ssize_t j1
) {
    // TODO
}

template <typename NT, typename TC, typename TM>
void type_A_dq_EL(
    const Diagram<NT, TM> &dgm,
    std::vector<SparseFact<TC>> &facts,
    std::vector<TM> &mats,
    ssize_t j0,
    ssize_t j1
) {
    // first check if j1 - j0 is small enough - then no recursion
    if (j1 - j0 < 4) {
        // TODO
        // sweep from left to right LEUP
        // sweep from right to left L-EL commutations

    }
    // else, we will recurse on two sub-quivers
    type_A_dq_common(dgm, facts, mats, j0, j1);

    // now update the right-hand side to be EL-type
    // start with Q term, then propagate rightward


    // at the end, we have:
    // a hanging L on the left
    // a hanging UP on the right

}

template <typename NT, typename TM>
auto barcode_form_divide_conquer(const Diagram<NT, TM> &dgm) {

    size_t m = dgm.nedge();

    using TC = typename TM::col_type;
    std::vector<SparseFact<TC>> facts(m);

    // copy matrices on edges
    std::vector<TM> mats = dgm.edata;

    type_A_dq_common(dgm, facts, mats, m-1, 0);
    // we don't do second sweep if we only want barcode form

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
    auto mats = barcode_form_leftright(dgm);
    auto bars = barcode_from_barcode_form(mats, dgm);
    return bars_to_pairs(bars, hdim);
}

template <typename NT, typename TM>
std::vector<PersistencePair<size_t>> barcode_sparse_rightleft(const Diagram<NT, TM> &dgm, size_t hdim) {
    auto mats = barcode_form_rightleft(dgm);
    auto bars = barcode_from_barcode_form(mats, dgm);
    return bars_to_pairs(bars, hdim);
}
