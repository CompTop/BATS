#pragma once

#include <vector>
#include <map>

#include <multigraph/diagram.hpp>
#include <linalg/col_matrix.hpp>
#include <linalg/sparse_fact.hpp>
#include <persistence/barcode.hpp>

#include <omp.h> // openMP header

namespace bats {

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
    facts[i].U = TM::identity(facts[i].U.nrow());
    facts[i].P = TM::identity(facts[i].P.nrow());
}

template <typename NT, typename TC, typename TM>
void pass_P_right(
    const Diagram<NT, TM> &dgm,
    std::vector<SparseFact<TC>> &facts,
    ssize_t i
) {
    if (is_left_arrow(dgm.elist[i+1])) {
        if (is_left_arrow(dgm.elist[i])) {
            facts[i+1].E = facts[i].P * facts[i+1].E;
        } else {
            // A = U^{-1} * P^{-1} * A
            facts[i+1].E = facts[i].P.T() * facts[i+1].E;
        }
    } else {
        if (is_left_arrow(dgm.elist[i])) {
             // A = A * P^{-1} * U^{-1}
             facts[i+1].E = facts[i+1].E * facts[i].P.T();
        } else {
             // A = A * P * U
             facts[i+1].E = facts[i+1].E * facts[i].P;
        }
    }
    // clear facts[i].P
    facts[i].P = TM::identity(facts[i].P.nrow());
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
    // clear facts[i].L
    facts[i].L = TM::identity(facts[i].L.nrow());
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

// pass facts[i].L, P to mats[i-1]
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
    facts[i].L = TM::identity(facts[i].L.nrow());
    facts[i].P = TM::identity(facts[i].P.nrow());
}

// pass facts[i].L, P to mats[i-1]
template <typename NT, typename TC, typename TM>
void pass_P_left(
    const Diagram<NT, TM> &dgm,
    std::vector<SparseFact<TC>> &facts,
    ssize_t i
) {
    if (is_left_arrow(dgm.elist[i-1])) {
        if (is_left_arrow(dgm.elist[i])) {
            // A = A * P * L
            facts[i-1].E = facts[i-1].E * facts[i].P;
        } else {
            // A = A * P^{-1} * L^{-1}
            facts[i-1].E = facts[i-1].E * facts[i].P.T();
        }
    } else {
        if (is_left_arrow(dgm.elist[i])) {
            // A = L^{-1} * P^{-1} * A
            facts[i-1].E = facts[i].P.T() * facts[i-1].E;
        } else {
            // A = L * P * A
            facts[i-1].E = facts[i].P * facts[i-1].E;
        }
    }
    facts[i].P = TM::identity(facts[i].P.nrow());
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
    facts[i].U = TM::identity(facts[i].U.nrow());
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
ssize_t type_A_dq_common(
    const Diagram<NT, TM> &dgm,
    std::vector<SparseFact<TC>> &facts,
    std::vector<TM> &mats,
    ssize_t j0,
    ssize_t j1
) {
    // phase 1 - divide into two sub-quivers
    size_t j2 = (j0 + j1) / 2; //edge at j2 will be LQU factorization
    // std::cout << "partitioning at " << j2 <<  std::endl;
    size_t j0b = j2 - 1;
    size_t j1a = j2 + 1;


    // left side will be EL-type, right-side will be EU-type
    #pragma omp task shared(dgm, facts, mats) firstprivate(j0, j0b)
    (void) type_A_dq_EL(dgm, facts, mats, j0, j0b);
    #pragma omp task shared(dgm, facts, mats) firstprivate(j1a, j1)
    (void) type_A_dq_EU(dgm, facts, mats, j1a, j1);
    #pragma omp taskwait

    // now do the LQU factorization in the middle
    // first pass the UP and PL terms from left and right
    pass_PL_left(dgm, facts, mats, j1a);
    pass_UP_right(dgm, facts, mats, j0b);

    // take LQU factorization
    facts[j2] = is_left_arrow(dgm.elist[j2]) ? LQU(mats[j2]) : UQL(mats[j2]);

    // pass L and U factors out
    pass_U_right(dgm, facts, j2);
    pass_L_left(dgm, facts, j2);

    // commute L and U factors out
    // TODO: could specialize for identity
    #pragma omp task shared(dgm, facts) firstprivate(j0, j0b)
    (void) type_A_leftright_sweep2(dgm, facts, j0, j0b);
    #pragma omp task shared(dgm, facts) firstprivate(j1a, j1)
    (void) type_A_rightleft_sweep2(dgm, facts, j1a, j1);
    #pragma omp taskwait

    // quiver now has pivot matrices everywhere
    // there is a dangling L on the left and a dangling U on the right.

    return j2; // return break location

}


template <typename NT, typename TC, typename TM>
void type_A_dq_EU(
    const Diagram<NT, TM> &dgm,
    std::vector<SparseFact<TC>> &facts,
    std::vector<TM> &mats,
    ssize_t j0,
    ssize_t j1
) {
    // first check if j1 - j0 is small enough - then no recursion
    // alternatively, if we shouldn't be spawining new levels of recursion
    // if (j1 - j0 < 4) {
    if (j1 - j0 < 4 ) {
        type_A_rightleft_sweep1(dgm, facts, mats, j0, j1);
        type_A_rightleft_sweep2(dgm, facts, j0, j1);
        return;
    }
    // else, we will recurse on two sub-quivers
    ssize_t j2 = type_A_dq_common(dgm, facts, mats, j0, j1);

    // now update the left-hand side to be EU-type
    // start with Q term (at j2), then propagate leftward
    for (ssize_t j = j2; j > j0; j--) {
        facts[j] = is_left_arrow(dgm.elist[j]) ? PLEU(facts[j].E) : UELP(facts[j].E);
        // pass P term to left (L term is identity)
        pass_P_left(dgm, facts, j);
        // pass_PL_left(dgm, facts, j);
    }
    // handle last node - has dangling L term.
    facts[j0] = is_left_arrow(dgm.elist[j0]) ? PLEU(facts[j0].L * facts[j0].E) : UELP(facts[j0].E * facts[j0].L);
    // now commute U term all the way back
    type_A_rightleft_sweep2(dgm, facts, j0, j1);

    // at the end, we have:
    // a hanging PL on the left
    // a hanging U on the right
    return;
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
        type_A_leftright_sweep1(dgm, facts, mats, j0, j1);
        type_A_leftright_sweep2(dgm, facts, j0, j1);
        return;
    }
    // else, we will recurse on two sub-quivers
    ssize_t j2 = type_A_dq_common(dgm, facts, mats, j0, j1);

    // now update the right-hand side to be EL-type
    // start with Q term (at j2), then propagate rightward
    for (ssize_t j = j2; j < j1; j++) {
        facts[j] = is_left_arrow(dgm.elist[j]) ? LEUP(facts[j].E) : PUEL(facts[j].E);
        // pass P term to right (U term is identity)
        pass_P_right(dgm, facts, j);
        // pass_UP_right(dgm, facts, j);
    }
    // handle last node - has dangling U term.
    facts[j1] = is_left_arrow(dgm.elist[j1]) ? LEUP(facts[j1].E * facts[j1].U) : PUEL(facts[j1].U * facts[j1].E);
    // now commute L term all the way back
    type_A_leftright_sweep2(dgm, facts, j0, j1);

    // at the end, we have:
    // a hanging L on the left
    // a hanging UP on the right
    return;

}

template <typename NT, typename TM>
auto barcode_form_divide_conquer(const Diagram<NT, TM> &dgm) {

    size_t m = dgm.nedge();

    // compute number of levels of recursion to use
    int nthreads = (int) omp_get_max_threads();
    int nlevels = 0;
    while (nthreads >>= 1) nlevels++;

    // check if we should just do sequential algorithm
    if (m < 5 || nlevels <= 1) { return barcode_form_leftright(dgm);}

    using TC = typename TM::col_type;
    std::vector<SparseFact<TC>> facts(m);

    // copy matrices on edges
    std::vector<TM> mats = dgm.edata;

    // begin parallel region
    #pragma omp parallel default(none) shared(dgm, facts, mats, m, nlevels)
    {
        #pragma omp single nowait
        {
            (void) type_A_dq_common(dgm, facts, mats, 0, m-1);
        }
    }

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

namespace flags {
/**
Flag to choose a divide and conquer algorithm
*/
struct divide_conquer {};

/**
Flag to choose rightward algorithm
*/
struct rightward {};

/**
Flag to choose leftward algorithm
*/
struct leftward {};

} // namespace flags

/**
Compute barcode from diagram of vector spaces and linear maps

Uses divide and conquer algorithm by default.
*/
template <typename NT, typename TM>
inline auto barcode(const Diagram<NT, TM> &dgm, size_t hdim) {
    auto mats = barcode_form_divide_conquer(dgm);
    auto bars = barcode_from_barcode_form(mats, dgm);
    return bars_to_pairs(bars, hdim);
}

/**
Compute barcode from diagram of vector spaces and linear maps

Uses divide and conquer algorithm
*/
template <typename NT, typename TM>
inline auto barcode(const Diagram<NT, TM> &dgm, size_t hdim, flags::divide_conquer) {
    auto mats = barcode_form_divide_conquer(dgm);
    auto bars = barcode_from_barcode_form(mats, dgm);
    return bars_to_pairs(bars, hdim);
}

/**
Compute barcode from diagram of vector spaces and linear maps

Uses rightward algorithm (sweeps right to left)
*/
template <typename NT, typename TM>
inline auto barcode(const Diagram<NT, TM> &dgm, size_t hdim, flags::leftward) {
    auto mats = barcode_form_rightleft(dgm);
    auto bars = barcode_from_barcode_form(mats, dgm);
    return bars_to_pairs(bars, hdim);
}

/**
Compute barcode from diagram of vector spaces and linear maps

Uses leftward algorithm (sweeps left to right)
*/
template <typename NT, typename TM>
inline auto barcode(const Diagram<NT, TM> &dgm, size_t hdim, flags::rightward) {
    auto mats = barcode_form_leftright(dgm);
    auto bars = barcode_from_barcode_form(mats, dgm);
    return bars_to_pairs(bars, hdim);
}


/**
extracts dimenion k from a diagram with stacked dimensions
*/
template <typename NT, typename TM>
auto extract_dimension(const Diagram<NT, std::vector<TM>>& D, size_t k) {
    size_t n = D.nnode();
	size_t m = D.nedge();

    Diagram<void*, TM> TD(n, m);
    #pragma omp parallel for
    for (size_t i = 0; i < m; ++i) {
        auto s = D.elist[i].src;
		auto t = D.elist[i].targ;
        TD.set_edge(i, s, t, D.edata[i][k]);
    }

    return TD;

}

/**
Compute barcode in all dimensions of diagram
*/
template <typename NT, typename TM, typename ...Args>
inline auto barcode(const Diagram<NT, std::vector<TM>>& dgm, Args ...args) {
    size_t maxdim = dgm.edata[0].size();
    std::vector<PersistencePair<size_t>> pairs;
    for (size_t k = 0; k < maxdim; ++k) {
        auto dgmk = extract_dimension(dgm, k);
        auto pairsk = barcode(dgmk, k, args...);
        pairs.insert(pairs.end(), pairsk.begin(), pairsk.end());
    }
    return pairs;
}


} // namespace bats
