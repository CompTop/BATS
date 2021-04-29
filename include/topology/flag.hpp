#pragma once

namespace bats {

// template over filtration type
// CpxT should be a
template <typename CpxT, typename NT>
void add_dimension_recursive_flag(
    CpxT &X,
    const NT &nbrs, // lists of neighbors
    const size_t d, // dimension
    const size_t maxd, // max dimension
    const std::vector<size_t> &iter_idxs,
    std::vector<size_t> &spx_idxs
) {
    // sorted simplices will end up here
    std::vector<size_t> spx_idxs2(spx_idxs.size() + 1);
    if (d == maxd) {
        // no recursion - we're adding maximal dimension cells
        for (auto k : iter_idxs) {
            // append k to spx_idxs, sort
            spx_idxs.push_back(k);
            bats::util::sort_into(spx_idxs, spx_idxs2);

            // add to X
            X.add(spx_idxs2);

            // pop k off spx_idxs
            spx_idxs.pop_back();
        }
    } else if (d < maxd) { // d < maxd
        // recursion
        std::vector<size_t> iter_idxs2; // indices for recursing on
        iter_idxs2.reserve(iter_idxs.size());
        for (auto k : iter_idxs) {
            // append k to spx_idxs, sort
            spx_idxs.push_back(k);
            bats::util::sort_into(spx_idxs, spx_idxs2);

            // add to X
            X.add(spx_idxs2);

            // recurse
            bats::util::intersect_sorted_lt(iter_idxs, nbrs[k], k, iter_idxs2);
            add_dimension_recursive_flag(X, nbrs, d+1, maxd, iter_idxs2, spx_idxs2);

            // pop k off spx_idxs
            spx_idxs.pop_back();
        }
    }
    // else do nothing
}


// Flag complex using list of edges
// (edges[2*k], edges[2*k+1]) = (i, j) is an edge
// n - number of vertices
// maxdim - maximum dimension of simplices
template <typename CpxT>
CpxT FlagComplex(
    const std::vector<size_t> &edges,
    const size_t n, // number of 0-cells
    const size_t maxdim
) {

    // check that dimensions agree
    size_t m = edges.size() / 2;
    if (!(edges.size() == 2 * m)) {
        throw std::logic_error("edge vector must have length multiple of 2!");
    }

    // X = SimplicialComplex(maxdim);
    // F = Filtration<T, SimplicialComplex>(X);
    // reset simplicial complex
    CpxT X(n, maxdim);

    // sets 0-cells
    std::vector<size_t> spx_idxs(1);
    for (size_t k = 0; k < n; k++) {
        spx_idxs[0] = k;
        X.add(spx_idxs);
    }

    std::vector<std::vector<size_t>> nbrs(n);

    spx_idxs.resize(2); // now time to add edges
    std::vector<size_t> iter_idxs;
    iter_idxs.reserve(n); // maximum size

    for (size_t k = 0; k < m; k++) {
        size_t i = edges[2*k];
        size_t j = edges[2*k + 1];
        spx_idxs[0] = i;
        spx_idxs[1] = j;
        X.add(spx_idxs);
        // std::cout << ret.second << std::endl;
        bats::util::intersect_sorted(nbrs[i], nbrs[j], iter_idxs);

        if (!iter_idxs.empty()) {
            add_dimension_recursive_flag(X, nbrs, 2, maxdim, iter_idxs, spx_idxs);
        }

        // TODO: use std::set for neighbors - insertion is log(n)
        // nbrs[i].emplace(j);
        // nbrs[j].emplace(i);
        // TODO: insertion sort
        nbrs[i].emplace_back(j);
        std::sort(nbrs[i].begin(), nbrs[i].end());
        nbrs[j].emplace_back(i);
        std::sort(nbrs[j].begin(), nbrs[j].end());
    }

    return X;
}


} // namespace bats
