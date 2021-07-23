#pragma once

#include <functional>
#include <vector>
#include <utility>
#include <complex/cubical_complex.hpp>

/*
extend filtrations
*/

template <typename CpxT, typename T>
auto extend_filtration(
    const CpxT& X, // simplicial complex type
    std::function<std::tuple<T,size_t>(const std::vector<size_t>&)>& f // function to apply to simplices
) {
    std::vector<std::vector<T>> vals(X.maxdim() + 1);
    std::vector<std::vector<size_t>> inds(X.maxdim() + 1);
    for (size_t k = 0; k < X.maxdim() + 1; k++) {
        vals[k].reserve(X.ncells(k));
        inds[k].reserve(X.ncells(k));
        for (auto& s : X.get_simplices(k)) {
            auto [v, i] = f(s);
            vals[k].emplace_back(v);
            inds[k].emplace_back(i);
        }
    }
    return std::make_tuple(vals, inds);
}

template <typename CpxT, typename T>
auto lower_star_filtration(
    const CpxT& X, // simplicial complex type
    const std::vector<T>& f0
) {
    std::function<std::tuple<double, size_t>(const std::vector<size_t>&)> filtfn = [&f0](
        const std::vector<size_t>& s
    ) -> std::tuple<double, size_t> {
        auto it = std::max_element(s.begin(), s.end(), [&f0](size_t i, size_t j) {return f0[i] < f0[j];});
        return std::make_tuple(f0[*it], *it);
    };
    return extend_filtration(X, filtfn);
}

/**
return maximum vertex value on a cube c
assume 2 dimensions
*/
template <typename T>
auto max_cube_val(
    const std::vector<size_t>& c,
    const std::vector<std::vector<T>>& f0
) {
    return std::max(
        std::max(
            f0[c[0]][c[2]],
            f0[c[0]][c[3]]),
        std::max(
            f0[c[1]][c[2]],
            f0[c[1]][c[3]])
    );
}

/**
return maximum vertex value on a cube c
assume 3 dimensions
*/
template <typename T>
auto max_cube_val(
    const std::vector<size_t>& c,
    const std::vector<std::vector<std::vector<T>>>& f0
) {
    T val = f0[c[0]][c[2]][c[4]];
    for (size_t i : {0, 1}) {
        for (size_t j : {2, 3}) {
            for (size_t k : {4, 5}) {
                T tmp = f0[c[i]][c[j]][c[k]];
                val = (tmp > val) ? tmp : val;
            }
        }
    }
    return val;
}

template <typename T>
auto lower_star_filtration(
    const bats::CubicalComplex& X,
    const std::vector<std::vector<T>>& f0
) {
    std::vector<std::vector<T>> vals(X.maxdim() + 1);
    // std::vector<std::vector<size_t>> inds(X.maxdim() + 1);
    for (size_t k = 0; k < X.maxdim() + 1; k++) {
        vals[k].reserve(X.ncells(k));
        // inds[k].reserve(X.ncells(k));
        for (auto& s : X.get_cubes(k)) {
            // auto [v, i] = f(s);
            auto v = max_cube_val(s, f0);
            vals[k].emplace_back(v);
            // inds[k].emplace_back(i);
        }
    }
    // return std::make_tuple(vals, inds);
    return vals;
}

template <typename T>
auto lower_star_filtration(
    const bats::CubicalComplex& X,
    const std::vector<std::vector<std::vector<T>>>& f0
) {
    std::vector<std::vector<T>> vals(X.maxdim() + 1);
    // std::vector<std::vector<size_t>> inds(X.maxdim() + 1);
    for (size_t k = 0; k < X.maxdim() + 1; k++) {
        vals[k].reserve(X.ncells(k));
        // inds[k].reserve(X.ncells(k));
        for (auto& s : X.get_cubes(k)) {
            // auto [v, i] = f(s);
            auto v = max_cube_val(s, f0);
            vals[k].emplace_back(v);
            // inds[k].emplace_back(i);
        }
    }
    // return std::make_tuple(vals, inds);
    return vals;
}
