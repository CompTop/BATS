#pragma once

#include <functional>
#include <vector>

/*
extend filtrations
*/

template <typename CpxT, typename T>
auto extend_filtration(
    const CpxT& X, // simplicial complex type
    std::function<T(const std::vector<size_t>&)>& f // function to apply to simplices
) {
    std::vector<std::vector<T>> vals;
    vals.resize(X.maxdim() + 1);
    for (size_t k = 0; k < X.maxdim() + 1; k++) {
        vals[k].reserve(X.ncells(k));
        for (auto& s : X.get_simplices(k)) {
            vals[k].emplace_back(f(s));
        }
    }
    return vals;
}

template <typename CpxT, typename T>
auto lower_star_filtration(
    const CpxT& X, // simplicial complex type
    const std::vector<T>& f0
) {
    std::function<double(const std::vector<size_t>&)> filtfn = [&f0](
        const std::vector<size_t>& s
    ) -> double {
        return f0[*std::max_element(s.begin(), s.end(), [&f0](size_t i, size_t j) {return f0[i] < f0[j];})];
    };
    return extend_filtration(X, filtfn);
}
