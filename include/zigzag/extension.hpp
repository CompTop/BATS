#pragma once
/*
Extension functions for right filtrations
*/
#include <vector>
#include <utility> // pair
#include <limits> // numeric limits
#include <algorithm> // std::min, std::max
#include "zigzag_filtration.hpp"

namespace bats {
namespace zigzag {

namespace detail {

template <typename T, typename Iterator>
std::pair<T,T> simplex_extrema(
    Iterator it,
    const Iterator& end,
    const std::vector<T>& f0
) {
    T cmax = std::numeric_limits<T>::min();
    T cmin = std::numeric_limits<T>::max();
    while (it != end) {
        cmin = std::min(cmin, f0[*it]);
        cmax = std::max(cmax, f0[*it]);
        ++it;
    }
    return std::make_pair(cmin, cmax);
}

} // namespace detail

/**
Extension of right filtration on 0-cells to right-filtration on
all cells in a complex.  Right filtration goes from low to high values

@param f0   function on 0-cells
@param X    complex representing topological space
@param eps  thickened levelset radius

A 0-cell x enters the filtration at parameter f0(x) - eps and
is removed from the filtration at parameter f0(x) + eps

Higher-dimensional cells are present for the interval that all faces are also present.
*/
template <typename T, typename CpxT>
auto extend_zigzag_filtration(
    const std::vector<T>& f0, // function on 0-cells
    const CpxT& X,
    const T eps
) {
    std::vector<std::vector<std::vector<std::pair<T,T>>>> val(X.maxdim() + 1);
    for (size_t k = 0; k < X.maxdim() + 1; ++k) {
        val[k].resize(X.ncells(k));
        for (size_t i = 0; i < X.ncells(k); ++i) {
            auto minmax = detail::simplex_extrema(
                X.simplex_begin(k, i), X.simplex_end(k, i), f0
            );
            T entry = minmax.second - eps;
            T removal = minmax.first + eps;
            if (entry > removal) {
                auto spx = X.get_simplex(k, i);
                for (auto s : spx) { std::cout << s << ","; }
                std::cout << std::endl;
                std::cout << minmax.first << ", " << minmax.second << ", eps = " << eps << std::endl;
                throw std::runtime_error("entry time after removal!");
            }
            val[k][i].emplace_back(std::make_pair(entry, removal));

        }
    }

    return ZigzagFiltration(X, val);
}


namespace detail {

// get value of function at i,j,k location
template <typename T>
inline T cube_val(
    const std::vector<T>& f0,
    const size_t i,
    const size_t j,
    const size_t k,
    const size_t n
) {
    return f0[k + (j + i * n) * n];
}

/**
get maximum and minimum value of a function on cube vertices

@param f0   function on vertices - in column-major order
@param cube representation of cube, can be degenerate
@param n    length of 3D cube
*/
template <typename T>
std::pair<T,T> cube_extrema(
    const std::vector<T>& f0,
    const std::vector<size_t>& cube,
    const size_t n
) {
    T cmax = std::numeric_limits<T>::min();
    T cmin = std::numeric_limits<T>::max();
    for (auto i : {0,1}) {
        for (auto j : {2,3}) {
            for (auto k : {4,5}) {
                T cval = cube_val(f0, cube[i], cube[j], cube[k], n);
                cmin = std::min(cmin, cval);
                cmax = std::max(cmax, cval);
            }
        }
    }
    return std::make_pair(cmin, cmax);
}

} // namespace detail

/**
Extension of zigzag filtration on vertices  to a zigzag filtration on a cubical
complex.  The filtration progresses from low to high values.

@param f0   function on vertices - stored in column-major format
@param X    CubicalComplex - on a vertex set of size n^3
@param eps  thickened levelset radius
@param n    length of cube

@return val vector of vector of pairs holding zigzag filtration values for each cube.
*/
template <typename T>
std::vector<std::vector<std::pair<T,T>>> extend_levelset(
    const std::vector<T>& f0, // function on vertices
    const CubicalComplex& X,
    const T eps,
    const size_t n // length of cube
) {

    // zigzag filtration values
    std::vector<std::vector<std::vector<std::pair<T,T>>>> val(X.maxdim() + 1);

    std::vector<size_t> cube;
    std::pair<T, T> minmax;
    for (size_t dim = 0; dim < X.maxdim() + 1; ++dim) {
        val[dim].resize(X.ncells(dim));
        for (size_t i = 0; i < X.ncells(dim); ++i) {
            X.get_cube(dim, i, cube);
            // get extremal function values on cube vertices
            minmax = detail::cube_extrema(f0, cube, n);
            // entry time is max - eps
            T entry = minmax.second - eps;
            // removal time is min + eps
            T removal = minmax.first + eps;
            if (entry > removal) {throw std::runtime_error("entry time after removal!");}
            val[dim][i].emplace_back(std::make_pair(entry, removal));
        }

    }

    return val;
}

template <typename T>
ZigzagFiltration<CubicalComplex, T> extend_zigzag_filtration(
    const std::vector<T>& f0, // function on vertices
    const CubicalComplex& X,
    const T eps,
    const size_t n // length of cube
) {

    ZigzagFiltration<CubicalComplex, T> F(3); // initialize filtration

    std::vector<size_t> cube;
    std::pair<T, T> minmax;
    for (size_t dim = 0; dim < X.maxdim() + 1; ++dim) {
        for (size_t i = 0; i < X.ncells(dim); ++i) {
            X.get_cube(dim, i, cube);
            // get extremal function values on cube vertices
            minmax = detail::cube_extrema(f0, cube, n);
            // entry time is max - eps
            T entry = minmax.second - eps;
            // removal time is min + eps
            T removal = minmax.first + eps;
            if (entry > removal) {continue;} // do not add
            F.add(entry, removal, cube); // cube is present in filtration
        }

    }
    return F;

}

} // namespace zigzag
} // namespace bats
