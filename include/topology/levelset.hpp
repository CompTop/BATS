#pragma once
/*
A variety of levelset constructions
*/
#include <limits>
#include <zigzag/zigzag_filtration.hpp>
#include "grid.hpp"

namespace bats {

/**
Create a zigzag filtration from a 3D image.

Extend zigzag filtration from toplexes i.e. the maximal cubes correspond to the
pixel grid.

Lower-dimensional cube filtration values are extended from cofaces

This means a n1 x n2 x n3 image will be on a vertex set of size
(n1 + 1) x (n2 + 1) x (n3 + 1)
*/
template <typename T>
zigzag::ZigzagFiltration<CubicalComplex, T> zigzag_toplex(
    const std::vector<std::vector<std::vector<T>>>& img
) {
    zigzag::ZigzagFiltration<CubicalComplex, T> X(3); // initialize filtration on 3-dim cubical complex

    // get grid dimensions
    size_t n1 = img.size() + 1;
    size_t n2 = img[0].size() + 1;
    size_t n3 = img[0][0].size() + 1;

    T b, d; // birth and death parameters

    // add 0-cubes
    for (size_t i = 0; i < n1; ++i) {
        for (size_t j = 0; j < n2; ++j) {
            for (size_t k = 0; k < n3; ++k) {
                // 8 cases to check ...
                b = std::numeric_limits<T>::max();
                d = std::numeric_limits<T>::lowest();
                if (i < n1-1 && j < n2-1 && k < n3-1) {
                    b = std::min(b, img[i][j][k]);
                    d = std::max(d, img[i][j][k]);
                }
                if (i > 0 && j < n2-1 && k < n3-1) {
                    b = std::min(b, img[i-1][j][k]);
                    d = std::max(d, img[i-1][j][k]);
                }
                if (i < n1-1 && j > 0 && k < n3-1) {
                    b = std::min(b, img[i][j-1][k]);
                    d = std::max(d, img[i][j-1][k]);
                }
                if (i < n1-1 && j < n2-1 && k > 0) {
                    b = std::min(b, img[i][j][k-1]);
                    d = std::max(d, img[i][j][k-1]);
                }
                if (i < n1-1 && j > 0 && k > 0) {
                    b = std::min(b, img[i][j-1][k-1]);
                    d = std::max(d, img[i][j-1][k-1]);
                }
                if (i > 0 && j < n2-1 && k > 0) {
                    b = std::min(b, img[i-1][j][k-1]);
                    d = std::max(d, img[i-1][j][k-1]);
                }
                if (i > 0 && j > 0 && k < n3-1) {
                    b = std::min(b, img[i-1][j-1][k]);
                    d = std::max(d, img[i-1][j-1][k]);
                }
                if (i > 0 && j > 0 && k > 0) {
                    b = std::min(b, img[i-1][j-1][k-1]);
                    d = std::max(d, img[i-1][j-1][k-1]);
                }
                X.add(b, d, {i, i, j, j, k, k});
            }
        }
    }

    // add 1-cubes
    for (size_t i = 0; i < n1-1; ++i) {
        for (size_t j = 0; j < n2; ++j) {
            for (size_t k = 0; k < n3; ++k) {
                // i direction
                b = std::numeric_limits<T>::max();
                d = std::numeric_limits<T>::lowest();
                if (j > 0 && k > 0) {
                    b = std::min(b, img[i][j-1][k-1]);
                    d = std::max(b, img[i][j-1][k-1]);
                }
                if (j < n2-1 && k < n3-1) {
                    b = std::min(b, img[i][j][k]);
                    d = std::max(b, img[i][j][k]);
                }
                if (j > 0 && k < n3-1) {
                    b = std::min(b, img[i][j-1][k]);
                    d = std::max(b, img[i][j-1][k]);
                }
                if (j < n2-1 && k > 0) {
                    b = std::min(b, img[i][j][k-1]);
                    d = std::max(b, img[i][j][k-1]);
                }
                X.add(b, d, {i, i+1, j, j, k, k});
            }
        }
    }
    for (size_t i = 0; i < n1; ++i) {
        for (size_t j = 0; j < n2-1; ++j) {
            for (size_t k = 0; k < n3; ++k) {
                // j direction
                b = std::numeric_limits<T>::max();
                d = std::numeric_limits<T>::lowest();
                if (i > 0 && k > 0) {
                    b = std::min(b, img[i-1][j][k-1]);
                    d = std::max(b, img[i-1][j][k-1]);
                }
                if (i < n1-1 && k < n3-1) {
                    b = std::min(b, img[i][j][k]);
                    d = std::max(b, img[i][j][k]);
                }
                if (i > 0 && k < n3-1) {
                    b = std::min(b, img[i-1][j][k]);
                    d = std::max(b, img[i-1][j][k]);
                }
                if (i < n1-1 && k > 0) {
                    b = std::min(b, img[i][j][k-1]);
                    d = std::max(b, img[i][j][k-1]);
                }
                X.add(b, d, {i, i, j, j+1, k, k});
            }
        }
    }
    for (size_t i = 0; i < n1; ++i) {
        for (size_t j = 0; j < n2; ++j) {
            for (size_t k = 0; k < n3-1; ++k) {
                // k direction
                b = std::numeric_limits<T>::max();
                d = std::numeric_limits<T>::lowest();
                if (i > 0 && j > 0) {
                    b = std::min(b, img[i-1][j-1][k]);
                    d = std::max(b, img[i-1][j-1][k]);
                }
                if (i < n1-1 && j < n2-1) {
                    b = std::min(b, img[i][j][k]);
                    d = std::max(b, img[i][j][k]);
                }
                if (i > 0 && j < n2-1) {
                    b = std::min(b, img[i-1][j][k]);
                    d = std::max(b, img[i-1][j][k]);
                }
                if (i < n1-1 && j > 0) {
                    b = std::min(b, img[i][j-1][k]);
                    d = std::max(b, img[i][j-1][k]);
                }
                X.add(b, d, {i, i, j, j, k, k+1});
            }
        }
    }

    // add 2-cubes
    for (size_t i = 0; i < n1-1; ++i) {
        for (size_t j = 0; j < n2-1; ++j) {
            for (size_t k = 0; k < n3; ++k) {
                // i-j plane
                b = std::numeric_limits<T>::max();
                d = std::numeric_limits<T>::lowest();
                if (k < n3-1) {
                    b = std::min(b, img[i][j][k]);
                    d = std::max(d, img[i][j][k]);
                }
                if (k > 0) {
                    b = std::min(b, img[i][j][k-1]);
                    d = std::max(d, img[i][j][k-1]);
                }
                X.add(b, d, {i, i+1, j, j+1, k, k});
            }
        }
    }
    for (size_t i = 0; i < n1; ++i) {
        for (size_t j = 0; j < n2-1; ++j) {
            for (size_t k = 0; k < n3-1; ++k) {
                // j-k plane
                b = std::numeric_limits<T>::max();
                d = std::numeric_limits<T>::lowest();
                if (i < n1-1) {
                    b = std::min(b, img[i][j][k]);
                    d = std::max(d, img[i][j][k]);
                }
                if (i > 0) {
                    b = std::min(b, img[i-1][j][k]);
                    d = std::max(d, img[i-1][j][k]);
                }
                X.add(b, d, {i, i, j, j+1, k, k+1});
            }
        }
    }
    for (size_t i = 0; i < n1-1; ++i) {
        for (size_t j = 0; j < n2; ++j) {
            for (size_t k = 0; k < n3-1; ++k) {
                // i-k plane
                b = std::numeric_limits<T>::max();
                d = std::numeric_limits<T>::lowest();
                if (j < n2-1) {
                    b = std::min(b, img[i][j][k]);
                    d = std::max(d, img[i][j][k]);
                }
                if (j > 0) {
                    b = std::min(b, img[i][j-1][k]);
                    d = std::max(d, img[i][j-1][k]);
                }
                X.add(b, d, {i, i+1, j, j, k, k+1});
            }
        }
    }

    // add 3-cubes
    for (size_t i = 0; i < n1-1; ++i) {
        for (size_t j = 0; j < n2-1; ++j) {
            for (size_t k = 0; k < n3-1; ++k) {
                X.add(
                    img[i][j][k],
                    img[i][j][k],
                    {i, i+1, j, j+1, k, k+1}
                );
            }
        }
    }

    return X;
}

/**
Function which creates a blocked levelset zigzag diagram from
a zigzag::ZigzagFiltration

f^{-1}([s_i,s_{i+1}]) \subseteq f^{-1}([s_i, s_{i+2}]) \supseteq f^{-1}([s_{i+1}, s_{i+2}])

where eps = s_{i+1} - s_i

@param X zigzag::ZigzagFiltration
@param eps window size
@param t0 lower bound on first window
@param t1 upper bound on last window
*/
template <typename CpxT, typename T>
auto zigzag_levelsets(
    const zigzag::ZigzagFiltration<CpxT, T>& X,
    T eps,
    T t0,
    T t1
) {

    // create sets
    std::vector<std::pair<T, T>> sets;
    T si = t0;
    sets.emplace_back(std::make_pair(si, si + eps)); // first set
    while (si < t1 - eps) {
        sets.emplace_back(std::make_pair(si, si + 2*eps));
        si += eps;
        sets.emplace_back(std::make_pair(si, si + eps));
    }

    size_t n = sets.size();

    // create diagram
    Diagram<CpxT, CellularMap> D(n, n-1);

    // fill in sets
    #pragma omp parallel for
    for (size_t i = 0; i < n; ++i) {
        D.set_node(i, X.levelset(sets[i].first, sets[i].second));
    }

    // fill in inclusion maps
    #pragma omp parallel for
    for (size_t i = 0; i < n-1; ++i) {
        size_t src = (i % 2 == 0) ? i : i+1;
        size_t targ = (i % 2 == 0) ? i+1 : i;
        D.set_edge(i, src, targ,
            CellularMap(D.node_data(src), D.node_data(targ))
        );
    }

    return D;

}



} // namespace bats
