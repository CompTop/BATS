#pragma once
/**
Constructions on a grid

Freudenthal triangulation
*/

namespace bats {

namespace rowmajor {
    /**
    Utility for translating to row-major index
    */
    template <typename T>
    inline T get_idx(T i, T j, T n) {return j + n * i;}

    template <typename T>
    inline T get_idx(T i, T j, T k, T n1, T n2) {
        return k + n2 * (j + i * n1);
    }
}



/**
Freudenthal triangulation of 2-dimensional grid on m x n vertices

@param m number of rows
@param n number of columns

The grid is indexed in row-major order.  The index for vertex `(i,j)` is `j + n * i`
*/
template <typename CpxT>
CpxT Freudenthal(size_t m, size_t n) {
    CpxT X(m*n, 2);

    for (size_t i = 0; i < m-1; ++i) {
        for (size_t j = 0; j < n-1; ++j) {
            auto k1 = rowmajor::get_idx(i,j,n);
            auto k2 = rowmajor::get_idx(i+1, j, n);
            auto k3 = rowmajor::get_idx(i, j+1, n);
            auto k4 = rowmajor::get_idx(i+1, j+1, n);
            X.add_recursive({k1, k2, k4});
            X.add_recursive({k1, k3, k4});
        }
    }

    return X;
}


/**
Freudenthal triangulation of 3-dimensional grid on n1 x n2 x n3 vertices

@param n1 grid size in 1st index
@param n2 grid size in 2nd index
@param n3 grid size in 3rd index

The grid is indexed in row-major order.  The index for vertex `(i,j,k)` is `k + n2 * (j + n1 * i)`
*/
template <typename CpxT>
CpxT Freudenthal(size_t n1, size_t n2, size_t n3) {
    CpxT X(n1*n2*n3, 3); // 3 is the max complex dimension

    for (size_t i = 0; i < n1-1; ++i) {
        for (size_t j = 0; j < n2-1; ++j) {
            for (size_t k = 0; k < n3-1; ++k) {
                auto s1 = rowmajor::get_idx(i,j,k, n1, n2);
                auto s2 = rowmajor::get_idx(i,j,k+1, n1, n2);
                auto s3 = rowmajor::get_idx(i,j+1,k, n1, n2);
                auto s4 = rowmajor::get_idx(i+1,j,k, n1, n2);
                auto s5 = rowmajor::get_idx(i,j+1,k+1, n1, n2);
                auto s6 = rowmajor::get_idx(i+1,j,k+1, n1, n2);
                auto s7 = rowmajor::get_idx(i+1,j+1,k, n1, n2);
                auto s8 = rowmajor::get_idx(i+1,j+1,k+1, n1, n2);
                X.add_recursive({s1,s2,s5,s8}); // k, j, i
                X.add_recursive({s1,s2,s6,s8}); // k, i, j
                X.add_recursive({s1,s4,s7,s8}); // i, j, k
                X.add_recursive({s1,s4,s6,s8}); // i, k, j
                X.add_recursive({s1,s3,s7,s8}); // j, i, k
                X.add_recursive({s1,s3,s5,s8}); // j, k, i
            }
        }
    }

    return X;
}

/**
Subset of Freudenthal triangulation of 3-dimensional grid on n1 x n2 x n3 vertices

@param n1 grid size in 1st index
@param n2 grid size in 2nd index
@param n3 grid size in 3rd index
@param i0 start of first index
@param i1 end of first index
@param j0 start of second index
@param j1 end of second index
@param k0 start of third index
@param k1 end of third index

The grid is indexed in row-major order.  The index for vertex `(i,j,k)` is `k + n2 * (j + n1 * i)`

This function is useful for constructing zigzags through a grid
*/
template <typename CpxT>
CpxT Freudenthal(
    size_t n1, size_t n2, size_t n3, // dimensions of grid
    size_t i0, size_t i1, // range of first index
    size_t j0, size_t j1, // range of second index
    size_t k0, size_t k1  // range of third index
) {
    CpxT X(n1*n2*n3, 3); // 3 is the max complex dimension

    // declare index helpers
    size_t s1, s2, s3, s4, s5, s6, s7, s8;

    // add 2-simplices on planes
    // do this in case we are on a 2-dimensional slice
    for (size_t i = i0; i < i1; ++i) {
        for (size_t j = j0; j < j1; ++j) {
            for (size_t k = k0; k < k1; ++k) {
                // i-j plane
                s1 = rowmajor::get_idx(i,   j,   k,   n1, n2);
                s2 = rowmajor::get_idx(i+1, j,   k,   n1, n2);
                s3 = rowmajor::get_idx(i,   j+1, k,   n1, n2);
                s4 = rowmajor::get_idx(i+1, j+1, k,   n1, n2);
                X.add_recursive({s1, s2, s4});
                X.add_recursive({s1, s3, s4});

                // j-k plane
                s1 = rowmajor::get_idx(i,   j,   k,   n1, n2);
                s2 = rowmajor::get_idx(i,   j+1, k,   n1, n2);
                s3 = rowmajor::get_idx(i,   j,   k+1, n1, n2);
                s4 = rowmajor::get_idx(i,   j+1, k+1, n1, n2);
                X.add_recursive({s1, s2, s4});
                X.add_recursive({s1, s3, s4});

                // i-k plane
                s1 = rowmajor::get_idx(i,   j,   k,   n1, n2);
                s2 = rowmajor::get_idx(i+1, j,   k,   n1, n2);
                s3 = rowmajor::get_idx(i,   j,   k+1, n1, n2);
                s4 = rowmajor::get_idx(i+1, j,   k+1, n1, n2);
                X.add_recursive({s1, s2, s4});
                X.add_recursive({s1, s3, s4});
            }
        }
    }

    // add 3-simplices
    for (size_t i = i0; i < i1-1; ++i) {
        for (size_t j = j0; j < j1-1; ++j) {
            for (size_t k = k0; k < k1-1; ++k) {
                s1 = rowmajor::get_idx(i,j,k, n1, n2);
                s2 = rowmajor::get_idx(i,j,k+1, n1, n2);
                s3 = rowmajor::get_idx(i,j+1,k, n1, n2);
                s4 = rowmajor::get_idx(i+1,j,k, n1, n2);
                s5 = rowmajor::get_idx(i,j+1,k+1, n1, n2);
                s6 = rowmajor::get_idx(i+1,j,k+1, n1, n2);
                s7 = rowmajor::get_idx(i+1,j+1,k, n1, n2);
                s8 = rowmajor::get_idx(i+1,j+1,k+1, n1, n2);
                X.add_recursive({s1,s2,s5,s8}); // k, j, i
                X.add_recursive({s1,s2,s6,s8}); // k, i, j
                X.add_recursive({s1,s4,s7,s8}); // i, j, k
                X.add_recursive({s1,s4,s6,s8}); // i, k, j
                X.add_recursive({s1,s3,s7,s8}); // j, i, k
                X.add_recursive({s1,s3,s5,s8}); // j, k, i
            }
        }
    }

    return X;
}

} // namespace bats
