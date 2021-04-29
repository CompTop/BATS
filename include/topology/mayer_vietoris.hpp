#pragma once

#include <complex/simplicial_complex.hpp>
#include <vector>

namespace bats {

/*
X = union(A, B)
Y = intersection(A, B)
return Mayer-Vietoris boundary C_{k+1}(X) \to C_k(Y)
a \mapsto boundary(a) |Y
b \mapsto 0
*/
template <typename CpxT> // SimplicialComplex
auto mayer_vietoris_boundary(
    const CpxT& A,
    const CpxT& X,
    const CpxT& Y,
    size_t k
) {
    using TC = SparseVector<int>;

    auto n = X.ncells(k+1);
    auto m = Y.ncells(k);


    std::vector<TC> col; // columns
    col.reserve(n);
    std::vector<size_t> face;
    face.reserve(k+1);
    std::vector<size_t> inds;
    std::vector<int> vals;

    for (auto& s : X.get_simplices(k+1)) {
        if (A.find_idx(s) != bats::NO_IND) {
            // compute boundary of s restricted to Y
            int c = -1;
            // loop over faces in lexicographical order
            size_t dim = k+1;
            size_t spx_len = dim + 1;
            inds.clear();
            vals.clear();
            for (size_t i = 0; i < spx_len; i++) {

                size_t i2 = dim-i; // index to skip
                face.clear();
                for (size_t j = 0; j < i2; j++) {
                    face.emplace_back(s[j]);
                }
                for (size_t j = i2+1; j < spx_len; j++) {
                    face.emplace_back(s[j]);
                }

                size_t ind = Y.find_idx(face);
                // if face is in Y, add it to boundary
                if (ind != bats::NO_IND) {
                    inds.emplace_back(ind);
                    vals.emplace_back(c);
                }
                c = -c;
            }
            col.emplace_back(TC(inds, vals));

        } else {
            // boundary is 0
            col.emplace_back(TC());
        }
    }
    return ColumnMatrix<TC>(m, n, col);
}

} // namespace bats
