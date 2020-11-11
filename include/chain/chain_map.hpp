#pragma once
// Chain map class
#include <vector>
#include <cstddef>
#include <complex/cell_map.hpp>
#include <util/permutation.hpp>

namespace bats {

template <typename TM>
struct ChainMap {
    std::vector<TM> map;

    ChainMap() {}

    ChainMap(const CellularMap &f) {
        size_t maxd = f.maxdim() + 1;
        map.resize(maxd);
        for (size_t k = 0; k < maxd; k++) {
            map[k] = TM(f[k]);
        }
    }

    inline TM& operator[](size_t k) { return map[k]; }
    inline const TM& operator[](size_t k) const { return map[k]; }

    inline void permute_row_basis(size_t k, const std::vector<size_t> &p) {
        map[k].permute_rows(bats::util::inv_perm(p));
    }

    inline void permute_column_basis(size_t k, const std::vector<size_t> &p) {
        map[k].permute_cols(p);
    }


};

} // namespace bats
