#pragma once

// differential graded map
#include <vector>
#include <cstddef>
#include <complex/cell_map.hpp>
#include <util/permutation.hpp>

namespace bats {

/**
A class for a map between two DGVectorSpaces

The map should be between DGVectorspaces of the same degree
and the map should commute with the differential.

If the DGVectorspaces are augmented, the map should be
augmentation preserving.
*/
template <typename TM>
struct DGLinearMap {
    std::vector<TM> map;

    DGLinearMap() {}

    DGLinearMap(const std::vector<TM>& map) : map(map) {}

    // chain map on d dimensions
    DGLinearMap(size_t d) : map(d+1) {}

    DGLinearMap(const CellularMap &f, int deg=-1) {
        size_t maxd = f.maxdim() + 1;
        map.resize(maxd);
        if (deg == -1) {
            // homology type
            for (size_t k = 0; k < maxd; k++) {
                map[k] = TM(f[k]);
            }
        } else if (deg == +1) {
            // cohomology type
            for (size_t k = 0; k < maxd; k++) {
                map[k] = TM(f[k]).T(); // transpose
            }
        } else {
            throw std::runtime_error("Degree of DGLinearMap should be +1 or -1");
        }

    }

    inline ssize_t maxdim() const { return map.size() - 1; }

    inline TM& operator[](ssize_t k) { return map[k]; }
    inline const TM& operator[](ssize_t k) const { return map[k]; }

    inline void permute_row_basis(ssize_t k, const std::vector<size_t> &p) {
        map[k].permute_rows(bats::util::inv_perm(p));
    }

    inline void permute_column_basis(ssize_t k, const std::vector<size_t> &p) {
        map[k].permute_cols(p);
    }


};

} // namespace bats
