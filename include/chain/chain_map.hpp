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

    ChainMap(const std::vector<TM>& map) : map(map) {}

    // chain map on d dimensions
    ChainMap(size_t d) : map(d+1) {}

    ChainMap(const CellularMap &f) {
        size_t maxd = f.maxdim() + 1;
        map.resize(maxd);
        for (size_t k = 0; k < maxd; k++) {
            map[k] = TM(f[k]);
        }
    }

    // relative chain map
    // A \subseteq X
    // B \subseteq Y
    // f: X \to Y
    // ASSUMES: f(A) \subseteq B
    template <typename CpxT>
    ChainMap(
        const CellularMap& f,
        const CpxT& X,
        const CpxT& A,
        const CpxT& Y,
        const CpxT& B
    ) {
        size_t maxd = f.maxdim() + 1;
        map.resize(maxd);
        auto Xinds = X.get_indices(A);
        auto Yinds = Y.get_indices(B);
        for (size_t k = 0; k < maxd; k++) {
            auto Yinds = Y.get_indices(B, k);
            std::sort(Yinds.begin(), Yinds.end());
            auto Xinds = X.get_indices(A, k);
            std::sort(Xinds.begin(), Xinds.end());
            map[k] = TM(f[k].submatrix(
                bats::util::sorted_complement(Yinds, Y.ncells(k)),
                bats::util::sorted_complement(Xinds, X.ncells(k))
            ));
        }
    }

    // relative chain map
    ChainMap relative_map(
        std::vector<std::vector<size_t>>& inds1,
        std::vector<std::vector<size_t>>& inds2
    ) const {
        std::vector<TM> newmap(map.size());
        for (size_t k = 0; k < map.size(); k++) {
            std::sort(inds1[k].begin(), inds1[k].end());
            std::sort(inds2[k].begin(), inds2[k].end());
            newmap[k] = map[k].submatrix(
                bats::util::sorted_complement(inds1[k], map[k].nrow()),
                bats::util::sorted_complement(inds2[k], map[k].ncol())
            );
        }
        return ChainMap(newmap);
    }

    inline size_t maxdim() const { return map.size() - 1; }

    inline TM& operator[](size_t k) { return map[k]; }
    inline const TM& operator[](size_t k) const { return map[k]; }

    inline void permute_row_basis(size_t k, const std::vector<size_t> &p) {
        map[k].permute_rows(bats::util::inv_perm(p));
    }

    inline void permute_column_basis(size_t k, const std::vector<size_t> &p) {
        map[k].permute_cols(p);
    }


};

// defualt return
template <typename T>
inline auto Chain(const CellularMap &F, T) {
	using VT = SparseVector<T, size_t>;
	using MT = ColumnMatrix<VT>;

	return ChainMap<MT>(F);
}

} // namespace bats
