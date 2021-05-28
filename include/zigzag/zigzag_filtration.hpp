#pragma once
/*
Right filtrations for complexes
*/

#include <vector>
#include <utility> // std::pair
#include <chain/chain_complex.hpp>
#include "reduction.hpp"
#include <type_traits>

namespace bats {



/*
@brief A class that wraps a complex with a right-filtration

A class that wraps a complex with a right-filtration. Cells can have entry times and removal times.
These intervals are stored as std::pair<T, T>
*/
template <typename CpxT, typename T=double>
class ZigzagFiltration
{
private:
    CpxT X; // underlying complex
    std::vector<std::vector<std::pair<T, T>>> val; // right filtration values

	inline void reserve(const size_t maxdim) {
        while(val.size() < maxdim+1) {val.emplace_back(std::vector<std::pair<T,T>>()); }
        return;
    }
    void reserve(const size_t dim, const size_t n) {
        reserve(dim);
        if (val[dim].size() < n) {
          val[dim].resize(n);
        }
    }

public:

    ZigzagFiltration() {}

    /**
    Initialization which passes arguments to initialize the underlying complex

    @param args...	passed to complex initialization
    */
    template <class ...Ts>
    ZigzagFiltration(const Ts (&...args)) : X(args...) {}

	/**
	Construct right filtration explicitly on a complex

	@param X	complex representing topological space
	@param val	right filtration values for each cell in X

	val should be a vector of vector of pairs.
	val[k][i] is the pair of entry times for cell i in dimension k in X.

	No checks are done to make sure the number of values matches the number of cells
	*/
	ZigzagFiltration(
		const CpxT& X,
		const std::vector<std::vector<std::pair<T, T>>>& val
	) : X(X), val(val) {}

	/**
	return const reference to underlying complex
	*/
    inline const CpxT& complex() const { return X; }

	/**
	return const reference to right filtration values
	*/
	inline const std::vector<std::vector<std::pair<T, T>>>& vals() const { return val;}

	/**
	return const reference to right filtration values

	@param k	dimension of values to return
	*/
	inline const std::vector<std::pair<T,T>>& vals(const size_t k) const { return val[k]; }

	/**
	return maximum dimension of cells
	*/
	inline size_t maxdim() const { return val.size() - 1; }

	/**
	return number of cells in specified dimension

	@param dim	dimension
	*/
	inline size_t ncells(const size_t dim) const { return val[dim].size(); }

	/**
	add cell to right filtration

	@param entry	entry parameter
	@param exit		exit parameter
	@param ...args - passed to add method of underlying complex
	*/
	template <class ...Ts>
	inline cell_ind add(const T entry, const T exit, Ts (&...args)) {
		cell_ind ret = X.add(args...);
		reserve(ret.dim, ret.ind+1);
		val[ret.dim][ret.ind] = std::make_pair(entry, exit);
		return ret;
	}

	/**
	Add recursively to right filtration.  Any cells added will take same entry
	and exit parameters

	@param entry	entry parameter
	@param exit		exit parameter
	@param ...args - passed to add_recursive method of underlying complex
	*/
	template <class ...Ts>
	std::vector<cell_ind> add_recursive(const T entry, const T exit, Ts (&...args)) {
		std::vector<cell_ind> ret = X.add_recursive(args...);
		for (auto &r : ret) {
			reserve(r.dim, r.ind+1);
			val[r.dim][r.ind] = std::make_pair(entry, exit);
		}
		return ret;
	}

};

// TODO: write a function to extract chain complex,
// permute complex to processing order for reduction
// and determine order in which to process everything
template <typename CpxT, typename T, typename FT>
auto prepare_ChainComplex(
	const ZigzagFiltration<CpxT, T>& F,
	FT // field for coefficients
) {
	// chain complex
	auto C = Chain(F.complex(), FT());
	std::cout << "valid complex: " << C.is_valid_complex() << std::endl;

	// compute permutations to order chain complex
	// order is by decreasing removal time
	std::vector<std::vector<size_t>> perm(C.maxdim() + 1);
	for (size_t k = 0; k < C.maxdim() + 1; ++k) {
		perm[k] = bats::util::sortperm(
			F.vals(k).begin(),
			F.vals(k).end(),
			[&](const std::pair<T,T>& a, const std::pair<T,T>& b) { return a.second > b.second; }
	 	);
	}
	C.permute_basis(perm);
	std::cout << "after perm valid complex: " << C.is_valid_complex() << std::endl;

	// compute order to process column
	// at same filtration value:
	// additions should be in increasing dimension order
	// removals should be in decreasing dimension order
	std::vector<rfilt_val<T>> filt_order;
	filt_order.reserve(F.complex().ncells() * 2);
	for (size_t k = 0; k < C.maxdim() + 1; ++k) {
		auto permk = perm[k]; //bats::util::inv_perm(perm[k]);
		for (size_t i = 0; i < F.ncells(k); ++i) {
			auto pair = F.vals(k)[i];
			filt_order.emplace_back(
				rfilt_val(
					k, // dimension
					permk[i], // index in permuted chain complex
					i, // index in original complex
					pair.first, // entry time
					true // entry
				)
			);
			filt_order.emplace_back(
				rfilt_val(
					k, // dimension
					permk[i], // index in permuted chain complex
					i, // index in original complex
					pair.second, // removal time
					false // removal
				)
			);
		}
	}
	std::sort(
		filt_order.begin(),
		filt_order.end(),
		[&](const rfilt_val<T>& a, const rfilt_val<T>& b) {
			return (a.val != b.val) ?
				a.val < b.val : // increasing value order
				(a.entry && b.entry) ? // else, if both entries
				a.dim < b.dim : // add smaller dimension first
				(a.entry ^ b.entry) ? // one entry and one removal
				a.entry : // process entry first - doesn't matter
				(a.dim == b.dim) ? // else, both are removals - check if dim is the same
				a.ind > b.ind : // process later index first
				a.dim > b.dim; // otherwise, we process larger dimension first for removal
			}
	);

	return std::make_tuple(C, filt_order);

}

template <typename CpxT, typename T, typename FT, typename opt_flag, typename reduction_flag>
auto barcode(
	const ZigzagFiltration<CpxT, T>& F,
	FT, // field for coefficients
	opt_flag,
	reduction_flag
) {
	auto [C, filt_order] = prepare_ChainComplex(F, FT());

	// if constexpr(std::is_same<pairs_flag, apparent_pairs_flag>::value) {
	// 	C.clear_compress_apparent_pairs();
	// }

	// return zigzag_barcode_reduction(C, filt_order);
	return zigzag_barcode_reduction(C, filt_order, opt_flag(), reduction_flag());
}


} // namespace bats
