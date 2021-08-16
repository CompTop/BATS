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
namespace zigzag {



/**
@brief A class that wraps a complex with a zigzag filtration

A class that wraps a complex with a right filtration. Cells can have entry times and removal times.
These intervals are stored as std::pair<T, T>
*/
template <typename CpxT, typename T=double>
class ZigzagFiltration
{
private:
    CpxT X; // underlying complex
    std::vector<std::vector<std::vector<std::pair<T, T>>>> val; // zigzag filtration values

	inline void reserve(const size_t maxdim) {
        while(val.size() < maxdim+1) {val.emplace_back(std::vector<std::vector<std::pair<T,T>>>()); }
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
		const std::vector<std::vector<std::vector<std::pair<T, T>>>>& val
	) : X(X), val(val) {}

	/**
	return const reference to underlying complex
	*/
    inline const CpxT& complex() const { return X; }

	/**
	return const reference to right filtration values
	*/
	inline const std::vector<std::vector<std::vector<std::pair<T, T>>>>& vals() const { return val;}

	/**
	return const reference to right filtration values

	@param k	dimension of values to return
	*/
	inline const std::vector<std::vector<std::pair<T,T>>>& vals(const size_t k) const { return val[k]; }

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
		val[ret.dim][ret.ind].emplace_back(std::make_pair(entry, exit));
		return ret;
	}
	inline cell_ind add(const T entry, const T exit, std::vector<size_t>&& s) {
		cell_ind ret = X.add(s);
		reserve(ret.dim, ret.ind+1);
		val[ret.dim][ret.ind].emplace_back(std::make_pair(entry, exit));
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
			val[r.dim][r.ind].emplace_back(std::make_pair(entry, exit));
		}
		return ret;
	}
	std::vector<cell_ind> add_recursive(const T entry, const T exit, std::vector<size_t>&& s) {
		std::vector<cell_ind> ret = X.add_recursive(s);
		for (auto &r : ret) {
			reserve(r.dim, r.ind+1);
			val[r.dim][r.ind].emplace_back(std::make_pair(entry, exit));
		}
		return ret;
	}

	/**
	Return thickened levelset f^{-1}([s0,s1])
	adds a cell to levelset if [b,d] \cap [s0,s1] is non-empty

	@param s0 lower bound of levelset
	@param s1 upper bound of levelset
	*/
	CpxT levelset(T s0, T s1) const {
		CpxT Y(X.ncells(0), X.maxdim());

		// loop over dimension
		for (size_t k = 0; k < X.maxdim() + 1; ++k) {
			for (size_t i = 0; i < X.ncells(k); ++i) {
				// assume only one interval
				if ((s0 < val[k][i][0].first && val[k][i][0].first < s1) || (s0 < val[k][i][0].second && val[k][i][0].second < s1)) {
					Y.add(X.get_cell(k, i));
				}
			}
		}

		return Y;
	}

};

/**
@brief A class that wraps a chain complex with a zigzag filtration

A class that wraps a chain complex with a zigzag filtration.
Unlike a Zigzag filtraion, every chain has a unique entry and removal time.
*/
template <typename MT, typename T=double>
struct ZigzagChainComplex
{

	ChainComplex<MT> C; // chain complex
	std::vector<std::vector<std::pair<T, T>>> val; // zigzag filtration values
	std::vector<std::vector<size_t>> cind; // maps back to original index in complex

	/**
	correct the indices in column j in dimension k

	assumes column j hasn't already been corrected
	assumes val[k][j] has been set, as well as val[k-1]
	extra_cells maps to duplicate cells
	*/
	void _correct_indices(
		size_t k,
		size_t j,
		const std::vector<std::vector<size_t>> extra_cells
	) {
		auto it = C[k][j].nzbegin();
		while (it != C[k][j].nzend()) {
			auto i = it->ind;
			// if this column enters after the original boundary exits, we need
			// to search for copy of boundary
			if (val[k][j].first > val[k-1][i].second) {
				for (size_t i2 = extra_cells[k-1][i]; i2 < extra_cells[k-1][i+1]; ++i2) {
					if (val[k][j].first < val[k-1][i2].second) {
						// we have found the correct copy
						it->ind = i2; // set new index
						break; // we're done with this index
					}
				}
			}
			++it;
		}
		C[k][j].sort(); // sort now that we've modified indices
	}

	ZigzagChainComplex() {}

	/**
	Construct a zigzag chain complex from a zigzag filtration

	constructs a distinct column for every time a cell enters
	*/
	template <typename CpxT>
	ZigzagChainComplex(
		const ZigzagFiltration<CpxT, T>& X
	) : C(X.complex()) {
		val.resize(C.maxdim() + 1);
		cind.resize(C.maxdim() + 1);

		// map to where extra cells are
		std::vector<std::vector<size_t>> extra_cells(C.maxdim()+1);

		// loop over each dimension
		// duplicate cells when there are multiple entry and exit times
		// correct the boundary for mulitple entry and exit times
		for (size_t k = 0; k < X.maxdim()+1; ++k) {
			auto& Xvalsk = X.vals(k);
			val[k].resize(X.ncells(k));
			cind[k].resize(X.ncells(k));
			extra_cells[k].resize(X.ncells(k)+1);
			// loops over all cells in dimension k, creates duplicates as needed.
			for (size_t i = 0; i < X.ncells(k); ++i) {
				auto it = Xvalsk[i].begin();
				val[k][i] = *it;
				cind[k][i] = i;
				++it;
				// now we append any additional copies of the cell to the end
				extra_cells[k][i] = val[k].size(); // this behaves like a pointer range
				// the extra cells for cell i are in the range [extra_cells[k][i], extra_cells[k][i+1])
				while (it != Xvalsk[i].end()) {
					C[k].append_column(C[k][i]); // duplicate column
					cind[k].emplace_back(i); // map to original index
					// add row one dimension up if appropriate
					if(k < X.maxdim()) {C[k+1].append_row();}
					// add to val[k]
					val[k].emplace_back(*it++);
					// modify indices of column to map to correct cells
					if (k > 0) {
						// correct the indices of the last column
						_correct_indices(k, C[k].ncol()-1, extra_cells);
					}
				}
			}
			extra_cells[k][X.ncells(k)] = val[k].size(); // total number of cells
		}
	}

	/**
	return maximum dimension of cells
	*/
	inline size_t maxdim() const { return C.maxdim(); }

	/**
	return number of cells in specified dimension

	@param dim	dimension
	*/
	inline size_t dim(const size_t k) const { return C.dim(k); }
	inline size_t dim() const { return C.dim(); }

	/**
	return const reference to right filtration values
	*/
	inline const std::vector<std::vector<std::pair<T, T>>>& vals() const { return val;}

	/**
	return const reference to right filtration values

	@param k	dimension of values to return
	*/
	inline const std::vector<std::pair<T,T>>& vals(const size_t k) const { return val[k]; }

};


// permute complex to processing order for reduction
// and determine order in which to process everything
template <typename CpxT, typename T, typename FT>
auto prepare_ChainComplex(
	const ZigzagFiltration<CpxT, T>& F,
	FT // field for coefficients
) {
	using VecT = SparseVector<FT>;
	using MatT = ColumnMatrix<VecT>;
	auto ZC = ZigzagChainComplex<MatT, T>(F);
	// chain complex
	// auto C = Chain(F.complex(), FT());

	// compute permutations to order chain complex
	// order is by decreasing removal time
	std::vector<std::vector<size_t>> perm(ZC.maxdim() + 1);
	for (size_t k = 0; k < ZC.maxdim() + 1; ++k) {

		perm[k] = bats::util::sortperm(
			ZC.vals(k).begin(),
			ZC.vals(k).end(),
			[&](const std::pair<T,T>& a, const std::pair<T,T>& b) { return a.second > b.second; }
	 	);

	}
	ZC.C.ipermute_basis(perm);

	// compute order to process column
	// at same filtration value:
	// additions should be in increasing dimension order
	// removals should be in decreasing dimension order
	std::vector<rfilt_val<T>> filt_order;
	filt_order.reserve(ZC.dim() * 2);
	for (size_t k = 0; k < ZC.maxdim() + 1; ++k) {
		// auto permk = perm[k];
		auto permk = bats::util::inv_perm(perm[k]);
		for (size_t i = 0; i < ZC.dim(k); ++i) {
			auto pair = ZC.vals(k)[i];
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

	return std::make_tuple(ZC.C, filt_order);

}

namespace detail {

/**
Compare simplices in lexicographical order
looking at largest vertex first
(v_0,...,v_p) < (w_0,...,w_q) if v_p < w_q

compare simplex i in dimension dimi with simplex j in dimension dimj
@param X simplicial complex
@param dimi dimension of first simplex
@param i index of first simplex
@param dimj dimension of second simplex
@param j index of second simplex
@return true if first simplex < second simplex, false otherwise
*/
bool lex_cmp(
	const bats::SimplicialComplex& X,
	size_t dimi, size_t i, size_t dimj, size_t j
) {
	// std::cout << "in lex_cmp" << std::endl;
	auto iti = X.simplex_end(dimi, i); --iti;
	auto itj = X.simplex_end(dimj, j); --itj;
	// make quick call
	if (*iti != *itj) {return *iti < *itj;}

	while (iti != X.simplex_begin(dimi, i) && itj != X.simplex_begin(dimj, j)) {
		--iti; --itj;
		if (*iti != *itj) { return *iti < *itj; } // we can make a comparison
		// else we continue
	}
	// if we make it this far, simplices argree for full postfix of shorter simplex
	// if first simplex is shorter, we return true.  Else, we return false
	return dimj > dimi; // true if second simplex is longer than first simplex
	// if both simplices are equal, then will return false
}

// reverse of above
bool rlex_cmp(
	const bats::SimplicialComplex& X,
	size_t dimi, size_t i, size_t dimj, size_t j
) {
	// std::cout << "in rlex_cmp" << std::endl;
	auto iti = X.simplex_end(dimi, i); --iti;
	auto itj = X.simplex_end(dimj, j); --itj;
	// make quick call
	if (*iti != *itj) {return *iti > *itj;}

	while (iti != X.simplex_begin(dimi, i) && itj != X.simplex_begin(dimj, j)) {
		--iti; --itj;
		if (*iti != *itj) { return *iti > *itj; } // we can make a comparison
		// else we continue
	}
	// if we make it this far, simplices argree for full postfix of shorter simplex
	// if first simplex is shorter, we return true.  Else, we return false
	return dimi > dimj; // true if second simplex is longer than first simplex
	// if both simplices are equal, then will return false
}

} // namespace detail

// permute complex to processing order for reduction
// and determine order in which to process everything
template <typename T, typename FT>
auto prepare_ChainComplex(
	const ZigzagFiltration<bats::SimplicialComplex, T>& F,
	FT // field for coefficients
) {

	using VecT = SparseVector<FT>;
	using MatT = ColumnMatrix<VecT>;
	auto ZC = ZigzagChainComplex<MatT, T>(F); // chain complex
	auto& X = F.complex(); // underlying complex

	// compute permutations to order chain complex
	// order is by decreasing removal time
	// TODO: need to map column indices back to indices in X
	std::vector<std::vector<size_t>> perm(ZC.maxdim() + 1);
	for (size_t k = 0; k < ZC.maxdim() + 1; ++k) {
		auto& valsk = ZC.vals(k);
		perm[k] = bats::util::sortperm(
			size_t(0),
			valsk.size(),
			[&](const size_t i, const size_t j) {
				return valsk[i].second == valsk[j].second ? detail::rlex_cmp(X, k, ZC.cind[k][j], k, ZC.cind[k][i]) : valsk[i].second > valsk[j].second;
			}
	 	);

	}
	ZC.C.ipermute_basis(perm);

	// compute order to process column
	// at same filtration value:
	// additions should be in increasing dimension order
	// removals should be in decreasing dimension order
	std::vector<rfilt_val<T>> filt_order;
	filt_order.reserve(ZC.dim() * 2);
	for (size_t k = 0; k < ZC.maxdim() + 1; ++k) {
		// auto permk = perm[k];
		auto permk = bats::util::inv_perm(perm[k]);
		for (size_t i = 0; i < ZC.dim(k); ++i) {
			auto pair = ZC.vals(k)[i];
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
				detail::lex_cmp(X, a.dim, ZC.cind[a.dim][a.cind], b.dim, ZC.cind[b.dim][b.cind]) : // use lexicographical comparison
				(a.entry ^ b.entry) ? // one entry and one removal
				a.entry : // process entry first - doesn't matter
				(a.dim == b.dim) ? // else, both are removals - check if dim is the same
				a.ind > b.ind : // process later index first
				// a.dim > b.dim; // otherwise, we process larger dimension first for removal
				detail::rlex_cmp(X, a.dim, ZC.cind[a.dim][a.cind], b.dim, ZC.cind[b.dim][b.cind]);
			}
	);

	return std::make_tuple(ZC.C, filt_order);

}

template <typename CpxT, typename T, typename FT, typename opt_flag, typename reduction_flag>
auto barcode(
	const ZigzagFiltration<CpxT, T>& F,
	ssize_t maxdim,
	FT, // field for coefficients
	opt_flag,
	reduction_flag
) {
	auto [C, filt_order] = prepare_ChainComplex(F, FT());

	return zigzag_barcode_reduction(C, filt_order, maxdim, opt_flag(), reduction_flag());
}

} // namespace zigzag
} // namespace bats
