#pragma once
#include <linalg/col_matrix.hpp>
#include <linalg/sparse_fact.hpp> // factorizations
#include <homology/reduction.hpp> // flags
// #include <type_traits> // std::is_same
#include <utility> // std::pair

/*
Zigzag reduction algorithm
*/
#include <vector>
//#include <persistence/barcode.hpp> // PersistencePair

namespace bats {
namespace zigzag {


// /**
// reduce column j past the pivot.
// Heuristic strategy to reduce fill-in during reduction
// */
// template <class VecT>
// void extra_col_reduction(
// 	const size_t j,
// 	ColumnMatrix<VecT> &M,
// 	const std::vector<size_t>& p2c,
// 	typename VecT::tmp_type& tmp
// ) {
//
// 	size_t end_offset = 1; // difference from end
// 	auto piv = M[j].nzend() - end_offset; // nonzero location
// 	while(piv - M[j].nzbegin() > 0) {
// 		// while the nonzero we are looking at is in the vector
// 		// piv is index-value nzpair
// 		size_t k = p2c[piv->ind];
// 		if (k != bats::NO_IND && k < j) {
// 			// eliminate pivot
// 			auto a = piv->val / M[k].lastnz().val;
// 			M[j].axpy(-a, M[k], tmp);
// 			piv = M[j].nzend() - end_offset; // next nonzero location
// 		} else {
// 			// we skip zeroing out this entry
// 			// need to increment offset
// 			end_offset++;
// 			piv--;
// 		}
// 	}
//
// 	return;
// }

/**
reduce a column of M.
Will eliminate pivots in column using columns to the left.
If a pivot is shared by column to right, we will continue reduction on
the column to the right.

@param j    column to reduce
@param M    partially reduced matrix
@param U	basis matrix
@param p2c  maps pivots to columns for reduction
@return     final column which was updated
*/
template <typename VecT, typename reduction_flag>
size_t reduce_column(
	size_t j,
	ColumnMatrix<VecT>& M,
	ColumnMatrix<VecT>& U,
	std::vector<size_t>& p2c,
	typename VecT::tmp_type& tmp,
	reduction_flag
) {

	while(M[j].nnz() > 0) {

		auto piv = M[j].lastnz();
		if (p2c[piv.ind] != bats::NO_IND) {
			size_t k = p2c[piv.ind];
			if (k > j) {
				// switch to reducing column to right
				p2c[piv.ind] = j;
				std::swap(k, j); // swap values of k and j
				piv = M[j].lastnz(); // update pivot after swap
			}
			// eliminate pivot in column j
			auto a = piv.val / M[k].lastnz().val;
			M[j].axpy(-a, M[k], tmp);
			U[j].axpy(-a, U[k], tmp); // apply same operation to basis U
		} else {
			// new pivot
			p2c[piv.ind] = j;
			break;
		}
	}

    return j;
}

/**
reduce column j past the pivot.
Heuristic strategy to reduce fill-in during reduction
*/
template <class VecT>
void extra_col_reduction(
	const size_t j,
	ColumnMatrix<VecT> &M,
	ColumnMatrix<VecT>& U,
	const std::vector<size_t>& p2c,
	typename VecT::tmp_type& tmp
) {

	size_t end_offset = 1; // difference from end
	auto piv = M[j].nzend() - end_offset; // nonzero location
	while(piv - M[j].nzbegin() > 0) {
		// while the nonzero we are looking at is in the vector
		// piv is index-value nzpair
		size_t k = p2c[piv->ind];
		if (k != bats::NO_IND && k < j) {
			// eliminate pivot
			auto a = piv->val / M[k].lastnz().val;
			M[j].axpy(-a, M[k], tmp);
			U[j].axpy(-a, U[k], tmp);
			piv = M[j].nzend() - end_offset; // next nonzero location
		} else {
			// we skip zeroing out this entry
			// need to increment offset
			end_offset++;
			piv--;
		}
	}

	return;
}

template <typename VecT>
size_t reduce_column(
	size_t j,
	ColumnMatrix<VecT>& M,
	ColumnMatrix<VecT>& U,
	std::vector<size_t>& p2c,
	typename VecT::tmp_type& tmp,
	bats::extra_reduction_flag
) {

	while(M[j].nnz() > 0) {

		extra_col_reduction(j, M, U, p2c, tmp);
		if (M[j].nnz() == 0) { break; } // zeroed out column
		auto piv = M[j].lastnz();
		size_t k = p2c[piv.ind];
		if (k != bats::NO_IND && k > j) {
			// continue reducing column to right which shares pivot
			p2c[piv.ind] = j;
			std::swap(k, j);
		} else {
			// new pivot
			p2c[piv.ind] = j;
			break;
		}
	}

    return j;
}

/**
A struct that packages information about entries and exits
in a right filtration
*/
template <typename T>
struct rfilt_val {
	size_t dim;   // dimension of cell
	size_t ind;   // index in permuted matrix
	size_t cind;  // index of cell in complex
	T val;        // value
	bool entry;   // true if entry, false if exit

	rfilt_val() {}
	rfilt_val(size_t dim, size_t ind, size_t cind, T val, bool entry)\
		: dim(dim), ind(ind), cind(cind), val(val), entry(entry) {}

	friend std::ostream& operator<<(
		std::ostream& os,
		const rfilt_val& v
	) {
	    os << "fval, dim="<<v.dim << ": ind=" << v.ind << ", cind=" << v.cind
			<< ", val=" << v.val << ", entry="<< v.entry;
	    return os;
	}

};

// store dimension, birth, death, and critical indices of pair
template <typename T>
struct ZigzagPair {
	size_t dim;
	size_t birth_ind;
	size_t death_ind;
	T birth;
	T death;
	bool birth_is_entry;
	bool death_is_entry;

	ZigzagPair() {}
	ZigzagPair(
		const size_t dim,
		const size_t birth_ind,
		const size_t death_ind,
		const T birth,
		const T death,
		const bool birth_is_entry,
		const bool death_is_entry
	) : dim(dim), birth_ind(birth_ind), death_ind(death_ind),
	 birth(birth), death(death),
	 birth_is_entry(birth_is_entry), death_is_entry(death_is_entry) {}

	// useful for checking equality of barcodes
	inline bool operator==(const ZigzagPair &other) const { return birth == other.birth && death == other.death && dim == other.dim; }
	inline bool operator!=(const ZigzagPair &other) const { return birth != other.birth || death != other.death || dim != other.dim; }

	inline size_t get_dim() const {return dim;}
	inline size_t get_birth_ind() const {return birth_ind;}
	inline size_t get_death_ind() const {return death_ind;}
	inline T get_birth() const {return birth;}
	inline T get_death() const {return death;}
	inline T length() const {return death - birth;}
	inline T mid() const {return (death + birth) / T(2);}


	friend std::ostream& operator<<(
		std::ostream& os,
		const ZigzagPair& p
	) {
	    os << p.dim << " : (" << p.birth << ',' << p.death << ") <"
			<< p.birth_ind << '(' << p.birth_is_entry << "),"
			<< (p.death_ind == bats::NO_IND ? (int) -1 : (int) p.death_ind) << '(' << p.death_is_entry << ")>";
	    return os;
	}

	std::string str() {
        std::ostringstream oss;
        oss << *this;
        return oss.str();
    }
};

namespace detail {
// create induced map by inserting new vector into HRB with index j
// I is HRB preferred rep set after j has been inserted
template <typename VecT>
ColumnMatrix<VecT> cycle_insertion_map(
	const std::vector<size_t>& I,
	const size_t j
) {
	size_t m = I.size();
	size_t n = m - 1;
	std::vector<VecT> col(n);
	for (size_t k = 0; k < n; ++k) {
		if (I[k] < j) {
			col[k] = VecT(k);
		} else {
			col[k] = VecT(k+1);
		}
	}
	return ColumnMatrix(m, n, col);
}

// create map induced by insertion of boundary that removed index i from HRB
template <typename VecT>
ColumnMatrix<VecT> boundary_insertion_map(
	const std::vector<size_t>& I,
	const size_t i,
	const VecT& v
) {
	size_t m = I.size();
	size_t n = m + 1;
	std::vector<VecT> col(n);
	bool seen_i = false;
	for (size_t k = 0; k < m; ++k) {
		if (I[k] < i) {
			col[k] = VecT(k);
		} else {
			if (!seen_i) {
				col[k] = v; // insert v
				seen_i = true;
			}
			col[k+1] = VecT(k);
		}
	}
	if (!seen_i) {
		// if last HRB is boundary pivot
		col[m] = v; // insert v
	}
	return ColumnMatrix(m, n, col);
}

// apply change of basis to matrix A
// direction bools:
// true  = -->
// false = <--
template <typename MT>
void apply_basis(MT& A, MT& L, MT& P, const bool prev_dir, const bool dir) {
	if (!prev_dir && !dir) {
		// * <-- * <-- *
		A = L * P * A;
	} else if (prev_dir && dir) {
		// * --> * --> *
		A = A * P * L;
	} else if (!prev_dir && dir) {
		// * <-- * --> *
		A = A * P.T() * l_inv(L);
	} else {
		// * --> * <-- *
		A = l_inv(L) * P.T() * A;
	}
	return;
}

/**
incrementally update barcode

see quiver/sparse.hpp barcode_from_barcode_form
*/
template <typename T, typename MT, typename Map>
void update_bars(
	std::vector<ZigzagPair<T>>& bars,
	const rfilt_val<T>& fval, // information about filtration
	const size_t hdim, // homology dimension we're updating
	MT& E, // E matrix in barcode form
	Map& piv_to_ind, // maps homology classes to birth bar
	Map& piv_to_ind2 // used to swap
) {

	piv_to_ind2.clear(); // updated pivots

	if (fval.entry) {
		// ->
		// extend bar by column
		// easier to work on transpose
		E = E.T();
		for (size_t j = 0; j < E.ncol(); j++) {
			auto piv = E[j].nzbegin();
			if (piv == E[j].nzend()) {
				// no extension to complete. start new bar
				// i.e. row was not in range of matrix
				piv_to_ind2[j] = bars.size();
				bars.emplace_back(ZigzagPair<T>(
					hdim, // homology dimension
					fval.cind, // birth - original cell index
					bats::NO_IND, // death - not set yet
					fval.val, // birth parameter
					std::numeric_limits<T>::max(), // death parameter - not set yet.
					fval.entry, // birth is entry
					true // death is entry - not set yet
				));
			} else {
				// extend bar that is in pivot row
				size_t pi = piv_to_ind[piv->ind];
				piv_to_ind2[j] = pi;
				piv_to_ind.erase(piv->ind);
			}
		}

	} else {
		// <-
		// extend bar by row
		for (size_t j = 0; j < E.ncol(); j++) {
			auto piv = E[j].nzbegin();
			if (piv == E[j].nzend()) {
				// no extension to complete. start new bar
				piv_to_ind2[j] = bars.size();
				bars.emplace_back(ZigzagPair<T>(
					hdim, // homology dimension
					fval.cind, // birth - original cell index
					bats::NO_IND, // death - not set yet
					fval.val, // birth parameter
					std::numeric_limits<T>::max(), // death parameter - not set yet.
					fval.entry, // birth is entry
					true // death is entry - not set yet
				));
			} else {
				// extend bar that is in pivot row
				size_t pi = piv_to_ind[piv->ind];
				piv_to_ind2[j] = pi;
				// erase in piv_to_ind
				piv_to_ind.erase(piv->ind);
			}
		}
	}
	// anything left in piv_to_ind was not erased
	// kill those bars
	for (auto it = piv_to_ind.begin(); it!= piv_to_ind.end(); ++it) {
		ZigzagPair<T>& pair = bars[it->second];
		pair.death_ind = fval.cind; // death mapped back to original cell
		pair.death = fval.val; // actual death value
		pair.death_is_entry = fval.entry; // actual death is entry
	}
	// swap roles of pivot lookups
	std::swap(piv_to_ind, piv_to_ind2);

}




} // namespace detail

/**
@brief computes zigzag barcode

Computes a zigzag barcode given a chain complex, entry times, and exit times

Computes reduced matrices directly, and assumes that they have already been
permuted into the correct order.

@param C	Chain complex (in zigzag reduction order)
@param filt_order	order in which zigzag filtration occurs
@param maxdim	maximum homology dimension to compute.
If set to -1, will default to computing all but top dimension.

TODO: don't compute top dimension homology - assume spurious
TODO: strategy for keeping homology dimensions small when adding a bunch
of cells at the same time.

*/
template <typename MT, typename T, typename opt_flag, typename reduction_flag>
auto zigzag_barcode_reduction(
	const ChainComplex<MT>& C,
	const std::vector<rfilt_val<T>>& filt_order,
	ssize_t maxdim,
	opt_flag,
	reduction_flag
) {

	using VecT = typename MT::col_type;
	typename VecT::tmp_type tmp;

	// handle maxdim
	if (maxdim == -1) {
		maxdim = C.maxdim() - 1;
	}

	ReducedChainComplex<MT> F; // maintain RU factorizations
	F.initialize(C); // initialize but do not reduce

	std::vector<std::vector<ZigzagPair<T>>> bars(maxdim + 1); // initialize barcode
	std::vector<std::unordered_map<size_t, size_t>> piv_to_ind(maxdim + 1), piv_to_ind2(maxdim + 1);
	std::vector<MT> prev_L(maxdim + 1); // used for applying running change of basis
	std::vector<MT> prev_P(maxdim + 1); // used for applying running change of basis
	std::vector<bool> prev_dir(maxdim + 1, false); // keep track of previous arrow direction


	for (auto& fval : filt_order) {
		size_t k = fval.dim; // dimension
		// std::cout << fval << std::endl;

		if (fval.entry) {
			// right arrow -->
			// PLEU factorization

			// add column by processing entry of boundary matrix
			size_t j = reduce_column(fval.ind, F.R[k], F.U[k], F.p2c[k], tmp, reduction_flag());
			// determine if homology was created or destroyed
			if ((k < maxdim + 2) && F.R[k][j].nnz() > 0) {
				// pivot in final column - homology destroyed one dimension down
				// look up pivot of this column
				auto piv = F.R[k][j].lastnz();
				// the added column kills homology generated by the pivot index
				size_t i = piv.ind;
				// delete from HRB basis in dimension k-1
				auto it = std::lower_bound(F.I[k-1].begin(), F.I[k-1].end(), i);
				if (i != *it) {throw std::runtime_error("pivot is not cycle!"); }
				F.I[k-1].erase(it);
				// determine induced map
				VecT v(i); // chain in HRB basis
				F.find_preferred_representative(v, k-1); // find preferred rep 1-dimension down
				v = v[F.I[k-1]]; // induced map
				MT induced_map = detail::boundary_insertion_map(F.I[k-1], i, v);
				// apply running change-of-basis
				detail::apply_basis(induced_map, prev_L[k-1], prev_P[k-1], prev_dir[k-1], true);
				// compute new factorization
				auto fact = PLEU(induced_map);
				// store L, P, and direction for next time
				prev_L[k-1] = fact.L;
				prev_P[k-1] = fact.P;
				prev_dir[k-1] = true; // --> arrow

				// need to look at fact.E to see what was killed
				detail::update_bars(bars[k-1], fval, k-1, fact.E, piv_to_ind[k-1], piv_to_ind2[k-1]);
				// delete bar from active_bars and transform initial indices of active_bars

			} else if (k < maxdim + 1) {
				// final column is cleared - homology created
				// determine where new vector inserted into HRB
				auto it = std::upper_bound(F.I[k].begin(), F.I[k].end(), j);
				// insert into list of indices
				it = F.I[k].insert(it, j);
				// induced map
				MT induced_map = detail::cycle_insertion_map<VecT>(F.I[k], j);
				// apply running change-of-basis
				detail::apply_basis(induced_map, prev_L[k], prev_P[k], prev_dir[k], true);
				// compute new factorization
				auto fact = PLEU(induced_map);
				// store L, P, and direction for next time
				prev_L[k] = fact.L;
				prev_P[k] = fact.P;
				prev_dir[k] = true; // --> arrow

				// need to transform indices using fact.E
				// add bar to active_bars and transform initial indices of active_bars
				detail::update_bars(bars[k], fval, k, fact.E, piv_to_ind[k], piv_to_ind2[k]);

			}
		} else {
			// left arrow <--
			// UELP factorization

			// we are deleting a column
			// determine if this will create or destroy homology
			if ((k < maxdim + 2) && F.R[k][fval.ind].nnz() > 0) {
				// homology will be created one dimension down since we remove column with pivot
				// look up pivot of this column
				auto piv = F.R[k][fval.ind].lastnz();
				size_t i = piv.ind;
				// the removed column creates homology generated by the pivot index
				VecT v(i); // chain in HRB basis
				F.find_preferred_representative(v, k-1); // find preferred rep 1-dimension down
				v = v[F.I[k-1]]; // induced map
				MT induced_map = detail::boundary_insertion_map(F.I[k-1], i, v);
				// add column i to HRB one dimension down
				auto it = std::upper_bound(F.I[k-1].begin(), F.I[k-1].end(), i);
				it = F.I[k-1].insert(it, i);

				// apply running change-of-basis
				detail::apply_basis(induced_map, prev_L[k-1], prev_P[k-1], prev_dir[k-1], false);
				// compute new factorization
				auto fact = UELP(induced_map);
				// store L, P, and direction for next time
				prev_L[k-1] = fact.L;
				prev_P[k-1] = fact.P;
				prev_dir[k-1] = false; // <-- arrow

				// need to transform indices using fact.E
				// add bar to active_bars and transform initial indices of active_bars
				detail::update_bars(bars[k-1], fval, k-1, fact.E, piv_to_ind[k-1], piv_to_ind2[k-1]);

				// remove i from list of pivots
				F.p2c[k][i] = bats::NO_IND;
			} else if (k < maxdim + 1) {
				// homology will be destroyed in this dimension since we remove zero column
				// induced map of insertion <--
				MT induced_map = detail::cycle_insertion_map<VecT>(F.I[k], fval.ind);
				// delete from HRB
				F.I[k].pop_back();

				// apply running change-of-basis
				detail::apply_basis(induced_map, prev_L[k], prev_P[k], prev_dir[k], false);
				// compute new factorization
				auto fact = UELP(induced_map);
				// store L, P, and direction for next time
				prev_L[k] = fact.L;
				prev_P[k] = fact.P;
				prev_dir[k] = false; // <-- arrow

				// need to look at fact.E to see what was killed
				detail::update_bars(bars[k], fval, k, fact.E, piv_to_ind[k], piv_to_ind2[k]);

				// this column should not have been recorded in p2c
			}
			F.R[fval.dim][fval.ind].clear_dealloc();
			F.U[fval.dim][fval.ind].clear_dealloc();
		}
	}

	return bars;

}


} // namespace zigzag
} // namespace bats
