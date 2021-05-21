#pragma once
#include <homology/reduction.hpp> // flags
// #include <type_traits> // std::is_same
/*
Zigzag reduction algorithm
*/
#include <vector>
//#include <persistence/barcode.hpp> // PersistencePair

namespace bats {


/**
reduce column j past the pivot.
Heuristic strategy to reduce fill-in during reduction
*/
template <class VecT>
void extra_col_reduction(
	const size_t j,
	ColumnMatrix<VecT> &M,
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

/**
reduce a column of M.
Will eliminate pivots in column using columns to the left.
If a pivot is shared by column to right, we will continue reduction on
the column to the right.

@param j    column to reduce
@param M    partially reduced matrix
@param p2c  maps pivots to columns for reduction
@return     final column which was updated
*/
template <typename VecT, typename reduction_flag>
size_t reduce_column(
	size_t j,
	ColumnMatrix<VecT>& M,
	std::vector<size_t>& p2c,
	typename VecT::tmp_type& tmp,
	reduction_flag
) {

	while(M[j].nnz() > 0) {

		auto piv = M[j].lastnz();
		if (p2c[piv.ind] != bats::NO_IND) {
			size_t k = p2c[piv.ind];
			// std::cout << j << ": "; M[j].print_row();
			// std::cout << k << ": "; M[k].print_row();
			if (k > j) {
				// std::cout << k << ',' << j << std::endl;
				// switch to reducing column to right
				p2c[piv.ind] = j;
				std::swap(k, j); // swap values of k and j
				// std::cout << k << ',' << j << std::endl;
				piv = M[j].lastnz(); // update pivot after swap
			}
			// eliminate pivot in column j
			auto a = piv.val / M[k].lastnz().val;
			M[j].axpy(-a, M[k], tmp);
		} else {
			// new pivot
			p2c[piv.ind] = j;
			break;
		}
	}

    return j;
}
template <typename VecT>
size_t reduce_column(
	size_t j,
	ColumnMatrix<VecT>& M,
	std::vector<size_t>& p2c,
	typename VecT::tmp_type& tmp,
	extra_reduction_flag
) {

	while(M[j].nnz() > 0) {

		extra_col_reduction(j, M, p2c, tmp);
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

/**
@brief computes zigzag barcode

Computes a zigzag barcode given a chain complex, entry times, and exit times

Computes reduced matrices directly, and assumes that they have already been
permuted into the correct order.

TODO: make struct ZigzagPair which can keep track of
whether births/deaths occured because of addition/removal
*/
template <typename MT, typename T, typename opt_flag, typename reduction_flag>
auto zigzag_barcode_reduction(
	ChainComplex<MT>& C,
	const std::vector<rfilt_val<T>>& filt_order,
	opt_flag,
	reduction_flag
) {

	using VecT = typename MT::col_type;
	typename VecT::tmp_type tmp;

	std::vector<std::vector<size_t>> p2c(C.maxdim() + 1);
	for (size_t k = 0; k < C.maxdim() + 1; k++) {
		p2c[k].resize(C[k].nrow(), bats::NO_IND);
	}

	std::vector<std::vector<ZigzagPair<T>>> finished_pairs(C.maxdim() + 1); // initialize finished pairs
	std::vector<std::unordered_map<size_t, ZigzagPair<T>>> unfinished_pairs(C.maxdim() + 1);

	for (auto& fval : filt_order) {
		if (fval.entry) {
			// std::cout << "adding col " << fval.ind << " in dim " << fval.dim <<std::endl;
			// add column by processing entry of boundary matrix
			size_t j = reduce_column(fval.ind, C[fval.dim], p2c[fval.dim], tmp, reduction_flag());
			// determine if homology was created or destroyed
			if (C[fval.dim][j].nnz() > 0) {
				// pivot in final column - homology destroyed
				// look up pivot of this column
				auto piv = C[fval.dim][j].lastnz();
				// the added column kills homology generated by the pivot index
				ZigzagPair<T>& pair = unfinished_pairs[fval.dim-1][piv.ind];
				pair.death_ind = fval.cind; // death mapped back to original cell
				pair.death = fval.val; // actual death value
				pair.death_is_entry = true; // actual death is entry
				// put in finished pairs
				finished_pairs[fval.dim-1].emplace_back(pair);
				// erase from unfinished_pairs
				unfinished_pairs[fval.dim-1].erase(piv.ind);
			} else {
				// final column is cleared - homology created
				unfinished_pairs[fval.dim][j] = ZigzagPair<T>(
					fval.dim, // homology dimension
					fval.cind, // birth - original cell index
					bats::NO_IND, // death - not set yet
					fval.val, // birth parameter
					T(0), // death parameter - not set yet.
					true, // birth is entry
					true // death is entry - not set yet
				);
			}
		} else {
			// std::cout << "deleting col " << fval.ind << " in dim " << fval.dim <<std::endl;
			// we are deleting a column
			// determine if this will create or destroy homology
			if (C[fval.dim][fval.ind].nnz() > 0) {
				// homology will be created since we remove column with pivot
				// look up pivot of this column
				auto piv = C[fval.dim][fval.ind].lastnz();
				// the removed column creates homology generated by the pivot index
				unfinished_pairs[fval.dim-1][piv.ind] = ZigzagPair<T>(
					fval.dim-1, // homology dimension
					fval.cind, // birth - original cell index. WARNING - this is not in the same dimension
					bats::NO_IND, // death - not set yet
					fval.val, // birth parameter
					T(0), // death parameter - not set yet.
					false, // birth_is_entry - false because birth from removal
					true // death is entry - not set yet
				);
				// we now want to clear the p2c information for this entry
				// before deletion
				p2c[fval.dim][piv.ind] = bats::NO_IND;
			} else {
				// homology will be destroyed since we remove zero column
				ZigzagPair<T>& pair = unfinished_pairs[fval.dim][fval.ind];
				pair.death_ind = fval.cind; // death mapped back to cell index
				pair.death = fval.val; // death value
				pair.death_is_entry = false; // death from removal
				// std::cout << "finished pair: " << pair.dim << ',' << pair.birth_ind << std::endl;
				// put in finished pairs
				finished_pairs[fval.dim].emplace_back(pair);
				// erase from unfinished pairs
				unfinished_pairs[fval.dim].erase(fval.ind);
				// this column should not have been recorded in p2c
			}
			C[fval.dim][fval.ind].clear();
		}
	}

	// finally, go through and put any unfinished pairs into finished pairs

	return finished_pairs;

}

// TODO: use some metaprogramming with bats::compression_flag to merge with above
template <typename MT, typename T, typename reduction_flag>
auto zigzag_barcode_reduction(
	ChainComplex<MT>& C,
	const std::vector<rfilt_val<T>>& filt_order,
	compression_flag,
	reduction_flag
) {

	using VecT = typename MT::col_type;
	typename VecT::tmp_type tmp;

	std::vector<std::vector<size_t>> p2c(C.maxdim() + 1);
	for (size_t k = 0; k < C.maxdim() + 1; ++k) {
		p2c[k].resize(C[k].nrow(), bats::NO_IND);
	}
	// compression indices
	std::vector<std::vector<bool>> compression_inds(C.maxdim() + 1);
	for (size_t k = 0; k < C.maxdim() + 1; ++k) {
		compression_inds[k].resize(C[k].nrow(), false); // start with no compression
	}


	std::vector<std::vector<ZigzagPair<T>>> finished_pairs(C.maxdim() + 1); // initialize finished pairs
	std::vector<std::unordered_map<size_t, ZigzagPair<T>>> unfinished_pairs(C.maxdim() + 1);

	for (auto& fval : filt_order) {
		if (fval.entry) {
			// std::cout << "adding col " << fval.ind << " in dim " << fval.dim <<std::endl;
			// use compression indices
			C[fval.dim][fval.ind].clear_inds(compression_inds[fval.dim], tmp);
			// add column by processing entry of boundary matrix
			size_t j = reduce_column(fval.ind, C[fval.dim], p2c[fval.dim], tmp, reduction_flag());
			// determine if homology was created or destroyed
			if (C[fval.dim][j].nnz() > 0) {
				// pivot in final column - homology destroyed
				// look up pivot of this column
				auto piv = C[fval.dim][j].lastnz();
				// add pivot index to compression_inds one dimension up
				if (fval.dim < C.maxdim())
					compression_inds[fval.dim+1][piv.ind] = true;

				// the added column kills homology generated by the pivot index
				ZigzagPair<T>& pair = unfinished_pairs[fval.dim-1][piv.ind];
				pair.death_ind = fval.cind; // death mapped back to original cell
				pair.death = fval.val; // actual death value
				pair.death_is_entry = true; // actual death is entry
				// put in finished pairs
				finished_pairs[fval.dim-1].emplace_back(pair);
				// erase from unfinished_pairs
				unfinished_pairs[fval.dim-1].erase(piv.ind);
			} else {
				// final column is cleared - homology created
				unfinished_pairs[fval.dim][j] = ZigzagPair<T>(
					fval.dim, // homology dimension
					fval.cind, // birth - original cell index
					bats::NO_IND, // death - not set yet
					fval.val, // birth parameter
					T(0), // death parameter - not set yet.
					true, // birth is entry
					true // death is entry - not set yet
				);
			}
		} else {
			// std::cout << "deleting col " << fval.ind << " in dim " << fval.dim <<std::endl;
			// we are deleting a column
			// determine if this will create or destroy homology
			if (C[fval.dim][fval.ind].nnz() > 0) {
				// homology will be created since we remove column with pivot
				// look up pivot of this column
				auto piv = C[fval.dim][fval.ind].lastnz();
				// the removed column creates homology generated by the pivot index
				unfinished_pairs[fval.dim-1][piv.ind] = ZigzagPair<T>(
					fval.dim-1, // homology dimension
					fval.cind, // birth - original cell index. WARNING - this is not in the same dimension
					bats::NO_IND, // death - not set yet
					fval.val, // birth parameter
					T(0), // death parameter - not set yet.
					false, // birth_is_entry - false because birth from removal
					true // death is entry - not set yet
				);
				// we now want to clear the p2c information for this entry
				// before deletion
				p2c[fval.dim][piv.ind] = bats::NO_IND;
			} else {
				// homology will be destroyed since we remove zero column
				ZigzagPair<T>& pair = unfinished_pairs[fval.dim][fval.ind];
				pair.death_ind = fval.cind; // death mapped back to cell index
				pair.death = fval.val; // death value
				pair.death_is_entry = false; // death from removal
				// std::cout << "finished pair: " << pair.dim << ',' << pair.birth_ind << std::endl;
				// put in finished pairs
				finished_pairs[fval.dim].emplace_back(pair);
				// erase from unfinished pairs
				unfinished_pairs[fval.dim].erase(fval.ind);
				// this column should not have been recorded in p2c
			}
			C[fval.dim][fval.ind].clear();
		}
	}

	// finally, go through and put any unfinished pairs into finished pairs

	return finished_pairs;

}

} // namespace bats
