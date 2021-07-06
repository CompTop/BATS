#pragma once

#include <vector>
#include <limits>
#include <chain/filtered_chain_complex.hpp>
#include <homology/basis.hpp>
#include <util/common.hpp>
#include "barcode.hpp"
#include <filtration/update_information.hpp>

namespace bats {

template <typename T, typename MT>
struct ReducedFilteredChainComplex {

	ReducedChainComplex<MT> RC; // reduced in permutation order
	std::vector<std::vector<T>> val; // stored in original order
	std::vector<std::vector<size_t>> perm; //from permutaiton order to original order

	ReducedFilteredChainComplex() {}

	// variadic template for passing arguments
	template <typename... Args>
	ReducedFilteredChainComplex(const FilteredChainComplex<T, MT>& C, Args ...args) :
		RC(C.complex(), args...),
		val(C.val),
		perm(C.perm) {}

	inline size_t maxdim() const { return RC.maxdim(); }
	inline size_t dim(const size_t k) const {return RC.dim(k);}
	inline size_t hdim(const size_t k) const {return RC.hdim(k);}

	/**
	persistence pairs in dimension k

	@param k	homology dimension
	@param permuted	set to true to return critical indices permuted by filtration parameter
	set to false to return with indices in original order.  Default: false
	*/
	std::vector<PersistencePair<T>> persistence_pairs(
		const size_t k,
		const bool permuted=false
	)const {
		std::vector<PersistencePair<T>> pairs;
		if (permuted) {
			for (size_t i =0; i < dim(k); i++) {
				if (RC.R[k][i].nnz() == 0) {
					// homology generated
					if (k == maxdim() || RC.p2c[k+1][i] == bats::NO_IND)  {
						// infinite bar
						pairs.emplace_back(
							PersistencePair(k, i, bats::NO_IND,
								val[k][perm[k][i]], std::numeric_limits<T>::infinity()
							)
						);
					} else {
						size_t j = RC.p2c[k+1][i];
						// finite bar
						pairs.emplace_back(
							PersistencePair(k, i, j,
								val[k][perm[k][i]], val[k+1][perm[k+1][j]]
							)
						);
					}
				}
			}
		} else { // undo filtration sort permutation
			for (size_t i =0; i < dim(k); i++) {
				if (RC.R[k][i].nnz() == 0) {
					// homology generated
					if (k == maxdim() || RC.p2c[k+1][i] == bats::NO_IND)  {
						// infinite bar
						pairs.emplace_back(
							PersistencePair(k, perm[k][i], bats::NO_IND,
								val[k][perm[k][i]], std::numeric_limits<T>::infinity()
							)
						);
					} else {
						size_t j = RC.p2c[k+1][i];
						// finite bar
						pairs.emplace_back(
							PersistencePair(k, perm[k][i], perm[k+1][j],
								val[k][perm[k][i]], val[k+1][perm[k+1][j]]
							)
						);
					}
				}
			}
		}

		return pairs;
	}

	/**
	returns representative for homology class corresponding to persistence pair

	@param p	persistence pair obtained from persistence_pairs
	@param permuted	set to true to return indices permuted by filtration parameter
	set to false to return with indices in original order.  Default: false

	if permuted is false, it is assumed the birth index is also not in permutation order.
	*/
	inline auto representative(const PersistencePair<T>& p, const bool permuted=false) const {
		auto v = RC.U[p.dim][p.birth_ind];
		if (!permuted) {
			auto it = std::find(perm[p.dim].begin(), perm[p.dim].end(), p.birth_ind);
			v = RC.U[p.dim][it - perm[p.dim].begin()];
			v.permute(perm[p.dim]); // undo permutation
		}
		return v;
	}

	// barcode w/out critical inds
	std::vector<T> barcode(const size_t k);

	// critical cells for barcode in dimension k
	std::vector<size_t> critical_cells(const size_t k);

	// update filtration
	void update_filtration(const std::vector<std::vector<T>> newval) {
		// step 1: determine permutation order for newval
		auto perms = filtration_sortperm(newval);

		// get inverse permutation for current permutation
		auto iperm = filtration_iperm(perm);

		// step 2: determine update to old permutation
		// iperm[k] will hold the updated permutation temporarily
		std::vector<size_t> tmp; // temporary vector
		for (size_t k = 0; k < iperm.size(); k++) {
			bats::util::apply_perm(iperm[k].data(), tmp, perms[k]);
		}

		// step 3: apply permutation updates to ReducedChainComplex RC
		RC.permute_basis(iperm);

		// step 4: store new permutation
		perm = perms;

		// store values as well
		val = newval;

	}

	// update filtration fast version
	template <typename Information_type, typename... Args>
	void update_filtration_general(
		const Information_type & updating_information,
		Args ...args
	) {
		// step 1: apply permutation updates to ReducedChainComplex RC
		RC.update_basis_general(updating_information, args...);

		// step 2: store new permutation
		perm = updating_information.F_Y_perms;

		// store values as well
		val = updating_information.F_Y_vals;

	}

	// get nnz of U
	std::vector<size_t> get_nnz_U(){
		std::vector<size_t> U_nnz;
		U_nnz.reserve(maxdim());
		for (size_t k = 0; k < maxdim() + 1; k++) {
			U_nnz.emplace_back(RC.U[k].nnz());
		}
		return U_nnz;
	}

	// get nnz of R
	std::vector<size_t> get_nnz_R(){
		std::vector<size_t> R_nnz;
		R_nnz.reserve(maxdim());
		for (size_t k = 0; k < maxdim() + 1; k++) {
			R_nnz.emplace_back(RC.R[k].nnz());
		}
		return R_nnz;
	}

	/**
	greedily introduce sparsity into basis
	*/
	inline void sparsify_basis() {
		RC.sparsify_basis();
	}

	/**
	remove extra cycles from U[k]
	*/
	inline void remove_extra_cycles() {
		RC.remove_extra_cycles();
	}

	// get subcomplex
	ReducedFilteredChainComplex get_subcomplex() const;

	void print_summary(bool print_nnz=false) const {
		std::cout << "ReducedFilteredChainComplex with " << maxdim() << " dimensions:" << std::endl;
		for (size_t k = 0; k < maxdim() + 1; k++) {
			std::cout << "\tdim " << k << ": " << dim(k)
			<< ", betti_" << k << ": " << hdim(k);
			if (print_nnz) {
				std::cout << " nnz(R): " << RC.R[k].nnz();
				if (RC.U.size() > k) { // handle if basis not computed
					std::cout << " nnz(U): " << RC.U[k].nnz();
				}
			}
			std::cout << "\n";
		}
	}

};

// check if two RFCC are the same
// (can be modified into new == operator of RFCC class)
template<typename T>
bool test_reduce_result(const T& RFCC2, const T& RFCC){

    auto R = RFCC.RC.R;
    auto R2 = RFCC2.RC.R;
    auto U = RFCC.RC.U;
    auto U2 = RFCC2.RC.U;


    // check persistent homology information
    bool test_result = true;
    for(size_t k =0; k < RFCC2.maxdim(); k++){
        auto R_ps_k = RFCC.persistence_pairs(k);
        auto R2_ps_k = RFCC2.persistence_pairs(k);
        bool test_result_k = barcode_equality(R2_ps_k, R_ps_k);
        if(!test_result_k){
            test_result = false;
            std::cout << "\ntwo RFCC are different!!" << std::endl;
        }
    }

    // check RU factorization
    for(size_t k = 0; k < RFCC2.maxdim(); k++){
        auto U_inv_k = u_inv(U[k]);
        auto U2_inv_k = u_inv(U2[k]);
        if(!(R[k] * U_inv_k == R2[k] * U2_inv_k)){
            test_result = false;
            std::cout << "\ntwo RFCC are different!!" << std::endl;
        }
    }

    return test_result;
}

// defualt return
template <typename FT, typename T, typename CpxT, typename... Args>
inline auto __ReducedFilteredChainComplex(const Filtration<FT, CpxT> &F, T, Args ...args) {
	using VT = SparseVector<T, size_t>;
	using MT = ColumnMatrix<VT>;

	return ReducedFilteredChainComplex(FilteredChainComplex<FT, MT>(F), args...);
}

template <typename FT, typename T, typename CpxT, typename... Args>
inline auto Reduce(const Filtration<FT, CpxT> &F, T, Args ...args) {
	return __ReducedFilteredChainComplex(F, T(), args...);
}

template <typename T, typename MT, typename... Args>
inline auto Reduce(const FilteredChainComplex<T, MT>& C, Args ...args) {
	return ReducedFilteredChainComplex(C, args...);
}


} // namespace bats
