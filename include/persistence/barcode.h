#pragma once

#include <string>
#include <iostream>

// store dimension, birth, death, and critical indices of pair
template <typename T>
struct PersistencePair {
	size_t dim;
	size_t birth_ind;
	size_t death_ind;
	T birth;
	T death;

	PersistencePair() {}
	PersistencePair(
		const size_t dim,
		const size_t birth_ind,
		const size_t death_ind,
		const T birth,
		const T death
	) : dim(dim), birth_ind(birth_ind), death_ind(death_ind), birth(birth), death(death) {}

	inline size_t get_dim() const {return dim;}
	inline size_t get_birth_ind() const {return birth_ind;}
	inline size_t get_death_ind() const {return death_ind;}
	inline T get_birth() const {return birth;}
	inline T get_death() const {return death;}

	std::string str() {
        std::ostringstream oss;
        oss << dim << " : (" << birth << ',' << death << ") <" << birth_ind << ',' << (death_ind == bats::NO_IND ? (int) -1 : (int) death_ind) << '>';
        return oss.str();
    }
};

// // template over filtration type, field type, complex type
// template <typename T, typename FT, typename CpxT>
// auto PersistencePairs(const Filtration<FT, CpxT> &F, FT, const size_t dim) {
// 	using VT = SparseVector<FT, size_t>;
// 	using MT = ColumnMatrix<VT>;
//
// 	auto FC = FilteredChainComplex<T, MT>(F);
// 	auto RFC = ReducedFilteredChainComplex(FC);
//
// 	return RFC.persistence_pairs(dim);
// }
