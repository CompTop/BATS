#pragma once

#include <string>
#include <iostream>
#include <vector>

namespace bats {

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

	// useful for checking equality of barcodes
	inline bool operator==(const PersistencePair &other) const { return birth == other.birth && death == other.death && dim == other.dim; }
	inline bool operator!=(const PersistencePair &other) const { return birth != other.birth || death != other.death || dim != other.dim; }

	inline size_t get_dim() const {return dim;}
	inline size_t get_birth_ind() const {return birth_ind;}
	inline size_t get_death_ind() const {return death_ind;}
	inline T get_birth() const {return birth;}
	inline T get_death() const {return death;}
	inline T length() const {return death - birth;}
	inline T mid() const {return (death + birth) / T(2);}

	std::string str() {
        std::ostringstream oss;
        oss << dim << " : (" << birth << ',' << death << ") <" << birth_ind << ',' << (death_ind == bats::NO_IND ? (int) -1 : (int) death_ind) << '>';
        return oss.str();
    }
};

template <typename T>
bool barcode_equality(
	const std::vector<PersistencePair<T>> &ps1,
	const std::vector<PersistencePair<T>> &ps2
) {
	// check for same number of pairs
	if (ps1.size() != ps2.size()) {
		return false;
	}

	size_t n = ps1.size();
	// vector for checking if a pair in ps2 has been matched
	std::vector<bool> matched(n, false);
	for (size_t i = 0; i < n; i++) {
		bool found_match = false;
		for (size_t j = 0; j < n; j++) {
			if (ps1[i] == ps2[j] && !matched[j]){
				found_match = true;
				matched[j] = true;
				break;
			}
		}
		if (!found_match) {
			return false; // no match for ps1[i]
		}
	}
	return true; // we found a match for every pair.

}

} // namespace bats
