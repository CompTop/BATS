#pragma once

#include <vector>

namespace bats {

/*
Helper functions for sort permutations
*/
// get sort permutation for filtration in dimension dim
template <typename T>
inline std::vector<size_t> filtration_sortperm(
	const std::vector<T> &v
) {
	return bats::util::stable_sortperm(v);
}

template <typename T>
std::vector<std::vector<size_t>> filtration_sortperm(
	const std::vector<std::vector<T>> &v
) {
	std::vector<std::vector<size_t>> perms;
	for (size_t dim = 0; dim < v.size(); dim++) {
		perms.emplace_back(filtration_sortperm(v[dim]));
	}
	return perms;
}

std::vector<std::vector<size_t>> filtration_iperm(
	const std::vector<std::vector<size_t>> &perms
) {
	std::vector<std::vector<size_t>> iperms;
	for (auto p : perms) {
		iperms.emplace_back(bats::util::inv_perm(p));
	}
	return iperms;
}


/**
@brief A filtration which can be used to wrap a simplicial/cubical/cell complex

A filtration class, templated over the type of the filtration parameter, and
the type of the underlying complex
*/
template <typename TF, class CpxT>
class Filtration
{
private:
	CpxT X; // complex type
    std::vector<std::vector<TF>> val; // filtration value for each cell

    inline void reserve(const size_t maxdim) {
        while(val.size() < maxdim+1) {val.emplace_back(std::vector<TF>()); }
        return;
    }
    void reserve(const size_t dim, const size_t n) {
        reserve(dim);
        if (val[dim].size() < n) {
          val[dim].resize(n);
        }
    }

public:


	Filtration() {}

    // initialize on complex
    Filtration(const CpxT &C) : X(C) {
        for (size_t dim = 0; dim < X.maxdim() + 1; dim++){
			// std::cout << "reserving " << P.ncells(dim) << std::endl;
			reserve(dim, X.ncells(dim));
		}
    }

	/**
	Initialization which passes arguments to initialize the underlying complex

	@param args...	passed to complex initialization
	*/
	template <class ...Ts>
	Filtration(const Ts (&...args)) : X(args...) {
		for (size_t dim = 0; dim < X.maxdim() + 1; dim++){
			// std::cout << "reserving " << P.ncells(dim) << std::endl;
			reserve(dim, X.ncells(dim));
		}
	}

	/**
	Initialization which passes arguments to initialize the underlying complex
	@param C		complex
	@param vals		filtration values for every cell in C
	*/
	Filtration(const CpxT& C, const std::vector<std::vector<TF>>& vals) : X(C), val(vals) {}



	/**
	return const reference to underlying complex
	*/
	inline const CpxT& complex() const { return X; }

	/**
	return const reference to filtration values
	*/
	inline const std::vector<std::vector<TF>>& vals() const { return val;}

	/**
	return const reference to filtration values

	@param k	dimension of values to return
	*/
	inline const std::vector<TF>& vals(const size_t k) const { return val[k]; }

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
	add cell to filtration

	@param t		filtration parameter
	@param ...args	passed to add method of underlying complex
	*/
	template <class ...Ts>
	inline cell_ind add(TF t, Ts (&...args)) {
		cell_ind ret = X.add(args...);
		reserve(ret.dim, ret.ind+1);
		val[ret.dim][ret.ind] = t;
		return ret;
	}


	/**
	Add recursively to filtration.  Any cells added will take filtration
	value t.

	@param t		filtration parameter
	@param ...args	passed to add_recursive method of underlying complex
	*/
	template <class ...Ts>
	inline std::vector<cell_ind> add_recursive(TF t, Ts (&...args)) {
		std::vector<cell_ind> ret = X.add_recursive(args...);
		for (auto &r : ret) {
			reserve(r.dim, r.ind+1);
			val[r.dim][r.ind] = t;
		}
		return ret;
	}

	// get sort permutation for filtration in dimension dim
	inline std::vector<size_t> sortperm(size_t dim) const {
		return filtration_sortperm(val[dim]);
	}

	// get sort permutation for filtration in all dimensions
	inline std::vector<std::vector<size_t>> sortperm() const {
		return filtration_sortperm(val);
	}

	/**
	Get sub-levelset of filtration (-inf, a]

	@param a upper bound of levelset
	*/
	CpxT sublevelset(const TF a) const {
		CpxT Y(X.ncells(0), X.maxdim()); // initialize new complex
		for (size_t k = 0; k < maxdim() + 1; ++k) {
			for (size_t i = 0; i < ncells(k); ++i) {
				if (val[k][i] <= a) {
					Y.add(X.get_cell(k, i));
				}
			}
		}
		return Y;
	}

};

} // namespace bats
