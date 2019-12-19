#pragma once

#include <vector>


// template over
//  TC - complex type
//  TF - filtration type
template <typename TF, class CpxT>
class Filtration
{
private:
    // MorsePairing<CpxT> P;
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

    // TODO: complete pairs in MorsePairing
    // TODO: extract pairs in MorseParing to get barcode
	Filtration() {};

    // initialize on complex
    Filtration(const CpxT &C) : X(C) {
        for (size_t dim = 0; dim < X.maxdim() + 1; dim++){
			// std::cout << "reserving " << P.ncells(dim) << std::endl;
			reserve(dim, X.ncells(dim));
		}
    };

	inline const CpxT& complex() const { return X; }
	inline const std::vector<std::vector<TF>>& vals() const { return val;}
	inline const std::vector<TF>& vals(const size_t k) const { return val[k]; }

	inline size_t maxdim() const { return val.size() - 1; }
	inline size_t ncells(const size_t dim) const { return val[dim].size(); }

	// add to underlying complex
	template <class ...Ts>
	inline cell_ind add(TF t, Ts (&...args)) {
		cell_ind ret = X.add(args...);
		reserve(ret.dim, ret.ind+1);
		val[ret.dim][ret.ind] = t;
		return ret;
	}

};
