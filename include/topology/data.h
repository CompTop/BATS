#pragma once

#include <linalg/naive_dense.h>

// data matrix - d x n
// where d is dimension, n is number of points
template<typename T>
using Matrix = A<Dense<T,ColMaj>>;

// dataset abstraction
// template over data type T
template<typename T>
struct DataSet {
	Matrix<T> data;

	DataSet(const Matrix<T> &d ) : data(d) {}

	inline size_t size() const { return data.ncol(); }
	inline size_t dim() const { return data.nrow(); }


	// subset of dataset
	// template over TI, which is assumed to have iterators
	template <typename TI>
	DataSet operator[](const TI &inds) const {
		size_t n = inds.size();
		Matrix<T> data2(dim(), n);
		size_t j = 0;
		for (auto it = inds.cbegin(); it != inds.cend(); ++it) {
			data2[j] = data[*it];
			j++;
		}
		return DataSet(data2);
	}

	inline VectorView<T> operator[](const int i) { return data[i]; }
	inline const VectorView<T> operator[](const int i) const { return data[i]; }
	inline VectorView<T> operator[](const size_t i) { return data[i]; }
	inline const VectorView<T> operator[](const size_t i) const { return data[i]; }

	T& operator()(const size_t i, const size_t j) { return data(i,j); }
	const T& operator()(const size_t i, const size_t j) const { return data(i,j); }


};
