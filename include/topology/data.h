#pragma once

// for I/O
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include <linalg/naive_dense.h>

// data matrix - d x n
// where d is dimension, n is number of points
template<typename T>
using Matrix = A<Dense<T,RowMaj>>;

// read point cloud from csv file.
// header=True indicates first line is a header
template<typename T=double>
Matrix<T> read_point_cloud(std::string &fname, bool header=false) {
	size_t m = 0;
	size_t n = 0;
	std::vector<T> data;
	std::ifstream file (fname, std::ios::in);
	if (file.is_open()) {
		std::string line;
		if (header) {getline(file, line);} // TODO: check validity of header?
		while (getline(file, line)) {
			m++;
			std::string token;
		    std::istringstream iss(line);
			size_t ncol = 0;
			while (getline(iss, token, ',')) {
				// std::cout << token << ',' << std::endl;
				if (token.size() > 0) {
					ncol++;
					data.emplace_back(T(std::stod(token)));
				}
			}
			if ( n== 0 ) {
				n = ncol;
			} else if (n != ncol) {
				std::cerr << "inconsistent number of columns at line " << m << std::endl;
				std::cerr << "saw " << ncol << " entries, expected " << n << std::endl;
			}
		}
		file.close();
	} else {
		std::cerr << "unable to open " << fname << std::endl;
	}
	T* newmat = new T[m*n];
	T* it = newmat;
	for (auto x : data) {
		*it++ = x;
	}
	return Matrix<T>(n, m, newmat); // swap rows and columns for col major format
}


// dataset abstraction
// template over data type T
template<typename T>
struct DataSet {
	Matrix<T> data;

	DataSet() {}
	DataSet(const Matrix<T> &d ) : data(d) {}

	inline size_t size() const { return data.ncol(); }
	inline size_t dim() const { return data.nrow(); }
	inline Matrix<T>& get_data() { return data; }


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
