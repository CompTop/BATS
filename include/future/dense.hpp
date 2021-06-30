#pragma once
/*
Dense row matrix implementation
*/
#include <string>
#include <iostream>
#include <utility>
#include <iterator>
#include <algorithm>
#include <vector>

namespace bats {
namespace future {

// forward declarations
template <typename T>
class Matrix;

template <typename T, typename I1, typename I2>
class MatrixView;

template <typename M1, typename M2>
auto gemm(const M1& A, const M2& B);

template <typename MT, typename VT>
auto gemv(const MT& A, const VT& x);


/*
A strided iterator type
*/
template <typename T>
class strided_iterator
	: public std::iterator<std::random_access_iterator_tag, T>
{
	typedef ssize_t difference_type;
	typedef T value_type;
	typedef T& reference;
	typedef T* pointer;
	typedef std::random_access_iterator_tag iterator_category;

	typedef strided_iterator<T> iterator;
	pointer pos_;
	ssize_t stride_;
public:
	strided_iterator() : pos_(nullptr), stride_(0) {}
	strided_iterator(T* v, ssize_t s) : pos_(v), stride_(s) {}
	~strided_iterator() {}

	iterator operator++(int) { return iterator(pos_+=stride_, stride_); } // postfix
	iterator& operator++() {pos_+=stride_; return *this; } // prefix
	iterator operator--(int) {return iterator(pos_-=stride_, stride_); } // postfix
	iterator& operator--() {pos_-=stride_; return *this; } // prefix
	reference operator*() const { return *pos_; }
	pointer operator->() const { return pos_; }
	bool operator==(const iterator& rhs) const { return pos_ == rhs.pos_; }
	bool operator!=(const iterator& rhs) const {return pos_ != rhs.pos_; }

	iterator operator+(size_t i) const {return iterator(pos_ + stride_*i, stride_);}
	reference operator[](size_t i) const {return *(pos_ + stride_*i); }

	// TODO: more operators
};
template <typename T>
class const_strided_iterator
	: public std::iterator<std::random_access_iterator_tag, T>
{
	typedef ssize_t difference_type;
	typedef T value_type;
	typedef const T& reference;
	typedef const T* pointer;
	typedef std::random_access_iterator_tag iterator_category;

	typedef const_strided_iterator<T> iterator;
	pointer pos_;
	ssize_t stride_;
public:
	const_strided_iterator() : pos_(nullptr), stride_(0) {}
	const_strided_iterator(const T* v, ssize_t s) : pos_(v), stride_(s) {}
	~const_strided_iterator() {}

	iterator operator++(int) { return iterator(pos_+=stride_, stride_); } // postfix
	iterator& operator++() {pos_+=stride_; return *this; } // prefix
	iterator operator--(int) {return iterator(pos_-=stride_, stride_); } // postfix
	iterator& operator--() {pos_-=stride_; return *this; } // prefix
	reference operator*() const { return *pos_; }
	pointer operator->() const { return pos_; }
	bool operator==(const iterator& rhs) const { return pos_ == rhs.pos_; }
	bool operator!=(const iterator& rhs) const {return pos_ != rhs.pos_; }

	iterator operator+(size_t i) const {return iterator(pos_ + stride_*i, stride_);}
	reference operator[](size_t i) const {return *(pos_ + stride_*i); }

	// TODO: more operators
};

template <typename T>
class range {
private:
	T begin_;
	T end_;
	T stride_;

public:
	range() : begin_(0), end_(0), stride_(0) {}
	range(T b, T e) : begin_(b), end_(e), stride_(1) {
		if (e < b) {
			begin_ = e;
			end_ = b;
			stride_ = T(-1);
		}
	}
	range(T b, T e, T s) : begin_(b), end_(e), stride_(s) {}

	T first() const { return begin_; }
	T last() const { return end_; }
	T stride() const { return stride_; }

	size_t size() const { return (end_ - begin_) / stride_; }

	T operator[](size_t i) const { return begin_ + i * stride_; }

	// TODO: iterator types for class
	class const_iterator
		: public std::iterator<std::random_access_iterator_tag, T>
	{
		typedef ssize_t difference_type;
		typedef T value_type;
		typedef const T& reference;
		typedef const T* pointer;
		typedef std::random_access_iterator_tag iterator_category;

		typedef const_iterator iterator;
		T val_;
		T stride_;

		// TODO: should implement this using a count to avoid floating point issues
	public:
		const_iterator() : val_(0), stride_(0) {}
		const_iterator(const T v, const T s) : val_(v), stride_(s) {}
		~const_iterator() {}

		iterator operator++(int) { return iterator(val_+=stride_, stride_); } // postfix
		iterator& operator++() {val_+=stride_; return *this; } // prefix
		iterator operator--(int) {return iterator(val_-=stride_, stride_); } // postfix
		iterator& operator--() {val_-=stride_; return *this; } // prefix
		reference operator*() const { return val_; }
		pointer operator->() const { return &val_; }
		bool operator==(const iterator& rhs) const { return val_ == rhs.val_; }
		bool operator!=(const iterator& rhs) const {return val_ != rhs.val_; }

		// TODO: more operators
	};
	typedef const_iterator iterator; // only use const iterator type

	const_iterator const begin() { return iterator(begin_, stride_); }
	const_iterator const end() { return iterator(end_, stride_); }

	const_iterator const cbegin() { return iterator(begin_, stride_); }
	const_iterator const cend() { return iterator(end_, stride_); }

};

// view type - doesn't own data
template <typename T>
class VectorView {
private:
	T* start_;
	T* end_;
	size_t stride_;

public:

	typedef strided_iterator<T> iterator;
	typedef const_strided_iterator<T> const_iterator;

	VectorView(T* s, T* e, size_t stride) : \
		start_(s), end_(e), stride_(stride) {}

	iterator begin() { return iterator(start_, stride_); }
	iterator end()   { return iterator(end_  , stride_); }

	const_iterator begin() const { return const_iterator(start_, stride_); }
	const_iterator end()   const { return const_iterator(end_, stride_); }

	const_iterator cbegin() const { return const_iterator(start_, stride_); }
	const_iterator cend()   const { return const_iterator(end_, stride_); }

	T& operator[](size_t i) {return begin()[i];}
	const T& operator[](size_t i) const {return cbegin()[i];}


	void print() const {
		T* it = start_;
		while (it != end_) {
			std::cout << *it << ' ';
			it += stride_;
		}
		std::cout << std::endl;
	}

	void axpy(T a, const VectorView &other) {
		auto it1 = begin();
	 	auto it2 = other.begin();
		while (it1 != end()) {
			*it1++ += a * (*it2++);
		}
	}
	void axpy(T a, const VectorView &&other) {
		auto it1 = begin();
	 	auto it2 = other.begin();
		while (it1 != end()) {
			*it1++ += a * (*it2++);
		}
	}

	void scale(T a) {
		auto it = begin();
		while (it != end()) {
			*it++ *= a;
		}
	}

};



struct RowMajor {};
struct ColumnMajor {};

// template over type
template <typename T>
class Matrix {
private:
	std::vector<T> data_;
	size_t dim[2]; // dimensions
	size_t stride[2]; // strides

public:

	typedef T value_type;

	size_t len() const { return dim[0]*dim[1]; }

	void fill(T a) {
		std::fill(data_.begin(), data_.end(), a);
	}

	// empty constructor
	Matrix() : data_(), dim{0,0}, stride{0,0} {}

	// construct empty row major matrix
	Matrix(size_t m, size_t n, RowMajor) : dim{m,n}, stride{n,1} {
		data_ = std::vector<T>(len());
	}

	Matrix(size_t m, size_t n, ColumnMajor) : dim{m,n}, stride{1,m} {
		data_ = std::vector<T>(len());
	}

	// construct an empty matrix
	Matrix(size_t m, size_t n) : dim{m, n}, stride{n, 1} {
		data_ = std::vector<T>(len());
	}

	// construct a matrix with a fill value a
	Matrix(size_t m, size_t n, T a) : Matrix(m, n) {
		fill(a);
	}
	template <typename Major>
	Matrix(size_t m, size_t n, T a, Major &ms): Matrix(m, n, ms) {
		fill(a);
	}

	Matrix(Matrix& other, ColumnMajor) : Matrix(other.dim[0], other.dim[1], ColumnMajor()) {
		for (size_t j = 0; j < dim[1]; j++) {
			for (size_t i = 0; i < dim[0]; i++) {
				operator()(i,j) = other(i,j);
			}
		}
	}

	T* data() { return data_.data(); }
	const T* data() const {return data_.data();}

	inline size_t nrow() const {return dim[0];}
	inline size_t ncol() const {return dim[1];}
	inline std::pair<size_t, size_t> dims() const { return std::make_pair(dim[0], dim[1]); }

	// access data
	inline T operator[](size_t k) const { return data_[k]; }
	inline T& operator[](size_t k) { return data_[k]; }

	inline T operator()(size_t k) const { return data_[k]; }
	inline T& operator()(size_t k) { return data_[k]; }

	template <typename I1, typename I2>
	inline MatrixView<T, I1, I2> view(I1 &rows, I2 &cols) {return MatrixView(this, rows, cols);}
	template <typename I1, typename I2>
	inline MatrixView<T, I1, I2> view(I1 &&rows, I2 &&cols) {return MatrixView(this, rows, cols);}
	// template <typename I1, typename I2>
	// inline MatrixView<T, I1, I2> operator()(I1 &rows, I2 &cols) {return view(rows, cols);}

	template <typename I1, typename I2>
	inline const MatrixView<T, I1, I2> view(I1 &rows, I2 &cols) const {return MatrixView(this, rows, cols);}
	template <typename I1, typename I2>
	inline const MatrixView<T, I1, I2> view(I1 &&rows, I2 &&cols) const {return MatrixView(this, rows, cols);}
	// template <typename I1, typename I2>
	// inline const MatrixView<T, I1, I2> operator()(I1 &rows, I2 &cols) const {return view(rows, cols);}
	// template <typename I1, typename I2>
	// inline const MatrixView<T, I1, I2> operator()(I1 &&rows, I2 &&cols) const {return view(rows, cols);}

	inline T operator()(size_t i, size_t j) const { return data_[stride[0]*i + stride[1]*j]; }
	inline T& operator()(size_t i, size_t j) { return data_[stride[0]*i + stride[1]*j]; }

	inline VectorView<T> row(size_t i) {return VectorView(data() + stride[0]*i, data() + stride[0]*i + stride[1]*dim[1], stride[1]);}
	inline VectorView<T> column(size_t j) {return VectorView(data() + stride[1]*j, data() + stride[1]*j + stride[0]*dim[0], stride[0]);}


	void print_info() const {
		std::cout << "[" << this << "] : " << nrow() << " x " << ncol() <<\
        " Matrix" << std::endl;
	}

	void print() const {
		print_info();
		for (size_t i = 0; i < nrow(); i++) {
			for (size_t j = 0; j < ncol(); j++) {
				std::cout << operator()(i, j) << ' ';
			}
			std::cout << std::endl;
		}
	}

	Matrix transpose() const {
		Matrix At(ncol(), nrow());
		for (size_t i = 0; i < nrow(); i++) {
			for (size_t j = 0; j < ncol(); j++) {
				At(j,i) = operator()(i,j);
			}
		}
		return At;
	}

	// elementary operations
	// swap rows of the matrix
	void swap_rows(size_t i1, size_t i2) {
		if (i1 == i2) {return;}
		for (size_t j = 0; j < ncol(); j++) {
			std::swap(operator()(i1,j), operator()(i2, j));
		}
	}

	void swap_columns(size_t j1, size_t j2) {
		if (j1 == j2) {return;}
		for (size_t i = 0; i < nrow(); i++) {
			std::swap(operator()(i,j1), operator()(i, j2));
		}
	}

	// append column - return true if successful.
	bool append_column(const T val=T(0)) {
		if (stride[1] == dim[0]) {
			// column major - just need to resize
			dim[1]++;
			data_.resize(len(), val);
			return true;
		} // TODO: implement row-major
		return false;
	}
	bool delete_column() {
		if (stride[1] == dim[0]) {
			// column major - just need to resize
			dim[1]--;
			data_.resize(len());
			return true;
		} // TODO: implement row-major
		return false;
	}

	// add a*row(i1) to row(i0)
	void add_row(T a, size_t i1, size_t i0) {
		row(i0).axpy(a, row(i1));
	}
	inline void add_row(size_t i1, size_t i0) { return add_row(T(1), i1, i0); }

	// add a * column(i1) to column(i0)
	void add_column(T a, size_t j1, size_t j0) {
		column(j0).axpy(a, column(j1));
	}
	inline void add_column(size_t j1, size_t j0) {return add_column(T(1), j1, j0); }

	void scale_row(T a, size_t i) {
		row(i).scale(a);
	}

	void scale_column(T a, size_t j) {
		column(j).scale(a);
	}

	// matrix-matrix multiplication
	template <typename TB>
	inline auto mm(const TB &B) const {
		return gemm(*this, B);
	}
	template <typename TB>
	inline Matrix operator*(const TB &B) const { return mm(B); }

	// matrix-vector multiplication
	template <typename Tx>
	inline auto mv(const Tx &x) const {
		return gemv(*this, x);
	}


	static Matrix identity(size_t n) {
		Matrix A(n, n, T(0));
		for (size_t k = 0; k < n; k++) {
			A(k,k) = T(1);
		}
		return A;
	}
};

// matrix-matrix multiplication
template <typename M1, typename M2>
auto gemm(const M1& A, const M2& B) {

	using T = typename M1::value_type;
	// TODO: static assert on value types

	// C = A * B
	size_t Cm = A.nrow();
	size_t Cn = B.ncol();
	Matrix C(Cm, Cn, T(0));
	for (size_t i = 0; i < Cm; i++) {
		for (size_t j = 0; j < Cn; j++) {
			for (size_t k = 0; k < A.ncol(); k++) {
				C(i,j) += A(i,k) * B(k,j);
			}
		}
	}
	return C;
}

// matrix-vector multiplication
// y <- A * x
template <typename MT, typename VT1, typename VT2>
VT2& gemv(const MT& A, const VT1& x, VT2&& y) {
	// clear y
	for (size_t i = 0; i < A.nrow(); i++) {
		y[i] = 0;
	}

	// perform matvec
	for (size_t j = 0; j < A.ncol(); j++) {
		for (size_t i = 0; i < A.nrow(); i++) {
			y[i] += A(i,j) * x[j];
		}
	}
	return y;
}

// matrix-vector multiplication
template <typename MT, typename VT>
auto gemv(const MT& A, const VT& x) {

	using T = typename MT::value_type;
	// TODO: static assert on value types

	// C = A * B
	size_t Cm = A.nrow();
	Matrix y(Cm, 1, T(0));
	return gemv(A, x, y);

}



// another view type
// template over types used to iterate over rows and columns.
template <typename T, typename I1, typename I2>
class MatrixView {
private:

	Matrix<T>* data;
	I1 rows;
	I2 cols;

public:

	typedef T val_type;


	MatrixView() : data(nullptr) {}

	// need to make copies of r and c anyways
	MatrixView(Matrix<T>* m, I1 r, I2 c) : data(m) {
		rows = std::move(r);
		cols = std::move(c);
	}
	MatrixView(Matrix<T>& m, I1 r, I2 c) : data(&m) {
		rows = std::move(r);
		cols = std::move(c);
	}

	inline T operator()(size_t i, size_t j) const { return data->operator()(rows[i], cols[j]); }
	inline T& operator()(size_t i, size_t j) { return data->operator()(rows[i], cols[j]); }

	inline size_t nrow() const { return rows.size(); }
	inline size_t ncol() const { return cols.size(); }

	void print_info() const {
		std::cout << "[" << this << "] : " << nrow() << " x " << ncol() <<\
        " MatrixView of [" << data << "]"  << std::endl;
	}

	void print() const {
		print_info();
		for (size_t i = 0; i < nrow(); i++) {
			for (size_t j = 0; j < ncol(); j++) {
				std::cout << operator()(i, j) << ' ';
			}
			std::cout << std::endl;
		}
	}
};

}
}
