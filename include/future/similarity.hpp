#pragma once

namespace bats {
namespace future {

// TODO: template over matrix type.
template <typename T>
struct SimilarityTransform {
	/*
	struct to hold similarity transformation of matrix A
	Maintains invariant A0 = S * A * Sinv
	*/
	Matrix<T> S;
	Matrix<T> A;
	Matrix<T> Sinv;


	// TODO: check that A is square
	SimilarityTransform(const Matrix<T> &A0) : \
	S(Matrix<T>::identity(A0.nrow())),
	A(A0),
	Sinv(Matrix<T>::identity(A0.nrow())) {}

	size_t size() const {return A.nrow();}

	// const access data
	inline T operator()(size_t i, size_t j) const {return A(i,j);}

	void print_info() const {
		std::cout << "[" << this << "] : " << \
        " SimilarityTransform" << std::endl;
	}
	void print() const {
		print_info();
		S.print();
		A.print();
		Sinv.print();
	}

	Matrix<T> prod() const { return ((S * A) * Sinv); }

	// swap rows/columns
	void swap_rows(size_t i0, size_t i1) {
		// have to do everything symmetric
		A.swap_rows(i0, i1);
		A.swap_columns(i0, i1);
		S.swap_columns(i0, i1);
		Sinv.swap_rows(i0, i1);
	}
	inline void swap_columns(size_t j0, size_t j1) { return swap_rows(j0, j1); }

	// add a * row(i1) to row(i0)
	void add_row(T a, size_t i1, size_t i0) {
		A.add_row(a, i1, i0);
		A.add_column(-a, i0, i1);
		S.add_column(-a, i0, i1);
		Sinv.add_row(a, i1, i0);
	}

	void scale_row(T a, size_t i) {
		T ai = a.inv();
		S.column(i).scale(ai);
		A.row(i).scale(a);
		A.column(i).scale(ai);
		Sinv.row(i).scale(a);
	}

	void scale_column(T a, size_t i) {
		T ai = a.inv();
		S.column(i).scale(a);
		A.row(i).scale(ai);
		A.column(i).scale(a);
		Sinv.row(i).scale(ai);
	}

};

}
}
