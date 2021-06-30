#pragma once

#include "dense.hpp"
#include "lu.hpp"
#include "util.hpp"

namespace bats {
namespace future {


/*
	Class which represents span
	TODO: make this behave somewhat like std::set
	http://www.cplusplus.com/reference/set/set/
	insert
	erase
	clear
	intersection  // how to do this? form projection matrix.  Not over finite fields
	union
*/
template <typename T>
class Span {
	typedef T value_type;
public:
	// use column-major matrix to store L
	Matrix<T> Pt; // TODO: swap with perumtation type
	Matrix<T> L;
	size_t _vdim;


public:

	Span() {}
	Span(T) {}
	Span(size_t n, T) {
		_vdim = n;
		Pt = Matrix<T>::identity(n);
		L = Matrix<T>(_vdim, 0, ColumnMajor());
	}

	// return dimension of span
	size_t dim() const { return L.ncol(); }

	// return dimension of vector space
	size_t vdim() const { return _vdim; }

	// add vector to span
	// returns false if vector is in span
	// returns true if vector is not in span
	// TV can be any random access iterator
	template <typename TV>
	bool add(const TV &v) {
		auto n = L.ncol();
		L.append_column();
		// put P^T*v in last column of L
		gemv(Pt, v, L.column(n));
		// apply LU
		l_residual(L.view(range<size_t>(0,L.nrow()), range<size_t>(0, n)), L.column(n));
		// permute first non-zero to front
		auto piv = find_pivot_high(L.column(n), n, L.nrow());
		if (piv == -1) {
			// vector already in span.
			// delete added column
			L.delete_column();

			return false;
		} else if (piv > n) {
			L.swap_rows(piv, n);
			Pt.swap_rows(piv, n);
		}

		// scale column to have unit diagonal
		L.scale_column(L(n,n).inv(), n);

		return true;
	}

	// check if vector is in span
	// returns true if vector is in span
	// returns false if vector is not in span
	template <typename TV>
	bool contains(const TV &v) const {
		std::vector<T> v_copy(dim());

		// put P^T*v in v_copy
		gemv(Pt, v, v_copy);

		l_residual(L, v_copy);

		auto piv = find_pivot_high(v, vdim(), dim());

		// if v_copy has been zeroed out, then in span.
		return piv == -1 ? true : false;
	}




};

}
}
