#pragma once
/*
Header file for univariate polynomial class with coefficients in T
*/
#include <vector>
#include <initializer_list>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include "col_matrix.hpp"
#include "sparse_vector.hpp"

template <typename T>
class UnivariatePolynomial
{
private:
    std::vector<T> c; // vector of coefficients

    // remove zero coeffients from back
    void clear_zeros() {
        while (!c.empty() && c.back() == 0) {
            c.pop_back();
        }
    }

public:

    // constructor from vector of coefficients
    UnivariatePolynomial(const std::vector<T>& c) : c(c) {clear_zeros();}

    // constructor from initializer list
    UnivariatePolynomial(std::initializer_list<T> l) : c(l) {clear_zeros();}

    // constructor of constant polynomial
    UnivariatePolynomial(const T& c0) : c{c0} {clear_zeros();}

	// empty constructor
	UnivariatePolynomial() {}

	// multipliciative identity
    static UnivariatePolynomial identity() {
		return UnivariatePolynomial(T(1));
	}

	// additive identity
    static UnivariatePolynomial zero() {
		return UnivariatePolynomial();
	}

	bool is_zero() const {
		return c.empty();
	}

	// monomial with degree d, scaled by scale
	static UnivariatePolynomial monomial(size_t d, T scale=T(1)) {
		std::vector<T> coeff(d + 1, T(0));
		coeff[d] = scale;
		return UnivariatePolynomial(coeff);
	}

    // access coefficients
    template <typename TI>
    inline const T operator[] (const TI i) const { return c[i]; }

    // update coefficient
    template <typename TI>
    inline T& operator[] (const TI i) { return c[i]; }


    // dimension of polynomial
    inline size_t dim() const { return c.size() - 1; }
	inline size_t degree() const { return dim(); }
    // size - useful for iteration
    inline size_t size() const { return c.size(); }

    // leading coefficient
    inline T leading_coeff() const {
        return c.back();
    }

	// apply polynomial
	T operator()(const T& v) const {
		if (is_zero()) {return T(0);}
		T sum = c[0];
		T vk = T(1);
		for (size_t k = 1; k < size(); k++) {
			vk = vk * v; // v^k
			sum += c[k] * vk;
		}
		return sum;
	}

	// apply polynomial of matrix to vector
	// p(A) * v
	template <typename TC>
	TC operator()(const ColumnMatrix<TC>& A, const TC& v) const {
		if (is_zero()){ return TC(); }
		TC sum; sum.axpy(c[0], v); // c[0]*v
		TC Akv(v); // A^k * v
		for (size_t k = 1; k < size(); k++) {
			Akv = A * Akv; // A^k * v
			sum.axpy(c[k], Akv);
		}
		return sum;
	}

    // negation
    UnivariatePolynomial operator-() const {
        std::vector<T> newc(c.size());
        for (size_t i = 0; i < c.size(); i++) {
            newc[i] = -c[i];
        }
        return UnivariatePolynomial(newc);
    }


    // addition
    UnivariatePolynomial operator+(const UnivariatePolynomial &other) const {
        std::vector<T> newc;
        newc.reserve(size());

        auto it1 = c.begin();
        auto it2 = other.c.begin();
        while (it1 != c.end() && it2 != other.c.end()) {
            newc.emplace_back(*it1++ + *it2++);
        }
        while (it1 != c.end()) {
            newc.emplace_back(*it1++);
        }
        while (it2 != other.c.end()) {
            newc.emplace_back(*it2++);
        }

        return UnivariatePolynomial(newc);
    }

	// addition
    UnivariatePolynomial& operator+=(const UnivariatePolynomial &other) {

        auto it1 = c.begin();
        auto it2 = other.c.begin();
        while (it1 != c.end() && it2 != other.c.end()) {
			*it1++ += *it2++;
        }
        while (it2 != other.c.end()) {
            c.emplace_back(*it2++);
        }

		clear_zeros();

        return *this;
    }

	// addition
    UnivariatePolynomial& operator-=(const UnivariatePolynomial &other) {

        auto it1 = c.begin();
        auto it2 = other.c.begin();
        while (it1 != c.end() && it2 != other.c.end()) {
			*it1++ -= *it2++;
        }
        while (it2 != other.c.end()) {
            c.emplace_back(-(*it2++));
        }

		clear_zeros();

        return *this;
    }

    // subtraction
    UnivariatePolynomial operator-(const UnivariatePolynomial &other) const {
        std::vector<T> newc;

        newc.reserve(size());

        auto it1 = c.begin();
        auto it2 = other.c.begin();
        while (it1 != c.end() && it2 != other.c.end()) {
            newc.emplace_back(*it1++ - *it2++);
        }
        while (it1 != c.end()) {
            newc.emplace_back(*it1++);
        }
        while (it2 != other.c.end()) {
            newc.emplace_back(-(*it2++));
        }


        return UnivariatePolynomial(newc);
    }

    // multiplication
    UnivariatePolynomial operator*(const UnivariatePolynomial &other) const {
		if (is_zero() || other.is_zero()) { return zero(); }
        std::vector<T> newc(size() + other.size() - 1, T(0));


        for (size_t i = 0; i < size(); i++) {
            for (size_t j = 0; j < other.size(); j++) {
                newc[i+j] += c[i] * other[j];
            }
        }

        return UnivariatePolynomial(newc);
    }

    // division + remainder
    std::tuple<UnivariatePolynomial, UnivariatePolynomial> divrem(const UnivariatePolynomial& other) const {
        UnivariatePolynomial quotient(0);
        UnivariatePolynomial remainder(*this);

        while (remainder.size() >= other.size()) {
            T coeff = remainder.leading_coeff() / other.leading_coeff();
			UnivariatePolynomial m = monomial(remainder.dim() - other.dim(), coeff);
			quotient += m;
			remainder -= (m * other);
        }

        return std::tie(quotient, remainder);
    }

    // division
    UnivariatePolynomial operator/(const UnivariatePolynomial& other) const {
        auto [q, r] = divrem(other);
        return q;
    }

    // remainder
    UnivariatePolynomial remainder(const UnivariatePolynomial& other) const {
        auto [q, r] = divrem(other);
        return r;
    }

	template <typename T2>
	bool operator==(const T2& other) const {
		return other == 0 ? is_zero() : (size() == 1 && c[0] == other);
	}

	bool operator==(const UnivariatePolynomial &other) const {
		if (size() != other.size()) return false;
		for (size_t i = 0; i < size(); i++) {
			if (c[i] != other[i]) return false;
		}
		return true;
	}

	template <typename T2>
	bool operator!=(const T2& other) const {
		return other == 0 ? !is_zero() : !(size() == 1 && c[0] == other);
	}

	inline bool operator!=(const UnivariatePolynomial &other) const {
		return !(*this == other);
	}


    // gcd using Euclidean algorithm
	UnivariatePolynomial gcd(const UnivariatePolynomial& other) const {
		UnivariatePolynomial a(*this), b(other), t;
		while (!b.is_zero()) {
			t = b;
			b = a.remainder(b);
			a = t;
		}
		return a;
	}

	inline bool is_monic() const {
		return leading_coeff() == T(1);
	}

    // companion matrix
	ColumnMatrix<SparseVector<T>> companion_matrix() const {
		if (!is_monic()) {throw std::runtime_error("Companion matrix requires monic polynomial!");}
		size_t n = dim(); // size of companion matrix
		std::vector<SparseVector<T>> cols(n);
		for (size_t j = 0; j < n-1; j++) {
			cols[j] = SparseVector<T>({j+1}, {T(1)});
		}
		std::vector<size_t> ind(n);
		std::vector<T> val(n);
		for (size_t i = 0; i < n; i++) {
			ind[i] = i;
			val[i] = -c[i];
		}
		cols[n-1] = SparseVector<T>(ind, val);

		return ColumnMatrix(n,n,cols);
	}


    // call as function to evaluate
    // input should be templated.  e.g. should be able to evaluate on matrix or field elt.

    // pretty printing
    friend std::ostream& operator<<( std::ostream& os, const UnivariatePolynomial &p) {
        for (size_t i = 0; i < p.c.size(); i++) {
            if (p.c[i] != 0) {
                os << p.c[i] << " x^" << i;
                if (i != p.c.size() - 1) os << " + ";
            }
        }
        return os;
    }

	void print() {
		std::cout << *this << std::endl;
	}


}; // end UnivariatePolynomial

// type trait
template <typename T>
struct is_UnivariatePolynomial : std::false_type {};

template <typename T>
struct is_UnivariatePolynomial<UnivariatePolynomial<T>> : std::true_type {};


// characteristic matrix
// xI - A
template <typename T>
auto characteristic_matrix(
	const ColumnMatrix<SparseVector<T>>& A
) {
	using PT = UnivariatePolynomial<T>;

	auto n = A.ncol();
	auto m = A.nrow();
	if (m != n) {throw std::runtime_error("Characteristic matrix must be square!");}

	std::vector<SparseVector<PT>> cols(n);
	std::vector<size_t> ind;
	std::vector<PT> val;
	for (size_t j = 0; j < n; j++) {
		ind.clear();
		val.clear();
		auto p = A[j].nzbegin();
		bool found_diagonal = false;
		while (p != A[j].nzend()) {
			if (p->ind == j) {
				ind.emplace_back(p->ind);
				val.emplace_back(PT({-p->val, T(1)}));
				found_diagonal=true;
				p++;
				break;
			} else if (p->ind < j) {
				ind.emplace_back(p->ind);
				val.emplace_back(PT(-p->val));
			} else {
				break;
			}
			p++;
		}
		if (!found_diagonal) {
			ind.emplace_back(j);
			val.emplace_back(PT({T(0), T(1)}));
		}
		// new loop without branch conditions
		// only loops over indices > j
		while (p != A[j].nzend()) {
			ind.emplace_back(p->ind);
			val.emplace_back(PT(-p->val));
			p++;
		}

		cols[j] = SparseVector<PT>(ind, val);
	}

	return ColumnMatrix<SparseVector<PT>>(n, n, cols);
}
