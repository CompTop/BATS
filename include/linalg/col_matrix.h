#pragma once

#include <vector>
#include <cstddef>
#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include "abstract_matrix.h"
#include "csc_matrix.h"
#include "abstract_vector.h"

// standard list of list implementation for sparse matrix

// template over column type
template <class TC>
class ColumnMatrix : public AbstractMatrix
{
private:
    size_t m = 0; // number of rows
    size_t n = 0; // number of columns
    std::vector<TC> col;
public:
    using val_type = typename TC::val_type;
	using col_type = TC;

    // default constructor
    ColumnMatrix() {}

    // construct empty matrix.  same as zero matrix
    ColumnMatrix(size_t m, size_t n) : m(m), n(n) {
        col.resize(n, TC());
    }

    ColumnMatrix(std::vector<TC> &col) : col(col) {
        n = col.size();
    }

    ColumnMatrix(size_t m, size_t n, std::vector<TC> &col) : m(m), n(n), col(col) {}

    // constructor from ColumnMatrix with integer type
    template <typename TC2>
    ColumnMatrix(const ColumnMatrix<TC2> &other) : m(other.nrow()), n(other.ncol()) {
        col.resize(n);
        for (size_t j = 0; j < n; j++) {
            col[j] = TC(other[j]);
        }
    }

    // constructor from CSCMatrix over the integers
    ColumnMatrix(const CSCMatrix<int, size_t> &A) : m(A.nrow()), n(A.ncol()) {
        col.reserve(n);
        auto colptr = A.get_colptr();
        auto rowind = A.get_rowind();
        auto val = A.get_val();
        auto rptr = rowind.cbegin();
        auto vptr = val.cbegin();
        for (size_t j = 0; j < n; j++) {
            // use iterator constructor
            col.emplace_back(TC(
                rptr + colptr[j],
                vptr + colptr[j],
                colptr[j+1] - colptr[j]
            ));
        }
    }

	void read(std::istream &io) {
		std::string line, token;
		getline(io, line); // TODO: should be "Sparse Matrix"
		getline(io, line); // get first line
		// get dimensions from first line
		std::istringstream iss(line);
		getline(iss, token, ',');
		m = std::stoi(token);
		getline(iss, token);
		n = std::stoi(token);
		col.resize(n);
		for (size_t j = 0; j < n; j++) {
			getline(io, line);
			col[j] = TC(line);
		}
	}

    ColumnMatrix(std::istream& io) {
		read(io);
    }

	ColumnMatrix(std::string &fname) {
		std::ifstream file (fname, std::ios::in);
        if (file.is_open()) {
			read(file);
            file.close();
        } else {
            std::cerr << "unable to read ColumnMatrix from " << fname << std::endl;
        }
	}


    inline size_t nrow() const { return m; }
    inline size_t ncol() const { return n; }
    inline std::vector<TC>& cols() { return col; }
    inline const std::vector<TC>& cols() const { return col; }

    auto getval(const size_t i, const size_t j) const {
        return col[j].getval(i);
    }

    // inline size_t width() {
    //     return col.size();
    // }

    inline TC& operator[](size_t index) { return col[index];}
    inline const TC& operator[](size_t index) const { return col[index];}

    inline auto operator()(size_t i, size_t j) const {
        return col[j][i];
    }

    // dump matrix into a dense array in column-major format
    // WARNING: allocates memory, which must be deleted
    val_type* dump_dense() const {
        val_type* arr  =  new val_type[m*n];
        for (size_t j = 0; j < n; j++) {
            for (size_t i = 0; i < m; i++) {
                arr[i + j*m] = col[j][i];
            }
        }
        return arr;
    }

    //
    // TC& constcol const (size_t index) {
    //   return col[index];
    // }

    // permutations: permute, permute_rows, permute_cols

    // permute columns in-place
    // TODO: evaluate if this is best method.
    inline void permute_cols(const std::vector<size_t> &colperm) {
        apply_perm(col, colperm);
    }

    // permute rows in-place
    void permute_rows(const std::vector<size_t> &rowperm) {
        // TODO: this is trivially parallelizable
        for (size_t i = 0; i < col.size(); i++) {
            col[i].permute(rowperm);
        }
    }

    // permute both rows and columns
    void permute(const std::vector<size_t> &rowperm,
        const std::vector<size_t> &colperm) {
        permute_cols(colperm);
        permute_rows(rowperm);
    }


    // addition, substraction, scalar multiplication

    // gemv
    TC gemv(const TC &x) const {
        TC y;  // zero initializer
        // loop over nonzero indices of x
        for (auto xit = x.nzbegin(); xit != x.nzend(); ++xit) {
            y.axpy((*xit).val, col[(*xit).ind]); // y <- x[j]*A[j]
        }
        return y;
    }

    inline TC operator*(const TC &x) const { return gemv(x); }

    // gemm C = self * B
    ColumnMatrix operator*(const ColumnMatrix &B) const {
        std::vector<TC> colC;
        colC.reserve(n);
        //ColumnMatrix C(m, B.n);
        for (size_t j = 0; j < B.n; j++) {
            colC.emplace_back(
                gemv(B.col[j])
            );
        }
        return ColumnMatrix(m, B.n, colC);
    }


    ColumnMatrix operator+(const ColumnMatrix &B) {
        ColumnMatrix C(B); // copy constructor
        for (size_t j = 0; j < n; j++) {
            C.col[j].axpy(1, col[j]);
        }
        return C;
    }

    ColumnMatrix operator-(const ColumnMatrix &B) {
        ColumnMatrix C(B); // copy constructor
        for (size_t j = 0; j < n; j++) {
            C.col[j].axpy(-1, col[j]);
        }
        return C;
    }

	ColumnMatrix transpose() const {
		std::vector<TC> tcol(m);
		// loop over columns
		for (size_t j = 0; j < n; j++) {
			for (auto it = col[j].nzbegin(); it != col[j].nzend(); it++) {
				tcol[it->ind].emplace_back(j, it->val);
			}
		}
		return ColumnMatrix(n, m, tcol);
	}

	inline ColumnMatrix T() const { return transpose(); }

	// apply J matrix from the right
	ColumnMatrix& J_right_inplace() {
		for (size_t j =0; j < n/2; j++) {
			std::swap(col[j], col[n-1-j]);
		}
		return *this;
	}

	inline ColumnMatrix J_right() const {
		ColumnMatrix other(*this);
		other.J_right_inplace();
		return other;
	}

	// apply J matrix from the left
	ColumnMatrix& J_left_inplace() {
		for (size_t j=0; j < n; j++) {
			col[j].J(m);
		}
		return *this;
	}

	inline ColumnMatrix J_left() const {
		ColumnMatrix other(*this);
		other.J_left_inplace();
		return other;
	}

	ColumnMatrix& J_conjugation_inplace() {
		J_left_inplace();
		J_right_inplace();
		return *this;
	}

	inline ColumnMatrix J_conjugation() const {
		ColumnMatrix other(*this);
		other.J_conjugation_inplace();
		return other;
	}

	// swap columns j1 and j2
	inline void swap_cols(size_t j1, size_t j2) {
		std::swap(col[j1], col[j2]);
	}

	// schur complement of i,j entry
	inline void schur_complement(size_t i, size_t j) {
		auto a11 = col[j][i];
		for (size_t j1 = j+1; j1 < n; j1++) {
			auto c = col[j1][i];
			if (c != 0) {
				col[j1].axpy(-c / a11, col[j], i+1, m);
			}
		}
	}

	// apply diagonal matrix with coeff vector along diagonal on the left
	ColumnMatrix& row_scale(const std::vector<val_type> &coeff) {
		for (size_t j = 0; j < n; j++) {
			col[j].scale_inplace(coeff);
		}
		return *this;
	}

    // triangular solve

    // schur complement friend

    // void print() {
    //     std::cout << m << " x " << n << " matrix. transpose: " << std::endl;
    //     for (size_t i = 0; i < col.size(); i++) {
    //         std::cout << i << " : ";
    //         col[i].print_row();
    //     }
    // }
    void print_size() const {
        std::cout << "[" << this << "] : " << m << " x " << n <<\
        " ColumnMatrix" << std::endl;
    }

    void print() const {
        print_size();
        // loop over rows
        for (size_t i = 0; i < m; i++) {
            for (size_t j = 0; j < n; j++) {
                std::cout << std::setw(3) << getval(i, j) << " ";
            }
            std::cout << std::endl;
        }
        return;
    }

    // expose as static member function?
    static ColumnMatrix identity(size_t n) {
        std::vector<TC> col(n);
        for (size_t j = 0; j < n; j++) {
            col[j] = TC(j);
        }
        return ColumnMatrix(n, n, col);
    }
    // dense matrix

    // template <typename IO>
    void write(std::ostream &io) const {
		io << "Sparse Matrix\n";
        io << nrow() << ',' << ncol() << '\n';
        for (size_t j = 0; j < ncol(); j++) {
            // write column as row
            col[j].write(io);
        }
    }

	std::string str() {
	  std::ostringstream oss;
	  write(oss);
	  return oss.str();
	}

    void save(std::string &fname) const {
        std::ofstream file (fname, std::ios::out);
        if (file.is_open()) {
            write(file);
            file.close();
        } else {
            std::cerr << "unable to write ColumnMatrix to " << fname << std::endl;
        }
    }
};



// matvec

// matmul

// return y = A*x
template <class TC>
TC gemv(ColumnMatrix<TC> &A, const TC &x) {
    return A.gemv(x);
}

// return y = U \ x
// solves x = U * y
// Assumes U is upper triangular, with unit diagonal
// pseudo-code
// for i = n:-1:1
//    y[i] = y[i] - U[i,j]x[j]
// this is not generic over TC!
// template <class TC>
// TC ut_solve(ColumnMatrix<TC> &U, const TC &x) {
//     //std::cout << "entering solve" << std::endl;
//     TC y(x);
//     auto yi = y.nzend() - 1;
//     //std::cout << yi->first << " : " << yi->second << std::endl;
//     while (yi >= y.nzbegin() && yi < y.nzend()) {
//         size_t i = yi - y.nzbegin();
//         size_t iind = y.nzind(i);
//         if (iind == 0) {
//             break;
//         }
//         y.axpy(-y.nzval(i), U[iind], 1); //y.axpy(-yp.second, U[i][:i-1])
//         // find next nonzero index
//         yi = y.find_last_nz(iind - 1);
//     }
//     return y;
// }
template <class TC>
TC u_solve(const ColumnMatrix<TC> &U, const TC &y) {
    TC x(y);
    if (x.nnz() == 0) { return x; }
    size_t n = U.ncol();
    //size_t m = U.nrow();
    size_t j = n;
    auto xit = x.upper_bound(j);

    while (xit != x.nzbegin()) {
        // exract j
        --xit;
        j = xit->ind;
        // x[j] = x[j] / U[j,j]
        auto Uj_it = U[j].lower_bound(j); // assume entry j exists and is invertible
        if (Uj_it == U[j].nzend()) {throw std::logic_error("diagonal doesn't exist");}
        auto a = (xit->val) / (Uj_it->val);
        x.replace(xit, a);
        if (j == 0) { break; } // we're done
        // x[0:j-1] -= x[j] * U[0:j-1, j]
        x.axpy(-a, U[j], 0, j);
        // get next nonzero
        xit = x.upper_bound(j-1);
    }
    return x;
}

template <class TC>
TC l_solve(const ColumnMatrix<TC> &L, const TC &y) {
    TC x = y;
    if (x.nnz() == 0) { return x; }
    //size_t n = L.ncol();
    size_t m = L.nrow();
    size_t j = 0;
    auto xit = x.lower_bound(j);
    while (xit != x.nzend()) {
        j = (*xit).ind;
        // x[j] = x[j] / L[j,j]
        auto Lj_it = L[j].lower_bound(j); // assume entry j exists and is invertible
        auto a = (*xit).val / (*Lj_it).val;
        x.replace(xit, a);
        //*xit = nzpair((*xit).ind, a);
        // x[j+1:] -= x[j] * L[j+1:m, j]
        x.axpy(-a, L[j], j+1, m);
        // get next nonzero
        xit = x.lower_bound(j+1);
    }
    return x;
}

// solve U \ A
template <class TC>
ColumnMatrix<TC> u_solve(const ColumnMatrix<TC> &U, const ColumnMatrix<TC> &A) {
    //std::cout << "entering solve" << std::endl;
    size_t m = A.nrow();
    size_t n = A.ncol();
    std::vector<TC> col;
    col.reserve(n);
    for (size_t j = 0; j < n; j++) {

        col.emplace_back(
            u_solve(U, A[j])
        );
    }
    return ColumnMatrix<TC>(m, n, col);
}

// factorization struct
