#pragma once

#include <vector>
#include <map>
#include <cstddef>
#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <sstream>
#include <random>
#include <chrono>
#include "abstract_matrix.hpp"
#include "csc_matrix.hpp"
#include "abstract_vector.hpp"

// standard list of list implementation for sparse matrix

// forward declarations
template <class TC>
class ColumnMatrix;

template <class TC>
TC l_solve(const ColumnMatrix<TC> &L, const TC &y);

template <class TC>
TC u_solve(const ColumnMatrix<TC> &U, const TC &y);

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
	using tmp_type = typename TC::tmp_type; // for use with axpy

    // default constructor
    ColumnMatrix() {}

    // construct empty matrix.  same as zero matrix
    ColumnMatrix(size_t m, size_t n) : m(m), n(n) {
        col.resize(n, TC());
    }

	/**
	construct column matrix filled with entry a
	*/
	ColumnMatrix(size_t _m, size_t _n, val_type a) : m(_m), n(_n) {
		col.resize(n, TC::dense(a, m));
	}

    ColumnMatrix(const std::vector<TC> &_col) : col(_col) {
        n = col.size();
    }

    ColumnMatrix(size_t _m, size_t _n, const std::vector<TC>& _col) : m(_m), n(_n), col(_col) {}

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
	inline void set_nrow(size_t mnew) { m = mnew; }

    auto getval(const size_t i, const size_t j) const {
        return col[j].getval(i);
    }

	// permute rows in-place
    void permute_rows(const std::vector<size_t> &rowperm) {
        // // TODO: this is trivially parallelizable
        for (size_t i = 0; i < col.size(); i++) {
            col[i].permute(rowperm);
        }
    }

	// permute rows in-place
    void ipermute_rows(const std::vector<size_t> &rowperm) {
        // // TODO: this is trivially parallelizable
        for (size_t i = 0; i < col.size(); i++) {
            col[i].ipermute(rowperm);
        }
    }

    // inline size_t width() {
    //     return col.size();
    // }

	/**
	append a column to the end of the matrix

	calls a column vector constructor on arguments passed in.
	*/
	template <class ...Ts>
	void append_column(Ts (&...args)) { col.emplace_back(TC(args...)); n++; }
	template <class ...Ts>
	void append_column(Ts (&&...args)) { col.emplace_back(TC(args...)); n++; }

	/**
	append an empty row to the end of the matrix
	*/
	inline void append_row() { m++; }

	/**
	Do not recommend,
	because if we want to insert multiple columns,
	the index will change after the first insertion.
	Thus the important thing we need to assume is
	the index list are in ascending order!
	*/
	template <class ...Ts>
	void insert_column(const size_t& index, Ts (&&...args)){
		// first append a row to the end
		append_column(args...);

		// second permute it
		std::vector<size_t> colperm;
		colperm.reserve(n);
		for (size_t i = 0; i < n ; i++){
			size_t j = i ;
			if (i == index) j = n-1;
			if (i > index) j = i-1;
			colperm.emplace_back(j);
		}
		ipermute_cols(colperm);
	}
	// append a zero column at the end
	void append_column(){ col.emplace_back(TC()); n++; }

	// insert a zero column at position index
	void insert_column(const size_t& index){
		// first append a row to the end
		append_column();
		// second permute it
		std::vector<size_t> colperm;
		colperm.reserve(n);
		for (size_t i = 0; i < n ; i++){
			size_t j = i ;
			if (i == index) j = n-1;
			if (i > index) j = i-1;
			colperm.emplace_back(j);
		}
		permute_cols(colperm);
	}

	/**
	insert list of columns at list of specified indices

	mutates input inserted columns
	*/
	void insert_columns(
		const std::vector<size_t>& c_inds,
		std::vector<TC>& insert_col
	) {
		n = n + c_inds.size();
		std::vector<TC> newcol(n);

		size_t i = 0;
		size_t newi = 0;
		size_t oldi = 0;
		while (oldi < col.size() && newi < insert_col.size()) {
			if (i == c_inds[newi]) {
				// insert next new column
				std::swap(newcol[i], insert_col[newi]);
				++i;
				++newi;
			} else {
				std::swap(newcol[i], col[oldi]);
				++i;
				++oldi;
			}
		}
		// at most one of the next two while-loops will execute
		// this one just puts the rest of the old columns in
		while (oldi < col.size()) {
			std::swap(newcol[i], col[oldi]);
			++i;
			++oldi;
		}
		// else, we just add rest of new columns to end
		while (newi < insert_col.size()) {
			std::swap(newcol[i], insert_col[newi]);
			++i;
			++newi;
		}

		// finally, swap these
		std::swap(col, newcol);

	}

	// append a row a specified value
	void append_row(const std::vector<val_type> & row) {
		assert (row.size() == ncol());
		++m;
		for (size_t j = 0; j < ncol(); j++) {
			if (row[j]!= val_type(0)) {
				col[j].emplace_back(m-1, row[j]);
			}
		}
	}

	// append a row a specified value
	void append_sparse_row(const col_type & new_row) {
		++m;
		for (auto it = new_row.nzbegin(); it != new_row.nzend(); it++){
			col[it->ind].emplace_back(m-1, it->val);
		}
	}

	// TODO: permute all rows at once
	void insert_row(const size_t& ind, const std::vector<val_type> & row){
		// first append a row to the end
		append_row(row);
		// second permute it
		std::vector<size_t> rowperm;
		rowperm.reserve(m);
		for (size_t i = 0; i < m ; i++){
			size_t j = i; // j is the new index, rows above ind keep the same
			if (i >= ind) j = i+1; // rows above ind will add up one
			if (i == m-1) j = ind; // the last one will be permuted the desired ind
			rowperm.emplace_back(j);
		}
		permute_rows(rowperm);
	}

	/**
	 * Insert sparse rows by permutation
	 * Input:
	 * 	r_inds: new row indices
	 *  new_rows: a vector of sparse row vectors
	*/
	void insert_sparse_rows(const std::vector<size_t>& r_inds,
	const std::vector<col_type> & new_rows){
		// first append rows to the end
		for (auto row: new_rows){
			append_sparse_row(row);
		}
		// second permute it
		auto length_new_rows = r_inds.size();
		std::vector<size_t> rowperm;
		rowperm.reserve(m);
		for (size_t i = 0; i < m-length_new_rows; i++){
			size_t j = i; // j is the new index, rows above ind keep the same
			// for(auto ele: r_inds){
			// 	if(j >= ele) j+=1; // rows above ele will add up one
			// }
			// or faster option by std::lower_bound
			auto it_r_inds = std::lower_bound(r_inds.begin(), r_inds.end(), j);
			j += std::distance(r_inds.begin(), it_r_inds);
			while (*it_r_inds <= j && it_r_inds != r_inds.end()) {
				++it_r_inds;
				j++;
			}
			rowperm.emplace_back(j);
		}
		// the last appended rows will be permuted the desired indices in r_inds
		rowperm.insert(std::end(rowperm), std::begin(r_inds), std::end(r_inds));
		permute_rows(rowperm);
	}


	// insert rows, e.g.,
	// r_inds = [1,3]
	// r_col_inds = [[2,4,6], [1,3,6]]
	// r_col_vals = [[1,2,1], [1,2,2]]
	// will add
	void insert_rows(
		const std::vector<size_t>& r_inds,
		const std::vector<std::vector<size_t>>& r_col_inds,
		const std::vector<std::vector<val_type>>& r_col_vals
	) {
		for (size_t i = 0; i < r_inds.size(); i++){
			for (size_t j= 0; j < r_col_inds.size(); j++){
				col[j].insert_rows({r_col_inds[i][j]}, {r_col_vals[i][j]});
			}
		}
	}

	// insert a zero row
	void insert_row(const size_t& ind){
		// first append a row to the end
		m++;
		// second permute it
		std::vector<size_t> rowperm;
		rowperm.reserve(m);
		for (size_t i = 0; i < m ; i++){
			size_t j = i; // j is the new index of all rows, rows above ind keep the same
			if (i >= ind) j = i+1; // rows above ind will add up one
			if (i == m-1) j = ind; // the last one will be permuted the desired ind
			rowperm.emplace_back(j);
		}
		permute_rows(rowperm);
	}

	/**
	insert zero rows at specified locations
	*/
	void insert_rows(const std::vector<size_t>& r_inds) {
		// figure out how many are inserted before first row
		size_t ninserted = 0;
		while (ninserted < r_inds.size() && r_inds[ninserted] == ninserted) { ++ninserted; }

		for (size_t j = 0; j < n; ++j) {
			col[j].insert_rows(r_inds, ninserted);
		}
		m = m + r_inds.size();
	}

	/**
	insert rows represented by sparse vectors
	*/
	void insert_rows(
		const std::vector<size_t>& r_inds,
		const std::vector<TC>& rows
	) {
		auto t0 = std::chrono::steady_clock::now();
		size_t iloc = m; // first insert location
		size_t m0 = m; // inital m
		m = m + r_inds.size(); // add zero rows to end
		// first figure out permutation
		std::vector<size_t> rperm(m);
		size_t i0 = 0; // original index
		size_t i1 = 0; // insertion index
		size_t i = 0; // permutation index
		while (i1 < r_inds.size() && i < m) {
			// insert rows until we get to an insertion index
			while (r_inds[i1] > i0) {
				rperm[i++] = i0++;
			}
			// insert index
			rperm[i++] = iloc++;
			++i1;
		}
		// insert any remaining rows
		while (i < m) {
			rperm[i++] = i0++;
		}
		auto t1 = std::chrono::steady_clock::now();
		std::cout << "\t\t\t  form perm : "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
			<< "ms" << std::endl;

		// now insert rows at end of columns
		t0 = std::chrono::steady_clock::now();
		for (i1 = 0; i1 < r_inds.size(); ++i1) {
			for (auto it = rows[i1].nzbegin(); it != rows[i1].nzend(); ++it) {
				col[it->ind].emplace_back(i1+m0, it->val);
			}
		}
		t1 = std::chrono::steady_clock::now();
		std::cout << "\t\t\t  append rows : "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
			<< "ms" << std::endl;

		// finally, apply permuation
		t0 = std::chrono::steady_clock::now();
		auto irperm = bats::util::sortperm(rperm);
		permute_rows(irperm);
		t1 = std::chrono::steady_clock::now();
		std::cout << "\t\t\t  apply perm : "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
			<< "ms" << std::endl;

	}

	/**
	insert zero rows at specified locations
	new_rows[i] at index r_inds[i]
	*/
	// void insert_rows(
	// 	const std::vector<size_t>& r_inds,
	// 	std::vector<TC>& new_rows
	// ) {
	// 	// new rows to columns
	// 	for (auto new_row: new_rows){

	// 	}
	// 	m = m + r_inds.size();
	// 	for (size_t j = 0; j < n; ++j) {
	// 		col[j].insert_rows(r_inds, new_rows[i]);
	// 	}
	// }

	void erase_column(const size_t& index){
    	col.erase(col.begin()+index);
		--n;
	}

	void erase_column(){
    	col.erase(--col.end());
		--n;
	}

	/**
	deletes ndelete columns at end of matrix
	*/
	void erase_final_columns(const size_t& ndelete){
    	// col.erase(col.end()-n_dele,col.end());//
		n -= ndelete;
		col.resize(n); // shrinks to new size
	}

	/**
	deletes ndelete columns at beginning of matrix
	*/
	void erase_initial_cols(const size_t ndelete) {
		if ( ndelete == 0 ) { return; } // do nothing
		// rotate first ndelete columns to end - left rotation by ndelete
		std::rotate(col.begin(), col.begin() + ndelete, col.end());
		n -= ndelete;
		col.resize(n);
	}

	/**
	erase specified row
	*/
	void erase_row(const size_t& index){
		for (size_t j = 0; j < ncol(); j++) {
			col[j].erase_for_matrix(index);
		}
		--m;
	}

	/**
	erase last row
	*/
	void erase_row(){
		for (size_t j = 0; j < ncol(); j++) {
			col[j].erase_last_row_of_matrix(m-1);
		}
		--m;
	}

	/**
	erase last a few rows
	*/
	void erase_rows(size_t n_delete_rows){
		assert(n_delete_rows <= nrow()); // the number of deleting rows should be less than total
		for (size_t j = 0; j < ncol(); j++) {
			// col[j].erase_last_rows_of_matrix(m-1, n_delete_rows);
			if (n_delete_rows < nrow()){
				col[j].erase_rows_after(m - n_delete_rows - 1);
			}else{
				col[j].erase_all();
			}
		}
		m -= n_delete_rows;
	}

	/**
	assumes that last row is zero
	so we only need to decrement number of rows
	*/
	void erase_row_unsafe(){
		--m;
	}

	/**
	assumes that last row is zero
	so we only need to decrement number of rows
	*/
	void erase_final_rows_unsafe(size_t n_delete_rows){
		assert(n_delete_rows <= nrow()); // the number of deleting rows should be less than total
		m -= n_delete_rows;
	}

	/**
	erases initial rows
	*/
	void erase_initial_rows(size_t ndelete) {
		for (auto& cj : col) {
			cj.erase_initial_rows(ndelete);
		}
		m -= ndelete;
	}

    inline TC& operator[](size_t index) { return col[index];}
    inline const TC& operator[](size_t index) const { return col[index];}

    inline auto operator()(size_t i, size_t j) const {
        return col[j][i];
    }

	bool operator==(const ColumnMatrix &other) const {
		if (m != other.m || n != other.n) {return false;}
		for (size_t j = 0; j < n; j++) {
			if (col[j] != other.col[j]) {return false;}
		}
		return true;
	}



	std::vector<std::vector<val_type>> to_row_array() const {
		std::vector<std::vector<val_type>> A(m);
		for (size_t i = 0; i < m; i++) {
			for (size_t j = 0; j < n; j++) {
				A[i].emplace_back(getval(i,j));
			}
		}
		return A;
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



	ColumnMatrix submatrix(
		const std::vector<size_t> &rind,
        const std::vector<size_t> &cind
    ) const {
		std::vector<TC> newcol;
		newcol.reserve(cind.size());

		auto prow = bats::util::partial_perm(rind, nrow());

        // loop in column permutation order
        for ( size_t j : cind) {
			newcol.emplace_back(col[j].subvector(prow));
        }
		return ColumnMatrix(rind.size(), cind.size(), newcol);
	}

	// get block A[i0:i1, j0:j1]
	// i1,j1 are not inclusive
	ColumnMatrix block(
		const size_t i0,
		const size_t i1,
		const size_t j0,
		const size_t j1
	) const {
		std::vector<TC> newcol;
		newcol.reserve(j1 - j0);

		for (size_t j = j0; j < j1; j++) {
			newcol.emplace_back(col[j].block(i0, i1));
		}

		return ColumnMatrix(i1 - i0, j1 - j0, newcol);
	}

	void set_block(
		const size_t i0,
		const size_t i1,
		const size_t j0,
		const size_t j1,
		const ColumnMatrix& B
	) {

		for (size_t j = j0, k=0; j < j1; j++, k++) {
			col[j].set_block(i0, i1, B[k]);
		}
	}

	/**
	clear rows i for which c[i] is true
	use vector of bools for quick lookup - vector of inds would require search
	*/
	void clear_rows(const std::vector<bool> &c) {
		if (c.size() != nrow()) {throw std::runtime_error("input vector does not match number of rows.");}
		for (size_t j = 0; j < ncol(); ++j) {
			col[j].clear_inds(c);
		}
	}

	/**
	clear columns j for which c[j] is true

	frees memory as well
	*/
	void clear_cols(const std::vector<bool> &c) {
		if (c.size() != ncol()) {throw std::runtime_error("input vector does not match number of columns.");}
		for (size_t j = 0; j < ncol(); ++j) {
			if (c[j]) {
				col[j].clear_dealloc();
			}
		}
	}

	// swap rows i, i2
	void swap_rows(const size_t i, const size_t i2) {
		for (size_t j = 0; j < ncol(); j++) {
			col[j].swap_rows(i,i2);
		}
	}

	void mix_rows(const size_t i, const size_t i2, const val_type& a, const val_type& b, const val_type& c, const val_type& d) {
		for (size_t j = 0; j < ncol(); j++) {
			col[j].mix_rows(i, i2, a, b, c, d);
		}
	}

	void add_rows(const size_t i, const val_type& c, const size_t i2) {
		for (size_t j = 0; j < ncol(); j++) {
			col[j].add_rows(i, c, i2);
		}
	}




    // addition, substraction, scalar multiplication

    // gemv
    TC gemv(const TC &x) const {
        TC y;  // zero initializer
		tmp_type tmp; // for use with axpy
        // loop over nonzero indices of x
        for (auto xit = x.nzbegin(); xit != x.nzend(); ++xit) {
            y.axpy((*xit).val, col[(*xit).ind], tmp); // y <- x[j]*A[j]
        }
        return y;
    }


	TC gemv(const TC &x, tmp_type& tmpx) const {
        TC y;  // zero initializer
        // loop over nonzero indices of x
        for (auto xit = x.nzbegin(); xit != x.nzend(); ++xit) {
            y.axpy((*xit).val, col[(*xit).ind], tmpx); // y <- x[j]*A[j]
        }
        return y;
    }

	// scalar multiplication
	ColumnMatrix operator*(const val_type a) const {
		std::vector<TC> newcol;
		newcol.reserve(n);
		for (size_t j = 0; j < n; j++) {
			newcol.emplace_back(col[j] * a);
		}
		return ColumnMatrix(m, n, newcol);
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

	/**
	multiplication x^T A
	*/
	friend TC operator*(const TC& x, const ColumnMatrix& A) {
		std::vector<size_t> ind;
		std::vector<typename TC::val_type> val;
		for (size_t j = 0; j < A.ncol(); ++j) {
			auto a = x * A[j];
			if (a != 0) {
				ind.emplace_back(j);
				val.emplace_back(a);
			}
		}
		return TC(ind, val);
	}


    ColumnMatrix operator+(const ColumnMatrix &B) const {
        ColumnMatrix C(B); // copy constructor
        for (size_t j = 0; j < n; j++) {
            C.col[j].axpy(1, col[j]);
        }
        return C;
    }

	ColumnMatrix& operator+=(const ColumnMatrix &B) {
        for (size_t j = 0; j < n; j++) {
            col[j].axpy(1, B.col[j]);
        }
        return *this;
    }

    ColumnMatrix operator-(const ColumnMatrix &B) const {
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
			for (auto it = col[j].nzbegin(); it != col[j].nzend(); ++it) {
				tcol[it->ind].emplace_back(j, it->val);
			}
		}
		return ColumnMatrix(n, m, tcol);
	}

	inline ColumnMatrix T() const { return transpose(); }

	// permutations: permute, permute_rows, permute_cols

    // permute columns in-place
    // TODO: evaluate if this is best method.
    inline void permute_cols(const std::vector<size_t> &colperm) {
        bats::util::apply_perm_swap(col, colperm);
    }

	inline void ipermute_cols(const std::vector<size_t> &colperm) {
        bats::util::apply_iperm_swap(col, colperm);
    }

	// vj <- a * vj + c * vk
	// vk <- b * vj + d * vk
	void mix_cols(const size_t j, const size_t k, const val_type& a, const val_type& b, const val_type& c, const val_type& d) {
		TC cj(col[j]);
		TC ck(col[k]);
		col[j].scale_inplace(a);
		col[j].axpy(c, ck);
		col[k].scale_inplace(d);
		col[k].axpy(b, cj);
	}


    // permute both rows and columns
    void permute(const std::vector<size_t> &rowperm,
        const std::vector<size_t> &colperm) {
        permute_cols(colperm);
        permute_rows(rowperm);
    }

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

	// apply inverse of diagonal matrix with coeff vector along diagonal on the left
	ColumnMatrix& col_inv_scale(const std::vector<val_type> &coeff) {
		for (size_t j = 0; j < n; j++) {
			col[j].scale_inplace(coeff[j].inv());
		}
		return *this;
	}

	// number of nonzeros
	size_t nnz() const {
		size_t ct = 0;
		for (size_t j = 0; j < n; j++) {
			ct += col[j].nnz();
		}
		return ct;
	}

	// tests to see if matrix is zero
	bool is_zero() const {
		for (size_t j = 0; j < n; j++) {
			if( col[j].nnz() > 0) {
				return false;
			}
		}
		return true;
	}

	// tests to see if matrix has structure
	bool is_upper() const {
		for (size_t j = 0; j < n; ++j) {
			auto iv = col[j].nzend();
			if (iv != col[j].nzbegin()) {
				--iv;
				if (iv->ind > j) {
					return false;
				}
			}
		}
		return true;
	}

	// tests to see if matrix is invertible and upper-triangular
	bool is_upper_invert() const {
		for (size_t j = 0; j < n; ++j) {
			auto iv = col[j].nzend();
			if (iv != col[j].nzbegin()) {
				--iv;
				if (iv->ind != j) {
					return false;
				}
			}
		}
		return true;
	}

	bool is_lower() const {
		for (size_t j = 0; j < n; j++) {
			auto iv = col[j].nzbegin();
			if (iv != col[j].nzend()) {
				if (iv->ind < j) {
					return false;
				}
			}
		}
		return true;
	}
	bool is_reduced() const {
		// check whether columns have unique lower index (if it exists)
		std::map<size_t, size_t> p2c;
		for (size_t j = 0; j < n; j++) {
			auto st = col[j].nzbegin();
			auto end = col[j].nzend();
			if (st != end) {
				end--;
				if (p2c.count(end->ind) > 0) {
					// pivot already exists - return false
					return false;
				} else {
					p2c[end->ind] = j;
				}
			}
		}
		return true;
	}
	bool is_pivot_matrix() const {
		for (size_t j = 0; j < n; j++) {
			if (col[j].nnz() > 1) {
				return false;
			}
		}
		// check that rows also have at most 1 nonzero
		ColumnMatrix row = transpose();
		for (size_t i = 0; i < m; i++) {
			if (row[i].nnz() > 1) {
				return false;
			}
		}
		return true;
	}
	bool is_EL() const {
		size_t j = 0;
		ssize_t i = -1;
		while (j < n) {
			auto st = col[j].nzbegin();
			auto end = col[j].nzend();
			// check if no pivots
			if (std::distance(st, end) == 0) { j++; break; }

			// check that pivot is strictly increasing
			if (std::distance(st, end) > 1 || ssize_t(st->ind) <= i) { return false; }

			i = st->ind;
			j++;
		}
		// enter this loop after we encounter a column with no zeros
		while (j < n) {
			if (col[j].nnz() == 0) {
				j++;
			} else {
				return false;
			}
		}
		return true;
	}
	inline bool is_EU() const { return transpose().is_EL(); }
	inline bool is_ELhat() const { return J_conjugation().is_EU(); }
	inline bool is_EUhat() const {return J_conjugation().is_EL(); }


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

    // template <typename IO>
    void write(std::ostream &io) const {
		io << "ColumnMatrix\n";
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

	// static methods
	static ColumnMatrix identity(size_t n) {
		std::vector<TC> col(n);
		for (size_t j = 0; j < n; j++) {
			col[j] = TC(j);
		}
		return ColumnMatrix(n, n, col);
	}

	static ColumnMatrix random(size_t m, size_t n, double p, int maxval, std::default_random_engine &generator) {
		std::vector<TC> col(n);
		for (size_t j = 0; j < n; j++) {
			col[j] = TC::random(m, p, maxval, generator);
		}
		return ColumnMatrix(m, n, col);
	}
	static ColumnMatrix random(size_t m, size_t n, double p, int maxval) {
		// obtain a seed from the system clock:
		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		std::default_random_engine generator(seed);
		return random(m, n, p, maxval, generator);
	}

	// tensor product A \otimes B
	// use Kronecker product ordering
	friend ColumnMatrix tensor_product(const ColumnMatrix& A, const ColumnMatrix& B) {
		size_t Am = A.nrow();
		size_t An = A.ncol();
		size_t Bm = B.nrow();
		size_t Bn = B.ncol();
		std::vector<TC> newcol;
		for (size_t i = 0; i < An; i++) {
			for (size_t j = 0; j < Bn; j++) {
				newcol.emplace_back(A[j].tensor(B[j]));
			}
		}
		return ColumnMatrix(Am * Bm, An * Bn, newcol);
	}

	// direct sum A \oplus B
	friend ColumnMatrix direct_sum(const ColumnMatrix& A, const ColumnMatrix& B) {
		size_t Am = A.nrow();
		size_t An = A.ncol();
		size_t Bm = B.nrow();
		size_t Bn = B.ncol();
		std::vector<TC> newcol;
		for (size_t j = 0; j < An; j++) {
			newcol.emplace_back(A[j]);
		}
		for (size_t j = 0; j < Bn; j++) {
			newcol.emplace_back(B[j].shift_inds(Am));
		}
		return ColumnMatrix(Am + Bm, An + Bn, newcol);
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
	typename TC::tmp_type tmp; // for use with axpy
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
        x.axpy(-a, U[j], 0, j, tmp);
        // get next nonzero
        xit = x.upper_bound(j-1);
    }
    return x;
}

template <class TC>
TC l_solve(const ColumnMatrix<TC> &L, const TC &y) {
    TC x = y;
    if (x.nnz() == 0) { return x; }
	typename TC::tmp_type tmp; // for use with axpy
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
        x.axpy(-a, L[j], j+1, m, tmp);
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

// solve L \ A
template <class TC>
ColumnMatrix<TC> l_solve(const ColumnMatrix<TC> &L, const ColumnMatrix<TC> &A) {
    //std::cout << "entering solve" << std::endl;
    size_t m = A.nrow();
    size_t n = A.ncol();
    std::vector<TC> col;
    col.reserve(n);
    for (size_t j = 0; j < n; j++) {

        col.emplace_back(
            l_solve(L, A[j])
        );
    }
    return ColumnMatrix<TC>(m, n, col);
}

// form inverse U^{-1}
template <class TC>
ColumnMatrix<TC> u_inv(const ColumnMatrix<TC> &U) {
    //std::cout << "entering solve" << std::endl;
    size_t m = U.nrow();
    // size_t n = U.ncol();
	return u_solve(U, ColumnMatrix<TC>::identity(m));
}

// form inverse L^{-1}
template <class TC>
ColumnMatrix<TC> l_inv(const ColumnMatrix<TC> &L) {
    //std::cout << "entering solve" << std::endl;
    size_t m = L.nrow();
    // size_t n = U.ncol();
	return l_solve(L, ColumnMatrix<TC>::identity(m));
}
