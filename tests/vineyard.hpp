#pragma once
#include <vector>
#include <iostream>
#include <bats.hpp>
#include <memory>
#include "linkded_list.hpp"

// forward declarations
// template <class TC>
// class VineyardMatrix;

// template over column type
template <typename TV>
class VineyardMatrix{
private:
    size_t m = 0; // number of rows
    size_t n = 0; // number of columns
    std::vector<LinkedList<TV>> cols; // vector of columns 
    std::vector<shared_ptr<size_t>> row_inds_ptr; // pointers of row indices shared by all columns
public:
	using col_type = LinkedList<TV>; // column type
    using TC = col_type;

    // default constructor
    VineyardMatrix() {}
    


    // construct empty matrix.  same as zero matrix
    VineyardMatrix(size_t m, size_t n) : m(m), n(n) {
        // col.reseve(n, TC());
        cols.resize(n);
        initial_row_ptrs(m);
    }

    // construct Vineyard matrix by passing vector of columns(sparse vectors)
    template <typename TF2> // TF2 : value type (of field) 
    VineyardMatrix(size_t m, size_t n, std::vector<SparseVector<TF2, size_t>> _cols) : m(m), n(n) {
        // assert (TF2 == val_type);
        cols.reserve(n);
        initial_row_ptrs(m);
        for (auto other_vec: _cols){ 
            cols.emplace_back(other_vec.nzbegin(),
                            other_vec.nzend(), row_inds_ptr.begin());
        }
        if (_cols.size() < n){ // if only given the first non-zero columns
            for (auto j = _cols.size(); j < n; j++){
                cols.emplace_back();
            }
        }
    }

    // constructor from CSCMatrix over the integers
    VineyardMatrix(const CSCMatrix<int, size_t> &A) : m(A.nrow()), n(A.ncol()) {
        // reserve capacity
        cols.reserve(n);
        initial_row_ptrs(m);

        auto colptr = A.get_colptr(); // column pointer of CSC
        auto rowind = A.get_rowind();
        auto val = A.get_val();
        auto rptr = rowind.cbegin();// row index pointer of CSC
        auto vptr = val.cbegin(); // value pointer of CSC
        
        for (size_t j = 0; j < n; j++) {
            size_t c_start_ind = colptr[j];
            size_t col_size =  colptr[j+1] - colptr[j];
            cols.emplace_back(rptr+c_start_ind, 
                                vptr+c_start_ind, 
                                row_inds_ptr.cbegin(),
                                col_size);
        }
    }

    // Initialize vector of row index pointers
    void initial_row_ptrs(size_t m){
        row_inds_ptr.reserve(m);
        for (size_t i = 0; i < m; i++){
            row_inds_ptr.emplace_back(std::make_shared<size_t>(i));
        }
    }

    inline col_type& operator[](size_t index) { return cols[index];}
    inline const col_type& operator[](size_t index) const { return cols[index];}

    auto getval(const size_t i, const size_t j) const {
        return cols[j].getval(i);
    }

    void print_size() const {
        std::cout << "[" << this << "] : " << m << " x " << n <<\
        " VineyardMatrix" << std::endl;
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
    void display(){

    }

    // permute rows in-place
    inline void permute_rows(const std::vector<size_t> &rowperm) {
        if (row_inds_ptr.size() != rowperm.size()){ // TODO: make it a if Debug 
            assert(row_inds_ptr.size() == rowperm.size());
        }
        auto it_row = row_inds_ptr.begin();
        auto it_perm = rowperm.begin();
        while (it_row < row_inds_ptr.end()){
            **it_row = *it_perm;
            it_perm++;
            it_row++;
        }
    }

    // permute columns in-place
    inline void permute_cols(const std::vector<size_t> &colperm) {
        bats::util::apply_perm_swap(cols, colperm);
    }

	inline void ipermute_cols(const std::vector<size_t> &colperm) {
        bats::util::apply_iperm_swap(cols, colperm);
    }

    // static methods
	static VineyardMatrix identity(size_t n) {
		std::vector<SparseVector<TV>> col(n);
		for (size_t j = 0; j < n; j++) {
			col[j] = SparseVector<TV>(j);
		}
		return VineyardMatrix(n, n, col);
	}

    // TODO: Simplicial complex to Vineyard directly without using CSC 
    
    void test_func(){
        cols[0].axpy(2, cols[1]);
    }

    inline std::vector<col_type>& get_cols() { return cols; }

};