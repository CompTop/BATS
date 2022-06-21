#pragma once
#include <vector>
#include <iostream>
#include <bats.hpp>
#include "linkded_list.hpp"

// forward declarations
// template <class TC>
// class VineyardMatrix;

// template over column type
// template <class TC>
class VineyardMatrix{
private:
    size_t m = 0; // number of rows
    size_t n = 0; // number of columns
    std::vector<LinkedList> cols; // vector of columns 
    std::vector<size_t *> row_inds_ptr; // pointers of row indices shared by all columns
public:
    // using val_type = typename TC::val_type;
	// using col_type = TC;
	// using tmp_type = typename TC::tmp_type; // for use with axpy

    // default constructor
    VineyardMatrix() {}
    
    ~VineyardMatrix() {
        clear();
    }
    // Delete the entire list
    void clear() {
        // delete all elements in a vector
        auto it = row_inds_ptr.begin(); 
        while (it != row_inds_ptr.end()) {
            it = row_inds_ptr.erase(it);
        }
    }


    // construct empty matrix.  same as zero matrix
    VineyardMatrix(size_t m, size_t n) : m(m), n(n) {
        // col.resize(n, TC());
        cols.resize(n);
        // row_inds_ptr.resize(m, *size_t);
    }

    // construct Vineyard matrix by passing vector of columns(sparse vectors)
    template <typename TF2> // TF2 : value type (of field) 
    VineyardMatrix(size_t m, size_t n, std::vector<SparseVector<TF2, size_t>> _cols) : m(m), n(n) {
        // assert (TF2 == val_type);
        // TODO: resize col and row_inds_ptr properly
        // col.resize(n, TC());
        // row_inds_ptr.resize(m, *size_t);
        for (size_t i = 0; i < m; i++){
            row_inds_ptr.emplace_back(new size_t(i));
        }
        
        for (auto other_vec: _cols){
            cols.emplace_back(LinkedList(other_vec.nzbegin(),
                                        other_vec.nzend(), 
                                        row_inds_ptr.cbegin()));
        }

    }

    // TODO: Row permutation (hard)

    // TODO: Column permutation (easy)

    // TODO: operator overload eg. indexing []

    // TODO: Simplicial complex to Vineyard directly without using CSC 

    // constructor from CSCMatrix over the integers
    // TODO: CSC to Vineyard 
    // ColumnMatrix(const CSCMatrix<int, size_t> &A) : m(A.nrow()), n(A.ncol()) {
    //     col.reserve(n);
    //     auto colptr = A.get_colptr();
    //     auto rowind = A.get_rowind();
    //     auto val = A.get_val();
    //     auto rptr = rowind.cbegin();
    //     auto vptr = val.cbegin();
    //     
    //     for (size_t j = 0; j < n; j++) {
    //         // use iterator constructor
    //         col.emplace_back(TC(
    //             rptr + colptr[j],
    //             vptr + colptr[j],
    //             colptr[j+1] - colptr[j]
    //         ));
    //     }
    // }


};