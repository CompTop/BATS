#include <iostream>
#include <bats.hpp>
#include "vineyard.hpp"

using namespace std;
#define F int
#define F3 ModP<int, 3>


int main() {
    // how to construct CSC matrix 
    // std::vector<size_t> rowind = {0,0,1,1,2};
    // std::vector<size_t> colptr = {0,1,2,3,5,5}; 
    // std::vector<int> val = {1,2,3,4,5};
    
    // auto M = CSCMatrix<int, size_t>(5,5, colptr, rowind, val);
    // M.print();

    // how to construct Vineyard matrix
    std::vector<SparseVector<F, size_t>> cols;
    std::vector<size_t> ind;
    std::vector<F> val;
    ind = {1,2};
    val = {1,2};
    SparseVector<F, size_t> a(ind, val);
    cols.emplace_back(a);

    ind = {3,4};
    val = {3,4};
    SparseVector<F, size_t> b(ind, val);
    cols.emplace_back(b);

    VineyardMatrix M (6, 6, cols);


    return 0;
}