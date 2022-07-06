#include <iostream>
#include <bats.hpp>
#include "vineyard.hpp"

using namespace std;
#define F ModP<int, 5>
#define F3 ModP<int, 3>


int main() {
    // Construct an Vineyard matrix
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

    VineyardMatrix<F> M (5, 5, cols);
    M.print();


    // std::cout << "after column permutation" << std::endl;
    // M.permute_cols({2,0,1,3,4});
    // M.print();
    // std::cout << "after row permutation" << std::endl;
    // M.permute_rows({0,2,1,3,4});
    // M.print();

    // construct a CSC matrix 
    std::vector<size_t> rowind = {0,0,1,1,2};
    std::vector<size_t> colptr = {0,1,2,3,5,5}; 
    std::vector<int> val2 = {1,2,3,4,5};
    
    auto M2 = CSCMatrix<int, size_t>(5,5, colptr, rowind, val2);
    std::cout << "A CSC Matrix" << std::endl;
    M2.print();

    std::cout << "CSC to Vineyard" << std::endl;
    VineyardMatrix<F> M3(M2);
    M3.print();
    std::cout << "Add cols[0] += 2 * cols[1] " << std::endl;
    M3.test_func();
    M3.print();
    return 0;
}