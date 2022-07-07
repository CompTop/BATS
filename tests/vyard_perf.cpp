#include <iostream>
#include <memory>
#include "vineyard.hpp"
#include <bats.hpp>
#include <chrono>
// C++ implementation of the approach
#include <bits/stdc++.h>


#define F ModP<int, 5>

// number of seeds for random checks
#define N_SEEDS 4


// Function to return the next random number
int getNum(vector<int>& v)
{
 
    // Size of the vector
    int n = v.size();
 
    // Generate a random number
    std::srand(time(NULL));
 
    // Make sure the number is within
    // the index range
    int index = rand() % n;
 
    // Get random number from the vector
    int num = v[index];
 
    // Remove the number from the vector
    std::swap(v[index], v[n - 1]);
    v.pop_back();
 
    // Return the removed number
    return num;
}
 
// Function to generate n non-repeating random numbers
std::vector<size_t> generateRandom(int n)
{
    std::vector<int> v(n);
 
    // Fill the vector with the values
    // 1, 2, 3, ..., n
    for (int i = 0; i < n; i++)
        v[i] = i;
    
    // While vector has elements
    // get a random number from the vector
    std::vector<size_t> perm;
    perm.reserve(n);
    while (v.size()) {
        perm.emplace_back(getNum(v));
    }
    return perm;
}

int main()
{   
    using MatT2 = VineyardMatrix<F>;
    using VT = SparseVector<F, size_t>;
	using MatT = ColumnMatrix<VT>;

    size_t n_rows = 1000;
    size_t n_cols = 1000;
    auto A = MatT::random(n_rows, n_cols, 0.4, 4);
    auto B = MatT2(n_rows, n_cols, A.cols());
    // A.print();
    // B.print();

    auto perm = generateRandom(n_rows);
    // bats::print_1D_vectors(perm);

    auto start = std::chrono::steady_clock::now();
    A.permute_rows(perm);
    auto end = std::chrono::steady_clock::now();
    std::cout << "Permute rows of standard sparse matrix: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
            << "ms" << std::endl;
    
    start = std::chrono::steady_clock::now();
    B.permute_rows(perm);
    end = std::chrono::steady_clock::now();
    std::cout << "Permute rows of Vineyard matrix: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
            << "ms" << std::endl;
    
    // A.print();
    // B.print();
    // bats::reduce_matrix(B, U);

    return 0;
}