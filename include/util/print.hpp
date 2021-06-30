#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>

namespace bats{
template <typename CpxT, typename T>
void print_summary_of_filtration(const CpxT& X, std::function<T(const std::vector<size_t>&)>& filtfn){
    std::cout << "\nLet's see the filtration value on each dimension" << std::endl;
    for (size_t k = 0; k < X.maxdim() + 1; k++) {
        std::cout << "For dimension "<< k << std::endl;
        for (auto& s : X.get_simplices(k)) {
            std::cout << "simplex with index "<< X.find_idx(s) << ", filtration value f(s) is " << filtfn(s) << std::endl; 
        }
    }
}


//two ways to display 2D vectors
template<typename T>
void print_2D_vectors (const T& perms) {
//perms represents permutations
    for (auto& valsk: perms) {
        for (auto& v: valsk) {
            std::cout << v << " ";
        }
        std::cout << "\n";
    }
}
template<typename T>
void print_1D_vectors (const T& perm) {
    std::cout << std::boolalpha; 
    for (size_t i = 0; i < perm.size(); i++)
    {
        std::cout << perm[i] << " ";
    }
    std::cout << "\n";
}

template<typename T>
void print_simplex (const T& perm) {
    std::cout << "{" ;
    for (auto v = perm.begin(); v < perm.end() - 1; v++ ) {
        std::cout << *v << ", ";
    }
    std::cout << *(perm.end()-1) << "} " ;
}

template <class Filtration>
void print_filtration_info(const Filtration& F){

    auto F_complex = F.complex();
    std::cout << "\nRips Filtration Values are" << std::endl;
    auto filtration_vals = F.vals();

    for (size_t i = 0; i <= F.maxdim(); i++) {
        std::cout << F.ncells(i) << " cells in dim " << i << ":"<< std::endl;
        auto simplices_i = F_complex.get_simplices(i);
        for (size_t j = 0; j < F.ncells(i); j++){
            print_simplex(simplices_i[j]);
            std::cout << ": "<< filtration_vals[i][j] << std::endl;   
        }
    }
}

}