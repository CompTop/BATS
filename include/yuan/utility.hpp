#pragma once
#include <bats.hpp>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>

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

// Reduction algorithm of matrix permutation version 
// of updating persistence by passing old RU factorzation 
// and row and column permutations 
template <class TC>
void update_reduction(ColumnMatrix<TC> &R_k, ColumnMatrix<TC> &U_k, const ColumnMatrix<TC> & row_perm, const ColumnMatrix<TC> &column_perm){
	
	using MatT = ColumnMatrix<TC>;
	MatT P_k_minus_1 = row_perm;
	MatT P_k = column_perm;
	    
	auto F = UQL(P_k.T() * U_k);                                                
    MatT U = F.U;        
    MatT L = F.L;
    MatT Q = F.E;

    MatT R_k_prime = P_k_minus_1 * R_k * l_inv(L) * Q.T();
    MatT U_double_prime = MatT::identity(R_k_prime.ncol());
    // bats::reduce_matrix(R_k_prime, U_double_prime);
    reduce_matrix(R_k_prime, U_double_prime);

    MatT U_k_prime = U * U_double_prime;

    R_k = R_k_prime;
    U_k = U_k_prime;    
}

// Reduction algorithm of list permutation version(maybe faster)
// for updating persistence by passing old RU factorzation 
// and row and column permutations (in the form of list permutation)
template <class TC>
void update_reduction(ColumnMatrix<TC> &R_k, ColumnMatrix<TC> &U_k, const std::vector<size_t> &row_perm, const std::vector<size_t> &column_perm){
	
	using MatT = ColumnMatrix<TC>;
    U_k.permute_rows(bats::util::inv_perm(column_perm));
	auto F = UQL(U_k);                                                
    MatT U = F.U;        
    MatT L = F.L;
    MatT Q = F.E;

    MatT R_k_prime = R_k * l_inv(L) * Q.T();
    R_k_prime.permute_rows(row_perm);
    MatT U_double_prime = MatT::identity(R_k_prime.ncol());
    // bats::reduce_matrix(R_k_prime, U_double_prime);
    reduce_matrix(R_k_prime, U_double_prime);

    MatT U_k_prime = U * U_double_prime;

    R_k = R_k_prime;
    U_k = U_k_prime;    
}

// Return an identity permutation, k is its length
std::vector<size_t> identity_perm(const size_t& k){
    std::vector<size_t> l(k);
    std::iota(l.begin(), l.end(), 0);
    return l;
}

/*
If we sort vector v, this function will return the indices of sorted values 
*/
template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

    // initialize original index locations
    std::vector<size_t> idx(v.size());
    // fill the vector with increasing numbers starting from 0
    std::iota(idx.begin(), idx.end(), 0); 

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values 
    std::stable_sort(idx.begin(), idx.end(),
        [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

    return idx;
}

/*
Make indices list to permutation, e.g., [5 0 2] to [2 0 1]

Given a list of correpsonding old indices of permutations, find the permutation
e.g., given [5 0 2], which is the corresponding old indices of an old vector v_old 
for a new vector v: 
v[0] = v_old[5]
v[1] = v_old[0]
v[2] = v_old[2]
need to find the permutation [2 0 1]
*/
template <typename T>
std::vector<size_t> find_perm_from_vector(const std::vector<T> &v) {

    // initialize original index locations
    std::vector<size_t> idx(v.size());
    // fill the vector with increasing numbers starting from 0
    std::iota(idx.begin(), idx.end(), 0); 

    // sort indexes based on comparing values in v
    std::stable_sort(idx.begin(), idx.end(),
        [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

    // get the desired permutation
    std::vector<size_t> perm(v.size());
    for(size_t i = 0; i < v.size(); i++){
        perm[idx[i]] = i;
    }
    return perm;
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



// find the indices of a boundary in the list of simplices of one dimensional lower
template <typename CpxT>
std::vector<size_t> get_bd_index(CpxT X, size_t dim, size_t index){
    auto [inds,vals] = X.boundary(dim, index);
    return inds;
}

// get the permutation for permuting i-th element to the end 
std::vector<size_t> perm_to_the_end(const size_t& index, const size_t& length){
    std::vector<size_t> v;
    v.reserve(length);

    for(size_t i = 0; i < length; i++){
        if(i!=index)v.emplace_back(i);
    }
    v.emplace_back(index);
    return v;
}

// get the permutation for permuting elements, with their indices in a list, to the end 
// the index_list should be sorted 
std::vector<size_t> perm_to_the_end(const std::vector<size_t>& index_list, const size_t& length){
    std::vector<size_t> v;
    v.reserve(length);
    auto it = index_list.begin();
    for(size_t i = 0; i < length; i++){
        if(i != *it){
            v.emplace_back(i);
        }else{
            it++;
        }
    }

    for(auto& element: index_list){
        v.emplace_back(element);
    }
    return v;
}

// extend a permutation to a desired lenght 
// with the elements appended unmoved, e.g.
// (2 0 1) to (2 0 1 3 4 5)
std::vector<size_t> extension_perm(const std::vector<size_t>& perm, const size_t& length){
    if(!perm.empty()){
        std::vector<size_t> v;
        v.reserve(length);
        size_t i = 0;
        for(auto& element: perm){
            v.emplace_back(element);
            i++;
        }
        while (i<length)
        {
            v.emplace_back(i);
            i++;
        }
        return v;
    }else
    {
        return identity_perm(length);
    }
    
}


// permute a general obeject(could be a vector) that could be dereferenced 
template <typename T>
void general_permute(const std::vector<size_t>  &perm, std::vector<T>& object) {
    bats::util::apply_perm_swap(object, perm);
}


void permutation_multi(std::vector<size_t>& perm1,const std::vector<size_t>& perm2){
    bats::util::apply_perm(perm1, perm2);
}

// check if two RFCC are the same 
template<typename T>
bool test_reduce_result(const T& RFCC2, const T& RFCC){

    auto R = RFCC.RC.R;
    auto R2 = RFCC2.RC.R;
    auto U = RFCC.RC.U;
    auto U2 = RFCC2.RC.U;


    // check homology dimension
    bool test_result = true;
    for(size_t k =0; k < RFCC2.maxdim(); k++){
        auto R_ps_k = RFCC.persistence_pairs(k);
        auto R2_ps_k = RFCC2.persistence_pairs(k);
        bool test_result_k = barcode_equality(R2_ps_k, R_ps_k);
        if(!test_result_k){
            test_result = false;
            std::cout << "\ntwo RFCC are different!!" << std::endl;
        }
    }
    
    // check RU factorization
    for(size_t k = 0; k < RFCC2.maxdim(); k++){
        auto U_inv_k = u_inv(U[k]); 
        auto U2_inv_k = u_inv(U2[k]); 
        if(!(R[k] * U_inv_k == R2[k] * U2_inv_k)){
            test_result = false;
            std::cout << "\ntwo RFCC are different!!" << std::endl;
        }
    }
    
    return test_result;
}


