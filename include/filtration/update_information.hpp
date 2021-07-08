#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include "util/print.hpp"
#include "util/permutation.hpp"
/**
We need to specify what is the conventional usage of permutation in bats!!
[2, 0, 1] apply to [2.0 , 1.0 ,5.0] 
a) is [5.0, 2.0, 1.0] in matrix permutation notation 
b) is [1.0, 5.0 ,2.0] in traditioanl notation(in book From Mathematics to Generic Programming)


We take the frist notation here, but notice that BATs, sometimes, will mix the two notations!!!
Check before use!!!
Since the above two notations are inverse to each other, inverse them if needed!
*/

namespace bats{

// two filtration dimensions are assumed to be the same!!
template <class FiltrationType>
struct Update_info{
    // addition information: 
    // a vector stores the indices of simplices in all dimensions
    // that will be added to the new filtration Y
    // (in ascending order of Y)
    std::vector<std::vector<size_t>> addition_indices;
    // alse need the indices of boundaries for k-simplex, 
    // i.e., the indices of (k-1)-simplice are boundaries of a k-simplex
    std::vector<std::vector<std::vector<size_t>>> boundary_indices; 

    // deletion information: 
    // a vector stores the indices of simplices in all dimensions
    // that will be deleted from the old filtration X
    // (in ascending order of X)
    std::vector<std::vector<size_t>> deletion_indices;

    // permutation information: 
    // a vector that will permute the intersected simplices in the new filration B
    // to the correct position
    std::vector<std::vector<size_t>> permutations;
    std::vector<std::vector<size_t>> intersection_indices_Y; // intersection indices in Y (sorted order)
    std::vector<std::vector<size_t>> intersection_indices_X; // their corresponding indices in X

    // max dimension of the new filtration = size of permutations -1 
    size_t max_dim;

    // Filtration values and permutations
    std::vector<std::vector<double>> F_X_vals;
    std::vector<std::vector<size_t>> F_X_perms;
    std::vector<std::vector<size_t>> perms_X_inv;
    std::vector<std::vector<double>> F_Y_vals;
    std::vector<std::vector<size_t>> F_Y_perms;
    std::vector<std::vector<size_t>> perms_Y_inv;

    FiltrationType F_old;
    FiltrationType F_new;
    bool filtered_boolean = false; 

    Update_info(const FiltrationType& F_X, const FiltrationType& F_Y){

        F_old = F_X;
        F_new = F_Y;

        max_dim = F_Y.maxdim();
        F_X_vals = F_X.vals();
        F_Y_vals = F_Y.vals();
        auto perms_Y = bats::filtration_sortperm(F_Y_vals);
        F_Y_perms = perms_Y;

        // step 2 loop over each cell in the new filtration Y
        // get simplicial complex for two filtrations
        auto smcplex_X = F_X.complex();
        auto smcplex_Y = F_Y.complex();

        for (size_t i = 0; i <= max_dim; i++){
            // addition information: 
            // a vector stores the indices of simplices in dimension i
            // that will be added to the new filtration Y
            std::vector<size_t> addition_in_i;
            std::vector<std::vector<size_t>> boundary_indices_in_i;

            // deletion information: 
            // a vector stores the indices of simplices in dimension i
            // that will be deleted from the old filtration X
            std::vector<size_t> deletion_in_i;

            // permutation information: 
            // a vector that will permute the intersected simplices in the new filration B
            // to the correct position
            std::vector<size_t> intersection_in_i_X;
            std::vector<size_t> intersection_in_i_Y;
            std::vector<size_t> perm_in_i;

            // find simplices in dimension i for both X and Y
            auto simplices_Y = smcplex_Y.get_simplices(i); // simplices in dimension i
            auto simplices_X = smcplex_X.get_simplices(i);

            // create a vector of boolean resutls of comparison
            // default values is false, which means a simplex in A is not in B 
            std::vector<bool> bool_results_i(simplices_X.size(), false);

            // loop over each simplex in Y
            for (size_t j = 0; j < simplices_Y.size(); j++){
                // find intersecting indices in X
                auto index = smcplex_X.find_idx(simplices_Y[j]); 
                // if the simplex can be found in X (intersection information)
                if(index != bats::NO_IND){
                    bool_results_i[index] = true;
                    intersection_in_i_X.emplace_back(index);
                    intersection_in_i_Y.emplace_back(j);
                }else{ // if not in X (addition information)
                    addition_in_i.emplace_back(j); //in ascending order!
                    // find the boundary indices of the simplex that will be added
                    if(i>0){ // dimension zero has boundary 0
                        auto [inds,vals] = smcplex_Y.boundary(i, j);
                        boundary_indices_in_i.emplace_back(inds);
                    }
                }
            }
            
            // find the permutation information
            perm_in_i = find_perm_from_vector(intersection_in_i_X);

            // find the deletion indices in X (in ascending order!)
            for(size_t k = 0; k < bool_results_i.size(); k++){
                if(!bool_results_i[k]){ // flase, need to deal with deletion information
                    deletion_in_i.emplace_back(k);
                }
            }

            deletion_indices.emplace_back(deletion_in_i);

            intersection_indices_X.emplace_back(intersection_in_i_X);
            intersection_indices_Y.emplace_back(intersection_in_i_Y);
            permutations.emplace_back(perm_in_i);

            addition_indices.emplace_back(addition_in_i);
            boundary_indices.emplace_back(boundary_indices_in_i);
        }
    }



    void filtered_info(const std::vector<std::vector<size_t>>&perms_X){
        // step 1: determine permutation order for new filtration values
        F_X_perms = perms_X;

        // Check if filtration is sorted when construction.
        // If so, we do not need this step!
        bool filtration_sort = true;
        for(auto&p: F_X_perms){
            if(p != identity_perm(p.size())) filtration_sort = false;
        }
        for(auto&p: F_Y_perms){
            if(p != identity_perm(p.size())) filtration_sort = false;
        }
        if(filtration_sort){
            std::cout << "\nfiltration is sorted when construction " << std::endl;
            return;
        }
        filtered_boolean = true;
        // compute inverse permutations
        perms_X_inv = bats::filtration_iperm(F_X_perms);
        perms_Y_inv = bats::filtration_iperm(F_Y_perms);

        // step 2: determine the deletion indices in X(old filtration) after permutation

        // Note, for the considertation of addtion and deletion convenience,
        // we need to sort the addition by ascending order and, 
        // deletion indices by descending order!!! 
        for (size_t i = 0; i < perms_X_inv.size(); i++)
        {
            // get the list of new indices
            auto perm_inv_X = perms_X_inv[i];

            // set new deletion indices
            for (size_t j = 0; j < deletion_indices[i].size(); j++){
                deletion_indices[i][j] = perm_inv_X[deletion_indices[i][j]];
            }
            
            // step 3: determine the addition indices in Y(new filtration) after permutation
            // get the list of new indices
            auto perm_inv_Y = perms_Y_inv[i];
            
            // set new indices
            for (size_t j = 0; j < addition_indices[i].size(); j++){
                addition_indices[i][j] = perm_inv_Y[addition_indices[i][j]];
            }
            
            // step 4: determine the permutation of intersection between X and Y (filtered)
            // first, get the list of new indices of intersection in X after filtration sort
            // "filtration sort": means to sort simplices in filtration by their filtration values 
            std::vector<size_t> ind_X_i = intersection_indices_X[i];

            std::vector<size_t> ind_X_i_sorted;
            ind_X_i_sorted.reserve(ind_X_i.size());
            for (auto& ind: ind_X_i){
                ind_X_i_sorted.emplace_back(perm_inv_X[ind]);
            }
            intersection_indices_X[i] = ind_X_i_sorted;

            // second, get the list of new indices of intersection in Y after filtration sort
            std::vector<size_t> ind_Y_i = intersection_indices_Y[i];

            std::vector<size_t> ind_Y_i_sorted;
            ind_Y_i_sorted.reserve(ind_Y_i.size());
            for (auto& ind: ind_Y_i){
                ind_Y_i_sorted.emplace_back(perm_inv_Y[ind]);
            }
            intersection_indices_Y[i] = ind_Y_i_sorted;

            // third, get the permutation needed for permute the new indices of X
            std::vector<size_t> p = bats::util::sortperm(ind_X_i_sorted);
            // forth, perform permutation on the new indices of Y
            bats::util::apply_perm_swap(ind_Y_i_sorted, bats::util::inv_perm(p));
            // finally, get the permutation for filtred filtration
            permutations[i] = bats::util::sortperm(ind_Y_i_sorted);

            // step 5, find the boundary indices
            if(i>0){
                for(auto& ind_list: boundary_indices[i]){
                    // update each index of the boundary of a simplex in dim i
                    for(size_t j = 0; j < ind_list.size(); j++){
                        ind_list[j] = perms_Y_inv[i-1][ind_list[j]];
                    }
                    std::sort(ind_list.begin(), ind_list.end()); // sort indices for addition of columns
                }
            }

            // step 6, for the convenience of adding and deleting,
            // we sort their indices
            
            // sort deletion indices by ascending order
            std::sort(deletion_indices[i].begin(), deletion_indices[i].end());

            // sort addition indices by ascending order
            // and sort the boundary indices correspondingly
            if(i ==0 ){ //no boundary indices in dim = 0
                std::sort(addition_indices[i].begin(), addition_indices[i].end());
            }
            else{ // sort addition indices and its boundary indices simutaneously
                auto sv = SparseVector(addition_indices[i], boundary_indices[i]);
                sv.sort();
                for (size_t j = 0; j < sv.nnz(); j++)
                {
                    addition_indices[i][j] = (sv.nzbegin()+j)->ind;
                    boundary_indices[i][j] = (sv.nzbegin()+j)->val;
                }
            }

            // p = sort_indexes(addition_indices[i]);
            // bats::util::apply_perm_swap(addition_indices[i], bats::util::inv_perm(p));
            // if(i>0)
            // bats::util::apply_perm_swap(boundary_indices[i], bats::util::inv_perm(p));
            
        }         
    }
    
    // return the permutation used to move deletion simplices in dimension i
    // to the end of the list of simplices
    std::vector<size_t> permutation_deletion_end(size_t i){

        if(!deletion_indices[i].empty()){
            std::vector<size_t> perm = perm_to_the_end(deletion_indices[i], F_X_vals[i].size());
            return perm;
        }else
        {
            return identity_perm(F_X_vals[i].size());
        }
        
    }

    // print summary of updating information
    void print_summary(){
        for(size_t i = 0 ; i <= max_dim; i++){
            std::cout << "\nFor dimension "<< i << std::endl;
            std::cout << "Deletion information:" << std::endl;
            print_1D_vectors(deletion_indices[i]);
            std::cout << "Addition information:" << std::endl;
            print_1D_vectors(addition_indices[i]);
            if(i>0){    
                std::cout << ",and the boundary indices are" << std::endl;
                print_2D_vectors(boundary_indices[i]);
            }
            std::cout << "Permutation information:" << std::endl;
            print_1D_vectors(permutations[i]);
        }
    }

    void print_detail(){
        auto smcplex_X = F_old.complex();
        auto smcplex_Y = F_new.complex();

        for(size_t i = 0 ; i <= max_dim; i++){
            // find simplices in dimension i for both X and Y
            auto simplices_Y = smcplex_Y.get_simplices(i); // simplices in dimension i
            auto simplices_X = smcplex_X.get_simplices(i);
            if(filtered_boolean){
                bats::util::apply_perm_swap(simplices_X, perms_X_inv[i]);
                bats::util::apply_perm_swap(simplices_Y, perms_Y_inv[i]);
            }
            std::cout << "\nFor dimension "<< i << std::endl;
            std::cout << "Deletion information:" << std::endl;
            for(size_t j = 0; j < deletion_indices[i].size(); j++){
                std::cout << "\t";
                print_simplex(simplices_X[deletion_indices[i][j]]);
                std::cout << "with index "<< deletion_indices[i][j] << std::endl;
            }
            
            std::cout << "Addition information:" << std::endl;
            for(size_t j = 0; j < addition_indices[i].size(); j++){
                std::cout << "\t";
                print_simplex(simplices_Y[addition_indices[i][j]]);
                std::cout << "with index "<< addition_indices[i][j] << std::endl;

                if(i>0){
                    std::cout << "\twith boundary indices: ";
                    for(auto& ind: boundary_indices[i][j]){
                        std::cout << ind <<" ";
                    }
                }
                std::cout << "\n";
            }

            // if(i>0){    
            //     std::cout << ",and the boundary indices are" << std::endl;
            //     print_2D_vectors(boundary_indices[i]);
            // }
            std::cout << "Intersection Information" << std::endl;
            std::cout << "\t(simplex, index in old and new filtration)" << std::endl;
            for(size_t j = 0; j < intersection_indices_X[i].size(); j++){
                std::cout << "\t";
                print_simplex(simplices_X[intersection_indices_X[i][j]]);
                std::cout << " "<< intersection_indices_X[i][j]; 
                std::cout << " "<< intersection_indices_Y[i][j];
                std::cout << "\n";
            }

            std::cout << "Transformed to Permutation:";
            print_1D_vectors(permutations[i]);
        }
    }
};

} // namespace bats