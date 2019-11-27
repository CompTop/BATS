#pragma once

#include <linalg/naive_dense.h>
#include <linalg/field.h>

using std::cout;

// type A quiver
template< typename F>
struct Type_A{
    // matrix types
    using DI = Dense<F,ColMaj>;
    using AD = A<DI>;
    
    //input data
    std::vector<AD> mats;
    std::vector<bool> arrow_dir;
    
    // number of edges
    size_t n;
    // dimension of each vector space
    std::vector<size_t> dims;
    
    // reducedm pivot matrices
    std::vector<EL<DI>> ELmats;
    std::vector<ELH<DI>> ELHmats;

    // factors from factorization
    std::vector<L<DI>> Lmats;
    std::vector<U<DI>> Umats;
    std::vector<P<DI>> Pmats;

    // basis and inverse basis matrices
    std::vector<AD> basis;
    std::vector<AD> inv_basis;
    
	bool fwd_sweep_complete = 0;
	bool bwd_sweep_complete = 0;

    /* constructor 
    *   ms - Vector of matrices
    *   as - Vector of bools indication arrow directions
    */
    Type_A(std::vector<AD>& ms,std::vector<bool>& as ) :mats(ms), arrow_dir(as) {
        n = ms.size();
        dims.resize(n+1);
        get_dims();
        
        // initialize basis matrices
        AD b1,b2;
        for(size_t i=0;i<n+1;i++){
            b1 = A<DI>::identity(dims[i]);
            basis.emplace_back(b1);
            b2 = A<DI>::identity(dims[i]);
            inv_basis.emplace_back(b2);
        }
    }
    
    
    /*
    *  Extract dimesnions for vector spaces from matrices
    */
    void get_dims(){
        for(size_t i=0;i<n;i++){
            if(arrow_dir[i]==0)
                dims[i] = mats[i].nrow();
            else
                dims[i] = mats[i].ncol();
        }
        if(arrow_dir[n-1]==0)
            dims[n] = mats[n-1].ncol();
        else
            dims[n] = mats[n-1].nrow();
    }
    
    /*
    *  Performs the Forward sweep of the quiver algorithm
    */
    void forward_sweep(){
        // variables to hold factors
        L<DI> Lmat;
        EL<DI> ELmat;
        P<DI> Pmat;
        U<DI> Umat;
        ELH<DI> ELHmat;
        for(size_t i=0;i<n;i++){
            if(arrow_dir[i]==0){
                // LEUP Factorization
                std::tie(Lmat, ELmat, Umat, Pmat) = LEUP_fact(mats[i]);

                auto Pinv = Pmat.Trp(); // inverse of transpose
                if(i<n-1){
                    if(arrow_dir[i+1]==0){
                        // pass P and U
                        apply_matmul_on_left(Pmat,mats[i+1]);
                        apply_matmul_on_left(Umat,mats[i+1]);

                    }else{
                        apply_matmul_on_right(mats[i+1],Pinv);
                        apply_inverse_on_right(mats[i+1],Umat);
                    }
                }
                //inv basis change
                apply_matmul_on_left(Pmat,inv_basis[i+1]);
                apply_matmul_on_left(Umat,inv_basis[i+1]);
                //basis change
                apply_matmul_on_right(basis[i+1],Pinv);
                apply_inverse_on_right(basis[i+1],Umat);

                // record all matrices
                ELmats.emplace_back(ELmat);
                ELHmats.emplace_back(ELH<DI>());
                Lmats.emplace_back(Lmat);
                Umats.emplace_back(Umat);
                Pmats.emplace_back(Pmat);
            }
            else{
                // PUEL Factorization
                std::tie(Pmat, Umat, ELHmat, Lmat) = PUEL_fact(mats[i]);
                auto Pinv = Pmat.Trp();
                if(i<n-1){
                    if(arrow_dir[i+1]==0){
                        // pass Pinv and Uinv
                        apply_matmul_on_left(Pinv,mats[i+1]);
                        apply_inverse_on_left(Umat,mats[i+1]);
                    }else{
                        // pass P and U
                        apply_matmul_on_right(mats[i+1],Pmat);
                        apply_matmul_on_right(mats[i+1],Umat);
                    }
                }
                //inv basis change
                apply_matmul_on_left(Pinv,inv_basis[i+1]);
                apply_inverse_on_left(Umat,inv_basis[i+1]);
                //basis change
                apply_matmul_on_right(basis[i+1],Pmat);
                apply_matmul_on_right(basis[i+1],Umat);

                // record all matrices
                ELmats.emplace_back(EL<DI>());
                ELHmats.emplace_back(ELHmat);
                Lmats.emplace_back(Lmat);
                Umats.emplace_back(Umat);
                Pmats.emplace_back(Pmat);
            }
        }
		fwd_sweep_complete = 1;
    }
    
    /*
    *  Computes inverse of a lower triangular matrix.
    *  Allocates new memory.
    */
    auto L_inverse(L<DI> La){
        auto Linv = L<DI>(La.m,La.n);
        make_diag_ones<F>(Linv);
        apply_inverse_on_left(La,Linv);
        return Linv;
    }
    
    
    /*
    *  Performs the Backward sweep of the quiver algorithm
    */
    void backward_sweep(){
		if( !fwd_sweep_complete ){
			throw std::runtime_error("Forward sweep not completed");
		}
        L<DI> Lt; // The L that gets passed
        for(size_t i=n-1; i<n; i--){ // i<n instead of i>=0 for unsignedint
            if(arrow_dir[i]==0){
                if (i==n-1){ // check if first step, then no incoming L
                    Lt = Lmats[i];
                }else{ // if not, commute incoming with  E and multiply with existing L
                    auto Lt2 = commute(ELmats[i],Lt); //allocates new mem 
                    Lt.free();
                    Lt = Lt2;
                    apply_matmul_on_left(Lmats[i],Lt);
                }
                // apply the passing L to basis mats
                apply_inverse_on_left(Lt,inv_basis[i]);
                apply_matmul_on_right(basis[i],Lt);  
                // If next arrow direction is opposite inverse the passing L
                if( i>0 && arrow_dir[i-1]!= 0 ){
                    auto Lt3 = L_inverse(Lt); //allocates new mem 
                    Lt.free();
                    Lt = Lt3;
                }

            }else{
                if (i==n-1){
                    Lt = Lmats[i];
                }else{ // commute the other way
                    auto Lt2 = commute(Lt,ELHmats[i]); //allocates new mem 
                    Lt.free();
                    Lt = Lt2;
                    apply_matmul_on_right(Lt,Lmats[i]);
                }

                // apply the passing L to basis mats
                apply_matmul_on_left(Lt,inv_basis[i]);
                apply_inverse_on_right(basis[i],Lt);  

                // If next arrow direction is opposite inverse the passing L
                if( i>0 && arrow_dir[i-1]!= 1 ){
                    auto Lt3 = L_inverse(Lt); //allocates new mem 
                    Lt.free();
                    Lt = Lt3;
                }

            }
        }
        bwd_sweep_complete = 1;
    }

	// copy for consistency check
	std::vector<AD> mats_copy;
	bool is_copy_saved = 0;
	/*
	* Create copy for consistency check
	*/
	void create_copy_of_mats(){
		AD a;
		for(size_t i=0;i<n;i++){
			a = mats[i].copy();
			mats_copy.emplace_back(a);
		}
		is_copy_saved = 1;
	}


	/*
	* Check consistency of factorization
	*/
	bool is_consistent(){
		if( !is_copy_saved )
			throw std::runtime_error("Copy of mats not saved");
		if( ! bwd_sweep_complete || !fwd_sweep_complete )
			throw std::runtime_error("Factorization not complete");
		std::vector<AD> mats_recon2;
		for(size_t i=0;i<n;i++){
			if(arrow_dir[i]==0){
				auto mtemp = matmul(matmul(
						basis[i],ELmats[i]),inv_basis[i+1]
					);
				mats_recon2.emplace_back(mtemp);
			}else{
				auto mtemp = matmul(matmul(
						basis[i+1],ELHmats[i]),inv_basis[i]
					);
				mats_recon2.emplace_back(mtemp);
			}
		}
		bool consistent = 1;
		for(size_t i=0;i<n;i++){
			consistent &=  (mats_copy[i]==mats_recon2[i]);
		}
		return consistent;
	}
};
