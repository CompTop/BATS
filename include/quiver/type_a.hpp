#pragma once

#include <linalg/naive_dense.hpp>
#include <linalg/field.hpp>
#include <persistence/barcode.hpp>
#include <string>

using std::cout;
using std::string;


namespace bats {


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

	// put barcode into vector of Persistence Pairs
	std::vector<PersistencePair<size_t>> barcode_pairs(size_t hdim=0) {

		std::vector<PersistencePair<size_t>> pairs;
		std::vector<int> active, active2;

		A<DI> Em;
		for(size_t i=0;i<dims[0];i++){
		    pairs.emplace_back(PersistencePair(hdim, size_t(0), size_t(0), size_t(0), size_t(0)));
		    active.emplace_back(i);
		}
		// loop over edges
		for(size_t kk=0; kk<ELmats.size(); kk++){
		    if(arrow_dir[kk]==0){
		        Em=ELmats[kk];
		        active2.clear();
		        active2.resize(Em.ncol());
		        for(size_t j=0; j<Em.ncol(); j++){
		            size_t i;
		            for(i=0; i<Em.nrow(); i++){
						// look for bar to continue
		                if(Em(i,j)!=F(0)){
		                    pairs[active[i]].death_ind = kk + 1;
							pairs[active[i]].death = kk + 1;

		                    active2[j]=active[i];
		                    break;
		                }
		            }

		            if(i==Em.nrow()){
						// we are starting a new bar
		                active2[j]=pairs.size();
		                pairs.emplace_back(PersistencePair(hdim, kk+1, kk+1, kk+1, kk+1));
		            }
		        }
		        std::swap(active,active2);
		    }else{
		        Em=ELHmats[kk];
		        active2.clear();
		        active2.resize(Em.nrow());
		        for(int i=Em.nrow()-1; i>=0; i--){
		            int j;
		            for(j=Em.ncol()-1; j>=0; j--){
		                if(Em(i,j)!=F(0)){
							// look for bar to continue
							pairs[active[i]].death_ind = kk + 1;
							pairs[active[i]].death = kk + 1;

		                    active2[i]=active[j];
		                    break;
		                }
		            }
		            if(j==-1){
						// we are starting a new bar
		                active2[i]=pairs.size();
		                pairs.emplace_back(PersistencePair(hdim, kk+1, kk+1, kk+1, kk+1));
		            }
		        }
		        std::swap(active,active2);
		    }
		}

		return pairs;
	}


	/*
	*   print barcodes
	*		- Each line has eaxactly one bar
	*/
	void print_barcodes(){
		std::vector<string> lines;
		std::vector<int> active,active2;
		A<DI> Em;
		for(size_t i=0;i<dims[0];i++){
		    lines.emplace_back("*");
		    active.emplace_back(i);
		}
		for(size_t kk=0;kk<ELmats.size();kk++){
		    if(arrow_dir[kk]==0){
		        Em=ELmats[kk];
		        active2.clear();
		        active2.resize(Em.ncol());
		        for(size_t j=0; j<Em.ncol(); j++){
		            size_t i;
		            for(i=0; i<Em.nrow(); i++){
		                if(Em(i,j)==F(1)){
		                    lines[active[i]]+="<--*";
		                    active2[j]=active[i];
		                    break;
		                }
		            }
		            if(i==Em.nrow()){
		                string plu ="";
		                for(size_t l=0;l<kk+1;l++)
		                    plu+="    ";
		                plu+="*";
		                active2[j]=lines.size();
		                lines.emplace_back(plu);
		            }
		        }
		        std::swap(active,active2);
		    }else{
		        Em=ELHmats[kk];
		        active2.clear();
		        active2.resize(Em.nrow());
		        for(int i=Em.nrow()-1; i>=0; i--){
		            int j;
		            for(j=Em.ncol()-1; j>=0; j--){
		                if(Em(i,j)==F(1)){
		                    lines[active[j]]+="-->*";
		                    active2[i]=active[j];
		                    break;
		                }
		            }
		            if(j==-1){
		                string plu ="";
		                for(size_t l=0;l<kk+1;l++)
		                    plu+="    ";
		                plu+="*";
		                active2[i]=lines.size();
		                lines.emplace_back(plu);
		            }
		        }
		        std::swap(active,active2);
		    }
		}
		for(size_t i=0;i<lines.size();i++)
		cout<<lines[i]<<"\n";
	}

	/*
	*   print barcodes in compressed format,
	*		- Each line can have multiple bars
	*/
	void print_barcodes_compressed(){
		std::vector<string> lines;
		std::vector<size_t> active,active2;
		A<DI> Em ;
		std::vector<int> lastup;

		for(size_t i=0;i<dims[0];i++){
			lines.emplace_back("*");
			active.emplace_back(i);
			lastup.emplace_back(0);
		}
		for(size_t kk=0;kk<ELmats.size();kk++){
			if(arrow_dir[kk]==0){
				Em=ELmats[kk];
				active2.clear();
				active2.resize(Em.ncol());
				for(size_t j=0; j<Em.ncol(); j++){
				    size_t i;
				    for(i=0; i<Em.nrow(); i++){
				        if(Em(i,j)==F(1)){
				            lines[active[i]]+="<--*";
				            lastup[active[i]]+=1;
				            active2[j]=active[i];
				            break;
				        }
				    }
				    if(i==Em.nrow()){
				        size_t l=0;
				        for(;l<lines.size();l++){
				            size_t jj=0;
				            bool bad=0;
				            for(;jj<active2.size() ;jj++)
				                if(active2[jj]==l ) bad=1;
				            if(!bad) break;
				        }
				        active2[j]=l;

				        if(l==lines.size()){
				            lines.emplace_back(" ");
				            lastup.emplace_back(0);
				        }

				        string plu = "";
				        for(size_t ll=lastup[l];ll<kk;ll++)
				            plu+="    ";
				        plu+="   *";

				        lines[l] += plu;
				        lastup[l] = kk+1;

				    }
				}
				std::swap(active,active2);
			}else{
				Em=ELHmats[kk];
				active2.clear();
				active2.resize(Em.nrow());
				for(int i=Em.nrow()-1; i>=0; i--){
				    int j;
				    for(j=Em.ncol()-1; j>=0; j--){
				        if(Em(i,j)==F(1)){
				            lines[active[j]]+="-->*";
				            lastup[active[j]]+=1;
				            active2[i]=active[j];
				            break;
				        }
				    }
				    if(j==-1){
				        size_t l=0;
				        for(;l<lines.size();l++){
				            size_t ii=0;
				            bool bad=0;
				            for(;ii<active2.size() ;ii++)
				                if(active2[ii]==l ) bad=1;
				            if(!bad) break;
				        }
				        active2[i]=l;

				        if(l==lines.size()){
				            lines.emplace_back(" ");
				            lastup.emplace_back(0);
				        }

				        string plu = "";
				        for(size_t ll=lastup[l];ll<kk;ll++)
				            plu+="    ";
				        plu+="   *";

				        lines[l] += plu;
				        lastup[l] = kk+1;
				    }
				}
				std::swap(active,active2);
			}
		}
		for(size_t i=0;i<lines.size();i++)
			cout<<lines[i]<<"\n";

	}
};

} // namespace bats
