

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <linalg/naive_dense.h>
#include <linalg/field.h>
#include <quiver/type_a.h>

#define F3 ModP<int, 3>

TEST_CASE_TEMPLATE("Matrix Multiplication", T, int, ModP<int, 2>, ModP<int,3>, ModP<int, 5>, Rational<int>) {

	using AD = A<Dense<T,ColMaj>>;

	T a1v[] = {
	2,3,4,
	1,2,3,
	8,5,2,
	};
	AD a1(3,3,a1v);

	T a2v[] = {
	39,32,25,
	28,22,16,
	37,44,51,
	};
	AD a2(3,3,a2v);

	CHECK( (matmul(a1,a1) == a2) );
}


TEST_CASE_TEMPLATE("VectorView, Transpose and JConjugate" ,
	Acc, ColMaj, RowMaj, RevColMaj, RevRowMaj 
) {
	using F = F3;
	using ADAcc = A<Dense<F,Acc>>;
	using AD = A<Dense<F,ColMaj>>;

	srand(1);

	ADAcc a1(5,3);
	fill_rand(a1);
	AD a2(5,5);
	fill_zeros(a2);
	AD a3(3,3);
	fill_zeros(a3);


	a2[2] = a1[1];
	a3[1] = a1.r(4);


	//------Compare----------
	F v1[] = {
	0,0,0,0,0,
	0,0,0,0,0,
	1,2,0,2,1,
	0,0,0,0,0,
	0,0,0,0,0,
	};
	AD cmat1(5,5,v1);
	CHECK( (a2==cmat1) );
	//-----------------------
	//------Compare----------
	F v2[] = {
	0,0,0,
	2,1,2,
	0,0,0,
	};
	AD cmat2(3,3,v2);
	CHECK( (a3==cmat2) );
	//-----------------------


	auto a1t = a1.Trp();
	a3[1] = a1t[2];
	a2[2] = a1t.r(1);

	//------Compare----------
	F v3[] = {
	1,1,0,
	1,2,1,
	1,0,0,
	1,2,1,
	2,1,2,
	};
	AD cmat3(3,5,v3);
	CHECK( (a1t==cmat3) );
	//-----------------------
	//------Compare----------
	F v4[] = {
	0,0,0,0,0,
	0,0,0,0,0,
	1,2,0,2,1,
	0,0,0,0,0,
	0,0,0,0,0,
	};
	AD cmat4(5,5,v4);
	CHECK( (a2==cmat4) );
	//-----------------------
	//------Compare----------
	F v5[] = {
	0,0,0,
	1,0,0,
	0,0,0,
	};
	AD cmat5(3,3,v5);
	CHECK( (a3==cmat5) );
	//-----------------------

	auto a1j = a1.JConj();
	a2[2] = a1j[2];
	a3[1] = a1j.r(2);

	//------Compare----------
	F v6[] = {
	2,1,0,1,0,
	1,2,0,2,1,
	2,1,1,1,1,
	};
	AD cmat6(5,3,v6);
	CHECK( (a1j==cmat6) );
	//-----------------------
	//------Compare----------
	F v7[] = {
	0,0,0,0,0,
	0,0,0,0,0,
	2,1,1,1,1,
	0,0,0,0,0,
	0,0,0,0,0,
	};
	AD cmat7(5,5,v7);
	CHECK( (a2==cmat7) );
	//-----------------------
	//------Compare----------
	F v8[] = {
	0,0,0,
	0,0,1,
	0,0,0,
	};
	AD cmat8(3,3,v8);
	CHECK( (a3==cmat8) );
	//-----------------------

	auto a1tj = a1.TJConj();
	a3[1] = a1tj[3];
	a2[2] = a1tj.r(1);


	//------Compare----------
	F v9[] = {
	2,1,2,
	1,2,1,
	0,0,1,
	1,2,1,
	0,1,1,
	};
	AD cmat9(3,5,v9);
	CHECK( (a1tj==cmat9) );
	//-----------------------
	//------Compare----------
	F v10[] = {
	0,0,0,0,0,
	0,0,0,0,0,
	1,2,0,2,1,
	0,0,0,0,0,
	0,0,0,0,0,
	};
	AD cmat10(5,5,v10);
	CHECK( (a2==cmat10) );
	//-----------------------
	//------Compare----------
	F v11[] = {
	0,0,0,
	1,2,1,
	0,0,0,
	};
	AD cmat11(3,3,v11);
	CHECK( (a3==cmat11) );
	//-----------------------

}




TEST_CASE("el_commute") {

	using DI = Dense<int,ColMaj>;
	int elmat1v[] = {1,0,0,0,
		          0,1,0,0,
		          0,0,0,1,
		          0,0,0,0};
	EL<DI> elmat1(4,4,elmat1v);

	int lmat1v[] = {2,2,3,4,
		          0,3,2,4,
		          0,0,4,4,
		          0,0,0,5};
	L<DI> lmat1(4,4,lmat1v);

	auto res = commute(elmat1,lmat1);

	CHECK( (matmul(res,elmat1)== matmul(elmat1,lmat1)) );
}

TEST_CASE_TEMPLATE("LEUP Factorization -self consistency",
	Acc, ColMaj, RowMaj, RevColMaj, RevRowMaj 
) {

	SUBCASE(""){ srand(0); }
	SUBCASE(""){ srand(1); }
	SUBCASE(""){ srand(2); }
	SUBCASE(""){ srand(3); }
	SUBCASE(""){ srand(4); }
	SUBCASE(""){ srand(5); }

	using DI = Dense<F3,Acc>;

	A<DI> mat3(4,5);
	A<DI> mat3copy;
	L<DI> Lmat;
	EL<DI> ELmat;
	P<DI> Pmat;
	U<DI> Umat;
	A<DI> leup;

	fill_rand(mat3);

	//reduce rank
	mat3.r(2)=mat3.r(3);
	mat3[1]=mat3[3];

	mat3copy = mat3.copy();

	std::tie(Lmat, ELmat, Umat, Pmat) = LEUP_fact(mat3);

	// check no change in mat3
	CHECK( (mat3==mat3copy) );

	leup = matmul(matmul(Lmat,ELmat),matmul(Umat,Pmat));

	CHECK( (leup==mat3) );

}

TEST_CASE_TEMPLATE("l_solve", F, ModP<int, 2>, ModP<int,3>, ModP<int, 5>, Rational<int>) {

	using DI = Dense<F,ColMaj>;
	
	auto I = A<DI>::identity(4);
	F lmat1v[] = {
	1,1,0,0,
	0,1,2,2,
	0,0,1,1,
	0,0,0,1,
	};
	L<DI> Lmat(4,4,lmat1v);
	auto tempI = I.copy();
	auto Linv = apply_inverse_on_left(Lmat,tempI);

	CHECK( (matmul(Lmat,Linv) == I) );
	CHECK( (matmul(Linv,Lmat) == I) );

}


TEST_CASE_TEMPLATE("Quiver Alg - self consistency", F, ModP<int,2>, ModP<int,3>, ModP<int, 5>) {

	using DI = Dense<F,ColMaj>;
	using AD = A<DI>;


	size_t n=10;
	std::vector<AD> mats;
	std::vector<bool> arrow_dir;
	std::vector<size_t> dims;
	
	// generate random type-A quiver with random dimensions
	AD a1;
	size_t dim1 = 1+rand()%10;
	dims.emplace_back(dim1);
	for(size_t i=0;i<n;i++){
		size_t dim_ip1 = 1+rand()%10;
		dims.emplace_back(dim_ip1);
		
		arrow_dir.emplace_back(rand()%2);
		
		if(arrow_dir[i]==0)
		    a1 = AD(dims[i],dims[i+1]);
		else
		    a1 = AD(dims[i+1],dims[i]);
		fill_rand<F>(a1);
		mats.emplace_back(a1);
//a1.print();
	}

	// create copy to compare later
	std::vector<AD> mats_copy;
	AD a3;
	for(size_t i=0;i<n;i++){
		a3 = mats[i].copy();
		mats_copy.emplace_back(a3);
	}

	// create object
	auto taq = Type_A<F>(mats,arrow_dir);
	// forward sweep
	taq.forward_sweep();

	// check consistency of forward sweep
	std::vector<AD> mats_recon;
	for(size_t i=0;i<taq.n;i++){
		if(taq.arrow_dir[i]==0){
		    auto mtemp = matmul(matmul(
							taq.basis[i],
							matmul(taq.Lmats[i],taq.ELmats[i])),
							taq.inv_basis[i+1]
						 );
		    mats_recon.emplace_back(mtemp);
		}else{
		    auto mtemp = matmul(matmul(
							taq.basis[i+1],
							matmul(taq.ELHmats[i],taq.Lmats[i])),
							taq.inv_basis[i]
						 );
		    mats_recon.emplace_back(mtemp);
		}
	}
	for(size_t i=0;i<n;i++){
		CHECK( (mats_copy[i]==mats_recon[i]) );
	}

	//Backward sweep
	taq.backward_sweep();

	// check quiver factorization consistency
	std::vector<AD> mats_recon2;
	for(size_t i=0;i<taq.n;i++){
		if(taq.arrow_dir[i]==0){
		    auto mtemp = matmul(matmul(
						taq.basis[i],taq.ELmats[i]),taq.inv_basis[i+1]
					);
		    mats_recon2.emplace_back(mtemp);
		}else{
		    auto mtemp = matmul(matmul(
						taq.basis[i+1],taq.ELHmats[i]),taq.inv_basis[i]
					);
		    mats_recon2.emplace_back(mtemp);
		}
	}

	for(size_t i=0;i<n;i++){
    	CHECK( (mats_copy[i]==mats_recon2[i]) );
	}
}
