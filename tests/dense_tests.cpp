

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <linalg/naive_dense.h>
#include <linalg/field.h>

#define F3 ModP<int, 3>

TEST_CASE_TEMPLATE("Matrix Multiplication", T, int, ModP<int, 2>, ModP<int,3>, ModP<int, 5>, Rational<int>) {

	T a1v[] = {
	2,3,4,
	1,2,3,
	8,5,2,
	};
	A<Dense<T>> a1(3,3,a1v);

	T a2v[] = {
	39,32,25,
	28,22,16,
	37,44,51,
	};
	A<Dense<T>> a2(3,3,a2v);

	CHECK( (matmul(a1,a1) == a2) );
}


TEST_CASE("ColView") {
	using AD = A<Dense<F3>>;
	srand(0);

	AD a1(4,4);
	AD a2(4,4);
	fill_rand(a2);
	a1[2]=a2[3];

	F3 resv[] = {
	0,0,0,0,
	0,0,0,0,
	1,0,1,1,
	0,0,0,0,
	};
	AD res(4,4,resv);

	CHECK( (a1==res) );
}


TEST_CASE("RowView") {
	using AD = A<Dense<F3>>;
	srand(0);

	AD a1(4,4);
	AD a2(4,4);
	fill_rand(a2);
	a1.r(2)=a2.r(3);

	F3 resv[] = {
	0,0,2,0,
	0,0,1,0,
	0,0,2,0,
	0,0,1,0,
	};
	AD res(4,4,resv);

	CHECK( (a1==res) );
}

TEST_CASE("el_commute") {
	int elmat1v[] = {1,0,0,0,
		          0,1,0,0,
		          0,0,0,1,
		          0,0,0,0};
	EL<Dense<int>> elmat1(4,4,elmat1v);

	int lmat1v[] = {2,2,3,4,
		          0,3,2,4,
		          0,0,4,4,
		          0,0,0,5};
	L<Dense<int>> lmat1(4,4,lmat1v);

	auto res = el_commute(elmat1,lmat1);

	CHECK( (matmul(res,elmat1)== matmul(elmat1,lmat1)) );
}

TEST_CASE("LEUP Factorization -self consistency") {

	A<Dense<F3>> mat3(4,5);
	A<Dense<F3>> mat3copy;
	L<Dense<F3>> Lmat;
	EL<Dense<F3>> ELmat;
	P<Dense<F3>> Pmat;
	U<Dense<F3>> Umat;
	A<Dense<F3>> leup;

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
	// A<Dense<F>> I(4,4);
	// make_diag_ones(I);
	auto I = A<Dense<F>>::identity(4);
	F lmat1v[] = {
	1,1,0,0,
	0,1,2,2,
	0,0,1,1,
	0,0,0,1,
	};
	L<Dense<F>> Lmat(4,4,lmat1v);

	auto Linv = l_solve(Lmat,I);

	CHECK( (matmul(Lmat,Linv) == I) );
	CHECK( (matmul(Linv,Lmat) == I) );

}
