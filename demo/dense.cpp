#include <linalg/field.h>
#include <linalg/naive_dense.h>
#include <assert.h>

#define F ModP<int, 3>

int main() {
	
	using std::cout;

   	A<Dense<F>> mat(5,5);
	A<Dense<F>> matcopy;
	
	L<Dense<F>> Lmat;
	EL<Dense<F>> ELmat;
	P<Dense<F>> Pmat;
	U<Dense<F>> Umat;
	
	A<Dense<F>> leup;


	fill_rand(mat);
	cout<<"A = "<<"\n";
	mat.print();
	cout<<"\n";

	// reduce rank
	mat.r(1)=mat.r(0);
	mat[1]=mat[3];

	std::tie(Lmat, ELmat, Umat, Pmat) = LEUP_fact(mat);

	leup = matmul(matmul(Lmat,ELmat),matmul(Umat,Pmat));

	cout<<"L = "<<"\n";
	Lmat.print();
	cout<<"\n";

	cout<<"EL = "<<"\n";
	ELmat.print();
	cout<<"\n";

	cout<<"U = "<<"\n";
	Umat.print();
	cout<<"\n";

	cout<<"P = "<<"\n";
	Pmat.print();
	cout<<"\n";

	cout<<"LEUP = "<<"\n";
	mat.print();
	cout<<"\n";

	cout<<"(A == LEUP) = "<<(mat==leup)<<"\n";
	assert(mat==leup);

	cout<<"----- EL L = L EL commutation -----\n";

	auto L_left = el_commute(ELmat,Lmat) ;

	cout<<"Commute check = ";
    cout<<( matmul(L_left,ELmat) == matmul(ELmat,Lmat) )<<"\n";
	assert( matmul(L_left,ELmat) == matmul(ELmat,Lmat) );

	cout<<"----- L solve -----\n";
	A<Dense<F>> I(5,5);
	make_diag_ones(I);
	auto Linv = l_solve(Lmat,I);
	cout<<"L solve check = "<< (matmul(Linv,Lmat)==I) <<"\n\n";
	assert(matmul(Linv,Lmat)==I) ;

    return 0;
}
