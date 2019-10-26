

#include <iostream>
#include <linalg/symbolic_implementation.h>
#include <quiver/quiver_a3.h>
#include <linalg/field.h>
#include <linalg/naive_dense.h>
#include <assert.h>

#define F ModP<int, 3>

int main() {
	using std::cout;

	A<SI> a1,a2;

	Quiver_A3<SI> q(a1,a2);
	auto ret = q.factorize();
		
	cout<<"Symbolic Test: Asserting that the types are correct...\n";
	bool check = (typeid(ret) == typeid(std::tuple<EL<SI>, EL<SI> >));
	assert(check);
	cout<<"Correct types = "<<check<<"\n";

	using DI = Dense<F>;

	A<DI> b1(5,5),b2(5,5);
	EL<DI> E1,E2;
	fill_rand(b1);
	fill_rand(b2);

	Quiver_A3<DI> q2(b1,b2);
	std::tie(E1,E2) = q2.factorize();
	
	cout<<"E1 = "<<"\n";
	E1.print();
	cout<<"\n";

	cout<<"E2 = "<<"\n";
	E2.print();
	cout<<"\n";

}
