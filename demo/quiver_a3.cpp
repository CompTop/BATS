

#include<iostream>
#include "linalg/symbolic_implementation.h"
#include "quiver/quiver_a3.h"


int main() {

	A<SI> a1,a2;
	Quiver_A3<SI> q(a1,a2);
	auto ret = q.factorize();
		
	std::cout<<"Asserting types are same..."<<std::endl;
	assert(typeid(ret) == typeid(std::tuple<EL<SI>, EL<SI> >));
	std::cout<<"Success"<<std::endl;

}
