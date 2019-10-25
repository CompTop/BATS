
#pragma once

#include "matrix_interface.h"
#include <tuple>
#include <assert.h>

struct SI{};

#define MM(A,B,C) \
C<SI> matmul(A<SI> a, B<SI> b){\
    C<SI> c;\
	(void) a;\
	(void) b;\
    return c;\
};


MM(A,A,A)
MM(D,D,D)
MM(L,L,L)
MM(U,U,U)
MM(L,U,A)

    
std::tuple<L<SI> , EL<SI> , U<SI> , P<SI>> lufact(A<SI> a){
    L<SI> l;
    EL<SI> el;
    U<SI> u;
    P<SI> p;
	(void)a;
    return std::make_tuple( l,el,u,p);
}



std::tuple< L<SI>,EL<SI>  > commute( L<SI> l1 ,EL<SI> el1, L<SI> l2){
    L<SI> l;
    EL<SI> el;
	(void) l1;
	(void) el1;
	(void) l2;
    return std::make_tuple( l,el); 
}
