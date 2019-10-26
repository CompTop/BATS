
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

    
std::tuple<L<SI> , EL<SI> , U<SI> , P<SI>> LEUP_fact(A<SI> a){
    L<SI> l;
    EL<SI> el;
    U<SI> u;
    P<SI> p;
	(void)a;
    return std::make_tuple( l,el,u,p);
}



L<SI> el_commute( EL<SI> el1, L<SI> l1 ){
    L<SI> l;
	(void) l1;
	(void) el1;
    return l;
}
