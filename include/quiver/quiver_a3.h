#include "../linalg/matrix_interface.h"
#include <tuple>

template <typename I>
class Quiver_A3{
    public:
    Quiver_A3(A<I> a1, A<I> a2) :A1(a1),A2(a2) {};
    A<I> A1, A2;
    auto factorize(){
        auto [l1, el1, u1, p1] = LEUP_fact(A1); //structured unbindings
        auto nA2 = matmul(matmul(u1,p1), A2);
        auto [l2, el2, u2, p2] = LEUP_fact(nA2);
        auto l3p = el_commute(el1,l2);
		auto l3 =  matmul(l1,l3p);
		//suppress warning
		(void) el2;
		(void) u2;
		(void) p2;
		(void) l3;
        return std::make_tuple(el2,el1);
    }
};
