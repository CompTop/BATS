#include <linalg/field.h>
#include <linalg/naive_dense.h>

#define F ModP<int, 3>

int main() {

    A<Dense<F>> Amat(2,2);
    Amat.print();



    return 0;
}
