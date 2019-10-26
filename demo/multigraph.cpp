#include <multigraph/multigraph.h>
#include <iostream>

#define NT Node<int>
#define ET Edge<int>
int main() {

    int a = 5;
    NT na(a);
    na.print();

    int b = 6;
    NT nb(b);
    nb.print();

    int c = 4;
    ET ec(na, na, c);
    ec.print();

    int d = 1;
    ET ed(na, nb, d);
    ed.print();

    na.print();
    nb.print();

    return 0;
}
