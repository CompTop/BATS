#include <bats.hpp>
#include <vector>

using VT = SparseVector<int, size_t>;
using MatT = ColumnMatrix<VT>;

int main() {

    MatT A = MatT::identity(3);
    A.print();

    smith_normal_form(A);
    A.print();

    std::vector col{
        VT({0,1,2},{2, -6, 10}),
        VT({0,1,2},{4, 6, 4}),
        VT({0,1,2},{4, 12, 16})
    };
    A = MatT(3, 3, col);
    A.print();

    // let's eliminate first non-zero in column 1
    size_t i(0), j(1);
    auto [x, y, g, w, z] = extended_euclidean(A(i,i), A(i,j));
    std::cout << x << ',' << y << ',' << g << ',' << w << ',' << z << std::endl;
    std::cout << A(i,i) * x + A(i,j) * y << std::endl;
    A.mix_cols(i, j, x, -z, y, w);

    j = 2;
    std::tie(x, y, g, w, z) = extended_euclidean(A(i,i), A(i,j));
    A.mix_cols(i, j, x, -z, y, w);

    A.print();

    // now get rows
    i = 0;
    j = 1;
    std::tie(x, y, g, w, z) = extended_euclidean(A(i,i), A(j,i));
    std::cout << x << ',' << y << ',' << g << ',' << w << ',' << z << std::endl;
    A.mix_rows(i, j, x, y, -z, w);
    A.print();

    j = 2;
    std::tie(x, y, g, w, z) = extended_euclidean(A(i,i), A(j,i));
    std::cout << x << ',' << y << ',' << g << ',' << w << ',' << z << std::endl;
    A.mix_rows(i, j, x, y, -z, w);
    A.print();

    // {
    //     int a = 8;
    //     int b = 20;
    //     auto [x, y, g, ax, by] = extended_euclidean(a, b);
    //     std::cout << x << ',' << y << ',' << g << ',' << ax << ',' << by << std::endl;
    // }
    // smith_normal_form(A);
    // A.print();

    A = MatT(3, 3, col);
    A.print();

    smith_normal_form(A);
    A.print();

    A = MatT(3,3,col);
    A.print();
    auto F = smith_factorization(A);
    std::cout << "R:  "; F.R.print();
    std::cout << "S:  "; F.S.print();
    std::cout << "C:  "; F.C.print();
    F.prod().print();


    return 0;
}
