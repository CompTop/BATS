#include <linalg/field.h>
#include <linalg/sparse_vector.h>
#include <linalg/set_vector.h>
#include <vector>
#include <chrono>

//#define F int
#define F ModP<int, 3>
//#define F Rational<int>

#define NITER 100000

template <typename T1, typename T2>
void time_axpy(T1 &x, T2 &y, size_t N) {
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < N; i++) {
        x.axpy( 1, y);
        x.axpy(-1, y);
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "microseconds: " << duration.count() << std::endl;
}

int main() {

    std::vector<size_t> ind;
    std::vector<int> val_int;
    std::vector<F> val;

    val_int = {-1,1,-1};
    // val = val_int;

    ind = {1,3,5};
    val = {-1,1,-1};

    // SetVector x(ind, val);
    SetVector<F, size_t> x(ind.cbegin(), val_int.cbegin(), ind.size());
    x.print_row();

    ind = {2, 4};
    val = {1, 1};

    SetVector y(ind, val);
    y.print_row();

    time_axpy(x, y, NITER);

    x.print_row();

    ind = {1,3,5};
    val = {-1,1,-1};
    // construct using iterator
    SparseVector<F, size_t> x2(ind.cbegin(), val_int.cbegin(), ind.size());

    x2.print_row();

    ind = {2, 4};
    val = {1, 1};

    SparseVector y2(ind, val);
    y2.print_row();

    time_axpy(x2, y2, NITER);

    x2.print_row();

    std::cout <<"SetVector, SetVector" << std::endl;
    time_axpy(x, y, NITER);
    time_axpy(x, y, NITER);
        std::cout <<"SparseVector, SparseVector" << std::endl;
    time_axpy(x2, y2, NITER);
    time_axpy(x2, y2, NITER);
        std::cout <<"SetVector, SparseVector" << std::endl;
    time_axpy(x, y2, NITER);
    time_axpy(x, y2, NITER);
        std::cout <<"SparseVector, SetVector" << std::endl;
    time_axpy(x2, y, NITER);
    time_axpy(x2, y, NITER);

    ind = {1,3};
    val = {-1,1};
    SparseVector<F, size_t> a(ind, val);
    ind = {0,3};
    val = {-1,1};
    SparseVector<F, size_t> b(ind, val);
    a.print_row();
    b.print_row();
    a.axpy(-1, b);
    a.print_row();

    // F v(2);
    // std::cout << v << std::endl;
    // std::cout << F(-1) * v << std::endl;
    //
    // nzpair p1(0, F(-1));
    // nzpair p2(p1);
    // std::cout <<p1 << p2 << std::endl;

    return 0;
}
