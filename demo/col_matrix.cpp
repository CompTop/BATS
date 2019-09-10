#include <linalg/field.h>
#include <linalg/set_vector.h>
#include <linalg/sparse_vector.h>
#include <linalg/col_matrix.h>
#include <vector>
#include <chrono>

//#define F int
#define F ModP<int, 3>
//#define F Rational<int>

#define NITER 10000

// template over column type
template <class CT>
void time_solve(size_t N) {
    // construct upper triangular matrix
    using MatT = ColumnMatrix<CT>;
    MatT I = MatT::identity(5);

    // construct vector to solve with
    std::vector<size_t> ind = {0,2,3};
    std::vector<F> val = {-1, 1, -1};
    CT y(ind, val);
    CT x;
    y.print_row();

    auto start = std::chrono::high_resolution_clock::now();
    for (size_t k = 0; k < N; k++) {
        x = solve_U(I, y);
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "Usolve ";
    std::cout << "microseconds: " << duration.count() << std::endl;
    x.print_row();

    start = std::chrono::high_resolution_clock::now();
    for (size_t k = 0; k < N; k++) {
        x = solve_L(I, y);
    }
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "Lsolve ";
    std::cout << "microseconds: " << duration.count() << std::endl;
    x.print_row();
    
    return;
}

int main() {

    std::cout << "SetVector: " << std::endl;
    time_solve<SetVector<F>>(NITER);
    time_solve<SetVector<F>>(NITER);

    std::cout << "SparseVector: " << std::endl;
    time_solve<SparseVector<F>>(NITER);
    time_solve<SparseVector<F>>(NITER);


    return 0;
}
